#!/usr/bin/env python3
"""Retrieve real source-month TerraClimate v1.1 1991-2025 at MOIS photo grids.

Only valid ORIGINAL CF packing metadata (.das scale_factor/add_offset) may
transform NCSS packed integer CSVs into physical units. NCSS point-NetCDF
also returns packed values without attributes; never trust its unit label.
Neither original photo label nor an unchecked flowering status can inform
climate acquisition. No 2026 climate interpolation or ecological effect.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,math,re,time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import date
from pathlib import Path
from typing import Any
import requests

VARS=("tmin","tmax","ppt","srad","def","vpd")
UNITS={"tmin":"degC","tmax":"degC","ppt":"mm","srad":"W/m2","def":"mm","vpd":"kPa"}
LIMITS={"tmin":(-100.,65.),"tmax":(-100.,75.),"ppt":(0.,3000.),
        "srad":(0.,450.),"def":(0.,3000.),"vpd":(0.,20.)}
ENDPOINT="https://thredds.northwestknowledge.net/thredds/ncss/agg_terraclimate_{var}_1950_CurrentYear_GLOBE.nc"
ATTR="https://thredds.northwestknowledge.net/thredds/dodsC/agg_terraclimate_{var}_1950_CurrentYear_GLOBE.nc.das"
HEADERS={"User-Agent":"Mois-2025-observed-flower-season-climate-CF-metadata/1.0"}
CLIMATE_COLS=("year_month","tc_lat_idx","tc_lon_idx","tc_lat_center","tc_lon_center",
              "n_photo_observations","tmin_C","tmax_C","tmean_C","ppt_mm","srad_W_m2","def_mm","vpd_kPa",
              "tmean_C_anomaly_1991_2020_same_month","ppt_mm_anomaly_1991_2020_same_month",
              "srad_W_m2_anomaly_1991_2020_same_month",
              "def_mm_anomaly_1991_2020_same_month","vpd_kPa_anomaly_1991_2020_same_month",
              "climate_is_actual_2025_not_a_historical_normal",
              "flowering_organ_is_botanically_unreviewed")
SEASON_COLS=("local_season","n_independent_original_observations","n_different_source_terraclimate_grid_months",
             "2025_mean_tmean_C","2025_mean_ppt_mm_per_month","2025_mean_srad_W_m2",
             "2025_mean_def_mm_per_month","2025_mean_vpd_kPa",
             "2025_tmean_anomaly_C_1991_2020","2025_ppt_anomaly_mm_1991_2020",
             "2025_srad_anomaly_W_m2_1991_2020",
             "2025_def_anomaly_mm_1991_2020","2025_vpd_anomaly_kPa_1991_2020",
             "not_a_white_flower_climate_effect")


def sha(raw:bytes)->str:return hashlib.sha256(raw).hexdigest()


def metadata(var:str,das_text:str)->dict:
    if var not in VARS:raise ValueError("Unexpected TerraClimate variable")
    if not re.search(r'\bversion\s+"V1\.1"\s*;',das_text,re.I):
        raise ValueError("TerraClimate v1.1 metadata missing; no v1.0 mixture")
    matched=re.search(r'(?m)^\s*'+re.escape(var)+r'\s*\{\s*(.*?)^\s*\}',das_text,re.S|re.M)
    if not matched:raise ValueError("Original TerraClimate variable CF metadata missing: "+var)
    body=matched.group(1)
    def val(key):
        m=re.search(r'(?mi)^\s*(?:Float64|Float32|Int16|Int32|Double|Short)\s+'
                    +re.escape(key)+r'\s+(-?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)\s*;',body)
        if not m:raise ValueError("Official original NetCDF CF "+key+" not documented for "+var)
        return float(m.group(1))
    scale,offset,fill=val("scale_factor"),val("add_offset"),val("_FillValue")
    if not math.isfinite(scale) or not 0<scale<=100 or not math.isfinite(offset) or not math.isfinite(fill):
        raise ValueError("Unsafe unverified TerraClimate packing metadata")
    unit=re.search(r'(?mi)^\s*String\s+units\s+"([^"]+)"\s*;',body)
    return {"var":var,"scale":scale,"offset":offset,"fill":fill,
            "source_unit":unit.group(1) if unit else None,
            "dataset_version":"V1.1","official_das_sha256":sha(das_text.encode("utf-8"))}


def unpack(var:str,raw:float,m:dict)->float:
    if not math.isfinite(raw) or raw==m["fill"]:raise ValueError("Unresolved TerraClimate source fill value")
    value=raw*m["scale"]+m["offset"]
    if not LIMITS[var][0]<=value<=LIMITS[var][1]:
        raise ValueError("Packed NCSS value not converted to physically plausible units: "+var)
    return value


def csv_series(var:str,raw:bytes,m:dict)->dict[str,float]:
    text=raw.decode("utf-8")
    with io.StringIO(text) as file:
        reader=csv.DictReader(file)
        fields=reader.fieldnames or []
        value=[x for x in fields if x==var or x.startswith(var+"[")]
        if not any(x.lower()=="time" for x in fields) or len(value)!=1:
            raise ValueError("Original provider NCSS CSV has no unique time + requested variable")
        out={}
        for row in reader:
            stamp=row["time"].strip()
            if len(stamp)<10:raise ValueError("Invalid original TerraClimate source monthly timestamp")
            ym=stamp[:7]
            date.fromisoformat(ym+"-01")
            if ym in out:raise ValueError("NCSS duplicate monthly climate timestamp")
            if not "1991-01"<=ym<="2025-12":
                raise ValueError("Unexpected year outside frozen physical month climate baseline")
            out[ym]=unpack(var,float(row[value[0]]),m)
    # Exactly 35 years, including the complete same-month 1991–2020 normal and 2025
    years=35
    if len(out)!=years*12 or any(f"{y:04d}-{mo:02d}" not in out
                                   for y in range(1991,2026) for mo in range(1,13)):
        raise ValueError(f"Incomplete 1991–2025 original monthly climate at site: {var}, {len(out)} rows")
    return out


def get(session:requests.Session,url:str,*,params:dict|None=None,timeout=(8,35))->bytes:
    error=None
    for attempt in range(3):
        try:
            response=session.get(url,params=params,timeout=timeout,headers=HEADERS)
            response.raise_for_status()
            if not response.content or len(response.content)>3_000_000:
                raise ValueError("Unexpected empty or overlarge NCSS point time series")
            return response.content
        except (requests.RequestException,ValueError) as exc:
            error=exc
            if attempt<2:time.sleep(attempt+1)
    raise RuntimeError("Official monthly TerraClimate source request failed") from error


def run(keys:Path,origin_receipt:Path,outdir:Path,*,workers:int=5,
        _fake_request=None)->dict:
    if outdir.exists() and any(outdir.iterdir()):raise ValueError("Never replace original verified climate outputs")
    doc=json.loads(origin_receipt.read_text(encoding="utf-8"))
    if (doc.get("status")!="REAL_DATE_AND_ROI_ATTEMPTS_NO_BOTANICAL_WHITE_CLIMATE_EFFECT" or
        doc.get("original_month_climate_grid_csv_sha256")!=sha(keys.read_bytes()) or
        doc.get("n_exact_2025_grid_and_month_climate_requests")!=40 or
        doc.get("n_original_source_grid_cells")!=21 or
        doc.get("open_flower_botanically_validated") is not False):
        raise ValueError("Mois climate keys have no authenticated original FCP source/ROI receipt")
    if not 1<=workers<=6:raise ValueError("Limit original provider concurrent requests")
    with keys.open(newline="",encoding="utf-8") as f:
        reader=csv.DictReader(f)
        wanted={"year_month","tc_lat_idx","tc_lon_idx","tc_lat_center","tc_lon_center",
                "n_photo_observations","anthesis_is_unverified"}
        if not wanted.issubset(reader.fieldnames or []):raise ValueError("Climate cell-month schema drift")
        source=list(reader)
    if len(source)!=40:raise ValueError("Missing original 2025 40 photo grid-month units")
    cells={}
    seen=set();counts_byseason=defaultdict(int)
    for r in source:
        ym=r["year_month"]
        if ym not in {f"2025-{m:02d}" for m in range(6,12)}:
            raise ValueError("Temporal source has non-frozen season-year")
        i,j=int(r["tc_lat_idx"]),int(r["tc_lon_idx"])
        key=(i,j,ym)
        if key in seen:raise ValueError("Repeated climate-month source key")
        seen.add(key)
        lat,lon=float(r["tc_lat_center"]),float(r["tc_lon_center"])
        if not (abs(lat-(90-(i+.5)/24))<=1e-6 and
                abs(lon-(-180+(j+.5)/24))<=1e-6 and
                r["anthesis_is_unverified"]=="True"):
            raise ValueError("Not original georeferenced photo-grid key")
        cells[(i,j)]=(lat,lon)
        counts_byseason["summer" if ym[-2:] in ("06","07","08") else "autumn"]+=int(r["n_photo_observations"])
    if len(cells)!=21 or dict(counts_byseason)!={"autumn":32,"summer":32}:
        raise ValueError("Climate year-photo season observer counts do not reconcile")
    outdir.mkdir(parents=True,exist_ok=True)
    cache=outdir/"source_ncss_cf_packing_and_request_receipts"
    cache.mkdir()
    attr={}
    with requests.Session() as session:
        for var in VARS:
            das=get(session,ATTR.format(var=var))
            attr[var]=metadata(var,das.decode("utf-8"))
            (cache/f"{var}_original_v11_cf_metadata.das").write_bytes(das)
    # Only the six provider CF metadata-scale definitions are accepted.
    if not math.isclose(attr["tmin"]["scale"],.01,rel_tol=1e-12) or not math.isclose(attr["tmin"]["offset"],-99.0,abs_tol=1e-12):
        raise ValueError("Official tmin v1.1 CF scale has changed; abort historical run")
    jobs=[(i,j,var,lat,lon) for (i,j),(lat,lon) in sorted(cells.items()) for var in VARS]
    series={}
    failures=[]
    def one(job):
        i,j,var,lat,lon=job
        params={"var":var,"latitude":f"{lat:.8f}","longitude":f"{lon:.8f}",
             "time_start":"1991-01-01T00:00:00Z",
             "time_end":"2025-12-31T23:59:59Z","accept":"csv"}
        with requests.Session() as s:
            raw=get(s,ENDPOINT.format(var=var),params=params)
        parsed=csv_series(var,raw,attr[var])
        return (i,j,var,raw,parsed)
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures={executor.submit(one,j):j for j in jobs}
        for n,future in enumerate(as_completed(futures),1):
            job=futures[future]
            try:
                i,j,var,raw,parsed=future.result()
                series[(i,j,var)]=parsed
                (cache/f"{var}_{i}_{j}_original_1991_2025_packed.csv").write_bytes(raw)
            except Exception as exc:
                failures.append({"job":job[:3],"reason":repr(exc)[:200]})
            if n%20==0 or n==len(jobs):
                print("MOIS_VERIFIED_TERRACLIMATE_RAW_CF_PROGRESS",n,"of",len(jobs),"failures",len(failures),flush=True)
    if failures:
        (outdir/"FAILURE_missing_actual_original_monthly_provider_values.json").write_text(
            json.dumps({"status":"NO_COMPLETE_CLIMATE_DATA_FIT","failures":failures[:40]},indent=2))
        raise ValueError("One or more authenticated TerraClimate v1.1 time series missing; no effects admitted")

    # Full physical monthly values + fixed-location 1991–2020 SAME-calendar-month anomalies.
    rows=[]
    for r in sorted(source,key=lambda r:(r["year_month"],int(r["tc_lat_idx"]),int(r["tc_lon_idx"]))):
        ym=r["year_month"];m=int(ym[-2:])
        i,j=int(r["tc_lat_idx"]),int(r["tc_lon_idx"])
        now={v:series[i,j,v][ym] for v in VARS}
        base={v:sum(series[i,j,v][f"{year:04d}-{m:02d}"] for year in range(1991,2021))/30 for v in VARS}
        if now["tmin"]>now["tmax"]:raise ValueError("Physical source tmin > tmax after CF decoding")
        tmean=(now["tmin"]+now["tmax"])/2
        btmean=(base["tmin"]+base["tmax"])/2
        row=dict(zip(CLIMATE_COLS,(
            ym,i,j,r["tc_lat_center"],r["tc_lon_center"],r["n_photo_observations"],
            *(f"{now[v]:.6f}" for v in VARS[:2]),f"{tmean:.6f}",
            *(f"{now[v]:.6f}" for v in VARS[2:]),
            f"{tmean-btmean:.6f}",f"{now['ppt']-base['ppt']:.6f}",
            f"{now['srad']-base['srad']:.6f}",f"{now['def']-base['def']:.6f}",
            f"{now['vpd']-base['vpd']:.6f}",True,True
        )))
        rows.append(row)
    observed=outdir/"authentic_2025_asclepias_terraclimate_physical_month_grid.csv"
    with observed.open("x",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=CLIMATE_COLS,lineterminator="\n");w.writeheader();w.writerows(rows)
    summaries=[]
    for season in ("summer","autumn"):
        relevant=[r for r in rows if ((season=="summer" and r["year_month"][-2:] in ("06","07","08")) or
                   (season=="autumn" and r["year_month"][-2:] in ("09","10","11")))]
        n=sum(int(r["n_photo_observations"]) for r in relevant)
        if n!=32:raise ValueError("Only real 32 original unique source observers/season valid")
        def avg(col):
            return sum(int(r["n_photo_observations"])*float(r[col]) for r in relevant)/n
        summaries.append(dict(zip(SEASON_COLS,(
            season,n,len(relevant),
            *[f"{avg(col):.6f}" for col in (
                "tmean_C","ppt_mm","srad_W_m2","def_mm","vpd_kPa",
                "tmean_C_anomaly_1991_2020_same_month",
                "ppt_mm_anomaly_1991_2020_same_month",
                "srad_W_m2_anomaly_1991_2020_same_month",
                "def_mm_anomaly_1991_2020_same_month",
                "vpd_kPa_anomaly_1991_2020_same_month")],
            True))))
    seasonal=outdir/"2025_source_observation_weighted_summer_autumn_climate_NOT_flower_effect.csv"
    with seasonal.open("x",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=SEASON_COLS,lineterminator="\n");w.writeheader();w.writerows(summaries)
    receipt={
      "schema":"mois_original_2025_terra_climate_v11_verified_cf_original_observation_month_v1",
      "status":"AUTHENTIC_2025_MONTHLY_CLIMATE_VALUES_NOT_WHITE_FLOWER_ASSOCIATION",
      "provider":"University of Idaho TerraClimate v1.1 THREDDS aggregated 1950-current NCSS / OPeNDAP attributes",
      "input_original_2025_photo_climate_request_sha256":sha(keys.read_bytes()),
      "input_original_2025_roi_source_receipt_sha256":sha(origin_receipt.read_bytes()),
      "baseline_period_years":[1991,2020],"observation_year":2025,
      "genuine_original_monthly_photo_observations":64,"independent_source_observers":64,
      "original_source_photo_grid_months":len(rows),"unique_original_source_grid_cells":len(cells),
      "provider_request_count":len(jobs),
      "physical_variables_and_units":UNITS,
      "official_cf_scaling":attr,
      "decoded_month_grid_climate_sha256":sha(observed.read_bytes()),
      "source_observation_weighted_seasonal_climate_sha256":sha(seasonal.read_bytes()),
      "monthly_climate_is_actual_2025_not_1991_2020_normal":True,
      "climatology_only_used_for_within_cell_same_month_anomalies":True,
      "never_assigned_future_month_or_year_from_a_different_year":True,
      "no_botanically_confirmed_open_flower_colour_data":True,
      "no_white_fraction_vs_climate_regression_or_causal_inference":True,
      "source_nc_sst_checksums":{
        p.name:sha(p.read_bytes()) for p in sorted(cache.iterdir()) if p.is_file()
      },
    }
    (outdir/"verified_2025_terraclimate_source_month_receipt.json").write_text(
       json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print("MOIS_VERIFIED_2025_REAL_TERRACLIMATE_OBSERVED_CLIMATE_ONLY",json.dumps({
        "n_grid_months":len(rows),"n_cells":len(cells),
        "summer":summaries[0],"autumn":summaries[1]}),flush=True)
    return receipt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-2025-terraclimate-grid-keys",type=Path,required=True)
    p.add_argument("--source-2025-roi-receipt",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    p.add_argument("--workers",type=int,default=5)
    a=p.parse_args()
    print(json.dumps(run(a.source_2025_terraclimate_grid_keys,
                         a.source_2025_roi_receipt,a.outdir,workers=a.workers),
                     sort_keys=True))


if __name__=="__main__":main()
