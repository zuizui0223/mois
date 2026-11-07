#!/usr/bin/env python3
"""True-year 2019,2021,2022,2023 VIIRS source sampling of worldwide MOIS land sites.

This extracts ONLY VIIRS, NOT population. Never treat a VIIRS-only row as
fully adjusted environmental support. Original source pixel/month metadata,
cloudfree band, and raw provider byte SHA are checked and stored.
"""
from __future__ import annotations
import argparse,calendar,concurrent.futures,csv,datetime as dt,hashlib,io,json,math,statistics
from pathlib import Path
import requests
from scripts.analysis.mois_global_multipoint_numeric_20261107 import (
    extract_samples,VIIRS,MONTHS)
SOURCE_SHA={
 2019:"01f4e2d2e1c8bc9a1392f4143f693563bfa8ea1a0f65ba800e1caa450803831c",
 2021:"f0d7cb134c09c74cccfb69199396ac30bc5b663d39dc5369d6e296d8c3a1cdca",
 2022:"594fe8e729cf54aa44859a841b5d3e49e6a66fb85b217da8aff0160044e14df5",
 2023:"e2825386f45a77da6e3dacdf714fc940cf9ab9cc5174671672b4d57466d04fba"}
SITES={2019:27020,2021:45966,2022:52075,2023:59862}

def require(ok,msg):
    if not bool(ok):raise ValueError("MOIS_YEAR_VIIRS: "+msg)

def digest(data):return hashlib.sha256(data).hexdigest()

def verify_source(rows,year):
    require(year in SITES and len(rows)==SITES[year],"Wrong original year/cohort")
    require(len({r["target_id"] for r in rows})==len(rows),"Duplicate original target ID")
    require(all(r["target_id"].startswith(str(year)+"|") for r in rows),
            "Original photo ID from a different observation year")
    require(all(r["land_class"] in ("island","mainland") for r in rows),
            "Unknown GSHHG landmass wrongly considered island/mainland")
    require(len({(float(r["query_longitude"]),float(r["query_latitude"])) for r in rows})==len(rows),
            "Duplicate original raster request coordinates")
    require(all(-180<=float(r["query_longitude"])<=180 and
                -90<=float(r["query_latitude"])<=90 for r in rows),
            "Invalid original geographic position")
    return rows

def plans(year,rows):
    require(year in SITES,"Unsupported source-year VIIRS")
    require(1<=len(rows)<=SITES[year],"Source subset too large")
    return [(month,start,min(512,len(rows)-start))
            for month in MONTHS for start in range(0,len(rows),512)]

def params(year,month,rows):
    require(year in SITES and month in MONTHS and 1<=len(rows)<=512,
            "Invalid exact year/month/point count")
    start=dt.datetime(year,month,1,tzinfo=dt.timezone.utc)
    end=dt.datetime(year,month,calendar.monthrange(year,month)[1],23,59,59,
                    tzinfo=dt.timezone.utc)
    geom={"points":[[float(r["query_longitude"]),float(r["query_latitude"])]
                     for r in rows],"spatialReference":{"wkid":4326}}
    return {"f":"json","geometryType":"esriGeometryMultipoint",
        "geometry":json.dumps(geom,separators=(',',':')),
        "returnFirstValueOnly":"true","outFields":"*","interpolation":"RSP_NearestNeighbor",
        "time":f"{int(start.timestamp()*1000)},{int(end.timestamp()*1000)}"}

def original_2022_october_provider_error(year,month,obj):
    """Never convert this provider 400 transformation error to a zero light pixel.

    The original source has 200 HTTP with per-provider raster error 400;
    its missing numerical values remain unavailable even when 3 other months work.
    """
    if year!=2022 or month!=10:return False
    err=obj.get("error") if isinstance(obj,dict) else None
    if not isinstance(err,dict):return False
    return (err.get("code")==400 and
            err.get("message")=="Invalid or missing input parameters." and
            all("Transformation is unavailable" in str(x) for x in (err.get("details") or [])))

def one(plan,year,rows,out,offline):
    month,start,n=plan
    subset=rows[start:start+n]
    filename=f"VIIRS_{year}_{month:02d}_{start:05d}_{n:03d}.response"
    path=out/"raw"/filename
    if offline or path.is_file():
        require(path.is_file(),"Offline replay lacks original VIIRS source bytes "+filename)
        raw=path.read_bytes()
        reused=True
    else:
        try:
            r=requests.get(VIIRS,params=params(year,month,subset),
                timeout=90,headers={"User-Agent":"MOIS Global Photo ALAN real-year source extraction /1.0"})
        except requests.RequestException as e:raise RuntimeError("VIIRS real source transport failure") from e
        require(r.status_code==200,"provider HTTP "+str(r.status_code)+" "+filename)
        raw=r.content
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_bytes(raw);reused=False
    try:obj=json.loads(raw)
    except ValueError as e:raise ValueError("Original provider bytes not JSON "+filename) from e
    if original_2022_october_provider_error(year,month,obj):
        vals=[{"status":"PROVIDER_RASTER_400","value":None,"cloud":None}
              for _ in subset]
    else:
        require(not obj.get("error"),"Source provider error "+filename)
        vals=extract_samples("viirs",month,subset,obj,year=year)
    require(len(vals)==n,"Source VIIRS dropped original locations")
    return plan,vals,{"file":filename,"raw_SHA256":digest(raw),
                     "new_request":not reused,"source_year":year,"month":month,
                     "source_n":len(subset),
                     "accepted_months":sum(x["status"]=="OK" for x in vals),
                     "cloudfree_failure":sum(x["status"]=="LOW_CLOUD_COUNT" for x in vals)}

def run(year,targets,out,offline=False,workers=4):
    require(1<=workers<=4,"More than four concurrent provider requests")
    p=Path(targets);b=p.read_bytes()
    require(digest(b)==SOURCE_SHA[year],"Original year source knownland manifest SHA changed")
    rows=verify_source(list(csv.DictReader(io.StringIO(b.decode()))),year)
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    results=[dict(x) for x in rows]
    audit=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        for plan,vs,meta in pool.map(lambda a:one(a,year,rows,out,offline),plans(year,rows)):
            month,start,n=plan
            for i,val in enumerate(vs,start=start):
                r=results[i];prefix=f"viirs_{month}"
                r[prefix+"_status"]=val["status"]
                r[prefix+"_value"]=val["value"]
                r[prefix+"_cloudfree"]=val["cloud"]
            audit.append(meta)
    for r in results:
        valid=[r[f"viirs_{m}_value"] for m in MONTHS
               if r[f"viirs_{m}_status"]=="OK"]
        r["viirs_valid_months"]=len(valid)
        r["viirs_median"]=statistics.median(valid) if len(valid)>=3 else None
        r["viirs_available"]=len(valid)>=3
    f=io.StringIO(newline='')
    writer=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for row in results for k in row)),
         lineterminator='\n')
    writer.writeheader();writer.writerows(results)
    data=f.getvalue().encode()
    (out/("replay_values.csv" if offline else "values.csv")).write_bytes(data)
    receipt={
       "schema":"MOIS_true_2019_2023_photo_year_VIIRS_raw_month_source_v1",
       "source_year":year,
       "source_original_photo_landmass_site_manifest_SHA256":digest(b),
       "all_original_year_photo_cells":len(rows),
       "source_island_sites":sum(r["land_class"]=="island" for r in results),
       "source_mainland_sites":sum(r["land_class"]=="mainland" for r in results),
       "VIIRS_months":[1,4,7,10],
       "source_VIIRS_at_least3months_sites":sum(r["viirs_available"] for r in results),
       "source_provider_errors_retained_as_missing":sum(
           r.get("viirs_10_status")=="PROVIDER_RASTER_400" for r in results),
       "2022_oct_source_error_not_equal_to_zero_or_official_monthly_data":True,
       "source_VIIRS_fewer3months_sites":sum(not r["viirs_available"] for r in results),
       "valid_monthly_VIIRS_values":sum(r[f"viirs_{m}_status"]=="OK"
                                      for r in results for m in MONTHS),
       "new_provider_requests":sum(x["new_request"] for x in audit),
       "original_provider_responses_reused":sum(not x["new_request"] for x in audit),
       "raw_original_response_SHA256":{x["file"]:x["raw_SHA256"] for x in audit},
       "values_SHA256":digest(data),
       "offline":offline,
       "population_source_not_fetched_or_imputed":True,
       "VIIRS_year_matches_original_photo_year":True,
       "VIIRS_four_month_raw_not_official_annual":True,
       "no_island_area_or_latitude_exclusions":True,
       "new_original_photo_downloads_or_ROI_reinferences":0,
       "colour_unknown_not_imputed":True,
       "no_ecological_COLOUR_ALAN_model_fitted":True}
    (out/("replay_receipt.json" if offline else "receipt.json")).write_text(
        json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print("MOIS_YEAR_SOURCE_VIIRS_ACTUAL",json.dumps({
      k:receipt[k] for k in ("source_year","all_original_year_photo_cells",
          "source_island_sites","source_mainland_sites",
          "source_VIIRS_at_least3months_sites","source_VIIRS_fewer3months_sites",
          "valid_monthly_VIIRS_values","new_provider_requests","original_provider_responses_reused")},
          sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--year",type=int,required=True,choices=sorted(SITES))
    p.add_argument("--targets",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    p.add_argument("--offline",action="store_true")
    p.add_argument("--workers",type=int,default=4)
    a=p.parse_args()
    run(a.year,a.targets,a.out,a.offline,a.workers)
