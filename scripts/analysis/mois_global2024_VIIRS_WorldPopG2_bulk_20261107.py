#!/usr/bin/env python3
"""2024 worldwide original VIIRS raw-month radiation + WorldPop Global2 population.

Global2 R2025A v1 raster measures PEOPLE per 30-arcsecond pixel, NOT
density. Convert each pixel using pyproj's WGS84 geodesic polygon area.
No Global1 2020 values, no floral selection, no minimum island area.
"""
from __future__ import annotations
import argparse,calendar,concurrent.futures,csv,datetime as dt,hashlib,io,json,math,statistics
from pathlib import Path
import numpy as np,pyproj,rasterio,requests
from scripts.analysis.mois_global_multipoint_numeric_20261107 import extract_samples,VIIRS,MONTHS

YEAR=2024
RASTER_NAME="global_pop_2024_CN_1km_R2025A_v1.tif"
RASTER_SHA="3539cbac76092c1e8c8d45426fc0dfb06067a1e3f28bfbee993dba21aa03194e"
RASTER_BYTES=289380599
SELECT_SHA="541b6134090b61a285ffa9c45ee3df8a382e66ae5628c394ca0c4d47af9207a0"
REST_SHA="a792602a05735c5049c960a80fc6a75f4a1241c0a3abaa9769f452330559bf2d"

def check(x,msg):
    if not bool(x):raise ValueError("MOIS_2024_WP_G2: "+msg)
def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(4<<20),b""):h.update(b)
    return h.hexdigest()
def time_range(m):
    a=dt.datetime(YEAR,m,1,tzinfo=dt.timezone.utc)
    b=dt.datetime(YEAR,m,calendar.monthrange(YEAR,m)[1],23,59,59,tzinfo=dt.timezone.utc)
    return f'{int(a.timestamp()*1000)},{int(b.timestamp()*1000)}'
def get_samples_params(rows,month):
    check(month in MONTHS and 1<=len(rows)<=512,"Invalid same-year VIIRS source month or batch")
    pts=[[float(r["query_longitude"]),float(r["query_latitude"])] for r in rows]
    geometry={"points":pts,"spatialReference":{"wkid":4326}}
    return {"f":"json","geometryType":"esriGeometryMultipoint",
       "geometry":json.dumps(geometry,separators=(",",":")),
       "returnFirstValueOnly":"true","outFields":"*",
       "interpolation":"RSP_NearestNeighbor","time":time_range(month)}
def validate_targets(rows,expected):
    check(len(rows)==expected and len({r["target_id"] for r in rows})==len(rows),"2024 ID cohort changed")
    check(all(str(r["target_id"]).startswith("2024|") and
              r["land_class"] in ("island","mainland") for r in rows),
          "Wrong photo year or mapped landmass")
    check(len({(float(r["query_longitude"]),float(r["query_latitude"])) for r in rows})==len(rows),
          "Duplicate source raster target coordinates")
    check(all(-180<=float(r["query_longitude"])<=180 and
              -90<=float(r["query_latitude"])<=90 for r in rows),"Invalid coordinates")
    return rows
def sample_population(raster,receipt,rows):
    path=Path(raster)
    check(path.name==RASTER_NAME and path.stat().st_size==RASTER_BYTES,
          "Not the official 2024 Global2 raster")
    check(digest(path)==RASTER_SHA,"Official 2024 Global2 full original bytes changed")
    obj=json.loads(Path(receipt).read_text())
    check(obj.get("original_file_SHA256")==RASTER_SHA and obj.get("source_year")==2024 and
          obj.get("release")=="R2025A v1",
          "Public original raster provenance not signed")
    geod=pyproj.Geod(ellps="WGS84")
    sampled=[];cached_area={}
    with rasterio.open(path) as ds:
        check(ds.crs.to_epsg()==4326 and ds.count==1 and
              ds.width==43200 and ds.height==17280,
              "Wrong release/global-grid EPSG or raster shape")
        check(abs(ds.res[0]-1/120)<1e-6 and abs(ds.res[1]-1/120)<1e-6,
              "WorldPop global raster 30 arcseconds grid changed")
        pts=[(float(r["query_longitude"]),float(r["query_latitude"])) for r in rows]
        for (lon,lat),data in zip(pts,ds.sample(pts,masked=True)):
            valid=(ds.bounds.left<=lon<ds.bounds.right and
                   ds.bounds.bottom<lat<=ds.bounds.top and
                   not np.ma.is_masked(data[0]))
            if not valid:
                sampled.append({"status":"RASTER_OUTSIDE_OR_NODATA",
                                "population_count":None,"density":None,"area_km2":None})
                continue
            value=float(data[0])
            if not math.isfinite(value) or value<0 or (ds.nodata is not None and value==ds.nodata):
                sampled.append({"status":"RASTER_OUTSIDE_OR_NODATA",
                                "population_count":None,"density":None,"area_km2":None})
                continue
            row,col=ds.index(lon,lat)
            check(0<=row<ds.height and 0<=col<ds.width,"Raster grid index out of bounds")
            if row not in cached_area:
                west,south,east,north=rasterio.windows.bounds(
                    rasterio.windows.Window(col,row,1,1),ds.transform)
                area,_=geod.polygon_area_perimeter(
                   [west,east,east,west],[south,south,north,north])
                cached_area[row]=abs(area)/1e6
            area=cached_area[row]
            check(math.isfinite(area) and area>0,"Invalid WGS84 raster pixel area")
            sampled.append({"status":"OK","population_count":value,
                            "density":value/area,"area_km2":area})
    check(len(sampled)==len(rows),"Lost original photo locations during population lookup")
    return sampled,{"source_release":"WorldPop Global2 R2025A v1",
       "source_raster_year":2024,"source_file_sha256":RASTER_SHA,
       "original_file_bytes":RASTER_BYTES,
       "source_units":"persons per 30 arcsecond pixel, NOT persons/km2",
       "density_conversion":"WGS84 Geod pixel polygon area km2",
       "number_unique_latitude_pixel_rows":len(cached_area)}
def viirs_task(plan,rows,out,offline):
    month,start,n=plan
    subset=rows[start:start+n]
    filename=f"VIIRS_2024_{month:02d}_{start:05d}_{n:03d}.response"
    path=out/"raw"/filename
    if offline or path.is_file():
        check(path.is_file(),"Original monthly source bytes missing on offline replay")
        raw=path.read_bytes();reused=True
    else:
        params=get_samples_params(subset,month)
        try:
            r=requests.get(VIIRS,params=params,timeout=90,
               headers={"User-Agent":"MOIS 2024 original all-land raw radiance source sampling"})
        except requests.RequestException as e:
            raise RuntimeError("VIIRS source transport failure") from e
        check(r.status_code==200,f"VIIRS source HTTP {r.status_code}")
        raw=r.content;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(raw)
        reused=False
    try:obj=json.loads(raw)
    except ValueError as e:raise ValueError("Monthly original raw VIIRS provider not JSON") from e
    check(not obj.get("error"),"VIIRS provider returned error "+str(obj.get("error")))
    values=extract_samples("viirs",month,subset,obj,year=YEAR)
    check(len(values)==n,"VIIRS provider dropped original photo-cell locations")
    return plan,values,{"filename":filename,"sha256":digest(path),"month":month,
        "number_sites":n,"previous_original_bytes_reused":reused,
        "cloudfree_and_year_month_verified":True}
def run(targets,raster,raster_receipt,out,expected_count=4096,
        source_sha=SELECT_SHA,offline=False,workers=4):
    check(1<=workers<=4,"Maximum four API requests in flight")
    source=Path(targets)
    check(digest(source)==source_sha,"Photo-selected target source SHA changed")
    rows=list(csv.DictReader(io.StringIO(source.read_text())))
    validate_targets(rows,expected_count)
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    pops,pop_meta=sample_population(raster,raster_receipt,rows)
    result=[dict(x) for x in rows]
    for rec,v in zip(result,pops):
        rec.update({"population_source_status":v["status"],
            "population_count_per_pixel":v["population_count"],
            "population_geodesic_pixel_km2":v["area_km2"],
            "population_density":v["density"]})
    plans=[(m,start,min(512,len(rows)-start))
        for m in MONTHS for start in range(0,len(rows),512)]
    logs=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        for plan,values,meta in executor.map(lambda x:viirs_task(x,rows,out,offline),plans):
            m,start,n=plan
            for i,v in enumerate(values,start=start):
                rec=result[i]
                rec[f"viirs_{m}_status"]=v["status"]
                rec[f"viirs_{m}_value"]=v["value"]
                rec[f"viirs_{m}_cloudfree"]=v["cloud"]
            logs.append(meta)
    for rec in result:
        valid=[rec[f"viirs_{m}_value"] for m in MONTHS
               if rec[f"viirs_{m}_status"]=="OK"]
        rec["viirs_valid_months"]=len(valid)
        rec["viirs_median"]=statistics.median(valid) if len(valid)>=3 else None
        rec["complete"]=len(valid)>=3 and rec["population_source_status"]=="OK"
    buf=io.StringIO(newline="")
    writer=csv.DictWriter(buf,fieldnames=list(dict.fromkeys(k for r in result for k in r)),
                          lineterminator="\n")
    writer.writeheader();writer.writerows(result)
    data=buf.getvalue().encode()
    (out/("replay_values.csv" if offline else "values.csv")).write_bytes(data)
    receipt={
        "schema":"MOIS_2024_real_sameyear_Global2_and_VIIRS_v1",
        "original_source_photo_grid_sha256":digest(source),
        "source_photo_year":2024,"VIIRS_year":2024,"population_raster_year":2024,
        "original_target_count":len(rows),
        "completed_targets":sum(r["complete"] for r in result),
        "incomplete_targets":sum(not r["complete"] for r in result),
        "attempted_by_land":{x:sum(r["land_class"]==x for r in result)
                             for x in ("island","mainland")},
        "complete_by_land":{x:sum(r["complete"] and r["land_class"]==x for r in result)
                            for x in ("island","mainland")},
        "valid_VIIRS_month_count":sum(
            r[f"viirs_{m}_status"]=="OK" for r in result for m in MONTHS),
        "valid_WorldPop_2024_population_pixels":sum(
            r["population_source_status"]=="OK" for r in result),
        "new_VIIRS_provider_requests":sum(not l["previous_original_bytes_reused"] for l in logs),
        "saved_original_VIIRS_replay_requests":sum(l["previous_original_bytes_reused"] for l in logs),
        "original_VIIRS_response_SHA256":{l["filename"]:l["sha256"] for l in logs},
        "original_WorldPop_source":pop_meta,
        "output_values_SHA256":hashlib.sha256(data).hexdigest(),
        "replayed_without_network":bool(offline),
        "population_2020_Global1_mixed_in":False,
        "no_island_area_or_latitude_cutoff":True,
        "zero_population_preserved_as_valid":True,
        "source_unknown_to_nonwhite_imputation":False,
        "new_photos_or_ROI_reinference":0,
        "VIIRS_fourmonth_raw_not_official_annual_product":True,
        "ecological_association_not_fitted":True}
    (out/("replay_receipt.json" if offline else "receipt.json")).write_text(
        json.dumps(receipt,sort_keys=True,indent=2)+"\n")
    print("MOIS_2024_ACTUAL_VIIRS_WORLDPOP_GLOBAL2_SAMEYEAR",
       json.dumps({k:receipt[k] for k in ("original_target_count","completed_targets",
           "incomplete_targets","complete_by_land","valid_VIIRS_month_count",
           "valid_WorldPop_2024_population_pixels",
           "new_VIIRS_provider_requests","saved_original_VIIRS_replay_requests")},
           sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--targets",type=Path,required=True)
    p.add_argument("--worldpop-tiff",type=Path,required=True)
    p.add_argument("--worldpop-receipt",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--expected-count",type=int,default=4096)
    p.add_argument("--source-sha",default=SELECT_SHA)
    p.add_argument("--offline",action="store_true")
    p.add_argument("--workers",type=int,default=4)
    args=p.parse_args()
    run(args.targets,args.worldpop_tiff,args.worldpop_receipt,args.out,args.expected_count,
        args.source_sha,args.offline,args.workers)
