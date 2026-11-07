#!/usr/bin/env python3
"""WorldPop 2024 R2024B UC: separate, same-year global population covariate.

NOT interchangeable with the already extracted WorldPop R2025A CN.
Do not convert existing CN-NODATA sites to zeros. Separate full-cohort
R2024B-only population estimates and published-source receipts.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from pathlib import Path
import numpy as np,pandas as pd,rasterio
from rasterio.windows import Window,bounds
from pyproj import Geod

ORIGINAL_2024_SHA="6e0a077d496e90848d05d043bffae57e4772519f4572a457221c5ef4da0fb813"
CN_SOURCE_SHA="3539cbac76092c1e8c8d45426fc0dfb06067a1e3f28bfbee993dba21aa03194e"
UC_FILENAME="global_pop_2024_UC_1km_R2024B_v1.tif"
UC_EXPECTED_BYTES=868021365
UC_URL="https://worldpop-public-data.soton.ac.uk/GIS/Population/Global_2015_2030/R2024B/2024/0_Mosaicked/v1/1km/unconstrained/"+UC_FILENAME

def require(ok,msg):
    if not bool(ok):raise ValueError("MOIS_2024_UC_SOURCE: "+msg)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(4<<20),b""):h.update(b)
    return h.hexdigest()

def parse_count(num,nodata):
    if not math.isfinite(float(num)):return None,"NODATA"
    if nodata is not None and float(num)==float(nodata):return None,"NODATA"
    if float(num)<0:return None,"INVALID_NEGATIVE"
    return float(num),"OK"

def area_km2_for_row(ds,row,geod):
    require(0<=int(row)<ds.height,"invalid source raster latitude row")
    left,bottom,right,top=bounds(Window(0,int(row),1,1),ds.transform)
    val,_=geod.polygon_area_perimeter([left,right,right,left],[bottom,bottom,top,top])
    result=abs(val)/1e6
    require(math.isfinite(result) and result>0,"zero/missing actual geodesic pixel area")
    return result

def raster_samples(ds,original):
    require(ds.count==1 and ds.crs.to_epsg()==4326,"expected original 2024 WGS84 single raster band")
    require(abs(ds.res[0]-1/120)<1e-6 and abs(ds.res[1]-1/120)<1e-6,
            "population source pixel grid changed")
    require(original.target_id.is_unique and original.target_id.str.startswith("2024|").all(),
            "not exact year-2024 original photograph sites")
    xy=original[["query_longitude","query_latitude"]].astype(float).to_numpy()
    rc=np.asarray([ds.index(float(x),float(y)) for x,y in xy],dtype=int)
    ok=(rc[:,0]>=0)&(rc[:,0]<ds.height)&(rc[:,1]>=0)&(rc[:,1]<ds.width)
    bh,bw=ds.block_shapes[0]
    tiles=(rc[:,0]//max(1,bh))*100000+(rc[:,1]//max(1,bw))
    order=np.argsort(tiles,kind="stable")
    ix=order[ok[order]]
    counts=np.full(len(original),np.nan)
    if len(ix):
        nums=ds.sample([tuple(r) for r in xy[ix]],indexes=1,masked=False)
        counts[ix]=np.fromiter((float(a[0]) for a in nums),dtype=np.float64,count=len(ix))
    area=np.full(len(original),np.nan)
    geod=Geod(ellps="WGS84")
    area_by_row={int(r):area_km2_for_row(ds,int(r),geod) for r in np.unique(rc[ok,0])}
    for i in np.flatnonzero(ok):area[i]=area_by_row[int(rc[i,0])]
    status=np.empty(len(original),dtype=object)
    for i in range(len(original)):
        if not ok[i]:status[i]="OUTSIDE_SOURCE_BOUNDS";continue
        value,state=parse_count(counts[i],ds.nodata)
        status[i]=state
        if state!="OK":counts[i]=np.nan
    density=np.divide(counts,area,out=np.full(len(original),np.nan),
                     where=(status=="OK")&np.isfinite(area)&(area>0))
    return pd.DataFrame({"target_id":original.target_id.astype(str),
             "population_UC_R2024B_status":status,
             "population_UC_R2024B_count_per_cell":counts,
             "population_UC_R2024B_pixel_area_km2":area,
             "population_UC_R2024B_density_per_km2":density}),{
             "outside_source_bounds":int((~ok).sum()),
             "source_nodata":int((status=="NODATA").sum()),
             "invalid_negative":int((status=="INVALID_NEGATIVE").sum()),
             "valid_population":int((status=="OK").sum()),
             "valid_zero_population":int((counts==0).sum()),
             "distinct_pixel_rows":int(len(area_by_row))}

def run(source,uc_file,output):
    output=Path(output)
    require(not output.exists(),"Refuse to replace a signed source output")
    source=Path(source);uc_file=Path(uc_file)
    require(sha(source)==ORIGINAL_2024_SHA,"Original 2024 66,456 photo-site values SHA drift")
    require(uc_file.name==UC_FILENAME and uc_file.stat().st_size==UC_EXPECTED_BYTES,
            "Not full original 2024 WorldPop UC R2024B source")
    d=pd.read_csv(source,dtype={"target_id":str},low_memory=False)
    require(len(d)==66456 and d.target_id.is_unique,"Original known land cohort changed")
    require(int(d.viirs_valid_months.ge(3).sum())==64646 and
            int(d.population_source_status.eq("OK").sum())==44667 and
            int(d.complete.sum())==43595,
            "Original R2025A constrained population/VIIRS source changed")
    raster_SHA=sha(uc_file)
    with rasterio.open(uc_file) as ds:
        vals,checks=raster_samples(ds,d)
        raster_receipt={"file_sha256":raster_SHA,"bytes":uc_file.stat().st_size,
            "crs":str(ds.crs),"width":ds.width,"height":ds.height,
            "bounds":list(ds.bounds),"nodata":ds.nodata,"resolution":list(ds.res)}
    merged=d[["target_id","land_class","viirs_valid_months","viirs_median",
        "population_source_status","population_density","complete","n_attempted",
        "n_strict9","source_subarms"]].merge(vals,on="target_id",how="left",validate="one_to_one")
    require(len(merged)==len(d) and merged.target_id.is_unique,"Population extraction changed source IDs")
    merged["UC_source_complete"]=merged.viirs_valid_months.ge(3)&merged.population_UC_R2024B_status.eq("OK")
    aggregate=merged.groupby("land_class").agg(
        targets=("target_id","size"),
        constrained_complete=("complete","sum"),
        unconstrained_complete=("UC_source_complete","sum"),
        unconstrained_valid_population=("population_UC_R2024B_status",lambda z:int(z.eq("OK").sum())),
        constrained_population_nodata=("population_source_status",lambda z:int(z.ne("OK").sum())),
        preexisting_strict9_success=("n_strict9","sum")).reset_index()
    output.mkdir(parents=True)
    merged.to_csv(output/"2024_66k_UC_R2024B_fullcohort_values.csv",index=False,float_format="%.12g")
    aggregate.to_csv(output/"2024_CN_R2025A_vs_UC_R2024B_population_coverage.csv",index=False)
    receipt={"schema":"MOIS_WorldPop_2024_UNCONSTRAINED_R2024B_separate_66k_v1",
        "original_2024_constrained_2024_values_SHA256":ORIGINAL_2024_SHA,
        "original_R2025A_CN_full_raster_SHA256":CN_SOURCE_SHA,
        "new_official_population_source_url":UC_URL,
        "UC_original_GeoTIFF":raster_receipt,
        "quality":checks,"source_cohort":len(merged),
        "2024_viirs3plus":int(d.viirs_valid_months.ge(3).sum()),
        "previous_CN_pop2024_valid":int(d.population_source_status.eq("OK").sum()),
        "previous_CN_complete":int(d.complete.sum()),
        "UC_pop2024_valid":int(merged.population_UC_R2024B_status.eq("OK").sum()),
        "UC_complete":int(merged.UC_source_complete.sum()),
        "UC_added_over_CN_missing":int((merged.UC_source_complete&
                             merged.population_source_status.ne("OK")).sum()),
        "UC_lost_where_CN_was_complete":int((merged.complete&
                             ~merged.UC_source_complete).sum()),
        "source_releases_NOT_mixed_in_one_covariate":True,
        "no_original_photo_downloads_or_ROI_reinference":True,
        "no_missing_population_imputed_to_zero":True,
        "no_new_ecological_2024_light_effect_fitted":True,
        "output_SHA256":{p.name:sha(p) for p in output.glob("*.csv")}}
    (output/"UC_R2024B_source_receipt.json").write_text(
        json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print("MOIS_2024_WORLDPOP_UC_REAL_SOURCE",json.dumps({
        k:receipt[k] for k in ("source_cohort","previous_CN_complete","UC_pop2024_valid",
          "UC_complete","UC_added_over_CN_missing","UC_lost_where_CN_was_complete")},sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",required=True,type=Path)
    p.add_argument("--uc-raster",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    run(a.source,a.uc_raster,a.output)
