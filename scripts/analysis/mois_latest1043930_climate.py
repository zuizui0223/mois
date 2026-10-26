#!/usr/bin/env python3
"""Actually append original WorldClim moisture/heat/light/elevation to ALL
latest 1,043,930 SHA-signed previously image-processed MOIS photo IDs.

Reuse existing exact 350,399 individual-photo climatology rows wherever
photo IDs intersect. For EVERY other genuine original image photo ID, sample
the same frozen WorldClim BIO5/BIO12/BIO14, monthly VAPR/SRAD and 2.5m DEM.
No photo inference, missing-colour imputation, "flowering" assumption or
metadata-only candidate ID masquerading as a measured ROI.
"""
from __future__ import annotations
import argparse,json,hashlib
from pathlib import Path
import numpy as np
import pandas as pd

from scripts.analysis.mois_sample_248312_worldclim_static_and_monthnormals import sample_tiff

HUES=("white","yellow_orange","red_pink","blue_purple")
FROZEN_SHA="9a4f3dae6fbcee967b71cb29576e71d52c48e1f1f790b7d3ea4750c4bcb0e6d4"
COLUMNS=("wc2_bio05_10m","wc2_bio12_10m","wc2_bio14_10m",
         "wc2_annual_vapr_10m","annual_solar_Wm2","wc2_elev_2p5m")
SCHEMA="mois_real_latest1043930_photo_id_aligned_worldclim_humidity_rain_altitude_v1"
def require(k,why):
    if not k:raise ValueError("MOIS_LATEST1043930_CLIMATE: "+why)
def sha(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as f:
        for v in iter(lambda:f.read(4<<20),b""):h.update(v)
    return h.hexdigest()
def sample_new(d,root):
    n=len(d)
    lon=pd.to_numeric(d.longitude,errors="coerce").to_numpy(float)
    lat=pd.to_numeric(d.latitude,errors="coerce").to_numpy(float)
    result=pd.DataFrame({"photo_id":d.photo_id.to_numpy()})
    for bio in (5,12,14):
        target=root/"bio"/("wc2.1_10m_bio_%d.tif"%bio)
        require(target.is_file(), "Original raster missing: "+str(target))
        result["wc2_bio%02d_10m"%bio]=sample_tiff(target,lon,lat)
        print("MOIS_LATEST1043930_CLIMATE_RASTER_SAMPLED BIO",bio,"new",n,flush=True)
    wet=np.zeros(n,dtype=float);sun=np.zeros(n,dtype=float)
    wcount=np.zeros(n,dtype=np.int8);scount=np.zeros(n,dtype=np.int8)
    for month in range(1,13):
        vap=root/"vapr"/("wc2.1_10m_vapr_%02d.tif"%month)
        sol=root/"srad"/("wc2.1_10m_srad_%02d.tif"%month)
        require(vap.is_file() and sol.is_file(),"Original monthly humidity/radiation tiff unavailable")
        x=sample_tiff(vap,lon,lat);y=sample_tiff(sol,lon,lat)
        v=np.isfinite(x);z=np.isfinite(y)
        wet[v]+=x[v];wcount[v]+=1
        sun[z]+=y[z]/86.4;scount[z]+=1
        if month%3==0:
            print("MOIS_LATEST1043930_CLIMATE_MONTH",month,flush=True)
    vap=np.full(n,np.nan,dtype=np.float32)
    sr=np.full(n,np.nan,dtype=np.float32)
    valid=wcount==12;vap[valid]=(wet[valid]/12).astype("float32")
    valid=scount==12;sr[valid]=(sun[valid]/12).astype("float32")
    result["wc2_annual_vapr_10m"]=vap
    result["annual_solar_Wm2"]=sr
    elev=root/"elev"/"wc2.1_2.5m_elev.tif"
    require(elev.is_file(),"DEM source tiff missing")
    result["wc2_elev_2p5m"]=sample_tiff(elev,lon,lat)
    require(len(result)==n and result.photo_id.is_unique,"Original new photo sampled once failed")
    return result[["photo_id",*COLUMNS]]

def run(photo,photo_receipt,old_env,old_receipt,worldclim,out):
    photo,photo_receipt,old_env,old_receipt,worldclim,out=map(Path,
        (photo,photo_receipt,old_env,old_receipt,worldclim,out))
    require(not out.exists(),"Do not overwrite previous source photo-climate result")
    meta=json.loads(photo_receipt.read_text(encoding="utf-8"))
    assert meta["source_pixel_attempts"]==1043930 and meta["source_four_hue_classified"]==459235
    require(meta["photo_parquet_sha256"]==FROZEN_SHA==sha(photo),"Not genuinely measured latest full photo source SHA")
    earlier=json.loads(old_receipt.read_text(encoding="utf-8"))
    source_env_name=old_env.name
    require(earlier["source_photo_attempts"]==350399
            and earlier["sha256"][source_env_name]==sha(old_env),
            "Previously published exact350399 environmental source changed")
    d=pd.read_parquet(photo,columns=["photo_id","latitude","longitude","morph"])
    require(len(d)==1043930 and d.photo_id.astype(str).is_unique,"Original source photo IDs not unique")
    d["photo_id"]=d.photo_id.astype(str)
    old=pd.read_parquet(old_env,columns=["photo_id",*COLUMNS])
    old["photo_id"]=old.photo_id.astype(str)
    require(len(old)==350399 and old.photo_id.is_unique,
            "Signed 350399 environment photo observations missing")
    merged=d[["photo_id","latitude","longitude"]].merge(
        old,on="photo_id",how="left",sort=False,validate="one_to_one",indicator=True)
    require(len(merged)==len(d) and merged.photo_id.tolist()==d.photo_id.tolist(),
            "Original photo order/identity lost by climate join")
    reused=merged["_merge"].eq("both")
    number_old=int(reused.sum())
    require(number_old>=300000,
            "Historical 350399 climate/photo overlaps unexpectedly sparse: "+
            str(number_old))
    sampled=sample_new(merged.loc[~reused,["photo_id","latitude","longitude"]],worldclim)
    require(len(sampled)+number_old==1043930 and sampled.photo_id.is_unique,
            "Original new source photo sites did not reconcile")
    outenv=merged.drop(columns=["_merge","latitude","longitude"]).copy()
    for col in COLUMNS:
        outenv.loc[~reused,col]=sampled[col].to_numpy()
        outenv[col]=outenv[col].astype("float32")
    outenv=outenv[["photo_id",*COLUMNS]]
    require(outenv.photo_id.is_unique and len(outenv)==1043930
            and outenv.photo_id.tolist()==d.photo_id.tolist(),
            "Annual climate must preserve all source photo IDs in exact original order")
    complete=np.isfinite(outenv[list(COLUMNS)]).all(axis=1).to_numpy()
    classified=d.morph.isin(HUES).to_numpy()
    # Re-used original older environmental values MUST be retained by original ID.
    rejoined=outenv.loc[reused,list(COLUMNS)].to_numpy(dtype=np.float32)
    before=old.set_index("photo_id").loc[
        outenv.loc[reused,"photo_id"].tolist(),list(COLUMNS)].to_numpy(dtype=np.float32)
    require(np.array_equal(np.isnan(rejoined),np.isnan(before))
            and np.allclose(rejoined,before,equal_nan=True,rtol=0,atol=0),
            "Previously SHA-frozen WorldClim 350399 values modified")
    out.mkdir(parents=True)
    target=out/"mois_latest1043930_original_photo_climate_six_columns.parquet"
    outenv.to_parquet(target,index=False,compression="zstd")
    receipt={
        "schema":SCHEMA,
        "original_exact_photo_attempts":1043930,
        "four_hue_classified":459235,
        "reused_old350399_same_photo_ID_annual_climate_rows":number_old,
        "new_original_photo_ID_climate_sites_sampled":len(sampled),
        "original_photo_climate_complete":int(complete.sum()),
        "four_hue_classified_climate_complete":int((complete&classified).sum()),
        "old350399_environment_sha256":sha(old_env),
        "new1043930_true_photo_source_sha256":FROZEN_SHA,
        "new1043930_climate_parquet_sha256":sha(target),
        "worldclim_zip_sha256":{
            group:sha(worldclim/group/("wc2.1_2.5m_elev.zip" if group=="elev"
                  else "wc2.1_10m_"+group+".zip"))
            for group in ("bio","srad","vapr","elev")},
        "climate_source":"WorldClim2.1 1970-2000 normals, 10arcmin BIO/SRAD/VAPR and 2.5arcmin DEM",
        "solar_unit":"W/m2 = WorldClim SRAD kJ/m2/day / 86.4, annual 12-month mean",
        "all_photo_ids_preserved":True,"unknowns_never_imputed_as_petal_colour":True,
        "additional_source_image_downloads":0,
        "biological_organ_anthesis_or_fitness_verified":False,
    }
    (out/"mois_latest1043930_photo_climate_receipt.json").write_text(
        json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print("MOIS_LATEST1043930_ACTUAL_PHOTO_CLIMATE_JOIN",json.dumps(receipt,sort_keys=True),flush=True)
    return receipt

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("photo","photo-receipt","old-env","old-receipt","worldclim","output"):
        p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args()
    run(a.photo,a.photo_receipt,a.old_env,a.old_receipt,a.worldclim,a.output)
if __name__=="__main__":
    main()
