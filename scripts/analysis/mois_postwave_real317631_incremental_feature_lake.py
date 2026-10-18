#!/usr/bin/env python3
"""On receipt of ALL 16 genuine fifth-wave ROI workers, append environment once.

313535->317631 actual photo RGB attempts, not new source-only IDs. Original
visible hue and observation date are immutable. Source-verified official
WorldClim 1970-2000 BIO19, 12 VAPR, 12 solar and DEM elevation/terrain are
sampled ONLY at the new 4096 image-photo coordinates. Original soil missing
at new sites is never set to zero. No genetic/fitness/50:50 equilibrium claim.
"""
from __future__ import annotations

import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd

from scripts.analysis import mois_global313535_water_light_heat_habitat_partition as base
from scripts.analysis import mois_original313535_full19_bio_vapr_terrain_lake as expanded
from scripts.analysis import mois_primary_global_white_balance_feature_lake as rates
from scripts.analysis.mois_sample_248312_worldclim_static_and_monthnormals import sample_tiff,raster_terrain

NEW=4096
TOTAL=317631
INPUT_SOURCE="mois_full317631_actual_original_photo_RGB_hues_and_dates.csv.gz"
INPUT_RECEIPT="mois_fifth4096_real_original_photo_source_reconciled_receipt.json"
PRIOR_SOURCE="mois_full313535_original_pixel_ROIs_and_observation_dates.csv.gz"
PRIOR_RECEIPT="mois_fourth_4096_actual_original_ROI_pixels_and_source_dates_decision.json"
PHOTO_OUT="primary317631_actual_photo_hue_RGB_attempts.parquet"
ANNUAL_OUT="primary317631_original_annual_bio_solar.parquet"
TERRAIN_OUT="primary317631_partial_soil_and_measured_terrain.parquet"
ALL19_OUT="primary317631_all19_BIO_vapr_solar_terrain.parquet"
SUMMARY="mois_317631_actual_new4096_original_photo_white_balance.json"
RECEIPT="mois_317631_new4096_incremental_climate_feature_lake_receipt.json"


def need(ok,message):
    if not ok:raise ValueError(message)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def load_verified(real,pixrec,oldpixel,oldrec,prior,extended_dir):
    receipts=json.loads(pixrec.read_text())
    previous=json.loads(oldrec.read_text())
    need(receipts["schema"]=="mois_fifth4096_original_real_ROI_photo_source_after313535_v1"
         and receipts["previous_real_original_photo_attempts"]==base.N_PHOTOS
         and receipts["new_original_RGB_image_attempts"]==NEW
         and receipts["total_source_photo_pixels_actually_processed"]==TOTAL
         and receipts["previous_313535_original_photo_rows_unchanged"] is True
         and receipts["previous_complete_original_source_SHA256"]==sha(oldpixel)
         and receipts["expanded_complete_original_source_SHA256"]==sha(real)
         and previous["current_fourth4096_source_SHA256"]==sha(oldpixel)
         and previous["expanded_genuine_original_photo_ids"]==base.N_PHOTOS
         and receipts["separately_prioritised_season4096_source_photo_IDs_excluded"] is True,
         "New genuine fifth original photo cohort absent or unverified")
    main,pr=base.load_source(prior)
    extrec=json.loads((extended_dir/expanded.OUT_REC).read_text())
    abio=extended_dir/expanded.OUT_ENV
    need(extrec["schema"]=="mois_actual313535_original_photo_worldclim_full19_vapr_dem_v1"
         and extrec["real_original_processed_photo_ids"]==base.N_PHOTOS
         and extrec["new_environment_parquet_sha256"]==sha(abio)
         and extrec["source_prior_release_three_Parquet_SHA256"]=={
            k:pr["files_sha256"][k] for k in
            (base.PHOTO,base.ANNUAL,base.STATIC)},
         "Prior original PHOTO ID/abiotic table bytes do not share source provenance")
    frozen=pd.read_csv(oldpixel,compression="gzip",dtype=str,
                       keep_default_na=False,low_memory=False)
    new=pd.read_csv(real,compression="gzip",dtype=str,
                    keep_default_na=False,low_memory=False)
    need(len(frozen)==base.N_PHOTOS and len(new)==TOTAL and
         new.photo_id.is_unique and frozen.equals(
            new.iloc[:base.N_PHOTOS].reset_index(drop=True)),
         "Earlier original RGB/colour/date source rows were overwritten")
    p=pd.read_parquet(prior/base.PHOTO)
    ann=pd.read_parquet(prior/base.ANNUAL)
    terrain=pd.read_parquet(prior/base.STATIC)
    abiotic=pd.read_parquet(abio)
    need(all(len(x)==base.N_PHOTOS for x in (p,ann,terrain,abiotic))
         and p.photo_id.tolist()==frozen.photo_id.tolist()
         and p.photo_id.tolist()==abiotic.photo_id.tolist()
         and p.photo_id.tolist()==ann.photo_id.tolist()
         and p.photo_id.tolist()==terrain.photo_id.tolist(),
         "Prior photo and all environment blocks are not original-photo aligned")
    for key in ("inat_taxon_id","species","morph","observation_id",
                "qualified_date","old_mois_or_new_fcp"):
        need(p[key].astype(str).tolist()==frozen[key].tolist(),
             "Previously source-certified photo colour/taxon/day has changed")
    extra=new.iloc[base.N_PHOTOS:].reset_index(drop=True)
    need(len(extra)==NEW and extra.photo_id.is_unique
         and not extra.photo_id.isin(p.photo_id).any()
         and extra.qualified_date.str.fullmatch(r"\d{4}-\d{2}-\d{2}").all()
         and int(extra.morph.isin(base.HUES).sum())==receipts[
             "new_original_classifiable_ROIs"],
         "Extra true RGB original photo colour outcomes have lost identity")
    return extra,p,ann,terrain,abiotic,receipts,extrec


def new_raster_variables(extra,raster):
    lon=pd.to_numeric(extra.longitude,errors="coerce").to_numpy(float)
    lat=pd.to_numeric(extra.latitude,errors="coerce").to_numpy(float)
    r=pd.DataFrame({"photo_id":extra.photo_id.to_numpy()})
    for i,key in enumerate(expanded.BIO,1):
        p=raster/"bio"/f"wc2.1_10m_bio_{i}.tif"
        need(p.is_file(),"WorldClim 19 BIO original raster incomplete")
        r[key]=sample_tiff(p,lon,lat)
    v=np.zeros(NEW);vc=np.zeros(NEW,dtype=np.int8)
    solar=np.zeros(NEW);sc=np.zeros(NEW,dtype=np.int8)
    for month in range(1,13):
        pv=raster/"vapr"/f"wc2.1_10m_vapr_{month:02d}.tif"
        ps=raster/"srad"/f"wc2.1_10m_srad_{month:02d}.tif"
        need(pv.is_file() and ps.is_file(),
             "Original WorldClim VAPR/SRAD month raster incomplete")
        a=sample_tiff(pv,lon,lat);b=sample_tiff(ps,lon,lat)
        va=np.isfinite(a);so=np.isfinite(b)
        v[va]+=a[va];vc[va]+=1
        solar[so]+=b[so]/86.4;sc[so]+=1
    vapr=np.full(NEW,np.nan,dtype=np.float32)
    rad=np.full(NEW,np.nan,dtype=np.float32)
    mask=vc==12;vapr[mask]=(v[mask]/12).astype(np.float32)
    mask=sc==12;rad[mask]=(solar[mask]/12).astype(np.float32)
    r[expanded.VAPR]=vapr
    r[base.SOLAR]=rad
    topo=raster_terrain(raster/"elev"/"wc2.1_2.5m_elev.tif",lon,lat)
    for name in expanded.TERRAIN:r[name]=topo[name]
    need(len(r)==NEW and r.photo_id.is_unique,
         "New original climate results lost one-to-one photo IDs")
    return r


def append(original,ext,out):
    need(len(ext)==NEW and original.photo_id.is_unique and ext.photo_id.is_unique
         and not ext.photo_id.isin(original.photo_id).any()
         and list(original.columns)==list(ext.columns),
         "Cannot append measured photo or environmental rows with differing columns")
    combined=pd.concat((original,ext),ignore_index=True)
    need(len(combined)==TOTAL and combined.photo_id.is_unique
         and combined.iloc[:base.N_PHOTOS].reset_index(drop=True).equals(original),
         "Earlier photo/abiotic source was not retained exactly")
    return combined


def run(real,pixelrec,oldpix,oldrec,prior,extended_dir,raster,out):
    need(not out.exists() or not any(out.iterdir()),
         "Cannot replace verified source-immutable publication")
    extra,p,annual,static,abiotic,source,abiorec=load_verified(
        real,pixelrec,oldpix,oldrec,prior,extended_dir)
    zip_sources={
        "bio":raster/"bio"/"wc2.1_10m_bio.zip",
        "srad":raster/"srad"/"wc2.1_10m_srad.zip",
        "vapr":raster/"vapr"/"wc2.1_10m_vapr.zip",
        "elev":raster/"elev"/"wc2.1_2.5m_elev.zip",
    }
    for key,path in zip_sources.items():
        need(path.is_file() and
             sha(path)==abiorec["source_climatology_zip_SHA256"].get(
                  key,pr_hash_compat(key,abiorec)),
             "Official original WorldClim raster archive unexpectedly changed")
    need(sha(zip_sources["srad"])==expanded.PIN_SRAD,
         "Original WorldClim SRAD source unit or raster archive changed")
    newenv=new_raster_variables(extra,raster)
    ep=extra[["photo_id","inat_taxon_id","species","morph",
              "observation_id","qualified_date","old_mois_or_new_fcp",
              "latitude","longitude"]].copy()
    for col in ("latitude","longitude"):
        ep[col]=pd.to_numeric(ep[col],errors="coerce")
    ep["classified"]=ep.morph.isin(base.HUES)
    ep["white"]=ep.morph.eq("white")
    P=append(p,ep[p.columns],PHOTO_OUT)
    A=append(annual,newenv[annual.columns],ANNUAL_OUT)
    fresh_static=pd.DataFrame({"photo_id":extra.photo_id.to_numpy()})
    for k in static.columns:
        if k=="photo_id":continue
        fresh_static[k]=(newenv[k] if k in newenv else np.nan)
    S=append(static,fresh_static[static.columns],TERRAIN_OUT)
    X=append(abiotic,newenv[abiotic.columns],ALL19_OUT)
    out.mkdir(parents=True,exist_ok=True)
    frames={PHOTO_OUT:P,ANNUAL_OUT:A,TERRAIN_OUT:S,ALL19_OUT:X}
    for k,df in frames.items():
        df.to_parquet(out/k,compression="zstd",index=False,engine="pyarrow")
        need(pd.read_parquet(out/k,columns=["photo_id"]).photo_id.tolist()==
             P.photo_id.tolist(),"Saved original photo IDs no longer aligned")
    rates_by_group=P.groupby("inat_taxon_id").white.mean()
    g=P.assign(genus=P.species.map(base.genus))
    genus_rates=(g[g.classified].groupby(["genus","inat_taxon_id"]).white.mean()
                  .groupby(level=0).mean())
    actual=P[P.classified]
    white=int(actual.white.sum())
    report={
        "schema":"mois_real317631_fifth4096_incremental_RGB_annual_bio19_terrain_v1",
        "previous_photo_RGB_original_attempts":base.N_PHOTOS,
        "new_processed_original_RGB_photo_attempts":NEW,
        "all_original_photo_attempts":TOTAL,
        "new_4hue_original_visible_ROIs":source["new_original_classifiable_ROIs"],
        "all_classified_four_hue_photo_ROIs":len(actual),
        "all_white_photo_ROIs":white,
        "white_fraction_classifiable_photo_ROIs":float(white/len(actual)),
        "all_chromatic_photo_ROIs":len(actual)-white,
        "total_original_observation_dates_with_known_colour":
            int((P.classified&P.qualified_date.astype(str).str.fullmatch(
                r"\d{4}-\d{2}-\d{2}")).sum()),
        "taxon_equal_photo_white_fraction":float(
            P[P.classified].groupby("inat_taxon_id").white.mean().mean()),
        "genus_equal_photo_white_fraction":float(genus_rates.mean()),
        "only_new_4096_photo_sites_extracted_this_execution":True,
        "original_photo_pixel_results_unchanged":True,
        "previous_313535_yearly_environment_data_unchanged":True,
        "annual_WorldClim_v2_1_normal_period":"1970-2000",
        "new_photo_soil_not_imputed":True,
        "source_supplementary_season_targeted_photo_ids_excluded":True,
        "new_source_fifth4096_photo_RGBA_original_SHA":sha(real),
        "source_original_raster_sha256":
            {k:sha(x) for k,x in zip_sources.items()},
        "photo_ID_aligned_Zstd_files_sha256":
            {k:sha(out/k) for k in frames},
        "not_an_unbiased_worldwide_flowering_individual_census":True,
        "not_a_maintained_50_50_colour_selection_or_pigment_fitness_test":True,
    }
    (out/RECEIPT).write_text(
        json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("MOIS_REAL317631_INCREMENTAL_PIXEL_CLIMATE_FEATURE_LAKE",
          json.dumps({k:report[k] for k in (
              "all_original_photo_attempts","all_classified_four_hue_photo_ROIs",
              "all_white_photo_ROIs","white_fraction_classifiable_photo_ROIs",
              "only_new_4096_photo_sites_extracted_this_execution")},
              sort_keys=True),flush=True)
    return report


def pr_hash_compat(k,d):
    # Previous BIO and SRAD sources use pinned SHA, while previous all19
    # origin manifest holds the official vapour and DEM source archives.
    if k=="srad":return expanded.PIN_SRAD
    if k=="bio":return expanded.PIN_BIO
    return ""


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--new317631-original-photo-ROI",type=Path,required=True)
    p.add_argument("--new317631-source-receipt",type=Path,required=True)
    p.add_argument("--original313535-RGB-source",type=Path,required=True)
    p.add_argument("--original313535-source-receipt",type=Path,required=True)
    p.add_argument("--original313535-feature-release",type=Path,required=True)
    p.add_argument("--original313535-full19-release",type=Path,required=True)
    p.add_argument("--real-WorldClim-raster-root",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    run(a.new317631_original_photo_ROI,a.new317631_source_receipt,
        a.original313535_RGB_source,a.original313535_source_receipt,
        a.original313535_feature_release,a.original313535_full19_release,
        a.real_WorldClim_raster_root,a.output)
