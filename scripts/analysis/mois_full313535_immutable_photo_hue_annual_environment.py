#!/usr/bin/env python3
"""MOIS actual full4096 new photographed ROI pixels -> annual abiotic Parquet v2.

REPLACE all prior 512/1792 partial fourth-wave snapshots with one full
immutable original 4096 photo cohort. Exclude separately season-enriched 256
photos from global means. The goal is descriptive near-half WHITE vs chromatic
ecology, not a globally maintained 50:50 adaptive equilibrium or pigment cost.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis import mois_primary_global_white_balance_feature_lake as old
from scripts.analysis.mois_sample_248312_worldclim_static_and_monthnormals import sample_tiff

TOTAL=313535
NEW=4096
N_CLASSES_NEW=1577
CLASSES=151239
DATES=151230
FULL_SOURCE="mois_full313535_original_pixel_ROIs_and_observation_dates.csv.gz"
FULL_RECEIPT="mois_fourth_4096_actual_original_ROI_pixels_and_source_dates_decision.json"
ORIGINAL_RELEASE="mois-primary-309439-featurelake-546266e2bbaa-v1"
OUT_PHOTO="primary313535_fullfourth4096_original_photo_hues.parquet"
OUT_ANNUAL="primary313535_fullfourth4096_annual_worldclim.parquet"
OUT_STATIC="primary313535_original_partial_soil_terrain.parquet"
OUT_BALANCE="primary313535_white_chromatic_balance.json"
OUT_DECILES="primary313535_annual_environment_deciles.csv"
OUT_LAT="primary313535_latitude_photo_white.csv"
OUT_MANIFEST="primary313535_fullfourth4096_featurelake_manifest.json"


def require(v, text):
    if not v: raise ValueError(text)


def file_sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):h.update(chunk)
    return h.hexdigest()


def verified_sources(full_path,receipt_path,release_root):
    source=json.loads(receipt_path.read_text())
    require(
        source["schema"]=="mois_actual_fourth4096_genuine_original_RGB_ROI_309439plus4096_v1"
        and source["before_real_measured_photo_ids"]==old.N
        and source["new_photo_ids_actually_processed"]==NEW
        and source["new_photo_four_hue_original_ROI_classified"]==N_CLASSES_NEW
        and source["expanded_genuine_original_photo_ids"]==TOTAL
        and source["expanded_original_four_hue_ROI_photo_ids"]==CLASSES
        and source["expanded_actual_date_and_four_hue_ROI_photo_ids"]==DATES
        and source["initial_measured_original_photo_GZIP_SHA256"]==old.SHA_PHOTO
        and source["no_old_source_hue_taxon_photo_or_observation_date_changed"] is True
        and source["botanically_confirmed_open_petal_per_photo"] is False
        and source["current_fourth4096_source_SHA256"]==file_sha(full_path),
        "New 4096 photographed original pixels not source-complete or authentic")
    full=pd.read_csv(full_path,compression="gzip",dtype=str,
                     keep_default_na=False,low_memory=False)
    require(len(full)==TOTAL and full.photo_id.is_unique and
            int(full.morph.isin(old.HUES).sum())==CLASSES and
            int((full.morph.isin(old.HUES)&
                 full.qualified_date.str.fullmatch(
                     r"\d{4}-\d{2}-\d{2}")).sum())==DATES,
            "Full original RGB-photo ID or classified/dated colour count drift")
    base_receipt=json.loads((release_root/old.MANIFEST).read_text())
    require(base_receipt["original_photo_SHA"]==old.SHA_PHOTO
            and base_receipt["original_annual_climate_SHA"]==old.SHA_CLIMATE
            and base_receipt["original_partial_248312_terrain_SHA"]==old.SHA_STATIC
            and base_receipt["original_measured_RGB_photo_IDs"]==old.N
            and base_receipt["known_four_hue_original_ROIs"]==old.CLASS,
            "Base released original309439 environmental snapshot not immutable")
    frames=[]
    for key in (old.P_OUT,old.E_OUT,old.S_OUT):
        path=release_root/key
        require(path.is_file() and
                file_sha(path)==base_receipt["files_SHA256"][key],
                "Old primary Parquet feature-lake source SHA not verified")
        frames.append(pd.read_parquet(path))
    orig,annual,static=frames
    require(len(orig)==len(annual)==len(static)==old.N and
            orig.photo_id.is_unique and
            orig.photo_id.tolist()==annual.photo_id.tolist()==static.photo_id.tolist(),
            "Original Parquet photo-ID rows or ordering changed")
    check=full.iloc[:old.N]
    require(check.photo_id.tolist()==orig.photo_id.tolist(),
            "Fourth4096 primary 309439 source prefix not retained")
    for key in ("inat_taxon_id","species","morph","observation_id",
                "qualified_date","old_mois_or_new_fcp"):
        require(check[key].tolist()==orig[key].tolist(),
                "Original photo taxon, hue, date or recruitment source changed")
    new=full.iloc[old.N:].reset_index(drop=True)
    require(len(new)==NEW and new.photo_id.is_unique and
            not new.photo_id.isin(orig.photo_id).any() and
            new.qualified_date.str.fullmatch(
                r"\d{4}-\d{2}-\d{2}").all(),
            "New fourth4096 photos are incomplete, overlapping or undated")
    require(int(new.morph.isin(old.HUES).sum())==N_CLASSES_NEW,
            "Completed new source image ROI hue outcomes changed")
    return full,orig,annual,static,source,base_receipt


def annual_new_site(new,raster_root):
    lat=pd.to_numeric(new.latitude,errors="coerce").to_numpy(float)
    lon=pd.to_numeric(new.longitude,errors="coerce").to_numpy(float)
    q=pd.DataFrame({"photo_id":new.photo_id.to_numpy()})
    for field,num in (
        ("wc2_bio05_10m",5),("wc2_bio12_10m",12),
        ("wc2_bio14_10m",14),("wc2_bio15_10m",15)):
        tif=raster_root/"bio"/f"wc2.1_10m_bio_{num}.tif"
        require(tif.is_file(),"Missing true original WorldClim BIO raster")
        q[field]=sample_tiff(tif,lon,lat)
    totals=np.zeros(len(q),dtype=np.float64)
    counts=np.zeros(len(q),dtype=np.int16)
    for mo in range(1,13):
        tif=raster_root/"srad"/f"wc2.1_10m_srad_{mo:02d}.tif"
        require(tif.is_file(),"Missing original monthly solar raster")
        vals=sample_tiff(tif,lon,lat)
        known=np.isfinite(vals)
        totals[known]+=vals[known]/86.4
        counts[known]+=1
    solar=np.full(len(q),np.nan,dtype=np.float32)
    good=counts==12
    solar[good]=(totals[good]/12).astype(np.float32)
    q["annual_solar_Wm2"]=solar
    return q


def assemble(full,orig,annual,static,raster_root):
    new=full.iloc[old.N:].copy().reset_index(drop=True)
    newp=new[["photo_id","inat_taxon_id","species","morph",
              "observation_id","qualified_date","old_mois_or_new_fcp",
              "latitude","longitude"]].copy()
    for key in ("latitude","longitude"):
        newp[key]=pd.to_numeric(newp[key],errors="coerce")
    newp["classified"]=newp.morph.isin(old.HUES)
    newp["white"]=newp.morph.eq("white")
    require(list(newp.columns)==list(orig.columns),
            "Previous primary original photo Parquet schema changed")
    p=pd.concat((orig,newp),ignore_index=True)
    extra=annual_new_site(new,raster_root)
    require(list(extra.columns)==list(annual.columns),
            "New WorldClim original site annual six-feature schema changed")
    env=pd.concat((annual,extra),ignore_index=True)
    newly_missing=pd.DataFrame({"photo_id":new.photo_id.to_numpy()})
    for key in static.columns:
        if key!="photo_id":newly_missing[key]=np.nan
    terrain=pd.concat((static,newly_missing[static.columns]),ignore_index=True)
    require(len(p)==len(env)==len(terrain)==TOTAL and
            p.photo_id.is_unique and
            p.photo_id.tolist()==env.photo_id.tolist()==terrain.photo_id.tolist(),
            "Cannot append complete authentic photographed pixel/environment rows")
    require(p.iloc[:old.N].reset_index(drop=True).equals(orig)
            and env.iloc[:old.N].reset_index(drop=True).equals(annual)
            and terrain.iloc[:old.N].reset_index(drop=True).equals(static),
            "Frozen previous photo/environment data were overwritten")
    return p,env,terrain


def run(fullfile,receiptfile,release_root,rasterroot,outdir):
    require(not outdir.exists() or not any(outdir.iterdir()),
            "Never overwrite a signed photo-feature lake release")
    full,p0,e0,s0,sourcerec,oldreceipt=verified_sources(
        fullfile,receiptfile,release_root)
    p,e,s=assemble(full,p0,e0,s0,rasterroot)
    stats=old.balances(p.rename(columns={"white":"white",
                                          "classified":"classified"}))
    stats.update({
      "schema":"mois_source313535_real_original_photo_white_chromatic_nearhalf_v1",
      "original_309439_classifiable":old.CLASS,
      "new_4096_classifiable":N_CLASSES_NEW,
      "original_photo_id_attempts":TOTAL,
      "four_hue_classifiable_photo_ids":CLASSES,
      "source_date_and_four_hue":DATES,
      "four_hue_counts":{h:int(p.morph.eq(h).sum()) for h in old.HUES},
      "taxon_genus_photo_depth_sensitivity":old.balance_sampling_sensitivity(p),
      "geographically_equal_10degree_source_cells":old.geographic_cell_sensitivity(p),
      "not_representative_true_white_flower_global_fraction":True,
      "not_a_maintained_half_half_evolutionary_equilibrium":True,
      "fcp_and_season_targeted256_NOT_pooled_as_unbiased_photo_source":True,
      "pigment_synthesis_cost_and_reproductive_fitness_NOT_measured":True,
      "environmental_extraction_new_4096_WorldClim_annual_1970_2000":True,
      "new_4096_soil_terrain_not_available_and_NOT_imputed":True,
    })
    outdir.mkdir(parents=True,exist_ok=True)
    for filename,table in ((OUT_PHOTO,p),(OUT_ANNUAL,e),(OUT_STATIC,s)):
        table.to_parquet(outdir/filename,compression="zstd",
                         engine="pyarrow",index=False)
        restored=pd.read_parquet(outdir/filename,columns=["photo_id"])
        require(restored.photo_id.tolist()==p.photo_id.tolist(),
                "Saved Parquet changed source photo original ID order")
    (outdir/OUT_BALANCE).write_text(json.dumps(
        stats,indent=2,sort_keys=True,allow_nan=False)+"\n")
    joined=p.merge(e,on="photo_id",validate="one_to_one")
    old.deciles(joined).to_csv(outdir/OUT_DECILES,
                              index=False,lineterminator="\n")
    old.latitude_table(p).to_csv(outdir/OUT_LAT,
                                 index=False,lineterminator="\n")
    product=(OUT_PHOTO,OUT_ANNUAL,OUT_STATIC,OUT_BALANCE,
             OUT_DECILES,OUT_LAT)
    report={
      "schema":"mois_313535_genuine_original_ROI_feature_lake_v2",
      "full_16_workers_original_pixel_source_SHA256":file_sha(fullfile),
      "source_first309439_photo_GZIP_SHA256":old.SHA_PHOTO,
      "previous_persisted_309439_feature_snapshot_verified":True,
      "new_original_4096_photo_IDs_not_source_only":True,
      "old512_1792_photos_not_double_counted":True,
      "new_photo_pixels":NEW,
      "new_classifiable_visible_hues":N_CLASSES_NEW,
      "actual_original_photo_attempts":TOTAL,
      "four_hue_labels":CLASSES,
      "dated_four_hue_labels":DATES,
      "new_soil_terrain_not_imputed":True,
      "new_4096_annual_climate_1970_2000_only":True,
      "actual_year_month_weather_measured":False,
      "season_targeted_256_source_excluded_from_primary_hue_balance":True,
      "50_50_evolutionary_equilibrium_tested":False,
      "floral_cost_fitness_measured":False,
      "files_sha256":{k:file_sha(outdir/k) for k in product},
    }
    (outdir/OUT_MANIFEST).write_text(json.dumps(
        report,indent=2,sort_keys=True)+"\n")
    print("MOIS_TRUE_FULL4096_SOURCE313535_BALANCE_AND_ENVIRONMENT",
          json.dumps({"photos":TOTAL,"classified":CLASSES,
                      "white":stats["photo_rates"]["white"],
                      "photo_white_fraction":
                        stats["near_half_rather_than_exact_50_50_fraction"]["photos"],
                      "taxon_equal_white":
                        stats["near_half_rather_than_exact_50_50_fraction"]["species_equal"],
                      "genus_equal_white":
                        stats["near_half_rather_than_exact_50_50_fraction"]["genus_equal"]},
                     sort_keys=True),flush=True)
    return report


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--original313535-pixels",type=Path,required=True)
    ap.add_argument("--original313535-source-receipt",type=Path,required=True)
    ap.add_argument("--previous309439-feature-release",type=Path,required=True)
    ap.add_argument("--worldclim-source",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    run(args.original313535_pixels,args.original313535_source_receipt,
        args.previous309439_feature_release,args.worldclim_source,args.output)
