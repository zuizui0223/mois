#!/usr/bin/env python3
"""MOIS original 100543 photo water method on ACTUAL new original ROI pixels.

Exactly the already frozen six macro climate anchors and genus-equal 999
shared resample/maxT used in the original 100543-photo study, applied to a
strict photo-ID enlargement from previously scored original pixels. Only
sample new WorldClim site values at the newly measured photo coordinates;
reuse immutable, independently verified 248312 original photo-site climate.
No intraspecific FCP D/W/polymorphism estimands enter this analysis.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis import mois_refit_original_100k_water_six_anchors_248k as old
from scripts.analysis.mois_sample_248312_worldclim_static_and_monthnormals import sample_tiff

HUES=("white","yellow_orange","red_pink","blue_purple")
BASE_PHOTOS=248312
BASE_HUES=120807
WAVE=4096
BOOT=999
INPUT=tuple(old.READ_COLS)


def sha(path):
    x=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):x.update(b)
    return x.hexdigest()


def need(ok,message):
    if not ok:raise ValueError(message)


def site_bio(new:pd.DataFrame,worldclim:Path)->pd.DataFrame:
    coords=new[["longitude","latitude"]].apply(pd.to_numeric,errors="coerce")
    lon=coords.longitude.to_numpy(float)
    lat=coords.latitude.to_numpy(float)
    for b in (5,12,14,15):
        path=worldclim/"bio"/f"wc2.1_10m_bio_{b}.tif"
        need(path.is_file(),f"Missing true WorldClim BIO{b} source")
        new[f"wc2_bio{b:02d}_10m"]=sample_tiff(path,lon,lat)
    sums=np.zeros(len(new),dtype=np.float64)
    counts=np.zeros(len(new),dtype=np.int16)
    for month in range(1,13):
        path=worldclim/"srad"/f"wc2.1_10m_srad_{month:02d}.tif"
        need(path.is_file(),f"Missing true WorldClim solar {month} normal")
        values=sample_tiff(path,lon,lat)
        valid=np.isfinite(values)
        sums[valid]+=values[valid]
        counts[valid]+=1
    solar=np.full(len(new),np.nan,dtype=np.float32)
    good=counts==12
    solar[good]=(sums[good]/12).astype(np.float32)
    new["wc2_annual_srad_normal_10m_Wm2"]=solar
    # WorldClim source radiation is reported in kJ m^-2 day^-1, despite
    # the misleading historical Wm2 column label. A rescaling does not
    # change the original standard-deviation regression estimands.
    return new


def fit(photo_union:Path,union_receipt:Path,old_site:Path,old_site_receipt:Path,
        worldclim:Path,vapr:Path,vapr_zip:Path,outdir:Path,
        *,expected_waves:int=1):
    need(expected_waves in (1,2),"Only 1 or 2 SHA-verified original photo waves supported")
    need(not outdir.exists() or not any(outdir.iterdir()),
         "Refuse overwriting genuine original 100k-method water source result")
    rec=json.loads(union_receipt.read_text())
    env=json.loads(old_site_receipt.read_text())
    total=BASE_PHOTOS+WAVE*expected_waves
    need((rec["all_original_unique_photos"],rec["new_original_source_photos"])==
         (total,WAVE*expected_waves) and rec["all_four_hue_classified_photos"]>BASE_HUES
         and rec["original_model_SHA_verified_for_every_photo"]
         and sha(photo_union)==rec["all_real_source_photo_hue_csv_sha256"],
         "New genuine original ROI-v4 model photo sources not all authenticated")
    need(env["n_original_photo_ids"]==BASE_PHOTOS and
         env["n_four_hue_classified"]==BASE_HUES and
         sha(old_site)==env["full_photo_static_monthly_normal_csv_sha256"],
         "Historical 248312-photo real environmental pixels changed")

    # The real ROI-shard receipts gate the appended source, while the original
    # WorldClim site products gate all old photo IDs and environmental normals.
    df=pd.read_csv(photo_union,compression="gzip",low_memory=False,
                   usecols=["photo_id","inat_taxon_id","species","morph",
                            "latitude","longitude","source_original_measured"],
                   dtype={"photo_id":str,"inat_taxon_id":str})
    old_df=pd.read_csv(old_site,compression="gzip",low_memory=False,
                       usecols=list(INPUT),dtype={"photo_id":str,"inat_taxon_id":str})
    need(len(df)==total and df.photo_id.nunique()==total,"Original enlarged photo IDs not unique")
    need(len(old_df)==BASE_PHOTOS and old_df.photo_id.nunique()==BASE_PHOTOS,
         "Earlier original 248k site photo IDs not unique")
    oldids=set(old_df.photo_id)
    still=df.loc[df.photo_id.isin(oldids),["photo_id","inat_taxon_id","species","morph"]]
    need(len(still)==BASE_PHOTOS,"Previously measured photo IDs were lost")
    check=old_df[["photo_id","inat_taxon_id","species","morph"]].merge(
        still,on="photo_id",validate="one_to_one",suffixes=("_original","_new"))
    for key in ("inat_taxon_id","species","morph"):
        need(check[key+"_original"].fillna("").astype(str).equals(
             check[key+"_new"].fillna("").astype(str)),
             "An original MOIS 248312 species/photo/hue was silently re-labelled")
    added=df.loc[~df.photo_id.isin(oldids)].copy()
    need(len(added)==WAVE*expected_waves and
         added.photo_id.nunique()==len(added),
         "New source lacks the exact verified original photo ID terminal count")
    hue_new=int(added.morph.isin(HUES).sum())
    need(hue_new==rec["new_four_hue_measured_photo_ids"],
         "Original new photo ROI classification denominators differ")
    added=added.rename(columns={"source_original_measured":"original_source"})
    added=site_bio(added,worldclim)
    for field in INPUT:
        if field not in added:
            added[field]=np.nan
    both=pd.concat([old_df[list(INPUT)],added[list(INPUT)]],ignore_index=True)
    need(len(both)==total and both.photo_id.nunique()==total and
         int(both.morph.isin(HUES).sum())==rec["all_four_hue_classified_photos"],
         "Cannot reconstruct all original photographic genus/morph outcomes")
    outdir.mkdir(parents=True,exist_ok=True)
    expanded=outdir/f"mois_real_{total}_source_photo_with_six_original_climate_anchors.csv.gz"
    both.to_csv(expanded,index=False,compression={"method":"gzip","mtime":0},
                float_format="%.8f",lineterminator="\n")
    source_digest=sha(expanded)
    source_rec=outdir/"mois_new_photo_plus_prior_site_source_receipt.json"
    source_rec.write_text(json.dumps({
        "n_original_photo_ids":total,
        "n_four_hue_classified":rec["all_four_hue_classified_photos"],
        "full_photo_static_monthly_normal_csv_sha256":source_digest,
        "new_model_verified_photo_hue_count":hue_new,
        "new_source_genuinely_measured_pixels":True,
        "not_new_photo_metadata_only":True,
        "original_bio_site_sha256":sha(worldclim/"bio"/"wc2.1_10m_bio.zip"),
        "original_srad_site_sha256":sha(worldclim/"srad"/"wc2.1_10m_srad.zip"),
    },indent=2,sort_keys=True)+"\n")

    # Only numeric cohort size changes; ORIGINAL species/genus/geo estimator
    # methods and six-anchor 999 bootstrap/maxT stay byte-for-byte unchanged.
    old.EXPECTED_PHOTOS=total
    old.EXPECTED_CLASSIFIED=rec["all_four_hue_classified_photos"]
    _,valid,vapr_meta,_=old.prepare_data(expanded,source_rec,vapr)
    norm=(*old.GEO,*old.ANCHORS)
    mu=valid[list(norm)].mean()
    sd=valid[list(norm)].std(ddof=0).where(lambda s:s>1e-12,1.0)
    clean_old=valid.loc[valid.original_source.eq("historical")]
    results={
       "historical_source_same_new_photo_environment_standardization":
           old.estimate(clean_old,mu,sd,kind="historical",nboot=BOOT),
       "newly_expanded_MOIS_global_photos":
           old.estimate(valid,mu,sd,kind="original_248312_plus_new_pix",nboot=BOOT),
    }
    summary={
       "schema":"mois_original_100543_six_water_climate_anchor_999_shared_genus_bootstrap_appended_real_pixels_v1",
       "status":"REAL_NEW_ORIGINAL_PHOTO_PIXEL_COLOURS_ORIGINAL_SIX_ANCHOR_WATER_MODEL_FITTED",
       "original_photo_ids":BASE_PHOTOS,"new_photo_ids":len(added),
       "new_classified_photo_hues":hue_new,
       "all_original_photo_ids":total,
       "all_four_hue_classified_photos":rec["all_four_hue_classified_photos"],
       "actually_classified_photos_with_six_real_climatic_normals":len(valid),
       "original_original_photo_environment_data_unchanged":True,
       "original_existing_248312_environment_sha256":sha(old_site),
       "expanded_original_photo_source_sha256":source_digest,
       "WorldClim_solar_raw_source_units":"kJ m^-2 day^-1 (legacy column name misleading)",
       "original_genus_bootstrap_draws":BOOT,
       "original_conservative_fourhue_times_six_anchors_maxT_family":24,
       "original_six_abiotic_climate_and_water_anchors":list(old.ANCHORS),
       "source_geographic_conditioning":"abs(lat), longitude sin/cos; additional genus×coarse spatial-cell sensitivity",
       "FCP_within_species_colour_polymorphism_estimands_reanalysed":False,
       "independent_photographic_source_replication":False,
       "original_genuine_photo_pixels_in_new_wave":True,
       "actual_photo_year_month_weather_from_1970_2000_climatology":False,
       "source_climatology_archives_sha256":{
          "bio":sha(worldclim/"bio"/"wc2.1_10m_bio.zip"),
          "srad":sha(worldclim/"srad"/"wc2.1_10m_srad.zip"),
          "vapr":sha(vapr_zip)},
       "source_month_vapr_raster_coverage":vapr_meta,
       "historical_global_original_genus_water_result_NOT_overwritten":True,
       "water_results":{
           k:old.extract_water(v) for k,v in results.items()},
       "models":results,
    }
    path=outdir/f"mois_{total}_real_original_photo_six_water_maxT999_receipt.json"
    path.write_text(json.dumps(summary,ensure_ascii=False,indent=2,
                               sort_keys=True,allow_nan=False)+"\n")
    rows=[]
    for dataset,water in summary["water_results"].items():
        for key,v in water.items():
            if v.get("status")=="HOLD":continue
            if "point_pp_per_original_photo_env_sd" in v:
                rows.append({
                   "photo_dataset":dataset,"comparison":key,
                   "white_source_species_SD_pp":v["point_pp_per_original_photo_env_sd"],
                   "maxT95_low":v["scope_anchor_hue_maxT_ci95_pp"][0],
                   "maxT95_high":v["scope_anchor_hue_maxT_ci95_pp"][1],
                   "maxT_family_excludes_zero":v["maxT_excludes_zero"],
                })
    with (outdir/f"mois_{total}_water_white_summary.csv").open(
        "x",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n")
        w.writeheader();w.writerows(rows)
    print("MOIS_SAME_ORIGINAL_100K_MOISTURE_PROTOCOL_NEW_REAL_PHOTOS",
          json.dumps({"all_photo_ids":total,"total_classified":rec["all_four_hue_classified_photos"],
                      "new_classified":hue_new,
                      "maxT_white_water":summary["water_results"]["newly_expanded_MOIS_global_photos"]},
                     sort_keys=True),flush=True)
    return summary


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--expanded-photo-union",type=Path,required=True)
    p.add_argument("--expanded-photo-receipt",type=Path,required=True)
    p.add_argument("--original-248312-photo-worldclim",type=Path,required=True)
    p.add_argument("--original-248312-worldclim-receipt",type=Path,required=True)
    p.add_argument("--worldclim-raster-dir",type=Path,required=True)
    p.add_argument("--worldclim-vapr-dir",type=Path,required=True)
    p.add_argument("--worldclim-vapr-zip",type=Path,required=True)
    p.add_argument("--number-of-new-4096-waves",type=int,default=1)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    fit(a.expanded_photo_union,a.expanded_photo_receipt,
        a.original_248312_photo_worldclim,a.original_248312_worldclim_receipt,
        a.worldclim_raster_dir,a.worldclim_vapr_dir,a.worldclim_vapr_zip,
        a.output,expected_waves=a.number_of_new_4096_waves)
