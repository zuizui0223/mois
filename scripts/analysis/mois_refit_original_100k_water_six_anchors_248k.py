#!/usr/bin/env python3
"""MOIS: refit the ORIGINAL 100543-photo global six-anchor moisture protocol.

Reuses the historical published genus-conditional 999-bootstrap/maxT source
algorithm, now on already classified 248312 distinct photo-ID opportunities
(120807 four-hue labels). FCP original photos are allowed; FCP's within-species
D/W and spatial polymorphism estimands are NEVER computed or imported.
This is a new data-scope analysis, not an independent biological replication.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio

from scripts.analysis import fcp_abiotic_hue_direction_macro_20261013 as original
from scripts.analysis.mois_sample_248312_worldclim_static_and_monthnormals import sample_tiff

HUES=original.CLASSES
GEO=original.GEO
ANCHORS=original.ANCHOR_CLIMATE
EXPECTED_PHOTOS=248312
EXPECTED_CLASSIFIED=120807
BOOTSTRAPS=999
VAPR_NAMES=tuple(f"wc2.1_10m_vapr_{i:02d}.tif" for i in range(1,13))
READ_COLS=(
    "photo_id","inat_taxon_id","species","morph","original_source",
    "latitude","longitude","wc2_bio05_10m","wc2_bio12_10m",
    "wc2_bio14_10m","wc2_bio15_10m","wc2_annual_srad_normal_10m_Wm2",
)


def sha(path:Path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for buf in iter(lambda:f.read(1024*1024),b""):
            h.update(buf)
    return h.hexdigest()


def require(ok,msg):
    if not ok: raise ValueError(msg)


def annual_vapor(folder:Path,lon:np.ndarray,lat:np.ndarray):
    totals=np.zeros(len(lon),dtype=np.float64)
    count=np.zeros(len(lon),dtype=np.int16)
    raster_info=[]
    for name in VAPR_NAMES:
        files=list(folder.rglob(name))
        require(len(files)==1,f"Missing or duplicated original WorldClim vapr raster: {name}")
        with rasterio.open(files[0]) as dataset:
            samples=sample_tiff(files[0],lon,lat)
            ok=np.isfinite(samples)
            totals[ok]+=samples[ok]
            count[ok]+=1
            raster_info.append({
                "name":name,"grid_res_degrees":round(float(abs(dataset.transform.a)),8),
                "valid_photo_sites":int(ok.sum()),"raster_no_data":dataset.nodata
            })
    out=np.full(len(lon),np.nan,dtype=np.float32)
    ok=count==12
    out[ok]=(totals[ok]/12).astype(np.float32)
    require(int(ok.sum())>=200000,
            "Photo-site annual vapour coverage too sparse for the historical original population")
    return out,raster_info


def prepare_data(source:Path,source_receipt:Path,vapr_folder:Path):
    rec=json.loads(source_receipt.read_text())
    require(sha(source)==rec["full_photo_static_monthly_normal_csv_sha256"],
            "Photo-level original location/colour + 19 BIO/solar source changed")
    require((rec["n_original_photo_ids"],rec["n_four_hue_classified"])==
            (EXPECTED_PHOTOS,EXPECTED_CLASSIFIED),
            "Original 248312/120807 cohort changed")
    frame=pd.read_csv(source,compression="gzip",usecols=list(READ_COLS),
                      dtype={"photo_id":str,"inat_taxon_id":str},low_memory=False)
    require(len(frame)==EXPECTED_PHOTOS and frame.photo_id.nunique()==len(frame),
            "Original image IDs missing or duplicated")
    require(int(frame.morph.isin(HUES).sum())==EXPECTED_CLASSIFIED,
            "Wrong original photographic colour denominator")
    lat=pd.to_numeric(frame.latitude,errors="coerce").to_numpy(float)
    lon=pd.to_numeric(frame.longitude,errors="coerce").to_numpy(float)
    vapr,details=annual_vapor(vapr_folder,lon,lat)
    for our,old in (
        ("wc2_bio05_10m","wc_bio5"),("wc2_bio12_10m","wc_bio12"),
        ("wc2_bio14_10m","wc_bio14"),("wc2_bio15_10m","wc_bio15"),
        ("wc2_annual_srad_normal_10m_Wm2","wc_srad_annual_mean")
    ):
        frame[old]=pd.to_numeric(frame[our],errors="coerce")
    frame["wc_vapr_annual_mean"]=vapr
    frame["latitude"]=lat
    frame["longitude"]=lon
    frame["abs_lat"]=np.abs(lat)
    frame["lon_sin"]=np.sin(np.deg2rad(lon))
    frame["lon_cos"]=np.cos(np.deg2rad(lon))
    frame["genus"]=frame.species.astype("string").str.split().str[0]
    photo_ok=frame.morph.isin(HUES)
    clean=frame.loc[photo_ok].dropna(subset=list((*GEO,*ANCHORS))).copy()
    # This retains the original 100k estimator's one-species-one-vote
    # aggregation and genus+geographic conditioning.
    require(len(clean)>100000 and clean.inat_taxon_id.nunique()>20000,
            "Historical source photo support disappeared after six-anchor join")
    return frame,clean,details,rec


def estimate(photo:pd.DataFrame,means:pd.Series,scales:pd.Series,
             *,kind:str,nboot=BOOTSTRAPS)->dict:
    rows=original.species_frame(photo,ANCHORS)
    groups,meta=original.cluster_suff(rows,ANCHORS,means,scales,geos=True)
    require(groups is not None,
            f"Global source diversity insufficient for original genus method {kind}: {meta}")
    results={}
    for weight in ("genus_equal","species_equal"):
        result=original.regress_boot(
            groups,ANCHORS,nboot,
            original.BOOT_SEED+len("world")+len("climate")*11,
            equal_genus=(weight=="genus_equal"))
        require(result["status"]=="OK",
                f"Historic original six-anchor bootstrap not estimable: {weight}, {result}")
        results[weight]=result
    spatial_groups,spatial_meta=original.cluster_suff(
        rows,ANCHORS,means,scales,geos=False,cell_conditioned=True)
    if spatial_groups is None:
        spatial={"status":"HOLD_UNSUPPORTED_CONGENERIC_COARSE_CELL","support":spatial_meta}
    else:
        spatial=original.regress_boot(
            spatial_groups,ANCHORS,nboot,original.BOOT_SEED+len("world")+900,
            equal_genus=True,cell_conditioned=True)
        spatial["support"]=spatial_meta
    return {
        "photo_population":kind,"photo_rows_with_six_complete_anchors":int(len(photo)),
        "source_species":int(rows.inat_taxon_id.nunique()),
        "qualified_genus_and_species_counts":meta,
        "original_method_FCP_macro_six_anchor_genus_equal":results["genus_equal"],
        "original_method_FCP_macro_six_anchor_species_equal":results["species_equal"],
        "same_genus_coarse_cell_original_method_sensitivity":spatial,
        "original_photo_source_observations_are_not_independent_provider_replication":True
    }


def extract_water(est:dict):
    out={}
    for stage,model in (
        ("genus_equal",est["original_method_FCP_macro_six_anchor_genus_equal"]),
        ("species_equal",est["original_method_FCP_macro_six_anchor_species_equal"]),
        ("coarse_cell",est["same_genus_coarse_cell_original_method_sensitivity"])
    ):
        for axis in ("wc_bio12","wc_vapr_annual_mean"):
            key=f"{stage}__{axis}"
            if model.get("status")!="OK":
                out[key]={"status":model.get("status","HOLD"),"not_reported_positive":True}
                continue
            found=[r for r in model["estimates"] if r["anchor"]==axis and r["hue"]=="white"]
            require(len(found)==1,"Original 100543 study white-water direction not uniquely identified")
            out[key]=found[0]
    return out


def run(source:Path,source_receipt:Path,vapr_folder:Path,vapr_zip:Path,
        outdir:Path,nboot:int=BOOTSTRAPS):
    require(nboot==BOOTSTRAPS,"Historical-style confidence comparison requires 999 genus bootstrap draws")
    require(not outdir.exists() or not any(outdir.iterdir()),
            "Do not overwrite original six-anchor photographic water receipt")
    full,clean,raster_meta,source_info=prepare_data(source,source_receipt,vapr_folder)
    historical=clean[clean.original_source.eq("historical")].copy()
    # Normalize from the SAME full-source set so comparisons of photo additions
    # are on one common photographic exposure-standard-deviation scale.
    norm=(*GEO,*ANCHORS)
    means=clean[list(norm)].mean()
    scales=clean[list(norm)].std(ddof=0).where(lambda s:s>1e-12,1.)
    fits={
        "historical_first100543_original_photos_reestimated":
            estimate(historical,means,scales,kind="historical_source_only",nboot=nboot),
        "expanded248312_original_photo_union":
            estimate(clean,means,scales,kind="historical_plus_FCP_photo_labels",nboot=nboot),
    }
    summary={
        "schema":"mois_full_248312_same_historical_six_anchor_global_moisture_999_genus_bootstrap_v1",
        "status":"REESTIMATED_SIX_ORIGINAL_ANCHORS_WITH_REAL_EXPANDED_248312_PHOTO_COLOURS",
        "original_historical_photos":100543,
        "source_original_photos":EXPECTED_PHOTOS,
        "photo_source_four_hue_classified":EXPECTED_CLASSIFIED,
        "full_source_six_anchor_classified_photos":len(clean),
        "historical_subset_six_anchor_classified_photos":len(historical),
        "six_original_named_anchors":list(ANCHORS),
        "original_genus_cluster_bootstrap_draws":nboot,
        "original_bootstrap_maxT_family_size":len(ANCHORS)*4,
        "same_original_macro_source_functions":[
            "fcp_abiotic_hue_direction_macro_20261013.cluster_suff",
            "fcp_abiotic_hue_direction_macro_20261013.regress_boot",
            "fcp_abiotic_hue_direction_macro_20261013.species_frame"
        ],
        "primary_white_water_outcomes":{key:extract_water(value) for key,value in fits.items()},
        "original_moisture_result_reference_BIO12_pp":2.254545,
        "new_model_is_not_exact_26covariate_original_sample":
          "Original six-anchor genus/maxT estimator reproduced; photo cohorts and WorldClim vapr-site coverage changed. Historical strict-photo 26-climate and soil quality subpopulations not reimputed.",
        "worldclim_2p1_original_vapr_zip_sha256":sha(vapr_zip),
        "worldclim_vapr_original_raster_coverage":raster_meta,
        "source_existing248312_worldclim_sha256":sha(source),
        "calendar_month_real_weather_not_in_this_annual_climate_model":True,
        "original_WorldClim_srad_metadata_warning":"Existing column labeled Wm2 but original WorldClim solar-radiation raster is in kJ m-2 day-1; source-photo standardized slopes are scale-invariant. Metadata must be corrected before claiming raw radiation units.",
        "original_existing_MOIS_results_NOT_overwritten":True,
        "FCP_within_species_polymorphism_D_W_50km_not_computed":True,
        "not_an_independent_photo_provider_replication_or_fitness_experiment":True,
        "models":fits,
    }
    outdir.mkdir(parents=True,exist_ok=True)
    target=outdir/"mois_248312_same_original_water_maxT999_genus_decision.json"
    target.write_text(json.dumps(summary,ensure_ascii=False,sort_keys=True,
                                indent=2,allow_nan=False)+"\n",encoding="utf-8")
    readable=[]
    for label,e in fits.items():
        for stage,val in extract_water(e).items():
            if val.get("point_pp_per_original_photo_env_sd") is not None:
                readable.append({
                    "population":label,"stage_and_water_variable":stage,
                    "point_pp_per_source_env_SD":val["point_pp_per_original_photo_env_sd"],
                    "maxT_family95_low":val["scope_anchor_hue_maxT_ci95_pp"][0],
                    "maxT_family95_high":val["scope_anchor_hue_maxT_ci95_pp"][1],
                    "familywise_maxT_excludes_zero":val["maxT_excludes_zero"],
                })
    import csv
    table=outdir/"mois_248312_original_six_anchor_water_white_summary.csv"
    with table.open("x",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(readable[0]),lineterminator="\n")
        w.writeheader();w.writerows(readable)
    (outdir/"mois_original_six_anchor_files_SHA.json").write_text(
        json.dumps({"decision_sha256":sha(target),"summary_sha256":sha(table)},indent=2)+"\n")
    print("MOIS_ORIGINAL_100K_METHOD_REAL_248K_WATER_SIX_ANCHORS_MAXT999",
          json.dumps({"available_original_photos":len(full),
                      "usable_classified_photos":len(clean),
                      "genus_equal_white_moisture":extract_water(
                          fits["expanded248312_original_photo_union"]
                      )},sort_keys=True),flush=True)
    return summary


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--real-248k-site-worldclim",type=Path,required=True)
    ap.add_argument("--site-worldclim-receipt",type=Path,required=True)
    ap.add_argument("--original-vapr-monthly-raster-dir",type=Path,required=True)
    ap.add_argument("--original-vapr-zip",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args()
    run(a.real_248k_site_worldclim,a.site_worldclim_receipt,
        a.original_vapr_monthly_raster_dir,a.original_vapr_zip,a.output)
