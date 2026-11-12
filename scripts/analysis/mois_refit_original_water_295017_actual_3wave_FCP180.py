#!/usr/bin/env python3
"""Refit unchanged MOIS six-anchor, four-hue genus-999/maxT on 295017 REAL photos.

Inputs contain already completed, SHA-verified original photo ROIs. Only the
8191 photo IDs absent from the separately measured 286826-photo FCP180 source
need their original site WorldClim normals sampled. The source photo, hue,
species and site identity must agree for every prior photograph. Same frozen
historical genus estimators and 999-bootstrap seeds; FCP W/D within-species
polymorphism never imported. This is NOT a seasonal weather or anthesis model.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from scripts.analysis import mois_refit_original_water_with_actual_new_original_photos as legacy

OLD_N, OLD_HUES = 286826, 139405
NEW_N, NEW_HUES = 8191, 3133
FULL_N, FULL_HUES = OLD_N+NEW_N,OLD_HUES+NEW_HUES
FULL_SHA = "44e3cb6b4c6174dfb1035de208343c7c6a2e6c321e2557a509224c4d9f2ac7b3"
OLD_ENV_SHA = "338fbb0d50837c9cc25aece2dfc15d0bbb02e8b192728b8adf030c392c0f41f0"
HUES=set(("white","yellow_orange","red_pink","blue_purple"))
PHOTO_NAME="mois_real295017_original_photo_pixels_three_waves_FCP180.csv.gz"
MODEL_NAME="mois_295017_real_original_photo_six_water_maxT999_receipt.json"


def sha_file(path:Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for x in iter(lambda:f.read(1<<20),b""):
            h.update(x)
    return h.hexdigest()


def require(ok:bool,message:str):
    if not ok:
        raise ValueError(message)


def verify(original:Path,photo_receipt:Path,prior:Path,prior_receipt:Path) -> dict:
    """Authenticate actual photo denominators and literal ID/hue base support."""
    require(sha_file(original)==FULL_SHA,"Latest 295017 real original-photo CSV SHA changed")
    require(sha_file(prior)==OLD_ENV_SHA,"Original 286826 WorldClim photo-level source SHA changed")
    r=json.loads(photo_receipt.read_text())
    prior_r=json.loads(prior_receipt.read_text())
    require(r["schema"]=="mois_295017_genuine_all_original_photos_three4096_waves_FCP180_date_enriched_v1"
            and r["original_measured_photographs"]==FULL_N
            and r["original_measured_photo_four_hue_labels"]==FULL_HUES
            and r["files_SHA256"][PHOTO_NAME]==FULL_SHA
            and r["source_shared_conflicts"]["inat_taxon_id_source_conflicts"]==0
            and r["source_shared_conflicts"]["morph_source_conflicts"]==0
            and r["source_shared_conflicts"]["observed_on_source_conflicts"]==0,
            "Original 295017-photo measured pixel/photo-ID metadata differs from signed source")
    require(prior_r["n_original_photo_ids"]==OLD_N
            and prior_r["n_four_hue_classified"]==OLD_HUES
            and prior_r["full_photo_static_monthly_normal_csv_sha256"]==OLD_ENV_SHA,
            "Previous 286826 original environment + four-hue source drift")
    d=pd.read_csv(original,compression="gzip",low_memory=False,
       usecols=["photo_id","inat_taxon_id","species","morph","source_original_measured",
                "latitude","longitude","observed_on"],
       dtype={"photo_id":str,"inat_taxon_id":str})
    e=pd.read_csv(prior,compression="gzip",low_memory=False,
       usecols=["photo_id","inat_taxon_id","species","morph"],
       dtype={"photo_id":str,"inat_taxon_id":str})
    require(len(d)==FULL_N and d.photo_id.is_unique and
            int(d.morph.isin(HUES).sum())==FULL_HUES and
            len(e)==OLD_N and e.photo_id.is_unique and
            int(e.morph.isin(HUES).sum())==OLD_HUES,
            "Real source photograph IDs and classifier denominators changed")
    matched=e.merge(d[["photo_id","inat_taxon_id","species","morph"]],
                    on="photo_id",how="inner",validate="one_to_one",
                    suffixes=("_old","_new"))
    require(len(matched)==OLD_N,"Original 286826 photographed IDs not retained in the full source")
    for field in ("inat_taxon_id","species","morph"):
        same=matched[field+"_old"].fillna("").astype(str).eq(
             matched[field+"_new"].fillna("").astype(str))
        require(same.all(),f"Previously actually measured original photograph {field} was overwritten")
    new=d.loc[~d.photo_id.isin(e.photo_id)]
    require(len(new)==NEW_N and new.photo_id.nunique()==NEW_N
            and int(new.morph.isin(HUES).sum())==NEW_HUES
            and new.observed_on.notna().all(),
            "Missing new original MOIS photo pixel/colour/date terminal")
    require(new.source_original_measured.isin([
                "mois_wave2_new_original_ROI_v4",
                "mois_wave3_new_original_ROI_v4"
            ]).all(),
            "Additional 8191 photos must come only from measured MOIS wave2/3 pixel sources")
    return {
      "all_original_unique_photos":FULL_N,
      "new_original_source_photos":NEW_N,
      "all_four_hue_classified_photos":FULL_HUES,
      "new_four_hue_measured_photo_ids":NEW_HUES,
      "original_model_SHA_verified_for_every_photo":True,
      "all_real_source_photo_hue_csv_sha256":FULL_SHA,
      "origin_286826_source_hue_photo_rows":OLD_N,
      "original_295017_recovered_dates_not_fitted_as_causal_month_weather":True,
      "original_FCP_within_species_W_D_or_50km_not_imported":True,
    }


def run(full_photo:Path,full_receipt:Path,old_env:Path,old_env_receipt:Path,
        raster_dir:Path,vapr_dir:Path,vapr_zip:Path,outdir:Path) -> dict:
    require(not outdir.exists() or not any(outdir.iterdir()),
            "No overwrite of source-locked original MOIS result")
    source=verify(full_photo,full_receipt,old_env,old_env_receipt)
    outdir.mkdir(parents=True,exist_ok=True)
    stage=outdir/"source295017_unmodified_published_ROI_adapter.json"
    stage.write_text(json.dumps(source,indent=2,sort_keys=True)+"\n")
    # Existing tested original estimator demands BASE + WAVE * expected_waves.
    # Apply those numbers only to this new one-shot Python PROCESS, with its
    # statistic code, original standardized regressions, BOOT seed and maxT
    # family unchanged. This is not a redefinition of the source original 248k.
    legacy.BASE_PHOTOS=OLD_N
    legacy.BASE_HUES=OLD_HUES
    legacy.WAVE=NEW_N
    legacy.BOOT=999
    stats=legacy.fit(
        full_photo,stage,old_env,old_env_receipt,
        raster_dir,vapr_dir,vapr_zip,outdir/"original_295017_maxT",
        expected_waves=1)
    require(stats["all_original_photo_ids"]==FULL_N and
            stats["all_four_hue_classified_photos"]==FULL_HUES and
            stats["original_genus_bootstrap_draws"]==999 and
            stats["original_conservative_fourhue_times_six_anchors_maxT_family"]==24 and
            not stats["FCP_within_species_colour_polymorphism_estimands_reanalysed"],
            "Original statistical protocol or input actual photo count changed")
    source_receipt=outdir/"original_295017_maxT"/MODEL_NAME
    require(source_receipt.is_file(),"Actual 295017 photo water estimator produced no original signed output")
    white=stats["water_results"]["newly_expanded_MOIS_global_photos"]
    def extract(stage,anchor):
        t=white[f"{stage}__{anchor}"]
        require(t.get("status")!="HOLD" and
                t.get("point_pp_per_original_photo_env_sd") is not None,
                "Original 295017 water model unusable in original genus/geo population")
        return t
    result={
      "schema":"mois_genuine_295017_photo_original_six_axis_water_maxT999_v1",
      "actual_measured_photo_ids":FULL_N,
      "automatically_assigned_four_hue_ROI_labels":FULL_HUES,
      "new_actual_original_photo_ids_after_FCP180":NEW_N,
      "new_actual_classifiable_ROIs_after_FCP180":NEW_HUES,
      "full_real_photo_source_SHA256":FULL_SHA,
      "original_286826_photo_site_source_SHA256":OLD_ENV_SHA,
      "unchanged_original_water_model_receipt_SHA256":sha_file(source_receipt),
      "full_295017_input_site_plus_real_new_climate_csv_sha256":
          stats["expanded_original_photo_source_sha256"],
      "genus_bootstrap_draws":999,
      "six_annual_worldclim_climate_anchors":stats["original_six_abiotic_climate_and_water_anchors"],
      "simultaneous_maxT_family_hue_x_climate":24,
      "global_BIO12_white_genus_equal":extract("genus_equal","wc_bio12"),
      "global_vapour_white_genus_equal":extract("genus_equal","wc_vapr_annual_mean"),
      "coarse_cell_BIO12_white":extract("coarse_cell","wc_bio12"),
      "coarse_cell_vapour_white":extract("coarse_cell","wc_vapr_annual_mean"),
      "observed_year_month_weather_model":False,
      "photo_specific_open_flower_petals_independently_verified":False,
      "photo_pixel_measurement_recomputed_by_this_job":False,
      "all_FCP_256_original_pixel_partitions_complete":False,
      "FCP_intraspecific_D_W_or_50km_polymorphism_imported":False,
      "causal_pigment_fitness_inference":False,
      "source_processing_completion_order_not_uniform_sampling":True,
      "historical_100543_genus_moisture_benchmark_not_overwritten":True,
    }
    (outdir/"mois_295017_actual_original_water_decision.json").write_text(
        json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2,allow_nan=False)+"\n")
    print("MOIS_295017_GENUINE_PHOTO_WATER_999_GENUS_MAXT",
          json.dumps({k:result[k] for k in (
              "actual_measured_photo_ids",
              "automatically_assigned_four_hue_ROI_labels",
              "global_BIO12_white_genus_equal",
              "coarse_cell_BIO12_white")},sort_keys=True),flush=True)
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--real-295017-original-photo",type=Path,required=True)
    p.add_argument("--real-295017-photo-receipt",type=Path,required=True)
    p.add_argument("--original-286826-worldclim",type=Path,required=True)
    p.add_argument("--original-286826-climate-receipt",type=Path,required=True)
    p.add_argument("--worldclim-raster-root",type=Path,required=True)
    p.add_argument("--official-vapr-root",type=Path,required=True)
    p.add_argument("--official-vapr-zip",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    run(args.real_295017_original_photo,args.real_295017_photo_receipt,
        args.original_286826_worldclim,args.original_286826_climate_receipt,
        args.worldclim_raster_root,args.official_vapr_root,args.official_vapr_zip,
        args.output)
