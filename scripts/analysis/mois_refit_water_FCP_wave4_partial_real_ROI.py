#!/usr/bin/env python3
"""Re-run original MOIS 100543 photo water maxT999 on actually finished FCP ROI subset.

The FCP 49,800 source is INCOMPLETE in this snapshot. The original genuine
photo hue class assignments are admitted partition-by-partition, with prior
248,312 photo ID and original site BIO untouched. This retains the same old
four-hue six-climate/water-anchor genus-equal 999 bootstrap and familywise
maxT; its result is explicitly provisional for execution-order sampling.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from scripts.analysis import mois_refit_original_water_with_actual_new_original_photos as water
from scripts.analysis.mois_reuse_fcp_wave4_completed_roi_partitions import (
    BASE_PHOTOS,HUES,digest,require
)

ORIGINAL_BASE=248312
ORIGINAL_FOUR_HUE=120807
ORIGINAL_NEW_4096=4096


def apply(partial_photo:Path,partial_receipt:Path,site:Path,site_receipt:Path,
          climatology:Path,vapr:Path,vapr_zip:Path,outdir:Path):
    require(not outdir.exists() or not any(outdir.iterdir()),
            "Do not overwrite the first historical water-macro-source analysis")
    ledger=json.loads(partial_receipt.read_text())
    require(ledger["schema"]=="mois_authentic_completed_original_FCP_wave4_ROI_PARTITIONS_source_v1"
            and ledger["FCP_completed_ROI_partitions"]>=16
            and ledger["FCP_completed_ROI_partitions"]<=256
            and ledger["original_MOIS_prior_photo_IDs"]==BASE_PHOTOS
            and ledger["original_MOIS_prior_classified_hues"]==122383
            and ledger["combined_unique_photo_IDs"]>=BASE_PHOTOS
            and ledger["combined_four_hue_classified_photo_IDs"]>=122383
            and digest(partial_photo)==ledger["source_photo_level_csv_sha256"]
            and ledger["FCP_model_source_frozen"] is True
            and ledger["FCP_genetic_or_within_species_D_W_50km_polymorphism_statistics_imported"] is False,
            "Only sealed FCP original RGB terminal IDs can reach the MOIS water analysis")
    total=ledger["combined_unique_photo_IDs"]
    classified=ledger["combined_four_hue_classified_photo_IDs"]
    delta=total-ORIGINAL_BASE
    true_new_classified=classified-ORIGINAL_FOUR_HUE
    require(delta>=ORIGINAL_NEW_4096 and 0<true_new_classified<=delta,
            "New photo ID/colour additions to original 248k are invalid")
    outdir.mkdir(parents=True,exist_ok=True)
    normalized=outdir/"admitted_genuine_pixels_original_water_source_adaptor.json"
    normalized.write_text(json.dumps({
       "all_original_unique_photos":total,
       "new_original_source_photos":delta,
       "all_four_hue_classified_photos":classified,
       "new_four_hue_measured_photo_ids":true_new_classified,
       "original_model_SHA_verified_for_every_photo":True,
       "all_real_source_photo_hue_csv_sha256":ledger["source_photo_level_csv_sha256"],
       "underlying_FCP_real_partition_sample_source_sha256":digest(partial_receipt),
       "provisional_execution_order_subset":ledger["status"]!="COMPLETE_ACTUAL_PIXELS",
       "complete_FCP_wave4_49800_source_not_implied":ledger["FCP_completed_ROI_partitions"]<256,
    },indent=2,sort_keys=True)+"\n")
    # The count change affects only the cohort size. The original six-anchor
    # functions, 999 bootstrap algorithm, colours, geographic adjustment, and
    # maxT family are exactly the previously tested frozen methods.
    original_wave_count=water.WAVE
    water.WAVE=delta
    nested=outdir/"original_six_anchor_water_refit"
    try:
        estimate=water.fit(
            partial_photo,normalized,site,site_receipt,climatology,vapr,vapr_zip,
            nested,expected_waves=1)
    finally:
        water.WAVE=original_wave_count
    require(estimate["all_original_photo_ids"]==total and
            estimate["all_four_hue_classified_photos"]==classified and
            estimate["original_genus_bootstrap_draws"]==999 and
            estimate["original_conservative_fourhue_times_six_anchors_maxT_family"]==24,
            "Original 100543-photo moisture estimator was not actually reproduced")
    result={
       "schema":"mois_partial_FCP_completed_original_photo_water_999genus_maxT24_v1",
       "status":"AUTHENTIC_FINISHED_FCP_IMAGE_PIXELS_ORIGINAL_MOIS_WATER_REFITTED",
       "original_FCP_source_full_photo_count":49800,
       "original_FCP_completed_pixel_partitions":ledger["FCP_completed_ROI_partitions"],
       "original_FCP_total_pixel_partitions":256,
       "FCP_original_real_photo_rows_admitted_from_parts":ledger[
            "additional_new_original_photo_IDs"],
       "FCP_original_real_four_hue_photo_rows_admitted_from_parts":ledger[
            "additional_new_real_four_hue_classified_photo_IDs"],
       "actual_union_original_photo_IDs":total,
       "actual_union_four_hue_classified_photo_IDs":classified,
       "water_model_actual_genus_bootstrap_draws":999,
       "water_model_four_hue_times_six_anchors_familywise_maxT":24,
       "source_based_image_state_not_field_verified_open_petals":True,
       "water_rain_vapor_coefficients":estimate[
           "water_results"]["newly_expanded_MOIS_global_photos"],
       "FCP_partial_worker_completion_order_is_not_outcome_neutral_proof":True,
       "FCP_within_species_D_W_local_spatial_polymorphism_reanalysis":False,
       "original_first100543_photo_water_results_unchanged":True,
       "original_248312_existing_water_results_unchanged":True,
       "original_FCP_partial_source_receipt_sha256":digest(partial_receipt),
       "original_Water999_source_receipt_sha256":digest(
           nested/f"mois_{total}_real_original_photo_six_water_maxT999_receipt.json"),
    }
    (outdir/"mois_FCP_wave4_partial_real_ROI_original_water_maxT999_decision.json").write_text(
        json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False,
                   allow_nan=False)+"\n")
    rain=result["water_rain_vapor_coefficients"]["genus_equal__wc_bio12"]
    vapour=result["water_rain_vapor_coefficients"]["genus_equal__wc_vapr_annual_mean"]
    print("MOIS_PARTIAL_REAL_FCP_RGB_PIXELS_SAME_ORIGINAL_999_WATER",json.dumps({
         "photos":total,"classified":classified,
         "water_rain_pp":rain["point_pp_per_original_photo_env_sd"],
         "water_rain_maxT95":rain["scope_anchor_hue_maxT_ci95_pp"],
         "water_vapor_pp":vapour["point_pp_per_original_photo_env_sd"],
         "water_vapor_maxT95":vapour["scope_anchor_hue_maxT_ci95_pp"],
         "source_completed_partitions":ledger["FCP_completed_ROI_partitions"]},
         sort_keys=True),flush=True)
    return result


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--authentic-partial-FCP-photo-source",type=Path,required=True)
    ap.add_argument("--authentic-partial-FCP-photo-receipt",type=Path,required=True)
    ap.add_argument("--source-original-248312-site-worldclim",type=Path,required=True)
    ap.add_argument("--source-original-248312-worldclim-receipt",type=Path,required=True)
    ap.add_argument("--worldclim-raster-root",type=Path,required=True)
    ap.add_argument("--official-monthly-vapr-raster-root",type=Path,required=True)
    ap.add_argument("--official-monthly-vapr-zip",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args()
    apply(a.authentic_partial_FCP_photo_source,
          a.authentic_partial_FCP_photo_receipt,
          a.source_original_248312_site_worldclim,
          a.source_original_248312_worldclim_receipt,
          a.worldclim_raster_root,
          a.official_monthly_vapr_raster_root,
          a.official_monthly_vapr_zip,a.output)
