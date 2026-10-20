#!/usr/bin/env python3
"""Append source-signed p0004 genuine image ROI outcomes to immutable FCP+p0002 photos.

No image downloads, colour re-inference, species nativeness or ecological models.
The input 93m photo index is metadata; ONLY real image ROI attempts are appended.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from scripts.analysis.mois_actual_photo_cohort_season_audit_20261019 import (
    HUES, require, sha256,
)
from scripts.acquisition.mois_union_FCP447500_and_first93m18462_actual_photos import (
    _check_pixels,
    RESULT_PARQUET as OLD_PARQUET,
    RESULT_RECEIPT as OLD_RECEIPT,
)

OLD_N = 612483
OLD_CLASSIFIED = 295801
NEW_N = 18119
TOTAL_N = 630602
FULL_PHOTO_INDEX = 93209351
OLD_SCHEMA = "mois_original_photo_ROI_FCP447500_plus_first_full93m_p0002_ID_union_v1"
SOURCE_SCHEMA = "mois_full_source_observation_sampling_cohort_v1"
MEASURED_SCHEMA = "mois_global93m_observation_baseline_real_RGB_photo_results_v1"
OUT_PHOTOS = "mois_FCP447500_plus93m_p0002_p0004_unique_actual_photo_ROIs.parquet"
OUT_RECEIPT = "mois_FCP447500_plus93m_p0004_source_signed_pixel_union_receipt.json"
OUT_COUNTS = "mois_FCP447500_plus93m_p0004_actual_hue_source_arm_counts.csv"


def verified_sources(old_path, old_receipt, new_path, new_receipt,
                     new_source, first_receipt):
    before = json.loads(Path(old_receipt).read_text(encoding="utf-8"))
    added = json.loads(Path(new_receipt).read_text(encoding="utf-8"))
    plan = json.loads(Path(new_source).read_text(encoding="utf-8"))
    first = json.loads(Path(first_receipt).read_text(encoding="utf-8"))
    require(before.get("schema") == OLD_SCHEMA
            and before.get("all_genuinely_attempted_unique_original_photo_IDs") == OLD_N
            and before.get("all_original_photo_four_hue_classified") == OLD_CLASSIFIED
            and before.get("duplicate_source_photo_IDs_in_final_union") == 0
            and sha256(old_path) == before["completed_output_sha256"][OLD_PARQUET],
            "Previously measured FCP+p0002 photo pixels must match original SHA release")
    require(added.get("schema") == MEASURED_SCHEMA
            and plan.get("schema") == SOURCE_SCHEMA
            and first.get("schema") == MEASURED_SCHEMA
            and plan.get("full_frame_photo_ids") == FULL_PHOTO_INDEX
            and added.get("full_original_upstream_photo_candidates") == FULL_PHOTO_INDEX
            and plan.get("selection_unit") == "observation_uuid"
            and plan.get("selection_probability") == 0.0004
            and added.get("observation_level_selection_probability") == 0.0004
            and added.get("attempted") == plan.get("new_photo_attempts_planned") == NEW_N
            and plan.get("sampled_photo_ids_before_prior_dedup") == 36758
            and plan.get("previously_measured_photo_ids_skipped") == 18639
            and added.get("matching_prior_measured_photos_reused") == 18639
            and plan.get("selected_new_photos_unsupported_by_frozen_loader") == 0
            and added.get("technical_unknowns_not_assumed_nonwhite") is True
            and added.get("ecological_analysis_run") is False
            and added.get("raw_licensed_photo_pixels_redistributed") is False
            and added.get("frozen_original_ROI_method") == first.get("frozen_original_ROI_method")
            and sha256(new_path) == added["measured_photo_parquet_sha256"]
            and sha256(new_source) == added["source_plan_sha256"],
            "p0004 actual-image ROI identity, probability, SHA or frozen method changed")
    old = pd.read_parquet(old_path)
    new = pd.read_parquet(new_path)
    return old, new, before, added, plan


def append(old, new, added, plan, *, old_n=OLD_N, new_n=NEW_N):
    _check_pixels(old, "immutable previously measured FCP+p0002 pixels")
    require(len(old) == old_n and int(old.classified.sum()) >= 0,
            "Unexpected previous original-photo pixels")
    required = {
        "photo_id", "inat_taxon_id", "observation_id", "species",
        "morph", "measurement_status", "classified_four_hue", "source_arm",
        "latitude", "longitude", "observed_on",
    }
    require(required <= set(new), f"p0004 source missing original image columns {required-set(new)}")
    ids = new.photo_id.astype(str)
    require(len(new) == new_n and ids.str.fullmatch(r"\d+").all()
            and ids.is_unique and not ids.isin(set(old.photo_id.astype(str))).any(),
            "p0004 reused an old original photo or duplicated a new photo ID")
    require(new.source_arm.eq("global93m_probability_sample_observations").all()
            and new.morph.isin((*HUES, "mixed_uncertain")).all(),
            "p0004 is not authentic four-state photo ROI or unknown image")
    known = new.morph.isin(HUES)
    status = new.measurement_status.fillna("").astype(str)
    failed = status.eq("image_acquisition_failed")
    require(known.eq(status.eq("classified_four_state_morph")).all()
            and known.map({True: "true", False: "false"}).eq(
                new.classified_four_hue.astype(str).str.lower()).all()
            and not (known & failed).any(),
            "Unknown acquisition failures or ROI image colour states were relabelled")
    require(int(known.sum()) == added["classifiable"]
            and int(failed.sum()) == added["acquisition_failed"]
            and int((~known & ~failed).sum()) == added["unclassifiable_acquired"]
            and len(new) == added["attempted"]
            and added["downloaded"] + added["acquisition_failed"] == len(new)
            and plan["new_photo_attempts_planned"] == len(new),
            "Signed actual RGB labels, downloaded images and failures disagree")

    next_rows = pd.DataFrame(index=new.index)
    for column in old.columns:
        next_rows[column] = new[column] if column in new else pd.NA
    next_rows["photo_id"] = ids
    next_rows["classified"] = known.to_numpy()
    next_rows["white"] = new.morph.eq("white").to_numpy()
    next_rows["technical_failure"] = failed.to_numpy()
    next_rows["arm"] = "global93m_probability_observation_p0004"
    if "photo_colour_provenance" in next_rows:
        next_rows["photo_colour_provenance"] = "actual_original_93m_observation_p0004_ROI_v4"
    if "source_original_fcp_morph" in next_rows:
        next_rows["source_original_fcp_morph"] = ""
    if "fcp_source_cohort" in next_rows:
        next_rows["fcp_source_cohort"] = ""
    if "source_fcp_global_classifiable" in next_rows:
        next_rows["source_fcp_global_classifiable"] = False

    result = pd.concat([old, next_rows.loc[:, old.columns]], ignore_index=True)
    _check_pixels(result, "FCP+p0002+p0004 genuine unique original photo ROI union")
    require(len(result) == old_n + new_n and result.photo_id.is_unique,
            "Original photo-ID count failed during authentic ROI append")
    for col in ("photo_id", "inat_taxon_id", "morph", "classified", "white",
                "observation_id", "measurement_status"):
        require(result[col].iloc[:old_n].fillna("").astype(str).reset_index(drop=True).equals(
                old[col].fillna("").astype(str).reset_index(drop=True)),
                "Existing old photo pixel or source identity changed: " + col)
    require(int(result.classified.sum()) == int(old.classified.sum()) + int(known.sum()),
            "Only genuinely measured new four-hue ROIs may increase classifiable count")
    old_taxa = set(old.loc[old.classified, "inat_taxon_id"].astype(str))
    new_taxa = set(new.loc[known, "inat_taxon_id"].astype(str))
    report = {
        "schema": "mois_source_signed_original_FCP_p0002_p0004_real_RGB_union_v1",
        "prior_unique_actual_original_photo_IDs": old_n,
        "p0004_new_unique_actual_original_photo_IDs": new_n,
        "combined_unique_actual_original_photo_IDs": len(result),
        "prior_actual_four_hue_classified": int(old.classified.sum()),
        "p0004_new_actual_four_hue_classified": int(known.sum()),
        "combined_actual_four_hue_classified": int(result.classified.sum()),
        "p0004_acquisition_failed": int(failed.sum()),
        "p0004_acquired_unclassifiable": int((~known & ~failed).sum()),
        "newly_classified_nominal_taxa": len(new_taxa - old_taxa),
        "overlap_with_previously_processed_original_photo_IDs": 0,
        "all_prior_original_ROI_labels_preserved": True,
        "full_source_index_photo_IDs_NOT_equated_to_image_RGB_measurements": FULL_PHOTO_INDEX,
        "no_original_image_download_or_new_colour_inference_in_union": True,
        "unknown_p0004_RGB_ROIs_not_assumed_chromatic": True,
        "no_native_island_climate_pollinator_or_ecological_analysis_run": True,
    }
    return result, report


def run(old_path, old_receipt, new_path, new_receipt, new_source,
        first_receipt, output):
    output = Path(output)
    require(not output.exists() or not any(output.iterdir()),
            "Previously published genuine pixels are immutable")
    old, new, _, receipt, source = verified_sources(
        old_path, old_receipt, new_path, new_receipt, new_source, first_receipt)
    union, report = append(old, new, receipt, source)
    require(report["combined_unique_actual_original_photo_IDs"] == TOTAL_N,
            "Photo-ID union must not drift from frozen source exact counts")
    output.mkdir(parents=True, exist_ok=True)
    union.to_parquet(output / OUT_PHOTOS, index=False, compression="zstd")
    counts = []
    for arm, frame in union.groupby("arm", observed=True, sort=True):
        classified = frame.loc[frame.classified]
        counts.append({
            "source_arm": arm, "true_image_photo_attempts": len(frame),
            "classifiable_four_hue_ROIs": len(classified),
            "acquisition_failed": int(frame.technical_failure.fillna(False).sum()),
            "unknown_or_failed_ROIs": int((~frame.classified).sum()),
            **{h: int(classified.morph.eq(h).sum()) for h in HUES},
        })
    pd.DataFrame(counts).to_csv(output / OUT_COUNTS, index=False)
    report["source_input_sha256"] = {
        OLD_PARQUET: sha256(old_path), OLD_RECEIPT: sha256(old_receipt),
        "p0004_measured_photos.parquet": sha256(new_path),
        "p0004_receipt.json": sha256(new_receipt),
        "p0004_source.json": sha256(new_source),
        "p0002_receipt.json": sha256(first_receipt),
    }
    report["output_sha256"] = {
        OUT_PHOTOS: sha256(output / OUT_PHOTOS),
        OUT_COUNTS: sha256(output / OUT_COUNTS),
    }
    (output / OUT_RECEIPT).write_text(json.dumps(
        report, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False) + "\n")
    print("MOIS_FCP_P0002_P0004_GENUINE_RGB_UNION",
          json.dumps({key: report[key] for key in (
              "prior_unique_actual_original_photo_IDs",
              "p0004_new_unique_actual_original_photo_IDs",
              "combined_unique_actual_original_photo_IDs",
              "prior_actual_four_hue_classified",
              "p0004_new_actual_four_hue_classified",
              "combined_actual_four_hue_classified",
              "newly_classified_nominal_taxa",
          )}, sort_keys=True), flush=True)
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("old-parquet", "old-receipt", "new-parquet", "new-receipt",
                "new-source", "first-receipt", "output"):
        p.add_argument("--" + key, type=Path, required=True)
    args = p.parse_args()
    run(args.old_parquet, args.old_receipt, args.new_parquet,
        args.new_receipt, args.new_source, args.first_receipt, args.output)


if __name__ == "__main__":
    main()
