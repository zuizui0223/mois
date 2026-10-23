#!/usr/bin/env python3
"""Seal one completed R2 buckets4..11 worker group as signed append-only Parquet.

This consumes only already-finished tier06/all16 second 8-of-512 image worker artifacts. It
does NOT download source photos, reinfer flower colour, or reread historical
pixel Parquets. Each group becomes its own immutable data asset; other groups
may still be running. The existing frozen FCP-v4 collector provides the
per-image original photo/observation/terminal and model SHA proof.
"""
from __future__ import annotations

import argparse
from collections import Counter
import math
import re
import csv
import json
import os
from pathlib import Path
import shutil
import tempfile

from scripts.acquisition.mois_completed_pixel_checkpoint import (
    HUES, collect_wave, sha, require,
)

SCHEMA = "mois_tier06_all16_wave2_append_only_original_RGB_group_v1"
ORIGINAL_SOURCE_SHA = "57592ce9f12c8abb0adb7ce2fda8252044ca9204caf295dbac73ecbf76a01516"
ORIGINAL_RUN = 39919846830
EXPECTED_SOURCE_ATTEMPTS = 8765
EXPECTED_UPSTREAM_PHOTOS = 93_209_351



def recorded_month(value):
    """A recorded observation month is NOT a verified flowering month."""
    text = str(value or "").strip()
    m = re.fullmatch(r"([12]\d{3})-(0[1-9]|1[0-2])-(0[1-9]|[12]\d|3[01])", text)
    return m.group(1) + "-" + m.group(2) if m else "unknown"


def latitude_band(value):
    """Geographic description only, not a native/wild status assertion."""
    try:
        v = float(value)
    except (ValueError, TypeError):
        return "unknown"
    if not math.isfinite(v) or not -90 <= v <= 90:
        return "unknown"
    if -23.44 <= v <= 23.44:
        return "tropical"
    return "south_extratropical" if v < -23.44 else "north_extratropical"


def original_photo_source_coverage(rows):
    """Source-authenticated descriptive photo, taxon, date and latitude QA."""
    months = Counter()
    latitudes = Counter()
    taxa_attempted, taxa_classified = set(), set()
    labelled_observations = set()
    taxa_by_month = {}
    for row in rows:
        morph = row["morph"]
        source_taxon = str(row.get("inat_taxon_id") or "").strip()
        obs = str(row["observation_id"])
        is_classified = morph in HUES
        months[(recorded_month(row.get("observed_on")), morph)] += 1
        latitudes[(latitude_band(row.get("latitude")), morph)] += 1
        if source_taxon:
            taxa_attempted.add(source_taxon)
        if is_classified:
            labelled_observations.add(obs)
            if source_taxon:
                taxa_classified.add(source_taxon)
                taxa_by_month.setdefault(recorded_month(row.get("observed_on")), set()).add(source_taxon)
    return {
        "schema": "mois_wave2_real_photo_source_taxonomy_time_latitude_coverage_v1",
        "original_photo_attempts": len(rows),
        "distinct_nominal_inat_taxon_ids_attempted": len(taxa_attempted),
        "distinct_nominal_inat_taxon_ids_four_hue": len(taxa_classified),
        "distinct_original_observations_four_hue": len(labelled_observations),
        "recorded_observation_month_by_ROI_state": {
            f"{month}|{hue}": n for (month, hue), n in sorted(months.items())
        },
        "latitude_band_by_ROI_state": {
            f"{band}|{hue}": n for (band, hue), n in sorted(latitudes.items())
        },
        "distinct_nominal_taxon_ids_with_four_hue_by_recorded_month": {
            month: len(taxa) for month, taxa in sorted(taxa_by_month.items())
        },
        "source_taxon_ids_are_not_species_validated": True,
        "individual_photo_open_flowers_independently_verified": False,
        "individual_native_wild_status_independently_verified": False,
        "technical_unclassified_images_are_not_other_petal_hues": True,
        "source_arm_partial_probability_not_global_flora_prevalence": True,
    }


def seal(source, workers, group, output, *,
         expected_source_sha=ORIGINAL_SOURCE_SHA,
         expected_run=ORIGINAL_RUN,
         expected_attempts=EXPECTED_SOURCE_ATTEMPTS,
         expected_upstream=EXPECTED_UPSTREAM_PHOTOS,
         frozen_method=None):
    """Build an idempotent full-photo terminal partition, never a hue guess."""
    import pandas as pd
    from scripts.analysis.mois_reassemble_original_new4096_measurements import METHOD

    source, workers, output = map(Path, (source, workers, output))
    if frozen_method is None:
        frozen_method = METHOD
    require(not output.exists(), "Do not overwrite immutable real-photo group output")
    require(0 <= group <= 100_000, "Invalid original-photo work group number")
    manifest_file = source / "manifest.json"
    require(sha(manifest_file) == expected_source_sha,
            "Frozen original photo source SHA changed")
    plan = json.loads(manifest_file.read_text(encoding="utf-8"))
    require(plan.get("schema") == "mois_full_source_observation_sampling_cohort_v1"
            and plan.get("full_frame_photo_ids") == expected_upstream
            and plan.get("new_photo_attempts_planned") == expected_attempts
            and plan.get("source_tier") == 6
            and plan.get("source_group") == "all16"
            and plan.get("complete_tier06_source_photo_IDs") == 597070
            and plan.get("compute_bucket_count") == 512
            and plan.get("admitted_bucket_count") == 8
            and plan.get("admitted_bucket_start") == 4
            and plan.get("admitted_buckets") == list(range(4,12))
            and plan.get("previous_broad_original_photo_attempts") == 4468
            and plan.get("source_arm") == "global93m_tier06_all16_buckets004to011_broad_warm_R2"
            and plan.get("colour_measured") is False
            and plan.get("not_probability_complete_global_source") is True
            and plan.get("native_photo_status_not_inferred_or_hard_filtered") is True,
            "Source is not the intended signed original tier06/all16 cohort")
    expected_units = list(range(group * 4, min((group + 1) * 4, plan["work_units"])))
    require(expected_units and group in plan["group_ids"]
            and plan["worker_ids"] == list(range(plan["work_units"])),
            "Group was not scheduled in immutable original photo source")
    data, audit = collect_wave(source, workers, frozen_method)
    require(audit["completed_units"] == expected_units
            and audit["attempted"] == len(data) > 0
            and audit["planned_photo_attempts"] == expected_attempts,
            "Original group was not fully sealed or unrelated work units present")
    require(set(audit["missing_units"]) == set(plan["worker_ids"]) - set(expected_units),
            "Original source missing-unit roster is inconsistent")
    ids = [str(row["photo_id"]) for row in data]
    obs = [str(row["observation_id"]) for row in data]
    require(len(set(ids)) == len(ids) and all(x.isdecimal() for x in ids)
            and all(obs), "Duplicate/missing signed original photo or observation ID")
    for row in data:
        known = row["morph"] in HUES
        require(known == (row["measurement_status"] == "classified_four_state_morph")
                and row["source_arm"] == plan["source_arm"]
                and float(row["selection_probability"]) == plan["selection_probability"],
                "Original frozen hue, source arm or source observation probability changed")
    classified = sum(row["morph"] in HUES for row in data)
    require(classified == audit["classifiable"]
            and sum([audit["classifiable"], audit["acquisition_failed"],
                     audit["unclassifiable_acquired"]]) == len(data),
            "Original full-image unknown/failure composition altered")

    output.parent.mkdir(parents=True, exist_ok=True)
    temp = Path(tempfile.mkdtemp(prefix="mois-pixel-group-", dir=output.parent))
    try:
        table = pd.DataFrame(data)
        table["photo_id"] = table["photo_id"].astype(str)
        table = table.sort_values("photo_id", key=lambda s: s.map(int)).reset_index(drop=True)
        data_path = temp / "original_photo_ROIs.parquet"
        table.to_parquet(data_path, index=False, compression="zstd")
        coverage = original_photo_source_coverage(data)
        require(coverage["original_photo_attempts"] == len(data)
                and sum(coverage["recorded_observation_month_by_ROI_state"].values()) == len(data)
                and sum(coverage["latitude_band_by_ROI_state"].values()) == len(data),
                "Original photo taxonomy/season/latitude source QA lost observations")
        coverage_file = temp / "original_photo_source_coverage.json"
        coverage_file.write_text(json.dumps(coverage, indent=2, sort_keys=True) + "\n",
                                 encoding="utf-8")
        by_hue = Counter(row["morph"] for row in data)
        hue_path = temp / "group_original_photo_hue_counts.csv"
        with hue_path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f, lineterminator="\n")
            writer.writerow(("morph", "original_photo_attempts"))
            for hue in (*HUES, "mixed_uncertain"):
                writer.writerow((hue, by_hue[hue]))
        receipt = {
            "schema": SCHEMA,
            "original_acquisition_workflow_run": expected_run,
            "source_original_photo_manifest_sha256": expected_source_sha,
            "source_photo_opportunity_frame": expected_upstream,
            "original_source_warm_group": group,
            "original_work_units": expected_units,
            "completed_original_work_unit_source_SHA_audit": audit["completed_unit_audit"],
            "complete_expected_source_photo_attempts": expected_attempts,
            "distinct_genuinely_completed_photo_ids": len(data),
            "distinct_original_observation_ids": len(set(obs)),
            "four_hue_model_classified": classified,
            "acquisition_failed": audit["acquisition_failed"],
            "acquired_but_unclassifiable": audit["unclassifiable_acquired"],
            "hue_counts": {hue: by_hue[hue] for hue in (*HUES, "mixed_uncertain")},
            "original_source_taxon_IDs_four_hue": coverage["distinct_nominal_inat_taxon_ids_four_hue"],
            "original_source_date_latitude_coverage_separately_QA_verified": True,
            "selection_probability_retained_for_each_original_photo": True,
            "frozen_original_FCP_model_and_output_verified": True,
            "all_source_photo_and_original_terminal_SHA_verified": True,
            "native_or_stage_screen_preexclusion_enabled": False,
            "image_bytes_reacquired_in_sealing": False,
            "ROI_or_hue_reinference_in_sealing": False,
            "globally_probability_complete": False,
            "photo_IDs_globally_deduplicated_with_other_partitions": False,
            "saved_file_sha256": {
                data_path.name: sha(data_path),
                hue_path.name: sha(hue_path),
                coverage_file.name: sha(coverage_file),
            },
        }
        (temp / "original_photo_group_receipt.json").write_text(
            json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        os.replace(temp, output)
        print("MOIS_APPEND_ONLY_ORIGINAL_REAL_RGB_GROUP_SEALED",
              json.dumps({k: receipt[k] for k in (
                  "original_source_warm_group",
                  "distinct_genuinely_completed_photo_ids",
                  "four_hue_model_classified", "acquisition_failed",
                  "acquired_but_unclassifiable")}, sort_keys=True), flush=True)
        return receipt
    except BaseException:
        shutil.rmtree(temp, ignore_errors=True)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--workers", type=Path, required=True)
    parser.add_argument("--group", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    seal(args.source, args.workers, args.group, args.output)


if __name__ == "__main__":
    main()
