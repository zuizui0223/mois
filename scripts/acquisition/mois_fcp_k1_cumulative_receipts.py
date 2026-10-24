#!/usr/bin/env python3
"""Aggregate already SHA-sealed K1 original-photo worker receipts without photos.

This data-only registry is a partial-compute progress record, NOT a globally
deduplicated cross-arm photo total or independently botanical flower census.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile

SCHEMA = "mois_fcp_K1_original_photo_cumulative_signed_group_ledger_v1"
SOURCE_SCHEMA = "mois_fcp_K1_93m_balanced_true_original_ROI_append_group_v1"
SOURCE_RUN = 40066261878
SOURCE_SHA = "4da66dc0d34f440269c2202f95b0f39ad369531fd625676e8ec69df4ae939f1c"
PARENT_SHA = "c0ba3a316a336816065590540e75034bd8c2b9b7362034163673702103c44977"
SOURCE_ARM = "mois_global93m_tier06_taxon_month_lat_probability_onephoto_K1_v1"
SOURCE_PHOTOS = 93_209_351
SOURCE_SELECTED_OBSERVATIONS = 80_035
REUSED_OLD_ORIGINALS = 3_280
NEW_TARGET = 76_753
NUM_GROUPS = 75
UNIT_SIZE = 256
GROUP_UNITS = 4
HUES = ("white", "yellow_orange", "red_pink", "blue_purple")
CSV_COLUMNS = (
    "original_worker_group", "new_original_photo_attempts",
    "four_hue_classified", "image_acquisition_failed",
    "acquired_but_unclassifiable", "original_wall_seconds",
    "original_group_receipt_SHA256",
)


def require(ok, msg):
    if not ok:
        raise ValueError(msg)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def group_expected(group):
    require(type(group) is int and 0 <= group < NUM_GROUPS,
            "Group index outside signed K1 original queue")
    return min(GROUP_UNITS * UNIT_SIZE, NEW_TARGET - group * GROUP_UNITS * UNIT_SIZE)


def audit_receipt(path):
    path = Path(path)
    d = json.loads(path.read_text(encoding="utf-8"))
    group = d.get("original_worker_group")
    require(type(group) is int and 0 <= group < NUM_GROUPS,
            "Invalid original K1 completed group")
    require(
        d.get("schema") == SOURCE_SCHEMA
        and d.get("frozen_original_source_run") == SOURCE_RUN
        and d.get("full_93m_photo_target_catalogue") == SOURCE_PHOTOS
        and d.get("original_source_manifest_SHA256") == SOURCE_SHA
        and d.get("original_precolour_probability_source_SHA256") == PARENT_SHA
        and d.get("source_arm") == SOURCE_ARM
        and d.get("previously_measured_old_original_photo_IDs_reused_outside_this_group") == REUSED_OLD_ORIGINALS
        and d.get("original_new_selected_photo_targets_all_groups") == NEW_TARGET
        and d.get("original_completed_source_units") == list(
            range(group * GROUP_UNITS,
                  min((group + 1) * GROUP_UNITS, math.ceil(NEW_TARGET / UNIT_SIZE))))
        and d.get("frozen_original_FCP_roi_v4_source_method_attested") is True
        and d.get("original_photo_identity_and_image_terminal_SHA_audited") is True
        and d.get("source_inclusion_probability_retained_at_each_photo") is True
        and d.get("original_warm_FCP_model_inits") == 1
        and d.get("image_downloads_OR_ROI_reinferences_in_this_sealer") == 0
        and d.get("old_historical_original_photo_parquets_rewritten") is False
        and d.get("technical_unknown_ROIs_not_turned_into_other_hues") is True
        and d.get("true_botanical_plant_flower_or_native_status_validated") is False
        and d.get("probability_weighted_global_flora_white_rate_estimated") is False,
        "K1 original image/source SHA, provenance or unknown contract invalid")
    photos = d.get("actual_new_original_photo_attempts")
    classified = d.get("actual_automatic_four_hue_ROI_classifications")
    failure = d.get("original_image_acquisition_failed")
    unknown = d.get("original_acquired_but_four_hue_unclassifiable")
    obs = d.get("original_distinct_observations")
    require(
        all(type(x) is int for x in (photos, classified, failure, unknown, obs))
        and photos == group_expected(group) == obs
        and classified >= 0 and failure >= 0 and unknown >= 0
        and classified + failure + unknown == photos,
        "Completed K1 original photo or technical nonresponse counts altered")
    hues = d.get("original_hue_counts")
    require(
        isinstance(hues, dict) and set(hues) == set(HUES) | {"mixed_uncertain"}
        and all(type(x) is int and x >= 0 for x in hues.values())
        and sum(hues[x] for x in HUES) == classified
        and hues["mixed_uncertain"] == failure + unknown,
        "K1 original four-hue and mixed-unknown outcomes do not reconcile")
    runtime = d.get("original_full_roi_wall_seconds")
    eval_time = d.get("original_full_roi_evaluator_seconds")
    prefetch = d.get("original_photo_prefetch")
    require(
        isinstance(runtime, (float, int)) and math.isfinite(runtime) and runtime > 0
        and isinstance(eval_time, (float, int)) and math.isfinite(eval_time)
        and 0 <= eval_time <= runtime
        and isinstance(prefetch, dict)
        and prefetch.get("consumed") == photos
        and prefetch.get("download_calls") == photos
        and prefetch.get("permanent_failures") == failure
        and type(d.get("distinct_conditional_joint_photo_probabilities")) is int
        and d["distinct_conditional_joint_photo_probabilities"] >= 1,
        "K1 original warm inference/prefetch runtime not independently attested")
    terminal_audit = d.get("original_photo_terminal_source_receipts_SHA256")
    require(
        isinstance(terminal_audit, list)
        and len(terminal_audit) == len(d["original_completed_source_units"])
        and [r.get("unit") for r in terminal_audit] == d["original_completed_source_units"]
        and all(
            isinstance(r, dict) and all(
                isinstance(r.get(key), str) and len(r[key]) == 64
                and all(c in "0123456789abcdef" for c in r[key])
                for key in ("receipt_sha256", "photos_sha256", "terminal_sha256"))
            for r in terminal_audit),
        "Original terminal/source work-unit SHA attestations missing")
    outputs = d.get("artifact_SHA256")
    require(
        isinstance(outputs, dict)
        and set(("original_photo_ROIs.parquet", "original_photo_source_coverage.json",
                 "original_photo_hue_counts.csv")) <= set(outputs)
        and all(isinstance(x, str) and len(x) == 64
                and all(c in "0123456789abcdef" for c in x) for x in outputs.values()),
        "Original group append artifact SHA inventory incomplete")
    return {
        "original_worker_group": group,
        "new_original_photo_attempts": photos,
        "four_hue_classified": classified,
        "image_acquisition_failed": failure,
        "acquired_but_unclassifiable": unknown,
        "original_wall_seconds": float(runtime),
        "original_evaluator_seconds": float(eval_time),
        "original_hue_counts": hues,
        "original_group_receipt_SHA256": sha256(path),
        "distinct_conditional_joint_photo_probabilities":
            d["distinct_conditional_joint_photo_probabilities"],
    }


def aggregate(receipts, output):
    output = Path(output)
    require(bool(receipts), "No genuine completed original K1 workers to aggregate")
    require(not output.exists(), "Refuse to rewrite immutable cumulative original-photo ledger")
    rows = [audit_receipt(Path(p)) for p in receipts]
    groups = [r["original_worker_group"] for r in rows]
    require(len(set(groups)) == len(groups), "Duplicate signed K1 original ROI group")
    rows.sort(key=lambda r: r["original_worker_group"])
    group_ids = [r["original_worker_group"] for r in rows]
    counts = Counter()
    hue_counts = Counter()
    for r in rows:
        for key in ("new_original_photo_attempts", "four_hue_classified",
                    "image_acquisition_failed", "acquired_but_unclassifiable"):
            counts[key] += r[key]
        hue_counts.update(r["original_hue_counts"])
    assert counts["new_original_photo_attempts"] == (
        counts["four_hue_classified"] + counts["image_acquisition_failed"]
        + counts["acquired_but_unclassifiable"])
    fingerprint = hashlib.sha256("\n".join(
        f'{r["original_worker_group"]}:{r["original_group_receipt_SHA256"]}'
        for r in rows).encode()).hexdigest()
    report = {
        "schema": SCHEMA,
        "source_run": SOURCE_RUN,
        "full_93m_original_candidate_photo_IDs_METADATA_ONLY": SOURCE_PHOTOS,
        "selected_independent_original_observations": SOURCE_SELECTED_OBSERVATIONS,
        "previously_measured_K1_photo_IDs_reused_elsewhere_NOT_NEW": REUSED_OLD_ORIGINALS,
        "new_original_photo_target_total": NEW_TARGET,
        "total_scheduled_original_worker_groups": NUM_GROUPS,
        "completed_and_SHA_sealed_original_worker_groups": len(rows),
        "completed_original_worker_group_IDs": group_ids,
        "not_yet_sealed_worker_groups":
            [x for x in range(NUM_GROUPS) if x not in set(groups)],
        "completed_real_new_original_photo_attempts": counts["new_original_photo_attempts"],
        "completed_automatic_four_hue_classifications": counts["four_hue_classified"],
        "completed_original_acquisition_failures": counts["image_acquisition_failed"],
        "completed_acquired_but_automatic_four_hue_unclassifiable":
            counts["acquired_but_unclassifiable"],
        "completed_hue_photo_counts_CONDITIONAL_ON_AUTOMATIC_CLASSIFICATION":
            {h: hue_counts[h] for h in HUES},
        "completed_mixed_unknown_photo_count": hue_counts["mixed_uncertain"],
        "sum_independent_original_worker_wall_seconds":
            sum(r["original_wall_seconds"] for r in rows),
        "sum_original_evaluator_seconds":
            sum(r["original_evaluator_seconds"] for r in rows),
        "source_original_group_receipt_aggregate_SHA256": fingerprint,
        "original_group_receipt_metadata_SHA_checked": True,
        "source_precolour_probability_photo_sampling_preserved": True,
        "within_K1_original_source_group_units_disjoint": True,
        "new_original_images_downloaded_or_colour_inferred_in_THIS_ledger": 0,
        "this_is_NOT_a_crossarm_photo_ID_deduplicated_global_total": True,
        "global_wild_flora_or_blooming_petal_status_not_independently_verified": True,
        "per_group": rows,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    parent = Path(tempfile.mkdtemp(prefix="mois-k1-cumulative-", dir=output.parent))
    try:
        (parent / "k1_cumulative_original_ROI_receipt.json").write_text(
            json.dumps(report, sort_keys=True, indent=2) + "\n", encoding="utf-8")
        with (parent / "k1_finished_original_group_counts.csv").open(
                "w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=CSV_COLUMNS, lineterminator="\n")
            writer.writeheader()
            for r in rows:
                writer.writerow({key: r[key] for key in CSV_COLUMNS})
        os.replace(parent, output)
    except BaseException:
        import shutil
        shutil.rmtree(parent, ignore_errors=True)
        raise
    print("MOIS_K1_TRUE_ORIGINAL_PIXEL_CUMULATIVE_PROGRESS", json.dumps({
        "verified_groups": len(rows),
        "new_actual_photos": counts["new_original_photo_attempts"],
        "classified_four_hue": counts["four_hue_classified"],
        "image_failed": counts["image_acquisition_failed"],
        "acquired_unclassifiable": counts["acquired_but_unclassifiable"],
        "receipt_set_sha256": fingerprint,
        "no_new_original_image_download_or_inference": True,
    }, sort_keys=True), flush=True)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receipt", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    aggregate(args.receipt, args.output)


if __name__ == "__main__":
    main()
