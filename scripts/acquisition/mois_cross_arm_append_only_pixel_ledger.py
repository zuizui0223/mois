#!/usr/bin/env python3
"""Append-only cross-arm index for signed real MOIS original-image colour results.

Read photo/observation IDs from each already-published Parquet exactly once.
SHA-pin prior real-pixel baselines and newly sealed worker-group evidence. Fail
closed on ANY duplicate photo or observation across arms. Do not download raw
images, infer colour, rewrite or concatenate historical Parquets or compute an
unweighted global flower-white frequency from mixed source designs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re

HUES = ("white", "yellow_orange", "red_pink", "blue_purple")
INDEX_SCHEMA = "mois_original_pixel_immutable_cross_arm_group_registry_v1"
PINNED_MAIN_RUN = 39884380759
MAIN_COUNT = 885552
MAIN_HUES = 400102
FIRST_TIER_COUNT = 4138
FIRST_TIER_HUES = 1573
BROAD_RUN = 39880330611
BROAD_SOURCE_SHA = "3e28ced8de85709ea4854cd3291aa00f56e98306b60a4c273c5c7bb47b68c858"
BROAD_PLANNED = 4468


def require(value, reason):
    if not value:
        raise ValueError(reason)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def inspect_photo_ids(path, columns=("photo_id", "observation_id")):
    import pandas as pd

    table = pd.read_parquet(path, columns=list(columns))
    require(all(col in table for col in columns), "Signed photo identity columns absent")
    photo_ids = table.photo_id.astype(str)
    require(photo_ids.str.fullmatch(r"[0-9]+").all() and photo_ids.is_unique,
            "Original photo IDs duplicated/non-numeric inside signed Parquet")
    obs_ids = table.observation_id.fillna("").astype(str).str.strip()
    # Baseline: some historical photos have not-yet-resolved source observations.
    require(not obs_ids.str.lower().isin(("nan", "none", "null")).any(),
            "Explicitly invalid photo observation IDs are not accepted")
    return table, set(photo_ids), {o for o in obs_ids if o}


def register(main_parquet, main_receipt, old_parquet, old_receipt,
             source_manifest, completed_groups, output,
             *, expected_main=(MAIN_COUNT, MAIN_HUES),
             expected_old=(FIRST_TIER_COUNT, FIRST_TIER_HUES),
             expected_source_sha=BROAD_SOURCE_SHA,
             expected_broad_total=BROAD_PLANNED):
    import pandas as pd

    output = Path(output)
    require(not output.exists(), "Immutable cross-arm original image index already exists")
    main = load_json(main_receipt)
    old = load_json(old_receipt)
    source = load_json(source_manifest)
    require(main.get("schema") == "mois_sealed_actual_pixel_checkpoint_v1"
            and (main.get("combined_unique_photo_attempts"),
                 main.get("combined_four_hue_classified")) == tuple(expected_main)
            and main["output_sha256"][Path(main_parquet).name] == sha(main_parquet),
            "Main genuine original photo baseline receipt or bytes do not match")
    require(old.get("schema") ==
            "mois_global93m_observation_baseline_real_RGB_photo_results_v1"
            and (old.get("attempted"), old.get("classifiable")) == tuple(expected_old)
            and old.get("measured_photo_parquet_sha256") == sha(old_parquet),
            "Old tier06 original photo pixel receipt/bytes do not match")
    require(sha(source_manifest) == expected_source_sha
            and source.get("schema") ==
            "mois_full_source_observation_sampling_cohort_v1"
            and source.get("source_tier") == 6
            and source.get("source_group") == "all16"
            and source.get("new_photo_attempts_planned") == expected_broad_total
            and source.get("full_frame_photo_ids") == 93209351
            and source.get("compute_bucket_count") == 512
            and source.get("admitted_bucket_count") == 4
            and source.get("colour_measured") is False
            and source.get("not_probability_complete_global_source") is True,
            "Broad first4/512 source is not original precolour SHA-locked input")
    require(len(completed_groups) > 0, "No newly completed original image worker Parquets")

    base, photo_ids, observations = inspect_photo_ids(main_parquet)
    require(len(base) == expected_main[0], "Main original image photo count changed")
    del base
    older, old_ids, old_obs = inspect_photo_ids(old_parquet)
    require(len(older) == expected_old[0], "Previous tier06 original photo count changed")
    del older
    require(not photo_ids.intersection(old_ids)
            and not observations.intersection(old_obs),
            "Main and preceding tier06 photo/observation IDs unexpectedly overlap")
    photo_ids.update(old_ids)
    observations.update(old_obs)

    group_results = []
    used_units = set()
    for group_dir in completed_groups:
        group_dir = Path(group_dir)
        receipt_path = group_dir / "original_photo_group_receipt.json"
        parquet_path = group_dir / "original_photo_ROIs.parquet"
        colors_path = group_dir / "group_original_photo_hue_counts.csv"
        receipt = load_json(receipt_path)
        group = receipt.get("original_source_warm_group")
        require(type(group) is int and group in source["group_ids"]
                and receipt.get("schema") ==
                    "mois_tier06_all16_append_only_original_RGB_group_v1"
                and receipt.get("original_acquisition_workflow_run") == BROAD_RUN
                and receipt.get("source_original_photo_manifest_sha256") ==
                    expected_source_sha
                and receipt.get("complete_expected_source_photo_attempts") ==
                    expected_broad_total
                and receipt.get("all_source_photo_and_original_terminal_SHA_verified") is True
                and receipt.get("frozen_original_FCP_model_and_output_verified") is True
                and receipt.get("image_bytes_reacquired_in_sealing") is False
                and receipt.get("ROI_or_hue_reinference_in_sealing") is False
                and receipt.get("globally_probability_complete") is False,
                "Unknown or altered original pixel worker-group evidence")
        expected_units = list(range(group*4, min((group+1)*4, source["work_units"])))
        units = receipt.get("original_work_units")
        require(units == expected_units and not used_units.intersection(units),
                "Original photo work units missing or duplicated across groups")
        used_units.update(units)
        fingerprints = receipt.get("saved_file_sha256", {})
        require(sha(parquet_path) == fingerprints.get(parquet_path.name)
                and sha(colors_path) == fingerprints.get(colors_path.name),
                "Immutable original group photo or hue count bytes changed")
        cols = ["photo_id", "observation_id", "morph", "measurement_status",
                "source_arm", "selection_probability"]
        table, ids, obs = inspect_photo_ids(parquet_path, cols)
        require(not photo_ids.intersection(ids) and not observations.intersection(obs),
                "Original photo or original observation reused in source arms")
        require(len(table) == receipt.get("distinct_genuinely_completed_photo_ids")
                and len(obs) == receipt.get("distinct_original_observation_ids")
                and table.observation_id.notna().all(),
                "Original image photo/observation counts mismatch")
        known = table.morph.isin(HUES)
        require(table.morph.isin((*HUES, "mixed_uncertain")).all()
                and known.eq(table.measurement_status.eq("classified_four_state_morph")).all()
                and int(known.sum()) == receipt.get("four_hue_model_classified")
                and (table.source_arm == source["source_arm"]).all(),
                "Original photo classification or source-arm changed")
        numeric = pd.to_numeric(table.selection_probability, errors="coerce")
        require(numeric.notna().all()
                and numeric.eq(source["selection_probability"]).all(),
                "Original observation sampling probability changed")
        requires = (
            int(known.sum())
            + int(receipt.get("acquisition_failed", -1))
            + int(receipt.get("acquired_but_unclassifiable", -1))
        )
        require(requires == len(table)
                and receipt.get("hue_counts", {}).get("mixed_uncertain")
                   == int((~known).sum()),
                "Original photo acquire/unknown/disposition counts disagree")
        photo_ids.update(ids)
        observations.update(obs)
        group_results.append({
            "group":group,
            "source_shards":expected_units,
            "photo_attempts":len(table),
            "four_hue_classified":int(known.sum()),
            "unclassifiable":int((~known).sum()),
            "acquisition_failed":int(receipt["acquisition_failed"]),
            "distinct_original_observations":len(obs),
            "original_worker_receipt_sha256":sha(receipt_path),
            "original_photo_parquet_sha256":sha(parquet_path),
            "original_hue_count_csv_sha256":sha(colors_path),
        })
    group_results.sort(key=lambda x:x["group"])
    require(len({r["group"] for r in group_results}) == len(group_results),
            "Original photo group processed more than once")
    completed_count = sum(r["photo_attempts"] for r in group_results)
    new_classified = sum(r["four_hue_classified"] for r in group_results)
    require(completed_count <= expected_broad_total
            and len(used_units) <= source["work_units"], "Original source photo work overstated")
    # This is a union of distinct photo IDs, not an unbiased worldwide
    # count of blooming flowers. No coloured-petal truth label is synthesized.
    audit = {
        "schema":INDEX_SCHEMA,
        "source_registry_photo_opportunities_not_measured":93209351,
        "original_broad_source_manifest_sha256":sha(source_manifest),
        "main_checkpoint_run":PINNED_MAIN_RUN,
        "main_checkpoint_photo_attempts":expected_main[0],
        "main_checkpoint_four_hue_classified":expected_main[1],
        "main_checkpoint_parquet_sha256":sha(main_parquet),
        "previous_tier06_source_photo_attempts":expected_old[0],
        "previous_tier06_four_hue_classified":expected_old[1],
        "previous_tier06_parquet_sha256":sha(old_parquet),
        "broad_source_original_photo_attempts_planned":expected_broad_total,
        "broad_sealed_groups":group_results,
        "broad_completed_photo_attempts":completed_count,
        "broad_completed_four_hue_classified":new_classified,
        "broad_pending_photo_attempts":expected_broad_total-completed_count,
        "original_photo_attempts_distinct_across_all_verified_arms":len(photo_ids),
        "four_hue_classifications_across_all_verified_arms":
            expected_main[1]+expected_old[1]+new_classified,
        "distinct_original_observation_ids_across_all_arms":len(observations),
        "original_photo_or_observation_ID_cross_arm_overlap":0,
        "all_original_photo_per_partition_SHA_receipts_verified":True,
        "original_FCP_roi_v4_model_attestations_reused_without_new_inference":True,
        "original_photo_images_fetched_by_this_registry":0,
        "old_pixel_Parquets_rewritten":False,
        "source_arms_have_distinct_observation_selection_probabilities":True,
        "global_flora_white_flower_rate_estimable_by_naive_pooled_rows":False,
        "independent_blooming_petals_confirmed":False,
        "photo_hue_labels_not_yet_species_level_or_wild_native_colour":True,
        "pending_photo_attempts_not_assumed_to_be_any_colour":True,
    }
    require(audit["original_photo_attempts_distinct_across_all_verified_arms"]
            == expected_main[0]+expected_old[0]+completed_count,
            "Cross-arm unique original-photo row accounting changed")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n",
                      encoding="utf-8")
    print("MOIS_SHA_VERIFIED_APPEND_ONLY_ORIGINAL_PHOTO_CROSS_ARM_LEDGER",
          json.dumps({k:audit[k] for k in (
              "broad_completed_photo_attempts",
              "broad_completed_four_hue_classified",
              "broad_pending_photo_attempts",
              "original_photo_attempts_distinct_across_all_verified_arms",
              "four_hue_classifications_across_all_verified_arms")},sort_keys=True),
          flush=True)
    return audit


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for label in ("main-parquet", "main-receipt", "old-parquet", "old-receipt",
                  "source-manifest", "output"):
        p.add_argument("--"+label, type=Path, required=True)
    p.add_argument("--group", type=Path, action="append", required=True)
    a=p.parse_args()
    register(a.main_parquet,a.main_receipt,a.old_parquet,a.old_receipt,
             a.source_manifest,a.group,a.output)


if __name__=="__main__":
    main()
