#!/usr/bin/env python3
"""Seal ONLY new SHA-pinned p0032 original-photo ROI groups after a real checkpoint.

Read one prior SHA-verified main-pixel checkpoint (photo IDs / hue projection),
the original source roster and ONLY newly completed 1024-photo group artifacts.
Do not download any original licensed photos, infer colour, scan the other five
waves' worker ZIPs or rewrite the previous ~0.9m original-photo Parquet.

A prior checkpoint snapshot is a strict, immutable ledger: every historical
p0032 group artifact must remain identical in the current provider manifest;
only newly sealed group artifacts may contribute a pixel delta. Future deltas
must advance the baseline to incorporate this snapshot, never replay the same
previously counted group.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import re

from scripts.acquisition.mois_completed_pixel_checkpoint import (
    WAVES, HUES, sha, collect_wave, require,
)
from scripts.analysis.mois_reassemble_original_new4096_measurements import METHOD

SOURCE_RUN, TARGET_PHOTOS, SOURCE_SHA = WAVES["p0032"]
PRIOR_CHECKPOINT_RUN = 39889436140
PRIOR_ATTEMPTS = 887600
PRIOR_FOUR_HUE = 400904
ARTIFACT_PREFIX = "mois-93m-p0032-RGB-group-"
DELTA_SCHEMA = "mois_p0032_new_sealed_original_RGB_photo_delta_from_signed_checkpoint_v1"
SHA_REGEX = re.compile(r"^sha256:[0-9a-f]{64}$")


def fingerprinted_diff(prior_receipt, current_artifacts, *,
                       baseline_groups=None, expected_run=SOURCE_RUN):
    """Compare source artifact ZIP IDs+digests, not mutable group names alone."""
    checkpoint = json.loads(Path(prior_receipt).read_text(encoding="utf-8"))
    require(checkpoint.get("schema") == "mois_sealed_actual_pixel_checkpoint_v1"
            and checkpoint.get("image_downloads_or_colour_reinference_performed") is False,
            "Previous real-pixel checkpoint isn't a source-validated immutable photo set")
    if baseline_groups is None:
        baseline_groups = checkpoint["source_artifact_snapshots"]["p0032"]
    require(baseline_groups.get("source_run_id") == expected_run,
            "Previous p0032 source worker run changed")
    old = baseline_groups.get("artifacts")
    live = current_artifacts.get("artifacts")
    require(current_artifacts.get("source_run_id") == expected_run
            and isinstance(old, list) and isinstance(live, list),
            "Snapshot absent or genuine original source run changed")

    def inventory(items):
        by_id, by_name, by_group = {}, {}, {}
        for asset in items:
            name = asset.get("name")
            require(isinstance(name, str) and name.startswith(ARTIFACT_PREFIX)
                    and name[len(ARTIFACT_PREFIX):].isdecimal()
                    and 0 <= int(name[len(ARTIFACT_PREFIX):]) < 145,
                    "Not an authenticated p0032 original RGB group artifact")
            group = int(name[len(ARTIFACT_PREFIX):])
            ident = asset.get("id")
            require(type(ident) is int and ident > 0
                    and isinstance(asset.get("digest"), str)
                    and SHA_REGEX.fullmatch(asset["digest"]) is not None
                    and type(asset.get("size_in_bytes")) is int
                    and asset["size_in_bytes"] > 0
                    and asset.get("expired") is not True,
                    "Invalid original-photo group artifact ID/digest/size or expired")
            # GitHub Release snapshots omit expired, live API provides it.
            require(ident not in by_id and name not in by_name and group not in by_group,
                    "Duplicated original-photo source worker artifact")
            canonical = {
                k: asset[k] for k in ("id", "name", "digest", "size_in_bytes")
            }
            by_id[ident] = canonical
            by_name[name] = canonical
            by_group[group] = canonical
        return by_id, by_name, by_group

    old_ids, old_names, old_groups = inventory(old)
    new_ids, new_names, new_groups = inventory(live)
    for ident, sealed in old_ids.items():
        require(ident in new_ids and new_ids[ident] == sealed,
                "Previously signed and classified original-photo worker changed")
    require(len(new_groups) >= len(old_groups), "p0032 previous true-photo workers disappeared")
    additions = [asset for group, asset in sorted(new_groups.items()) if group not in old_groups]
    require(all(item["id"] not in old_ids for item in additions),
            "Previously sealed RGB group was silently renamed or reused")
    require(len(additions) > 0, "No genuinely new original-photo group artifacts")
    return checkpoint, additions


def verify_new_pixel_delta(prior_parquet, prior_receipt, source, workers,
                           frozen_artifacts, output, *,
                           prior_run=PRIOR_CHECKPOINT_RUN,
                           prior_attempts=PRIOR_ATTEMPTS,
                           prior_hues=PRIOR_FOUR_HUE,
                           source_hash=SOURCE_SHA,
                           expected_photo_population=TARGET_PHOTOS):
    import pandas as pd

    output=Path(output)
    require(not output.exists(), "Refuse to overwrite an immutable original photo delta")
    prior_receipt, source = map(Path, (prior_receipt, source))
    snapshots=json.loads(Path(frozen_artifacts).read_text(encoding="utf-8"))
    previous, delta_assets = fingerprinted_diff(prior_receipt, snapshots)
    require((previous.get("combined_unique_photo_attempts"),
             previous.get("combined_four_hue_classified")) == (prior_attempts, prior_hues)
            and previous["output_sha256"][Path(prior_parquet).name] == sha(prior_parquet),
            "Main signed original-photo checkpoint counts or SHA changed")
    require(sha(source/"manifest.json") == source_hash,
            "Frozen complete p0032 original-photo source SHA changed")
    manifest=json.loads((source/"manifest.json").read_text(encoding="utf-8"))
    require(manifest["schema"] == "mois_full_source_observation_sampling_cohort_v1"
            and manifest["new_photo_attempts_planned"] == expected_photo_population
            and manifest["full_frame_photo_ids"] == 93209351
            and manifest["worker_ids"] == list(range(manifest["work_units"]))
            and len(manifest["group_ids"]) == 145,
            "Wrong signed original 93m image opportunity/source wave")
    prior=pd.read_parquet(prior_parquet, columns=["photo_id", "morph", "checkpoint_wave", "observation_id"])
    photo_ids=prior.photo_id.astype(str)
    require(len(prior) == prior_attempts and photo_ids.is_unique
            and photo_ids.str.fullmatch(r"\d+").all(),
            "Old classified original image IDs duplicated or missing")
    prior_set=set(photo_ids)
    del prior
    data, work_audit=collect_wave(source,workers,METHOD)
    expected_groups=[int(a["name"][len(ARTIFACT_PREFIX):]) for a in delta_assets]
    new_unit_indices=sorted({
        u for group in expected_groups
        for u in range(4*group,min(4*(group+1),manifest["work_units"]))
    })
    require(work_audit["completed_units"] == new_unit_indices
            and work_audit["planned_photo_attempts"] == expected_photo_population
            and work_audit["attempted"] == len(data)>0,
            "New actual FCP original photo unit set differs from artifact snapshot")
    newids={str(row["photo_id"]) for row in data}
    require(len(newids)==len(data)
            and not prior_set.intersection(newids)
            and all(p.isdecimal() for p in newids),
            "True original photo ID was processed twice across checkpoint and delta")
    require(len(data) <= 1024*len(expected_groups), "Original group contains more than 4x256 photo attempts")
    flags=Counter()
    probability_by_observation={}
    for row in data:
        known=row["morph"] in HUES
        failed=row["measurement_status"]=="image_acquisition_failed"
        require(row["morph"] in (*HUES,"mixed_uncertain")
                and known == (row["measurement_status"]=="classified_four_state_morph")
                and (not known or not failed),
                "Technical unknowns were falsely classified as flower petal colours")
        try:
            p=float(row["selection_probability"])
        except (TypeError,ValueError):
            raise ValueError("Original observation inclusion probability absent")
        require(0 < p <= 1,"Original observation inclusion probability outside valid range")
        obs=row["observation_id"]
        require(bool(obs),"Missing original observation ID in genuine photos")
        if obs in probability_by_observation:
            require(probability_by_observation[obs] == p,
                    "Within-observation original-photo sampling probability differs")
        probability_by_observation[obs]=p
        flags["attempted"]+=1
        flags["classified"]+=known
        flags["acquisition_failed"]+=failed
        flags["acquired_unclassifiable"]+=not known and not failed
        flags[row["morph"]]+=1
    require(flags["attempted"]==work_audit["attempted"]
            and flags["classified"]==work_audit["classifiable"]
            and flags["acquisition_failed"]==work_audit["acquisition_failed"]
            and flags["acquired_unclassifiable"]==work_audit["unclassifiable_acquired"],
            "Original source worker classification/failure receipt count changed")
    output.mkdir(parents=True, exist_ok=False)
    datafile=output/"new_p0032_original_photo_ROIs.parquet"
    frame=pd.DataFrame(data).sort_values("photo_id",key=lambda s:s.map(int)).reset_index(drop=True)
    frame.to_parquet(datafile,index=False,compression="zstd")
    countsfile=output/"new_p0032_original_photo_hue_counts.csv"
    with countsfile.open("w",newline="",encoding="utf-8") as f:
        writer=csv.writer(f,lineterminator="\n")
        writer.writerow(("morph","true_original_photo_attempts"))
        for key in (*HUES,"mixed_uncertain"):
            writer.writerow((key,flags[key]))
    new_snapshot={
        "source_run_id":SOURCE_RUN,
        "artifacts":delta_assets,
    }
    receipt={
        "schema":DELTA_SCHEMA,
        "prior_immutable_main_checkpoint_run":prior_run,
        "prior_immutable_main_photo_attempts":prior_attempts,
        "prior_immutable_main_four_hue_classified":prior_hues,
        "prior_checkpoint_sha256":sha(prior_receipt),
        "prior_checkpoint_parquet_sha256":sha(prior_parquet),
        "p0032_original_photo_source_SHA256":source_hash,
        "p0032_source_planned_photo_opportunities":expected_photo_population,
        "new_sealed_p0032_groups":expected_groups,
        "new_immutable_group_artifacts":new_snapshot,
        "source_worker_original_photo_Sha_audit":work_audit["completed_unit_audit"],
        "completed_true_original_photo_attempts":flags["attempted"],
        "new_true_original_four_hue_classified":flags["classified"],
        "new_true_original_acquisition_failed":flags["acquisition_failed"],
        "new_true_original_acquired_but_unclassifiable":flags["acquired_unclassifiable"],
        "new_original_photo_hue_counts":{h:flags[h] for h in (*HUES,"mixed_uncertain")},
        "combined_disjoint_main_plus_delta_photo_attempts":prior_attempts+flags["attempted"],
        "combined_disjoint_main_plus_delta_four_hue_classified":prior_hues+flags["classified"],
        "new_unique_original_observations":len(probability_by_observation),
        "source_arm":manifest.get("source_arm"),
        "unchanged_original_roi_model_and_image_terminal_byte_SHA_verified":True,
        "previous_photo_ids_and_new_pixel_ids_exactly_disjoint":True,
        "previous_photo_dataset_rewritten":False,
        "other_original_wave_worker_archives_redownloaded":0,
        "new_original_photos_downloaded_from_S3":0,
        "new_roi_or_hue_inference_performed":False,
        "previous_checkpoint_not_explicitly_includes_separate_tier06_source_arms":True,
        "sampling_partial_and_not_representative_global_flora":True,
        "immutable_output_sha256":{
            datafile.name:sha(datafile),
            countsfile.name:sha(countsfile),
        },
    }
    (output/"p0032_signed_incremental_original_photo_receipt.json").write_text(
        json.dumps(receipt,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("MOIS_P0032_REAL_ORIGINAL_PIXEL_DELTA_WITHOUT_OLD_WAVE_REDOWNLOAD",
          json.dumps({k:receipt[k] for k in (
              "new_sealed_p0032_groups",
              "completed_true_original_photo_attempts",
              "new_true_original_four_hue_classified",
              "combined_disjoint_main_plus_delta_photo_attempts",
              "combined_disjoint_main_plus_delta_four_hue_classified")},sort_keys=True),flush=True)
    return receipt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("prior-parquet","prior-receipt","source","workers","frozen-artifacts","output"):
        p.add_argument("--"+key,type=Path,required=True)
    args=p.parse_args()
    verify_new_pixel_delta(args.prior_parquet,args.prior_receipt,args.source,
                           args.workers,args.frozen_artifacts,args.output)


if __name__=="__main__":
    main()
