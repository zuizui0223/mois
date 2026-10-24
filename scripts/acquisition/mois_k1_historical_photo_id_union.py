#!/usr/bin/env python3
"""Exact-photo-ID union of SHA-authenticated historical MOIS and completed K1 pixels.

NO source photos downloaded, YOLO/EfficientSAM inference, or unweighted global
native-flower prevalence estimates. Reuse the same historical 929138-photo
source auditor that K1 used when freezing its new-only original image queue.
Read completed K1 full-resolution terminal Parquets, not unmeasured source IDs.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import math
import os
import tempfile

from scripts.acquisition.mois_completed_pixel_checkpoint import HUES, require, sha
from scripts.acquisition.mois_fcp_k1_cumulative_receipts import (
    SCHEMA as CUMULATIVE_SCHEMA,
    SOURCE_RUN, SOURCE_ARM, SOURCE_SHA, NEW_TARGET,
    NUM_GROUPS, REUSED_OLD_ORIGINALS, audit_receipt,
)
from scripts.acquisition.mois_93m_balanced_prior_ROI_reuse_bridge import (
    PRIOR_MAIN, PRIOR_OLD, PRIOR_R1, PRIOR_R2, prior_original_photo_index,
    prior_observation_unknown,
)

SCHEMA = "mois_genuine_original_929138_historical_and_K1_photo_ID_unique_union_v1"
OLD_COUNT = 929138
PRIOR_COUNTS = {"main":PRIOR_MAIN, "prior_group00":PRIOR_OLD,
                "prior_R1":PRIOR_R1, "prior_R2":PRIOR_R2}
SOURCE_SHAS = ("original_photo_ROIs.parquet", "original_photo_hue_counts.csv",
               "original_photo_source_coverage.json")


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def merge_verified(previous, original_obs, historical_audits, overlap, cumulative,
                   group_dirs, output, *, expected_prior=OLD_COUNT):
    """Read only genuine signed saved ROI rows; reject any wrong photo identity."""
    import pandas as pd

    output = Path(output)
    require(not output.exists(), "Refuse to overwrite source-certified global photo union")
    require(len(previous) == expected_prior == overlap["unique_old_photo_IDs"]
            and overlap["signed_arm_photo_rows"] == sum(x[0] for x in PRIOR_COUNTS.values())
            and not overlap.get("old_cross_arm_same_photo_evidence_conflicts", 0),
            "Historical global original-photo frame was not the signed unique photo census")
    require(cumulative.get("schema") == CUMULATIVE_SCHEMA
            and cumulative.get("source_run") == SOURCE_RUN
            and cumulative.get("new_original_photo_target_total") == NEW_TARGET
            and cumulative.get("total_scheduled_original_worker_groups") == NUM_GROUPS
            and cumulative.get("previously_measured_K1_photo_IDs_reused_elsewhere_NOT_NEW")
                == REUSED_OLD_ORIGINALS
            and cumulative.get("new_original_images_downloaded_or_colour_inferred_in_THIS_ledger") == 0
            and cumulative.get("this_is_NOT_a_crossarm_photo_ID_deduplicated_global_total") is True,
            "Source K1 signed previous-only cumulative ledger changed or is not genuine")
    attest = cumulative.get("per_group")
    require(isinstance(attest, list) and attest
            and cumulative.get("completed_and_SHA_sealed_original_worker_groups") == len(attest)
            and len(group_dirs) == len(attest),
            "Not all actually finished original K1 image groups are available")
    by_group = {}
    for a in attest:
        group = a.get("original_worker_group")
        require(type(group) is int and 0 <= group < NUM_GROUPS and group not in by_group,
                "Duplicate or malformed K1 cumulative original group")
        by_group[group] = a
    require(sorted(by_group) == cumulative["completed_original_worker_group_IDs"],
            "Signed K1 per-group roster changed")
    fingerprint=hashlib.sha256("\n".join(
        f'{g}:{by_group[g]["original_group_receipt_SHA256"]}'
        for g in sorted(by_group)).encode()).hexdigest()
    require(fingerprint==cumulative["source_original_group_receipt_aggregate_SHA256"],
            "Cumulative K1 source receipt set SHA does not match published selection")

    seen=set()
    obs_seen=set()
    sum_state=Counter()
    each=[]
    overlap_observation_ids=set()
    for folder in group_dirs:
        folder=Path(folder)
        recpath=folder/"original_K1_group_receipt.json"
        rec=audit_receipt(recpath)
        group=rec["original_worker_group"]
        require(group in by_group and group not in seen
                and rec["original_group_receipt_SHA256"]==by_group[group]["original_group_receipt_SHA256"]
                and rec["new_original_photo_attempts"]==by_group[group]["new_original_photo_attempts"]
                and rec["four_hue_classified"]==by_group[group]["four_hue_classified"],
                "Completed K1 original-photo group differs from frozen cumulative source")
        seen.add(group)
        original=read_json(recpath)
        for name in SOURCE_SHAS:
            f=folder/name
            require(f.is_file() and sha(f)==original["artifact_SHA256"][name],
                    "K1 saved original ROI/coverage/photo-hue artifact SHA changed: "+name)
        cover=read_json(folder/"original_photo_source_coverage.json")
        require(cover.get("schema")==
                "mois_fcp_K1_signed_original_taxon_recorded_month_latitude_coverage_v1"
                and cover.get("sampled_original_photo_inclusion_probability_varies_per_photo") is True
                and cover.get("one_selected_photo_per_original_observation") is True
                and cover.get("individual_photo_open_flowers_independently_verified") is False,
                "Completed original-photo coverage incorrectly claims botanical flower verification")
        df=pd.read_parquet(folder/"original_photo_ROIs.parquet",
            columns=["photo_id","observation_id","morph","measurement_status",
                     "image_sha256","source_arm","selection_probability"])
        photo=df.photo_id.fillna("").astype(str)
        observation=df.observation_id.fillna("").astype(str)
        morph=df.morph.fillna("").astype(str)
        status=df.measurement_status.fillna("").astype(str)
        file_sha=df.image_sha256.fillna("").astype(str)
        classified=morph.isin(HUES)
        failed=status.eq("image_acquisition_failed")
        prob=pd.to_numeric(df.selection_probability,errors="coerce")
        require(
            len(df)==rec["new_original_photo_attempts"]
            and photo.str.fullmatch(r"[0-9]+").all() and photo.is_unique
            and observation.str.strip().ne("").all() and observation.is_unique
            and set(photo).isdisjoint(previous) and set(photo).isdisjoint(sum_state.get("_photo_ids",set()))
            and set(observation).isdisjoint(obs_seen)
            and morph.isin((*HUES,"mixed_uncertain")).all()
            and status.eq("classified_four_state_morph").eq(classified).all()
            and ((file_sha.str.fullmatch(r"[0-9a-f]{64}")) | failed).all()
            and df.source_arm.eq(SOURCE_ARM).all()
            and prob.notna().all() and prob.gt(0).all() and prob.le(1).all()
            and int(classified.sum())==rec["four_hue_classified"]
            and int(failed.sum())==rec["image_acquisition_failed"]
            and int((~classified & ~failed).sum())==rec["acquired_but_unclassifiable"],
            "K1 actual photo IDs repeat, original terminal changed, or original photo joint probability invalid")
        with (folder/"original_photo_hue_counts.csv").open(
                newline="", encoding="utf-8") as handle:
            counted={r["morph"]:int(r["original_photo_attempts"])
                     for r in csv.DictReader(handle)}
        require(counted==original["original_hue_counts"]
                and all(int((morph==h).sum())==counted[h]
                        for h in (*HUES,"mixed_uncertain")),
                "K1 genuine original photo hue counts differ from signed data")
        # Observation identity is retained for dependence diagnostics; photo ID
        # is the actual deduplication key. Distinct photos from one observation
        # do NOT imply independent botanical observations.
        for obs in observation:
            if obs in original_obs and not prior_observation_unknown(obs):
                overlap_observation_ids.add(obs)
        obs_seen.update(observation)
        sum_state["_photo_ids"]=sum_state.get("_photo_ids",set()) | set(photo)
        for key in ("new_original_photo_attempts","four_hue_classified",
                    "image_acquisition_failed","acquired_but_unclassifiable"):
            sum_state[key]+=rec[key]
        each.append({"group":group,"new_photo_IDs":len(df),
                     "classified_hues":int(classified.sum()),
                     "saved_parquet_sha256":sha(folder/"original_photo_ROIs.parquet"),
                     "group_original_receipt_sha256":sha(recpath)})
    require(seen==set(by_group),"One or more genuine K1 groups omitted from photo ID union")
    require(sum_state["new_original_photo_attempts"]==
                cumulative["completed_real_new_original_photo_attempts"]
            and sum_state["four_hue_classified"]==
                cumulative["completed_automatic_four_hue_classifications"]
            and sum_state["image_acquisition_failed"]==
                cumulative["completed_original_acquisition_failures"]
            and sum_state["acquired_but_unclassifiable"]==
                cumulative["completed_acquired_but_automatic_four_hue_unclassifiable"],
            "Completed K1 original pixel identity/hue totals differ from group receipt census")

    historical_hues=sum(1 for item in previous.values() if item[2] in HUES)
    require(historical_hues == sum(x[1] for x in PRIOR_COUNTS.values()),
            "Historical original ROI four-hue total has unexpected label conflicts")
    global_attempts = len(previous) + sum_state["new_original_photo_attempts"]
    global_labels = historical_hues + sum_state["four_hue_classified"]
    report = {
        "schema":SCHEMA,
        "historical_actual_unique_original_photo_IDs":len(previous),
        "historical_automatic_four_hue_photo_ROIs":historical_hues,
        "new_K1_sha_verified_original_photo_IDs_disjoint_from_history":
            sum_state["new_original_photo_attempts"],
        "new_K1_frozen_four_hue_ROIs":sum_state["four_hue_classified"],
        "new_K1_original_image_acquisition_failures":sum_state["image_acquisition_failed"],
        "new_K1_acquired_but_colour_unclassifiable":sum_state["acquired_but_unclassifiable"],
        "exact_global_historical_plus_K1_distinct_original_photo_IDs":global_attempts,
        "exact_global_historical_plus_K1_automatic_four_hue_photo_ROIs":global_labels,
        "K1_source_finished_original_worker_groups":sorted(seen),
        "K1_source_pending_groups":[i for i in range(NUM_GROUPS) if i not in seen],
        "K1_vs_old_crossarm_exact_photo_ID_overlaps":0,
        "K1_vs_old_crossarm_literal_original_observation_ID_overlaps":
            len(overlap_observation_ids),
        "K1_original_observation_IDs_independent_within_K1":True,
        "K1_vs_historical_observation_level_independence_NOT_assumed":True,
        "historical_photo_source_arms_SHA_verified":historical_audits,
        "historical_crossarm_photo_overlap_summary":overlap,
        "signed_cumulative_original_K1_receipt_set_SHA256":
            cumulative["source_original_group_receipt_aggregate_SHA256"],
        "new_original_image_downloads_or_colour_reinference_in_this_union":0,
        "unweighted_globally_representative_white_flower_rate_estimated":False,
        "flower_blooming_wild_native_stage_independently_verified":False,
        "methodological_scope":"exact original photo-ID union only; partial K1 workers",
        "new_K1_group_ID_and_real_photo_SHA_attestations":sorted(each,key=lambda x:x["group"]),
    }
    output.parent.mkdir(parents=True,exist_ok=True)
    tmp=Path(tempfile.mkdtemp(prefix="mois-global-K1-exact-photos-",dir=output.parent))
    try:
        result=tmp/"historical_plus_K1_source_signed_original_photo_ID_union.json"
        result.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n",encoding="utf-8")
        os.replace(tmp,output)
    except BaseException:
        import shutil
        shutil.rmtree(tmp,ignore_errors=True)
        raise
    print("MOIS_HISTORICAL_AND_K1_TRUE_PHOTO_ID_DEDUP_UNION",json.dumps({
        "historical_photos":len(previous),"K1_new_photos":sum_state["new_original_photo_attempts"],
        "global_unique_photos":global_attempts,"global_classified":global_labels,
        "K1_groups":len(seen),"image_downloads_or_inference":0
    },sort_keys=True),flush=True)
    return report


def run(old, cumulative, groups, output):
    sources=[
        ("main", old/"main"/"combined_actual_photo_checkpoint.parquet",
         old/"main"/"checkpoint_receipt.json"),
        ("prior_group00", old/"group00"/"measured_photos.parquet",
         old/"group00"/"receipt.json"),
        ("prior_R1", old/"r1"/"measured_photos.parquet",old/"r1"/"receipt.json"),
        ("prior_R2", old/"r2"/"measured_photos.parquet",old/"r2"/"receipt.json"),
    ]
    prior,obs,source_audit,overlap,_,_=prior_original_photo_index(sources,PRIOR_COUNTS)
    return merge_verified(prior,obs,source_audit,overlap,read_json(cumulative),groups,output)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--old",type=Path,required=True)
    parser.add_argument("--cumulative",type=Path,required=True)
    parser.add_argument("--group",type=Path,action="append",required=True)
    parser.add_argument("--output",type=Path,required=True)
    p=parser.parse_args()
    run(p.old,p.cumulative,p.group,p.output)


if __name__=="__main__":
    main()
