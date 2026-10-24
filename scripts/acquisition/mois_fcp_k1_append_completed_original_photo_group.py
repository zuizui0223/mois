#!/usr/bin/env python3
"""Immediately seal a completed FCP-style K1 original-photo group as Parquet.

This is a DATA-ONLY transformation of genuine previously completed FCP ROI-v4
YOLO11/EfficientSAM workers from GitHub run 40066261878. It never starts models,
fetches original S3 pixels, or rewrites the 929138 prior genuine ROI photo rows.
Unlike the older tier06 source arm, individual K1 photo inclusion propensities
vary by original observation photo multiplicity and taxon x month x latitude
stratum. Preserve them on EACH image; a scalar global probability is invalid.

The SHA-pinned 80035-observation one-photo balanced cohort reused 3280 actual
old original photo results and only scheduled 76753 not-previously-downloaded
source photographs. Those 3280 reuse outcomes remain in their prior receipts,
not counted as new outcomes by this per-worker sealer.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import json
import math
import os
from pathlib import Path
import shutil
import tempfile

from scripts.acquisition.mois_completed_pixel_checkpoint import (
    HUES, collect_wave, require, sha,
)
from scripts.acquisition.mois_warm_original_roi import FROZEN_SOURCE_BLOB
from scripts.acquisition.mois_append_only_second8_original_worker_group import (
    original_photo_source_coverage,
)

SCHEMA = "mois_fcp_K1_93m_balanced_true_original_ROI_append_group_v1"
K1_SOURCE_SHA = "4da66dc0d34f440269c2202f95b0f39ad369531fd625676e8ec69df4ae939f1c"
PARENT_SOURCE_SHA = "c0ba3a316a336816065590540e75034bd8c2b9b7362034163673702103c44977"
SOURCE_ARM = "mois_global93m_tier06_taxon_month_lat_probability_onephoto_K1_v1"
ORIGINAL_RUN = 40066261878
SOURCE_PLANNED = 76753
SELECTED_SOURCE_OBSERVATIONS = 80035
REUSED_OLD_ORIGINALS = 3280


def seal(source, workers, group, output, *,
         source_sha=K1_SOURCE_SHA, run=ORIGINAL_RUN,
         expected_new=SOURCE_PLANNED,
         expected_selected=SELECTED_SOURCE_OBSERVATIONS,
         expected_reused=REUSED_OLD_ORIGINALS):
    """Verify one sealed 1024-photo source group without any raw-photo network."""
    import pandas as pd
    from scripts.analysis.mois_reassemble_original_new4096_measurements import METHOD

    source, workers, output = map(Path, (source, workers, output))
    require(not output.exists(), "Do not overwrite genuine immutable ROI evidence")
    require(sha(source/"manifest.json")==source_sha,
            "FCP 93m original K1 photo source manifest SHA changed")
    plan=json.loads((source/"manifest.json").read_text(encoding="utf-8"))
    require(plan.get("schema")=="mois_full_source_observation_sampling_cohort_v1"
            and plan.get("source_arm")==SOURCE_ARM
            and plan.get("full_frame_photo_ids")==93209351
            and plan.get("source_precolour_balanced_manifest_SHA256")==PARENT_SOURCE_SHA
            and plan.get("new_photo_attempts_planned")==expected_new
            and plan.get("all_balanced_selected_original_observations")==expected_selected
            and plan.get("reused_original_photo_ids")==expected_reused
            and plan.get("original_photo_source_probabilities_unchanged_by_prior_pixel_reuse") is True
            and plan.get("images_downloaded_by_this_source_bridge")==0
            and plan.get("colour_measured") is False
            and plan.get("work_units")==len(plan.get("worker_ids",[]))
            and plan.get("worker_ids")==list(range(plan["work_units"]))
            and plan.get("group_ids")==list(range((plan["work_units"]+3)//4)),
            "Wrong signed FCP-style first-image balanced K1 real source")
    require(type(group) is int and group in plan["group_ids"],
            "Completed K1 original photo group was not in the frozen source")
    expected_units=list(range(group*4,min((group+1)*4,plan["work_units"])))
    rows,audit=collect_wave(source,workers,METHOD)
    require(audit["completed_units"]==expected_units
            and audit["attempted"]==len(rows)>0
            and audit["planned_photo_attempts"]==expected_new
            and set(audit["missing_units"])==set(plan["worker_ids"])-set(expected_units),
            "Genuine 256-image FCP photo work units are missing or repeated")
    photos=[x["photo_id"] for x in rows]
    obs=[x["observation_id"] for x in rows]
    require(len(photos)==len(set(photos)) and all(p.isdecimal() for p in photos)
            and len(obs)==len(set(obs)) and all(obs),
            "K1 must have one original photo per independent selected observation")
    success=0
    known_hues=Counter()
    different_probabilities=set()
    for row in rows:
        known=row["morph"] in HUES
        require(known==(row["measurement_status"]=="classified_four_state_morph")
                and row["morph"] in (*HUES,"mixed_uncertain")
                and row["source_arm"]==SOURCE_ARM,
                "Original image ROI hue/missing status or source arm altered")
        try:p=float(row["selection_probability"])
        except (TypeError,ValueError):raise ValueError("Missing original-photo joint probability")
        require(math.isfinite(p) and 0<p<=1,
                "Original K1 photo/stratum joint inclusion probability invalid")
        different_probabilities.add(p)
        known_hues[row["morph"]]+=1
        success+=known
    require(success==audit["classifiable"]
            and success+audit["acquisition_failed"]+
                audit["unclassifiable_acquired"]==len(rows),
            "K1 actual image photo terminal and nonresponse counters changed")
    diagnostics=list(workers.rglob(f"warm_roi_group_{group:03d}_diagnostic.json"))
    require(len(diagnostics)==1,"One warm FCP original ROI timing attestation required")
    diag=json.loads(diagnostics[0].read_text(encoding="utf-8"))
    require(diag.get("group")==group
            and diag.get("separate_original_main_calls")==len(expected_units)
            and diag.get("frozen_yolo_and_onnx_initializations")==1
            and diag.get("frozen_source_blob")==FROZEN_SOURCE_BLOB
            and diag.get("new_photos_downloaded")==len(rows)
            and diag.get("method_or_thresholds_modified") is False
            and diag.get("old_photo_RGB_workers_modified") is False,
            "Previously sealed original ROI worker altered its frozen FCP model or image count")

    output.parent.mkdir(parents=True,exist_ok=True)
    tmp=Path(tempfile.mkdtemp(prefix="mois-fcp-K1-original-pixel-",dir=output.parent))
    try:
        frame=pd.DataFrame(rows)
        frame["photo_id"]=frame.photo_id.astype(str)
        frame=frame.sort_values("photo_id",key=lambda x:x.map(int)).reset_index(drop=True)
        parquet=tmp/"original_photo_ROIs.parquet"
        frame.to_parquet(parquet,index=False,compression="zstd")
        coverage=original_photo_source_coverage(rows)
        coverage["schema"]="mois_fcp_K1_signed_original_taxon_recorded_month_latitude_coverage_v1"
        coverage["sampled_original_photo_inclusion_probability_varies_per_photo"]=True
        coverage["one_selected_photo_per_original_observation"]=True
        cov=tmp/"original_photo_source_coverage.json"
        cov.write_text(json.dumps(coverage,indent=2,sort_keys=True)+"\n")
        huefile=tmp/"original_photo_hue_counts.csv"
        with huefile.open("w",encoding="utf-8",newline="") as f:
            w=csv.writer(f,lineterminator="\n")
            w.writerow(["morph","original_photo_attempts"])
            for hue in (*HUES,"mixed_uncertain"):w.writerow([hue,known_hues[hue]])
        receipt={
            "schema":SCHEMA,
            "frozen_original_source_run":run,
            "full_93m_photo_target_catalogue":93209351,
            "original_source_manifest_SHA256":source_sha,
            "original_precolour_probability_source_SHA256":PARENT_SOURCE_SHA,
            "source_arm":SOURCE_ARM,
            "previously_measured_old_original_photo_IDs_reused_outside_this_group":expected_reused,
            "original_new_selected_photo_targets_all_groups":expected_new,
            "original_worker_group":group,
            "original_completed_source_units":expected_units,
            "original_photo_terminal_source_receipts_SHA256":audit["completed_unit_audit"],
            "actual_new_original_photo_attempts":len(rows),
            "actual_automatic_four_hue_ROI_classifications":success,
            "original_image_acquisition_failed":audit["acquisition_failed"],
            "original_acquired_but_four_hue_unclassifiable":audit["unclassifiable_acquired"],
            "original_distinct_observations":len(set(obs)),
            "distinct_conditional_joint_photo_probabilities":len(different_probabilities),
            "original_hue_counts":{h:known_hues[h] for h in (*HUES,"mixed_uncertain")},
            "frozen_original_FCP_roi_v4_source_method_attested":True,
            "original_photo_identity_and_image_terminal_SHA_audited":True,
            "source_inclusion_probability_retained_at_each_photo":True,
            "original_warm_FCP_model_inits":1,
            "original_full_roi_wall_seconds":diag["group_wall_seconds"],
            "original_full_roi_evaluator_seconds":diag["total_frozen_evaluator_seconds"],
            "original_photo_prefetch":diag.get("photo_prefetch",{}),
            "image_downloads_OR_ROI_reinferences_in_this_sealer":0,
            "old_historical_original_photo_parquets_rewritten":False,
            "technical_unknown_ROIs_not_turned_into_other_hues":True,
            "true_botanical_plant_flower_or_native_status_validated":False,
            "probability_weighted_global_flora_white_rate_estimated":False,
            "artifact_SHA256":{
                parquet.name:sha(parquet),
                cov.name:sha(cov),
                huefile.name:sha(huefile),
                diagnostics[0].name:sha(diagnostics[0])
            },
        }
        (tmp/"original_K1_group_receipt.json").write_text(
            json.dumps(receipt,indent=2,sort_keys=True)+"\n")
        os.replace(tmp,output)
        print("MOIS_FCP_K1_GENUINE_ORIGINAL_PHOTO_ROI_PARTITION",
              json.dumps({k:receipt[k] for k in (
                  "original_worker_group","actual_new_original_photo_attempts",
                  "actual_automatic_four_hue_ROI_classifications",
                  "original_image_acquisition_failed","original_full_roi_wall_seconds")},
              sort_keys=True),flush=True)
        return receipt
    except BaseException:
        shutil.rmtree(tmp,ignore_errors=True)
        raise


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("source","workers","output"):
        p.add_argument("--"+key,type=Path,required=True)
    p.add_argument("--group",type=int,required=True)
    args=p.parse_args()
    seal(args.source,args.workers,args.group,args.output)


if __name__=="__main__":
    main()
