#!/usr/bin/env python3
"""Data-only all-group collector of NEW FCP depth1 authentic original-photo ROI.

Only the NEW 3,141 original images which were NEVER present in signed FCP/MOIS
existing photo-ID ledgers or the 80,035 in-flight K1 photo source are measured.
The other 17,101 source candidates remain prior pixel reuse or reserved K1
and are NOT added here as new colour outcomes.

The existing production collect_wave performs independent original source
photo/observation and terminal/image/model SHA comparisons per each FCP
256-photo work unit, with source-identity and unclassifiable/failed outcomes.
Never fetch licensed original images or run detector/segmentation at collection.
"""
from __future__ import annotations
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import math
import pandas as pd

from scripts.acquisition.mois_completed_pixel_checkpoint import (
    HUES,collect_wave,sha,require,
)
from scripts.analysis.mois_reassemble_original_new4096_measurements import METHOD

SCHEMA="mois_fcp_U100_progressive_prefix1_new_real_original_photo_ROI_v1"
SOURCE_SHA="fe49344241fbc051eebbbcbc8eb627a040198dc2b3ee56c89401876b58d28a6"
EXPECTED=3_141
SOURCE_N=20_242
PRIOR_REUSED=15_483
K1_RESERVED=1_618
MODEL_BLOB="b1d7091db44682559bbe4b48a9f3b8db648501a6"


def seal(source,workers,output,*,expected_source_sha=SOURCE_SHA,
         expected=EXPECTED,expected_source_n=SOURCE_N,
         prior_reused=PRIOR_REUSED,k1_reserved=K1_RESERVED):
    source,workers,output=map(Path,(source,workers,output))
    require(not output.exists(),"Do not overwrite authentic original-photo ROI outcome release")
    path=source/"manifest.json"
    require(sha(path)==expected_source_sha,
            "Signed FCP U100 first-photo/K1 prior reuse original source changed")
    plan=json.loads(path.read_text())
    require(plan.get("schema")=="mois_full_source_observation_sampling_cohort_v1"
            and plan.get("fcp_prefix_source_schema")==
                  "mois_FCP_U100_full93m_progressive_prefix1_prior_and_active_K1_photo_reuse_v1"
            and plan.get("selected_FCP_U100_eligible_source_taxa")==expected_source_n
            and plan.get("full_frame_photo_ids")==93_209_351
            and plan.get("new_photo_attempts_planned")==expected
            and plan.get("previously_measured_original_photo_ids_reused")==prior_reused
            and plan.get("inflight_K1_photo_IDs_reserved_without_redownload")==k1_reserved
            and expected_source_n==prior_reused+k1_reserved+
                plan.get("technical_pilot_photo_IDs_not_refetched",0)+
                plan.get("unsupported_original_photo_format_selected_as_unknown",0)+expected
            and plan.get("photo_already_reserved_by_K1_never_fetched_twice") is True
            and plan.get("FCP_top100_photo_prefix_inclusion_probability_corrected_to_1_over_N_M") is True
            and plan.get("colour_measured") is False,
            "Original FCP source/work schedule changed or previous pixels were double counted")
    rows,audit=collect_wave(source,workers,METHOD)
    require(audit["wave_complete"] and audit["attempted"]==expected
            and len(rows)==expected and audit["planned_photo_attempts"]==expected,
            "Original FCP photo source not all fully sealed; do not publish partial as complete")
    ids=[r["photo_id"] for r in rows]
    require(len(set(ids))==expected and all(pid.isdecimal() for pid in ids),
            "New original-photo identity duplicated or malformed in complete colour result")
    flags=Counter()
    obs=set()
    source_prob_by_obs={}
    for r in rows:
        known=r["morph"] in HUES
        failed=r["measurement_status"]=="image_acquisition_failed"
        require(r["morph"] in (*HUES,"mixed_uncertain")
                and known==(r["measurement_status"]=="classified_four_state_morph"),
                "Original photo image missing/uncertain was relabeled as visible flower petals")
        p=float(r["selection_probability"])
        require(math.isfinite(p) and 0<p<=1 and
                r["source_arm"]==plan["source_arm"],
                "Original FCP U100 prefix1 inclusion probability/source arm changed")
        if r["observation_id"] in source_prob_by_obs:
            require(source_prob_by_obs[r["observation_id"]]==p,
                    "Same source observation assigned multiple probability designs")
        source_prob_by_obs[r["observation_id"]]=p
        obs.add(r["observation_id"])
        flags["attempted"]+=1
        flags["classified"]+=known
        flags["failed"]+=failed
        flags["unclassified_acquired"]+=not known and not failed
        flags[r["morph"]]+=1
    require(flags["classified"]==audit["classifiable"]
            and flags["failed"]==audit["acquisition_failed"]
            and flags["unclassified_acquired"]==audit["unclassifiable_acquired"]
            and len(obs)==expected,
            "Original FCP U100 single-observation source/classification terminal accounting inconsistent")
    # All groups have exactly one old frozen model session once, unless a
    # recovery run reuses all completed source units from an older job.
    group_signatures=[]
    group_diagnostics={}
    for diagnostic in workers.rglob("warm_roi_group_*_diagnostic.json"):
        d=json.loads(diagnostic.read_text())
        group=d["group"]
        require(group in plan["group_ids"] and group not in group_diagnostics
                and d["frozen_source_blob"]==MODEL_BLOB
                and d["method_or_thresholds_modified"] is False
                and d["artificial_flower_labels_created"] is False
                and d["frozen_yolo_and_onnx_initializations"] in (0,1),
                "FCP original warm source ROI method or group diagnostic malformed")
        group_diagnostics[group]=d
        group_signatures.append({
            "group":group,
            "source_warm_diagnostic_SHA256":sha(diagnostic),
            "group_wall_seconds":d["group_wall_seconds"],
            "original_photo_fetch_attempts_in_this_worker":d["new_photos_downloaded"],
            "fixed_model_session_initializations":d["frozen_yolo_and_onnx_initializations"],
            "photo_loader_seconds":d["photo_prefetch"]["loader_seconds"],
            "photo_consumer_wait_seconds":d["photo_prefetch"]["consumer_wait_seconds"],
        })
    require(set(group_diagnostics)==set(plan["group_ids"]),
            "Original warm model source group diagnostic absent from completed photo workers")
    output.mkdir(parents=True,exist_ok=False)
    frame=pd.DataFrame(rows)
    frame["photo_id"]=frame["photo_id"].astype(str)
    frame=frame.sort_values("photo_id",key=lambda x:x.map(int)).reset_index(drop=True)
    photo_file=output/"new_FCP_U100_prefix1_original_photo_ROIs.parquet"
    frame.to_parquet(photo_file,index=False,compression="zstd")
    hue_file=output/"new_FCP_U100_prefix1_hue_counts.csv"
    with hue_file.open("w",newline="") as f:
        w=csv.writer(f,lineterminator="\n")
        w.writerow(["morph","original_photo_terminal_attempts"])
        for color in (*HUES,"mixed_uncertain"):
            w.writerow([color,flags[color]])
    timing_file=output/"original_warm_worker_throughput.json"
    timing_file.write_text(json.dumps(sorted(group_signatures,key=lambda x:x["group"]),
                                      indent=2,sort_keys=True)+"\n")
    receipt={
        "schema":SCHEMA,
        "source_full_93m_original_photo_ids":93_209_351,
        "FCP_U100_first_observation_taxon_IDs_in_source":expected_source_n,
        "prior_real_original_photo_ROI_reuse_without_new_inference":prior_reused,
        "overlapping_inflight_K1_photos_not_refetched":k1_reserved,
        "new_unique_original_photo_attempts":expected,
        "new_automatically_four_hue_ROI_classified":flags["classified"],
        "new_image_acquisition_failed":flags["failed"],
        "new_acquired_unclassifiable":flags["unclassified_acquired"],
        "new_unique_original_observations":len(obs),
        "original_frozen_model_fcp_roi_v4_method":METHOD,
        "new_original_photo_hue_counts":{h:flags[h] for h in (*HUES,"mixed_uncertain")},
        "original_photo_source_plan_sha256":sha(path),
        "independent_terminal_and_photo_source_ROI_SHA_verified":True,
        "model_photo_worker_groups_source_and_DIAG_sha_verified":True,
        "photo_group_diagnostics":sorted(group_signatures,key=lambda x:x["group"]),
        "per_original_photo_stage1_sample_inclusion_probabilities_retained":True,
        "historical_old_original_photo_ROIs_rewritten":False,
        "new_photo_fetch_or_ROI_inference_in_this_collector":0,
        "individual_photo_bloom_native_wild_status_confirmed":False,
        "FCP_sample_not_unbiased_global_native_flora_white_prevalence":True,
        "output_sha256":{
            photo_file.name:sha(photo_file),
            hue_file.name:sha(hue_file),
            timing_file.name:sha(timing_file),
        },
    }
    (output/"new_FCP_U100_prefix1_original_photo_receipt.json").write_text(
        json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print("MOIS_FCP_PREFIX1_NEW_FROZEN_ORIGINAL_PHOTO_ROIS_SEALED",
          json.dumps({k:receipt[k] for k in (
              "new_unique_original_photo_attempts","new_automatically_four_hue_ROI_classified",
              "new_image_acquisition_failed","new_acquired_unclassifiable")},sort_keys=True),flush=True)
    return receipt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("source","workers","output"):
        p.add_argument("--"+key,type=Path,required=True)
    args=p.parse_args()
    seal(args.source,args.workers,args.output)


if __name__=="__main__":
    main()
