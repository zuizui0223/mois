#!/usr/bin/env python3
"""Data-only partial FCP5 original-photo cumulative photo/hue receipt.

Genuine SHA-sealed per-group ROI image terminals are the sole input. No
attempt to count the 34,898 uncompleted original-image source opportunities
as pixels. The worldwide 1,009,032-photo baseline was independently
photo-ID checked against the precolour nested5 all-ID source before work.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import tempfile

from scripts.acquisition.mois_completed_pixel_checkpoint import HUES,require,sha
from scripts.acquisition.mois_fcp_u100_nested5_seal_finished_group import (
    SCHEMA as GROUP_SCHEMA, MODEL_RUN, SOURCE_MANIFEST_SHA,
    NEW, WORK_UNITS, WORK_GROUPS, PREVIOUS_ACTUAL,
    REUSED, CANDIDATES, ELIGIBLE_TAXA, SOURCE_PHOTOS,
)
from scripts.acquisition.mois_fcp_u100_k1_full_global_union import (
    SCHEMA as PRIOR_SCHEME,
)

SCHEMA="mois_FCP_nested5_genuine_partial_original_photo_roi_global_cumulative_v1"
PRIOR_GLOBAL_SHA="c924cb6d1a873428db197daecc51c67046c89a95193b66713789dcb0ce77a59c"
PREVIOUS_HUES=445805


def receipt(path, *, planned=NEW, units=WORK_UNITS, groups=WORK_GROUPS,
            source_sha=SOURCE_MANIFEST_SHA, previous=PREVIOUS_ACTUAL,
            reused=REUSED, candidates=CANDIDATES, eligible=ELIGIBLE_TAXA,
            first5=SOURCE_PHOTOS, model_run=MODEL_RUN):
    p=Path(path)
    row=json.loads(p.read_text(encoding="utf-8"))
    g=row.get("original_worker_group")
    require(type(g) is int and 0<=g<groups,
            "Signed FCP nested5 original photo group ID invalid")
    required_units=list(range(g*4,min((g+1)*4,units)))
    actual=min(1024,planned-g*1024)
    require(
        row.get("schema")==GROUP_SCHEMA
        and row.get("frozen_original_model_run")==model_run
        and row.get("source_SHA256")==source_sha
        and row.get("source_original_FCP_U100_first5_nominal_taxa")==eligible
        and row.get("source_101210_first5_original_photos")==first5
        and row.get("source_80968_new_nested_rank2to5_photo_candidates")==candidates
        and row.get("prior_complete_original_photo_IDs_reused_without_inference")==reused
        and row.get("prior_global_actual_photo_IDs_SHA_audited_before_new_selection")==previous
        and row.get("all_source_genuinely_new_original_photo_attempts")==planned
        and row.get("original_completed_source_work_units")==required_units
        and row.get("frozen_original_FCP_yolo_efficientSAM_roi_v4_method_attested") is True
        and row.get("per_original_photo_source_inclusion_probability_preserved") is True
        and row.get("old_photo_ID_reuses_not_relabelled") is True
        and row.get("unknown_technical_ROIs_never_cast_to_white") is True
        and row.get("original_image_bytes_downloaded_OR_colour_ROI_inferred_by_sealer")==0
        and row.get("botanical_wild_native_open_flower_status_independently_verified") is False
        and row.get("global_photo_observation_probability_complete_cohort") is False
        and row.get("worldwide_unbiased_white_flora_frequency_estimated") is False,
        "Completed original FCP5 source/ROI/probability/unknown attestation invalid")
    classifiable=row.get("new_actual_four_hue_classified")
    failed=row.get("new_original_image_access_failed")
    unknown=row.get("new_original_acquired_hue_unclassifiable")
    obs=row.get("new_distinct_original_observations")
    probs=row.get("distinct_original_per_photo_joint_probabilities")
    require(all(type(v) is int for v in (classifiable,failed,unknown,obs,probs))
            and actual>0 and row.get("new_actual_original_photo_ID_attempts")==actual
            and obs==actual and probs>=1
            and classifiable>=0 and failed>=0 and unknown>=0
            and classifiable+failed+unknown==actual,
            "Completed real photo groups must conserve attempts and unknown outcomes")
    hues=row.get("new_original_hue_counts")
    require(isinstance(hues,dict) and set(hues)==set(HUES)|{"mixed_uncertain"}
            and all(type(v) is int and v>=0 for v in hues.values())
            and sum(hues[h] for h in HUES)==classifiable
            and hues["mixed_uncertain"]==failed+unknown,
            "Zero-occurrence genuine four-hue groups must be zero, unknown unchanged")
    audits=row.get("original_terminal_source_receipts_SHA256")
    require(isinstance(audits,list) and len(audits)==len(required_units)
            and [x.get("unit") for x in audits]==required_units
            and all(isinstance(x,dict) and all(
                isinstance(x.get(key),str) and len(x[key])==64
                and all(c in "0123456789abcdef" for c in x[key])
                for key in ("receipt_sha256","photos_sha256","terminal_sha256"))
                for x in audits),
            "Original photo/image-terminal SHA work-unit evidence absent")
    saved=row.get("saved_file_SHA256")
    require(isinstance(saved,dict) and
            {"original_photo_ROIs.parquet","original_photo_hue_counts.csv",
             "original_photo_source_coverage.json"}<=set(saved)
            and all(isinstance(v,str) and len(v)==64
                    and all(x in "0123456789abcdef" for x in v)
                    for v in saved.values()),
            "Published FCP5 original-photo evidence file SHA missing")
    timing=row.get("model_group_wall_seconds")
    perf=row.get("original_photo_prefetch_stats")
    require(isinstance(timing,(float,int)) and math.isfinite(timing) and timing>0
            and isinstance(perf,dict) and perf.get("download_calls")==actual
            and perf.get("consumed")==actual and perf.get("permanent_failures")==failed,
            "Genuine FCP5 original image photo transfer/processing times inconsistent")
    return {"group":g,"attempted":actual,
            "four_hue_classified":classifiable,
            "image_acquisition_failed":failed,
            "acquired_unclassifiable":unknown,
            "hues":hues,"source_group_receipt_SHA256":sha(p),
            "wall_seconds":float(timing)}


def aggregate(group_receipts,previous_json,output, *,
              previous_sha=PRIOR_GLOBAL_SHA,
              prior_attempts=PREVIOUS_ACTUAL,prior_hues=PREVIOUS_HUES,
              planned=NEW,units=WORK_UNITS,groups=WORK_GROUPS,
              reused=REUSED,candidates=CANDIDATES,
              eligible=ELIGIBLE_TAXA,first5=SOURCE_PHOTOS):
    output=Path(output)
    require(not output.exists(),"Signed original-photo stage5 cumulative release immutable")
    require(sha(previous_json)==previous_sha,"Previous measured original 1009032 photo-ID exact union SHA changed")
    old=json.loads(Path(previous_json).read_text())
    require(old.get("schema")==PRIOR_SCHEME
            and old.get("exact_global_unique_original_photo_ID_attempts")==prior_attempts
            and old.get("exact_global_automatic_four_hue_photo_ROIs")==prior_hues,
            "Previously measured K1 plus FCP3141 true photo ROIs not authenticated")
    require(bool(group_receipts),"No completed actual new FCP5 image photo group to count")
    group_rows=[receipt(p,planned=planned,units=units,groups=groups,
                        previous=prior_attempts,reused=reused,
                        candidates=candidates,eligible=eligible,first5=first5)
                for p in group_receipts]
    group_ids=[x["group"] for x in group_rows]
    require(len(group_ids)==len(set(group_ids)),
            "Duplicate signed FCP5 original image-group receipt")
    group_rows.sort(key=lambda x:x["group"])
    sums=Counter()
    morphs=Counter()
    for row in group_rows:
        for field in ("attempted","four_hue_classified","image_acquisition_failed",
                      "acquired_unclassifiable"):
            sums[field]+=row[field]
        morphs.update(row["hues"])
    require(sums["attempted"]==sums["four_hue_classified"]+
                sums["image_acquisition_failed"]+sums["acquired_unclassifiable"]
            and sums["attempted"]<=planned,
            "Cannot convert unprocessed FCP source original photo IDs into finished hue")
    fingerprint=hashlib.sha256("\n".join(
        f'{r["group"]}:{r["source_group_receipt_SHA256"]}'
        for r in group_rows).encode()).hexdigest()
    result={
        "schema":SCHEMA,
        "source_original_candidate_photo_IDs_METADATA_ONLY":93209351,
        "original_FCP_U100_first5_nominal_eligible_taxa":eligible,
        "nested_first5_original_photo_source_IDs":first5,
        "nested_new_rank2to5_source_candidates_NOT_MEASURED":candidates,
        "previous_actual_independently_verified_original_photo_attempts":prior_attempts,
        "previous_actual_independently_verified_four_hue_ROIs":prior_hues,
        "previous_nested_source_exact_original_photo_ID_reuses":reused,
        "new_distinct_original_photo_target_IDs":planned,
        "new_actual_finished_original_photo_groups":len(group_rows),
        "new_finished_original_photo_group_IDs":[x["group"] for x in group_rows],
        "new_pending_model_worker_groups":[x for x in range(groups)
                                            if x not in set(group_ids)],
        "new_genuine_original_photo_pixel_attempts":sums["attempted"],
        "new_genuine_automatic_four_hue_ROIs":sums["four_hue_classified"],
        "new_original_photo_acquisition_failed":sums["image_acquisition_failed"],
        "new_acquired_photo_hue_unclassifiable":sums["acquired_unclassifiable"],
        "new_original_photo_morph_counts":{h:morphs[h] for h in (*HUES,"mixed_uncertain")},
        "exact_precolour_disjoint_source_global_original_photo_attempts":
            prior_attempts+sums["attempted"],
        "exact_precolour_disjoint_source_global_four_hue_ROIs":
            prior_hues+sums["four_hue_classified"],
        "source_sample_K1_FCP_nested5_photo_ID_previous_release_SHA256":previous_sha,
        "new_source_34898_photo_IDs_disjoint_from_1009032_prior_by_frozen_manifest":True,
        "original_group_receipt_set_SHA256":fingerprint,
        "original_group_receipts_recomputed_from_real_terminal_SHA_proofs":True,
        "counts_are_partial_model_completions_not_a_complete_probability_cohort":True,
        "old_original_photo_parquet_downloads_and_colour_inferences":0,
        "raw_original_source_photo_downloads_this_cumulative":0,
        "wild_native_open_flower_status_independently_reviewed":False,
        "global_flora_unweighted_white_hue_prevalence_estimated":False,
        "group_receipts":group_rows,
    }
    output.parent.mkdir(parents=True,exist_ok=True)
    staged=Path(tempfile.mkdtemp(prefix="mois-FCP5-new-genuine-cumulative-",dir=output.parent))
    try:
        (staged/"mois_FCP5_partial_genuine_original_pixel_cumulative.json").write_text(
            json.dumps(result,sort_keys=True,indent=2)+"\n",encoding="utf-8")
        os.replace(staged,output)
    except BaseException:
        shutil.rmtree(staged,ignore_errors=True)
        raise
    print("MOIS_FCP5_GENUINE_ORIGINAL_PIXELS_INCREMENTAL_SOURCE_DISJOINT_CUMULATIVE",json.dumps({
        "groups_sealed":len(group_rows),
        "new_genuine_original_photos":sums["attempted"],
        "new_classified":sums["four_hue_classified"],
        "global_original_photo_attempts":result["exact_precolour_disjoint_source_global_original_photo_attempts"],
        "global_automatic_hue_ROIs":result["exact_precolour_disjoint_source_global_four_hue_ROIs"],
        "new_original_image_downloads_or_colour_inferences":0,
        "source_group_set_SHA256":fingerprint,
    },sort_keys=True),flush=True)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--group-receipt",type=Path,action="append",required=True)
    parser.add_argument("--previous",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    x=parser.parse_args()
    aggregate(x.group_receipt,x.previous,x.output)


if __name__=="__main__":
    main()
