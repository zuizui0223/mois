#!/usr/bin/env python3
"""Extend original MOIS cross-arm unique-photo ledger with SHA-signed p0032 delta.

Only read the *new* p0032 original-photo Parquet and two small, independently
sealed tier06 arm Parquets. Prior 887,600-photo base was already fully verified
against those two arms in the immutable old cross-arm registry, and this delta
was independently photo-ID-deduplicated against that base by the p0032 collector.
Therefore do not re-download historic six-wave worker ZIPs or rewrite 0.9m rows.

All actual photo/observation IDs remain distinct ACROSS sampling arms. Repeated
original observations INSIDE the same source arm are allowed if photos differ.
No ecological white-flower prevalence from mixed, incomplete source arms.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re

from scripts.acquisition.mois_completed_pixel_checkpoint import HUES, sha, require

SCHEMA="mois_original_RGB_crossarm_plus_new_p0032_signed_delta_v1"
OLD_UNION_TAG="mois-original-RGB-cross-arm-append-index-main39889436140-groups-0-1-2-3-4-v1"
BASE_RUN=39889436140
BASE_ATTEMPTS=887600
BASE_HUES=400904
OLD_FIRST_ATTEMPTS=4138
OLD_FIRST_HUES=1573
OLD_BROAD_ATTEMPTS=4468
OLD_BROAD_HUES=1760
OLD_UNION_ATTEMPTS=896206
OLD_UNION_HUES=404237
DELTA_ATTEMPTS=12288
DELTA_HUES=4760
DELTA_GROUPS=list(range(121,133))
DELTA_SCHEMA="mois_p0032_new_sealed_original_RGB_photo_delta_from_signed_checkpoint_v1"
EXPECTED_NEW_ATTEMPTS=908494
EXPECTED_NEW_HUES=408997


def read_original_ids(path,receipt,*,n,hue,previous=False):
    import pandas as pd
    p=Path(path)
    d=json.loads(Path(receipt).read_text(encoding="utf-8"))
    require(d.get("schema")=="mois_global93m_observation_baseline_real_RGB_photo_results_v1"
            and d.get("attempted")==n
            and d.get("classifiable")==hue
            and d.get("measured_photo_parquet_sha256")==sha(p),
            "Historical independent tier06 original ROI photo rows not authenticated")
    df=pd.read_parquet(p,columns=["photo_id","observation_id"])
    ids=df.photo_id.astype(str)
    obs=df.observation_id.fillna("").astype(str)
    require(len(ids)==n and ids.is_unique and ids.str.fullmatch(r"[0-9]+").all()
            and obs.str.strip().ne("").all(),
            "Historical independently signed original photo/observation identities are missing")
    return set(ids),set(obs)


def merge_cross_arm(old_registry, delta_parquet, delta_receipt,
                    old_tier06_parquet, old_tier06_receipt,
                    broad_parquet, broad_receipt, output,
                    *,old_total=OLD_UNION_ATTEMPTS,old_hues=OLD_UNION_HUES,
                    delta_n=DELTA_ATTEMPTS,delta_hues=DELTA_HUES,
                    delta_groups=DELTA_GROUPS, expected_new=EXPECTED_NEW_ATTEMPTS,
                    expected_new_hues=EXPECTED_NEW_HUES,
                    main_run=BASE_RUN,base_n=BASE_ATTEMPTS,base_hues=BASE_HUES,
                    old_tier06_n=OLD_FIRST_ATTEMPTS,old_tier06_hues=OLD_FIRST_HUES,
                    broad_n=OLD_BROAD_ATTEMPTS,broad_hues=OLD_BROAD_HUES):
    import pandas as pd

    out=Path(output)
    require(not out.exists(),"Original photo cross-arm ledger cannot be rewritten")
    old=json.loads(Path(old_registry).read_text(encoding="utf-8"))
    new=json.loads(Path(delta_receipt).read_text(encoding="utf-8"))
    require(old.get("schema")=="mois_original_pixel_immutable_cross_arm_group_registry_v1"
            and old.get("main_checkpoint_run")==main_run
            and old.get("main_checkpoint_photo_attempts")==base_n
            and old.get("main_checkpoint_four_hue_classified")==base_hues
            and old.get("previous_tier06_source_photo_attempts")==old_tier06_n
            and old.get("previous_tier06_four_hue_classified")==old_tier06_hues
            and old.get("broad_completed_photo_attempts")==broad_n
            and old.get("broad_completed_four_hue_classified")==broad_hues
            and old.get("broad_pending_photo_attempts")==0
            and old.get("original_photo_attempts_distinct_across_all_verified_arms")==old_total
            and old.get("four_hue_classifications_across_all_verified_arms")==old_hues
            and old.get("original_photo_or_observation_ID_cross_arm_overlap")==0
            and old.get("all_original_photo_per_partition_SHA_receipts_verified") is True,
            "Original cross-arm full photo/photo-source arm SHA union is unverified")
    require(len(old["broad_sealed_groups"])==5,
            "Previous independently sealed all16 broad photo group audits incomplete")
    require(new.get("schema")==DELTA_SCHEMA
            and new.get("prior_immutable_main_checkpoint_run")==main_run
            and new.get("prior_immutable_main_photo_attempts")==base_n
            and new.get("prior_immutable_main_four_hue_classified")==base_hues
            and new.get("completed_true_original_photo_attempts")==delta_n
            and new.get("new_true_original_four_hue_classified")==delta_hues
            and new.get("new_sealed_p0032_groups")==delta_groups
            and new.get("previous_photo_ids_and_new_pixel_ids_exactly_disjoint") is True
            and new.get("unchanged_original_roi_model_and_image_terminal_byte_SHA_verified") is True
            and new.get("new_original_photos_downloaded_from_S3")==0
            and new.get("new_roi_or_hue_inference_performed") is False
            and new.get("previous_photo_dataset_rewritten") is False
            and new["combined_disjoint_main_plus_delta_photo_attempts"]==base_n+delta_n
            and new["combined_disjoint_main_plus_delta_four_hue_classified"]==base_hues+delta_hues,
            "New real original p0032 image pixel delta not SHA certified against same main")
    p=Path(delta_parquet)
    require(new["immutable_output_sha256"][p.name]==sha(p),
            "Sealed p0032 new original image photo Parquet hash changed")
    columns=["photo_id","observation_id","morph","measurement_status",
             "selection_probability","source_arm","image_sha256"]
    data=pd.read_parquet(p,columns=columns)
    ids=data.photo_id.astype(str)
    obs=data.observation_id.fillna("").astype(str)
    classified=data.morph.isin(HUES)
    require(len(ids)==delta_n and ids.is_unique
            and ids.str.fullmatch(r"[0-9]+").all()
            and obs.str.strip().ne("").all()
            and int(classified.sum())==delta_hues
            and data.morph.isin((*HUES,"mixed_uncertain")).all()
            and classified.eq(data.measurement_status.eq("classified_four_state_morph")).all(),
            "New original p0032 photo-ID/observation/hue terminal status mismatch")
    # The original nested p0032 source can contain multiple source-arm labels,
    # and older manifests need not advertise a single scalar source_arm.
    # collect_wave() ALREADY checked every observed per-photo source_arm and
    # observation propensity byte-for-byte against the SHA-pinned source CSV.
    # Never coerce actual per-photo heterogeneous propensities to one number.
    probabilities=pd.to_numeric(data.selection_probability,errors="coerce")
    require(probabilities.notna().all() and
            probabilities.between(0,1,inclusive="right").all(),
            "New original p0032 per-photo observation inclusion probabilities lost")
    actual_arms=data.source_arm.fillna("").astype(str).str.strip()
    advertised_arm=new.get("source_arm")
    if advertised_arm is not None and str(advertised_arm).strip():
        require(actual_arms.eq(str(advertised_arm)).all(),
                "New original p0032 source-arm receipt contradicted original worker source")
    source_arm_count=len(set(actual_arms))
    source_probabilities_have_multiple_values=probabilities.nunique()>1
    deltap=set(ids);deltaobs=set(obs)
    earlyp,earlyobs=read_original_ids(
        old_tier06_parquet,old_tier06_receipt,n=old_tier06_n,hue=old_tier06_hues)
    broadp,broadobs=read_original_ids(
        broad_parquet,broad_receipt,n=broad_n,hue=broad_hues)
    require(not deltap.intersection(earlyp|broadp)
            and not deltaobs.intersection(earlyobs|broadobs),
            "New genuine p0032 photos or observations overlap another existing sampling arm")
    require(not earlyp.intersection(broadp)
            and not earlyobs.intersection(broadobs),
            "Previously accepted independent tier06 sources have unexpectedly overlapping identities")
    require(old_total+delta_n==expected_new and old_hues+delta_hues==expected_new_hues,
            "Genuine published original photo and hue count addition changed")
    report={
        "schema":SCHEMA,
        "historical_immutable_cross_arm_SHA256":sha(old_registry),
        "new_p0032_original_photo_delta_parquet_SHA256":sha(delta_parquet),
        "new_p0032_original_photo_delta_receipt_SHA256":sha(delta_receipt),
        "signed_old_tier06_first_group_parquet_SHA256":sha(old_tier06_parquet),
        "signed_old_tier06_all16_first4_parquet_SHA256":sha(broad_parquet),
        "prior_historical_main_checkpoint_run":main_run,
        "prior_historical_main_photo_attempts":base_n,
        "prior_historical_main_four_hue_classified":base_hues,
        "prior_historical_cross_arm_distinct_photo_attempts":old_total,
        "prior_historical_cross_arm_four_hue_classified":old_hues,
        "newly_completed_p0032_worker_groups":delta_groups,
        "newly_completed_original_p0032_photo_attempts":delta_n,
        "newly_completed_original_p0032_four_hue_classified":delta_hues,
        "newly_completed_original_p0032_distinct_observations":len(deltaobs),
        "newly_completed_p0032_source_arm_labels_present":source_arm_count,
        "newly_completed_p0032_original_per_photo_probability_values_vary":
            bool(source_probabilities_have_multiple_values),
        "p0032_original_source_arm_and_probability_verified_by_original_full_source_worker_SHAs":True,
        "photo_ID_overlap_new_p0032_vs_other_original_sampling_arms":0,
        "observation_ID_overlap_new_p0032_vs_other_original_sampling_arms":0,
        "full_authenticated_distinct_original_photo_attempts":old_total+delta_n,
        "full_authenticated_four_hue_original_image_classifications":old_hues+delta_hues,
        "historical_six_arm_parquet_rewritten":False,
        "historical_other_waves_worker_artifacts_fetched":0,
        "new_original_images_downloaded_or_hues_reinferred_here":0,
        "separately_source_arm_audited":True,
        "native_wild_or_anthesis_image_status_confirmed":False,
        "unweighted_global_flora_white_rate_inferred":False,
        "source_tier06_scheduling_photo_counts_are_not_independent_blooming_species":True,
    }
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("MOIS_CROSS_ARM_REAL_PIXELS_INCREMENTAL_SOURCE_SHA_VERIFIED",
          json.dumps({k:report[k] for k in (
              "full_authenticated_distinct_original_photo_attempts",
              "full_authenticated_four_hue_original_image_classifications",
              "newly_completed_original_p0032_photo_attempts",
              "newly_completed_original_p0032_four_hue_classified")},sort_keys=True),flush=True)
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ("old-registry","delta-parquet","delta-receipt",
                "old-tier06-parquet","old-tier06-receipt",
                "broad-parquet","broad-receipt","output"):
        parser.add_argument("--"+key,type=Path,required=True)
    args=parser.parse_args()
    merge_cross_arm(args.old_registry,args.delta_parquet,args.delta_receipt,
                    args.old_tier06_parquet,args.old_tier06_receipt,
                    args.broad_parquet,args.broad_receipt,args.output)


if __name__=="__main__":
    main()
