#!/usr/bin/env python3
"""Verify and publish ONLY finished immutable original tier06 broad all16 RGB evidence.

All 4,468 photos have already been downloaded and measured by unchanged FCP
ROI-v4 in GitHub run 39880330611. Its original collector created a valid full
measured Parquet and source receipt, but the *release shell* had a missing fi.
This repair compares the collector photo-by-photo against five independently
SHA-sealed append-only worker Parquets, an immutable 93m precolour source
manifest and a prior published cross-arm photo+observation uniqueness ledger.

Zero photo downloads, zero model start-ups, zero colour inference. The 4,468
photo subset is NOT a probability-complete worldwide flora sample. Descriptive
season/realm coverage is not proof of local flowering stage or nativeness.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path
import re

from scripts.acquisition.mois_completed_pixel_checkpoint import HUES, sha, require
from scripts.analysis.mois_reassemble_original_new4096_measurements import METHOD

PROTOCOL = "mois_broad_tier06_4468_signed_original_RGB_completed_recovery_v1"
RUN = 39880330611
SOURCE_SHA = "3e28ced8de85709ea4854cd3291aa00f56e98306b60a4c273c5c7bb47b68c858"
N_ORIGINAL = 4468
N_CLASSIFIED = 1760
N_ACQ_FAILED = 7
N_UNCLASSIFIABLE = 2701
MAIN_CHECKPOINT = 39889436140
MAIN_PHOTOS = 887600
MAIN_CLASSIFIED = 400904
OLDER_TIER06_PHOTOS = 4138
OLDER_TIER06_CLASSIFIED = 1573
CROSS_ARM_TOTAL = 896206
CROSS_ARM_HUES = 404237
PHOTO_FIELDS = (
    "photo_id", "measurement_id", "observation_id", "inat_taxon_id",
    "source_arm", "selection_probability", "morph", "measurement_status",
    "image_sha256",
)


def month_str(value: object) -> str:
    raw = str(value or "").strip()
    match = re.match(r"^([12][0-9]{3})-(0[1-9]|1[012])-(0[1-9]|[12][0-9]|3[01])$", raw)
    return match.group(2) if match else "unknown"


def lat_band(value: object) -> str:
    try:
        z = float(value)
    except (ValueError, TypeError):
        return "unknown"
    if not math.isfinite(z) or not -90 <= z <= 90:
        return "unknown"
    if -23.44 <= z <= 23.44:
        return "tropical"
    return "south_extratropical" if z < -23.44 else "north_extratropical"


def photo_identity_and_hues(frame, label):
    import pandas as pd
    require(set(PHOTO_FIELDS) <= set(frame.columns),label+": expected original photo/ROI fields absent")
    photos = frame["photo_id"].astype(str)
    require(photos.str.fullmatch(r"[0-9]+").all() and photos.is_unique,
            label+": duplicated/malformed original-photo ID")
    observations = frame["observation_id"].fillna("").astype(str).str.strip()
    require(observations.ne("").all(),label+": original observation identity is missing")
    classified = frame["morph"].isin(HUES)
    require(frame["morph"].isin((*HUES, "mixed_uncertain")).all()
            and classified.eq(frame["measurement_status"].eq("classified_four_state_morph")).all(),
            label+": unknown/white/colour status was changed")
    acquisition_fail = frame["measurement_status"].eq("image_acquisition_failed")
    require(not (classified & acquisition_fail).any(),
            label+": missing source image assigned flower colour")
    for row in frame.loc[~acquisition_fail, "image_sha256"].fillna("").astype(str):
        require(bool(re.fullmatch(r"[0-9a-f]{64}",row)), label+": acquired original image SHA missing")
    return set(photos),set(observations),int(classified.sum()),int(acquisition_fail.sum())


def validate_collector_against_append_groups(
    collector, source_manifest, group_files, cross_arm,
    *, source_sha=SOURCE_SHA, expected_original=N_ORIGINAL,
    expected_classified=N_CLASSIFIED, expected_acq_failed=N_ACQ_FAILED,
    expected_unknown=N_UNCLASSIFIABLE, expected_groups=5,
    current_run=RUN, cross_arm_total=CROSS_ARM_TOTAL,
    cross_arm_hues=CROSS_ARM_HUES, main_run=MAIN_CHECKPOINT,
    main_photos=MAIN_PHOTOS, main_classified=MAIN_CLASSIFIED,
    older_tier06_photos=OLDER_TIER06_PHOTOS,
    older_tier06_classified=OLDER_TIER06_CLASSIFIED,
):
    import pandas as pd

    collector=Path(collector);source_manifest=Path(source_manifest)
    source=json.loads(source_manifest.read_text(encoding="utf-8"))
    require(sha(source_manifest)==source_sha
            and source.get("schema")=="mois_full_source_observation_sampling_cohort_v1"
            and source.get("source_group")=="all16"
            and source.get("source_tier")==6
            and source.get("admitted_bucket_count")==4
            and source.get("new_photo_attempts_planned")==expected_original
            and source.get("not_probability_complete_global_source") is True,
            "Complete original photo source plan changed or is not outcome-blind all16")
    receipt_path=collector/"receipt.json"
    receipt=json.loads(receipt_path.read_text(encoding="utf-8"))
    full_parquet=collector/"measured_photos.parquet"
    require(receipt.get("schema")=="mois_global93m_observation_baseline_real_RGB_photo_results_v1"
            and receipt.get("attempted")==expected_original
            and receipt.get("classifiable")==expected_classified
            and receipt.get("acquisition_failed")==expected_acq_failed
            and receipt.get("unclassifiable_acquired")==expected_unknown
            and receipt.get("downloaded")==expected_original-expected_acq_failed
            and receipt.get("measured_photo_parquet_sha256")==sha(full_parquet)
            and receipt.get("source_plan_sha256")==source_sha
            and receipt.get("frozen_original_ROI_method")==METHOD
            and receipt.get("raw_licensed_photo_pixels_redistributed") is False,
            "Genuine original photo full collector/ROI source receipts fail validation")
    require((collector/"source.json").read_bytes()==source_manifest.read_bytes(),
            "Collector source metadata does not equal immutable original source")
    full=pd.read_parquet(full_parquet)
    ids,obs,classified,failed=photo_identity_and_hues(full,"full original collector")
    require(len(ids)==expected_original and classified==expected_classified and failed==expected_acq_failed,
            "Original full photo pixel terminal counts changed")
    require(full["source_arm"].eq(source["source_arm"]).all(),
            "An unknown geographic/taxon source arm entered collector")
    prob=pd.to_numeric(full["selection_probability"],errors="coerce")
    require(prob.notna().all() and prob.eq(source["selection_probability"]).all(),
            "Source observation inclusion probability changed")

    audit=json.loads(Path(cross_arm).read_text(encoding="utf-8"))
    require(audit.get("schema")=="mois_original_pixel_immutable_cross_arm_group_registry_v1"
            and audit.get("main_checkpoint_run")==main_run
            and audit.get("main_checkpoint_photo_attempts")==main_photos
            and audit.get("main_checkpoint_four_hue_classified")==main_classified
            and audit.get("previous_tier06_source_photo_attempts")==older_tier06_photos
            and audit.get("previous_tier06_four_hue_classified")==older_tier06_classified
            and audit.get("broad_completed_photo_attempts")==expected_original
            and audit.get("broad_completed_four_hue_classified")==expected_classified
            and audit.get("broad_pending_photo_attempts")==0
            and audit.get("original_photo_attempts_distinct_across_all_verified_arms")==cross_arm_total
            and audit.get("four_hue_classifications_across_all_verified_arms")==cross_arm_hues
            and audit.get("original_photo_or_observation_ID_cross_arm_overlap")==0
            and audit.get("all_original_photo_per_partition_SHA_receipts_verified") is True,
            "Independent completed cross-arm photo/observation nonoverlap proof missing")
    require(len(group_files)==expected_groups
            and len(audit["broad_sealed_groups"])==expected_groups,
            "Not every original model worker has an immutable source/terminal group release")
    signed_by_group={int(x["group"]):x for x in audit["broad_sealed_groups"]}
    require(set(signed_by_group)==set(range(expected_groups)),
            "Cross-arm group receipt set differs from fully completed source")
    union_parts=[]
    total=Counter()
    for folder in group_files:
        folder=Path(folder)
        group_record=json.loads((folder/"original_photo_group_receipt.json").read_text())
        group=group_record.get("original_source_warm_group")
        require(type(group) is int and group in signed_by_group
                and group_record.get("schema")=="mois_tier06_all16_append_only_original_RGB_group_v1"
                and group_record.get("original_acquisition_workflow_run")==current_run
                and group_record.get("source_original_photo_manifest_sha256")==source_sha
                and group_record.get("all_source_photo_and_original_terminal_SHA_verified") is True
                and group_record.get("frozen_original_FCP_model_and_output_verified") is True
                and group_record.get("image_bytes_reacquired_in_sealing") is False,
                "Immutable original per-group ROI source and model attestation missing")
        g=signed_by_group[group]
        parcel=folder/"original_photo_ROIs.parquet"
        hues=folder/"group_original_photo_hue_counts.csv"
        record=folder/"original_photo_group_receipt.json"
        require(g["original_worker_receipt_sha256"]==sha(record)
                and g["original_photo_parquet_sha256"]==sha(parcel)
                and g["original_hue_count_csv_sha256"]==sha(hues)
                and group_record["saved_file_sha256"][parcel.name]==sha(parcel)
                and group_record["saved_file_sha256"][hues.name]==sha(hues),
                "Immutable original per-group Parquet/receipt bytes were altered")
        frame=pd.read_parquet(parcel)
        _,_,n,h=photo_identity_and_hues(frame,"original group")
        require(len(frame)==int(g["photo_attempts"])
                and n==int(g["four_hue_classified"])
                and h==int(g["acquisition_failed"]),
                "Original group morph/acquisition failure count differs from cross-arm audit")
        total["attempted"]+=len(frame)
        total["classified"]+=n
        total["failed"]+=h
        union_parts.append(frame)
    combined=pd.concat(union_parts,ignore_index=True)
    xids,_,xclass,xfail=photo_identity_and_hues(combined,"all append-only ROI groups")
    require(ids==xids and len(combined)==expected_original
            and xclass==classified and xfail==failed
            and total==Counter({"attempted":expected_original,
                                 "classified":expected_classified,
                                 "failed":expected_acq_failed}),
            "Five independently sealed original photo groups do not account for all collector rows")

    # Byte/content parity across BOTH independently archived terminal paths.
    compare_cols=list(PHOTO_FIELDS)
    compare_cols.sort()
    original=full.loc[:,compare_cols].fillna("").astype(str)
    appended=combined.loc[:,compare_cols].fillna("").astype(str)
    original=original.sort_values("photo_id",key=lambda a:a.map(int)).reset_index(drop=True)
    appended=appended.sort_values("photo_id",key=lambda a:a.map(int)).reset_index(drop=True)
    pd.testing.assert_frame_equal(original,appended,check_dtype=True,check_exact=True)

    # This report describes where actual sampled photos came from, not
    # whether a photo independently contains botanical anthesis or native taxa.
    status=Counter(full["measurement_status"].astype(str))
    by_hue=Counter(full["morph"].astype(str))
    by_lat=defaultdict(lambda:Counter())
    by_month=defaultdict(lambda:Counter())
    latitudes=full["latitude"] if "latitude" in full else [""]*len(full)
    dates=full["observed_on"] if "observed_on" in full else [""]*len(full)
    for morph,latitude,date in zip(full["morph"],latitudes,dates):
        k=str(morph)
        by_lat[lat_band(latitude)][k]+=1
        by_month[month_str(date)][k]+=1
    taxa=full["inat_taxon_id"].fillna("").astype(str)
    labelled_taxa=full.loc[full["morph"].isin(HUES),"inat_taxon_id"].fillna("").astype(str)
    descriptive={
        "schema":"mois_broad_all16_4468_actual_original_photo_hue_and_source_coverage_v1",
        "actual_original_photo_attempts":len(full),
        "distinct_original_observation_ids":len(obs),
        "source_taxon_ids_attempted":taxa[taxa!=""].nunique(),
        "source_taxon_ids_with_automatic_four_hue":labelled_taxa[labelled_taxa!=""].nunique(),
        "automatically_four_hue_classified":classified,
        "automatic_four_hue_distribution":{x:int(by_hue[x]) for x in HUES},
        "technical_mixed_uncertain_photo_ids":int(by_hue["mixed_uncertain"]),
        "technical_acquisition_failed_photo_ids":failed,
        "measurement_status":dict(sorted(status.items())),
        "latitude_band_by_automatic_ROI_morph":{
            band:dict(sorted(counts.items())) for band,counts in sorted(by_lat.items())
        },
        "observed_month_by_automatic_ROI_morph":{
            month:dict(sorted(counts.items())) for month,counts in sorted(by_month.items())
        },
        "unknown_photo_colour_is_not_nonwhite":True,
        "flowering_anthesis_confirmed":False,
        "individual_wild_native_flower_photo_identified":False,
        "source_arms_mixed_probability_not_worldwide_white_rate":True,
    }
    record={
        "schema":PROTOCOL,
        "original_completed_acquisition_run":current_run,
        "original_photo_attempts":expected_original,
        "genuine_four_hue_classifications":expected_classified,
        "acquisition_failed":expected_acq_failed,
        "acquired_unclassifiable":expected_unknown,
        "same_source_five_original_fcp_worker_groups_photo_by_photo_verified":True,
        "source_original_93m_frame_sha256":source_sha,
        "original_full_collector_receipt_sha256":sha(receipt_path),
        "original_full_collector_parquet_sha256":sha(full_parquet),
        "cross_arm_zero_overlap_registry_sha256":sha(cross_arm),
        "cross_arm_original_photo_attempts_distinct":cross_arm_total,
        "cross_arm_four_hue_classified_distinct":cross_arm_hues,
        "original_image_redownloads_for_this_release":0,
        "colour_roi_reinferences_for_this_release":0,
        "not_a_complete_world_flora_probability_cohort":True,
        "independent_photo_specific_flower_and_native_review_complete":False,
    }
    return record,descriptive


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--finished",type=Path,required=True)
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--cross-arm",type=Path,required=True)
    p.add_argument("--group",type=Path,action="append",required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=False)
    proof,coverage=validate_collector_against_append_groups(
        a.finished,a.source,a.group,a.cross_arm)
    (a.output/"complete_original_photo_release_evidence.json").write_text(
        json.dumps(proof,indent=2,sort_keys=True)+"\n")
    (a.output/"source_taxa_season_latitude_coverage.json").write_text(
        json.dumps(coverage,indent=2,sort_keys=True)+"\n")
    print("MOIS_FULL_BROAD_TIER06_REAL_ORIGINAL_PIXELS_RECOVERED_DATA_ONLY",
          json.dumps({"attempts":proof["original_photo_attempts"],
                      "classifiable":proof["genuine_four_hue_classifications"],
                      "unique_photo_ids_cross_arms":proof["cross_arm_original_photo_attempts_distinct"],
                      "taxon_ids_with_four_hue":coverage["source_taxon_ids_with_automatic_four_hue"]},
                     sort_keys=True),flush=True)


if __name__=="__main__":
    main()
