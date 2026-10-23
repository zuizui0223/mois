#!/usr/bin/env python3
"""Cross-arm ledger for ONLY genuinely completed tier06 wave2 original-photo groups.

The signed 908494-photo prior registry already proves that original p0032
pixel delta and original tier06 group00/R1 were disjoint from its 887600-photo
historical source. Wave2 was independently selected BEFORE colour from new
observation hash buckets4..11, after SHA-validating and excluding every photo
AND observation in a newer 897840-photo historic checkpoint plus the two
independent prior tier06 arms. Verify new wave2 original photos again against
both small tier06 arm tables and the full newest p0032 delta before counting.

DO NOT fetch the old historical ~0.9m-photo Parquet, any image bytes, or
re-run frozen ROI-v4. A record is an automatic ROI-hue classification, not
independent confirmation of wild/native/anthesis or global white-flower rate.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path

from scripts.acquisition.mois_completed_pixel_checkpoint import sha, require, HUES

SCHEMA="mois_crossarm_r1_r2_p0032_original_photo_source_arm_incremental_v1"
PRIOR_TOTAL=908_494
PRIOR_HUES=408_997
PRIOR_SCHEMA="mois_original_RGB_crossarm_plus_new_p0032_signed_delta_v1"
PRIOR_RUN=39889436140
R2_RUN=39919846830
R2_SOURCE_SHA="57592ce9f12c8abb0adb7ce2fda8252044ca9204caf295dbac73ecbf76a01516"
R2_EXPECTED_PHOTOS=8765
R2_SCHEMA="mois_tier06_all16_wave2_append_only_original_RGB_group_v1"
PRIOR_TIER06=(4_138,1_573)
PRIOR_BROAD=(4_468,1_760)
PRIOR_P0032_DELTA=(12_288,4_760)
SOURCE_ARM="global93m_tier06_all16_buckets004to011_broad_warm_R2"


def source_identity(path, photo_column="photo_id", observation_column="observation_id"):
    import pandas as pd
    rows=pd.read_parquet(path,columns=[photo_column,observation_column])
    pid=rows[photo_column].fillna("").astype(str)
    obs=rows[observation_column].fillna("").astype(str)
    require(pid.str.fullmatch(r"[0-9]+").all() and pid.is_unique
            and obs.str.strip().ne("").all(),
            "Prior genuine source photo/observation IDs are not unique or nonmissing")
    return set(pid),set(obs)


def validate(prior_registry, source_manifest,
             early_parquet, broad_parquet, p0032_delta_parquet,
             group_folders, output,
             *, prior_total=PRIOR_TOTAL, prior_hues=PRIOR_HUES,
             expected_source_sha=R2_SOURCE_SHA, expected_new=R2_EXPECTED_PHOTOS,
             prior_early=PRIOR_TIER06, prior_broad=PRIOR_BROAD,
             prior_delta=PRIOR_P0032_DELTA):
    import pandas as pd
    output=Path(output)
    require(not output.exists(),"Never mutate completed original photo cross-arm evidence")
    old=json.loads(Path(prior_registry).read_text(encoding="utf-8"))
    source=json.loads(Path(source_manifest).read_text(encoding="utf-8"))
    require(old.get("schema")==PRIOR_SCHEMA
            and old.get("prior_historical_main_checkpoint_run")==PRIOR_RUN
            and old.get("full_authenticated_distinct_original_photo_attempts")==prior_total
            and old.get("full_authenticated_four_hue_original_image_classifications")==prior_hues
            and old.get("photo_ID_overlap_new_p0032_vs_other_original_sampling_arms")==0
            and old.get("observation_ID_overlap_new_p0032_vs_other_original_sampling_arms")==0
            and old.get("separately_source_arm_audited") is True,
            "Prior complete source-arm original photo/observation audit invalid")
    require(sha(early_parquet)==old["signed_old_tier06_first_group_parquet_SHA256"]
            and sha(broad_parquet)==old["signed_old_tier06_all16_first4_parquet_SHA256"]
            and sha(p0032_delta_parquet)==old["new_p0032_original_photo_delta_parquet_SHA256"],
            "Previously SHA-audited cross-arm original photo Parquet changed")
    require(sha(source_manifest)==expected_source_sha
            and source.get("schema")=="mois_full_source_observation_sampling_cohort_v1"
            and source.get("full_frame_photo_ids")==93_209_351
            and source.get("complete_tier06_source_photo_IDs")==597_070
            and source.get("source_group")=="all16"
            and source.get("source_tier")==6
            and source.get("source_arm")==SOURCE_ARM
            and source.get("new_photo_attempts_planned")==expected_new
            and source.get("previous_broad_original_photo_attempts")==prior_broad[0]
            and source.get("previous_broad_source_parquet_sha256")==sha(broad_parquet)
            and source.get("metadata_original_photos_in_admitted_buckets")==9306
            and source.get("previously_measured_photo_ids_skipped")==541
            and source.get("compute_bucket_count")==512
            and source.get("admitted_buckets")==list(range(4,12))
            and source.get("prior_photo_and_observation_exclusion_verified") is True
            and source.get("not_probability_complete_global_source") is True
            and source.get("colour_measured") is False,
            "Wave2 was not selected before hue with prior source-arm IDs excluded")
    require(len(group_folders)>0,"No genuinely sealed original wave2 worker groups")

    early_id,early_obs=source_identity(early_parquet)
    broad_id,broad_obs=source_identity(broad_parquet)
    p32_id,p32_obs=source_identity(p0032_delta_parquet)
    require(len(early_id)==prior_early[0] and len(broad_id)==prior_broad[0]
            and len(p32_id)==prior_delta[0],
            "Prior separately measured source-arm photo count/identity changed")
    require(not (early_id&broad_id or early_id&p32_id or broad_id&p32_id)
            and not (early_obs&broad_obs or early_obs&p32_obs or broad_obs&p32_obs),
            "Previous original sampling arms have an unexpected repeated observation")
    old_photos=early_id|broad_id|p32_id
    old_obs=early_obs|broad_obs|p32_obs

    groups=set();used_units=set();new_photos=set();new_obs=set()
    records=[];previous_per_arm=Counter();counts=Counter()
    all_hues=set((*HUES,"mixed_uncertain"))
    for path in group_folders:
        folder=Path(path)
        r_path=folder/"original_photo_group_receipt.json"
        a=json.loads(r_path.read_text(encoding="utf-8"))
        group=a.get("original_source_warm_group")
        require(type(group) is int and group in source["group_ids"]
                and group not in groups
                and a.get("schema")==R2_SCHEMA
                and a.get("original_acquisition_workflow_run")==R2_RUN
                and a.get("source_original_photo_manifest_sha256")==expected_source_sha
                and a.get("complete_expected_source_photo_attempts")==expected_new
                and a.get("all_source_photo_and_original_terminal_SHA_verified") is True
                and a.get("frozen_original_FCP_model_and_output_verified") is True
                and a.get("image_bytes_reacquired_in_sealing") is False
                and a.get("ROI_or_hue_reinference_in_sealing") is False
                and a.get("globally_probability_complete") is False,
                "An original photo group receipt is not authenticated/unique")
        groups.add(group)
        expected_units=list(range(group*4,min((group+1)*4,source["work_units"])))
        units=a.get("original_work_units")
        require(units==expected_units and not used_units.intersection(units),
                "Original photo shards are missing/repeated across worker groups")
        used_units.update(units)
        f=folder/"original_photo_ROIs.parquet"
        receipt_hashes=a.get("saved_file_sha256",{})
        cov=folder/"original_photo_source_coverage.json"
        require(sha(f)==receipt_hashes.get(f.name)
                and sha(cov)==receipt_hashes.get(cov.name),
                "Original wave2 photo/coverage payload byte SHA changed")
        coverage=json.loads(cov.read_text(encoding="utf-8"))
        require(coverage.get("schema")==
                  "mois_wave2_real_photo_source_taxonomy_time_latitude_coverage_v1"
                and coverage.get("source_taxon_ids_are_not_species_validated") is True
                and coverage.get("individual_photo_open_flowers_independently_verified") is False,
                "Original photo source taxonomy/season QA missing")
        columns=["photo_id","observation_id","morph","measurement_status",
                 "selection_probability","source_arm","inat_taxon_id"]
        tab=pd.read_parquet(f,columns=columns)
        ids=tab.photo_id.astype(str)
        obs=tab.observation_id.fillna("").astype(str)
        labelled=tab.morph.isin(HUES)
        require(ids.is_unique and ids.str.fullmatch(r"[0-9]+").all()
                and obs.str.strip().ne("").all()
                and not new_photos.intersection(ids)
                and not old_photos.intersection(ids)
                and not old_obs.intersection(obs)
                and tab.morph.isin(all_hues).all()
                and labelled.eq(tab.measurement_status.eq("classified_four_state_morph")).all()
                and tab.source_arm.eq(SOURCE_ARM).all()
                and pd.to_numeric(tab.selection_probability,errors="coerce").eq(
                    source["selection_probability"]).all(),
                "Wave2 original photo/observation collision, propensity or frozen hue changed")
        require(len(tab)==a.get("distinct_genuinely_completed_photo_ids")
                and int(labelled.sum())==a.get("four_hue_model_classified")
                and len(set(obs))==a.get("distinct_original_observation_ids")
                and int(coverage["original_photo_attempts"])==len(tab)
                and int(coverage["distinct_nominal_inat_taxon_ids_four_hue"])==
                    tab.loc[labelled,"inat_taxon_id"].fillna("").astype(str).replace("",pd.NA).dropna().nunique(),
                "Wave2 photo totals, quality/coverage and original terminal labels disagree")
        new_photos.update(ids)
        same_arm_obs_spanning=len(new_obs.intersection(obs))
        new_obs.update(obs)
        counts["same_original_observations_across_wave2_compute_groups"]+=same_arm_obs_spanning
        counts["photo_attempts"]+=len(tab)
        counts["classified"]+=int(labelled.sum())
        counts["acquisition_failed"]+=a["acquisition_failed"]
        counts["acquired_unknown"]+=a["acquired_but_unclassifiable"]
        records.append({
            "group":group,
            "original_photo_attempts":len(tab),
            "four_hue_classified":int(labelled.sum()),
            "nominal_source_taxa_classified":coverage["distinct_nominal_inat_taxon_ids_four_hue"],
            "distinct_observations":len(set(obs)),
            "group_source_SHA256":sha(r_path),
            "group_photo_Parquet_SHA256":sha(f),
            "group_coverage_SHA256":sha(cov),
            "original_source_photo_units":units,
        })
    require(counts["photo_attempts"]<=expected_new
            and counts["classified"]+counts["acquisition_failed"]+
                    counts["acquired_unknown"]==counts["photo_attempts"]
            and len(used_units)<=source["work_units"],
            "New real original-photo counts exceed signed precolour source opportunity")
    records.sort(key=lambda x:x["group"])
    result={
        "schema":SCHEMA,
        "old_signed_908494_original_photo_registry_SHA256":sha(prior_registry),
        "wave2_precolour_source_SHA256":sha(source_manifest),
        "prior_base_908494_real_original_photo_attempts":prior_total,
        "prior_base_408997_four_hue_classified":prior_hues,
        "wave2_groups_completed":[x["group"] for x in records],
        "wave2_original_photo_attempts":counts["photo_attempts"],
        "wave2_original_four_hue_classified":counts["classified"],
        "wave2_original_photo_acquisition_failed":counts["acquisition_failed"],
        "wave2_original_photo_acquired_but_unclassifiable":counts["acquired_unknown"],
        "wave2_original_distinct_observations":len(new_obs),
        "wave2_same_original_observations_spanning_compute_groups":
            counts["same_original_observations_across_wave2_compute_groups"],
        "wave2_remaining_original_photo_attempts":expected_new-counts["photo_attempts"],
        "full_distinct_original_photo_attempts":
            prior_total+counts["photo_attempts"],
        "full_automatically_classified_four_hue_original_photo_ROIs":
            prior_hues+counts["classified"],
        "prior_separate_arms_and_new_wave2_no_duplicate_photo_or_observation":True,
        "source_tier06_buckets4to11_disjoint_from_first4_by_design":True,
        "historical_897840_photo_previous_selection_SHA_verified_by_r2_precolour_receipt":True,
        "post_897840_new_p0032_photo_IDs_and_observations_directly_compared":True,
        "all_wave2_group_photo_and_source_coverage_SHA_verified":True,
        "original_images_fetched_and_ROI_inferred_during_this_ledger":0,
        "historical_original_full_photo_Parquets_rewritten":False,
        "photo_hues_not_botanically_verified_as_flowering":True,
        "source_nominal_taxa_not_accepted_species_counts":True,
        "different_source_arm_sampling_probabilities_not_pooled_into_global_white_prevalence":True,
        "per_group_original_photo_source_audit":records,
    }
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print("MOIS_DISTINCT_NEW_ORIGINAL_ROI_CROSS_ARM_WAVE2_LEDGER",
          json.dumps({k:result[k] for k in (
              "wave2_groups_completed","wave2_original_photo_attempts",
              "wave2_original_four_hue_classified",
              "full_distinct_original_photo_attempts",
              "full_automatically_classified_four_hue_original_photo_ROIs")},sort_keys=True),flush=True)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("prior-registry","source-manifest","early-parquet",
                "broad-parquet","p0032-delta-parquet","output"):
        p.add_argument("--"+key,type=Path,required=True)
    p.add_argument("--group",type=Path,action="append",required=True)
    a=p.parse_args()
    validate(a.prior_registry,a.source_manifest,
             a.early_parquet,a.broad_parquet,a.p0032_delta_parquet,
             a.group,a.output)


if __name__=="__main__":
    main()
