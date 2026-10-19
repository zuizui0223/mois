#!/usr/bin/env python3
"""Append a completed authentic first 93m observation-probability RGB wave to MOIS.

Acquisition-only work: no native-status filtering, ecological model, or botanical
flower-petal inference. Every original source photo ID, source recruitment arm,
date and technical colour unknown remains represented. All source-byte hashes
and the exact prior 350,399 photo-ID set must match before publication.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis.mois_actual_photo_cohort_season_audit_20261019 import (
    HUES, assemble, load_pixels, require, sha256, summary,
)

OLD_PHOTOS=350399
FIRST_NEW=18462
FIRST_ATTEMPTED_SOURCE=18533
FIRST_PREVIOUSLY_MEASURED=71
EXPECTED_RESULT_SCHEMA="mois_global93m_observation_baseline_real_RGB_photo_results_v1"
EXPECTED_QUEUE_SCHEMA="mois_full_source_observation_sampling_cohort_v1"
OUTPUT_PHOTOS="mois368861_authentic_original_photo_ROI_source_union.parquet"
OUTPUT_REPORT="mois368861_original_real_photo_first93m_rgb_append_receipt.json"
OUTPUT_COHORT="mois368861_old_and_93m_new_photo_colour_coverage.csv"
COLUMNS=(
    "photo_id","observation_id","inat_taxon_id","species","latitude",
    "longitude","observed_on","morph","classified","white","arm",
    "measurement_status","technical_failure",
)


def source_photo_set(frame):
    require(frame.photo_id.is_unique, "Repeated original real-photo IDs")
    return set(frame.photo_id)


def normalise_source(frame):
    """A stable image-result view; raw source-Parquet Releases remain immutable."""
    z=frame.copy()
    if "measurement_status" not in z:
        z["measurement_status"]=""
    z["measurement_status"]=z["measurement_status"].fillna("").astype(str)
    z["technical_failure"]=z.measurement_status.eq("image_acquisition_failed")
    require(z.classified.eq(z.morph.isin(HUES)).all(),
            "Photo classifiability diverged from original four-hue ROI labels")
    for key in COLUMNS:
        require(key in z,"Original photo column missing: "+key)
    return z.loc[:,COLUMNS].copy()


def join(core,broad,wave2,cached8192,tail16384,
         new_parquet,new_receipt,source_manifest,
         *,expected_old=OLD_PHOTOS,expected_new=FIRST_NEW,
         expected_candidate=FIRST_ATTEMPTED_SOURCE,
         expected_reused=FIRST_PREVIOUSLY_MEASURED):
    old,arms=assemble(core,broad,wave2,cached8192,tail16384,strict=True)
    require(len(old)==expected_old,"Old authentic pixel photo baseline changed")
    receipt=json.loads(Path(new_receipt).read_text(encoding="utf-8"))
    plan=json.loads(Path(source_manifest).read_text(encoding="utf-8"))
    require(receipt["schema"]==EXPECTED_RESULT_SCHEMA
            and plan["schema"]==EXPECTED_QUEUE_SCHEMA,
            "New first 93m cohort must be source-signed genuine RGB photos")
    require(sha256(new_parquet)==receipt["measured_photo_parquet_sha256"]
            and sha256(source_manifest)==receipt["source_plan_sha256"],
            "Original 93m photo ROI file or SHA-pinned selection plan changed")
    require(plan["full_frame_photo_ids"]==93209351
            and plan["selection_unit"]=="observation_uuid"
            and plan["selection_probability"]==0.0002
            and receipt["observation_level_selection_probability"]==0.0002
            and plan["sampled_photo_ids_before_prior_dedup"]==expected_candidate
            and plan["previously_measured_photo_ids_skipped"]==expected_reused
            and plan["new_photo_attempts_planned"]==expected_new
            and plan["selected_new_photos_unsupported_by_frozen_loader"]==0
            and plan["new_photo_attempts_planned"]+expected_reused==expected_candidate
            and receipt["attempted"]==expected_new,
            "Not the immutable first observation-level p=0.0002 source photo batch")
    require(receipt["full_original_upstream_photo_candidates"]==93209351
            and receipt["matching_prior_measured_photos_reused"]==expected_reused
            and receipt["technical_unknowns_not_assumed_nonwhite"] is True
            and receipt["ecological_analysis_run"] is False,
            "Unverified original image ROI provenance or unsupported biological claim")
    new=load_pixels(new_parquet,"global93m_probability_observation")
    require(len(new)==expected_new
            and new["photo_id"].is_unique
            and (new["source_arm"]=="global93m_probability_sample_observations").all()
            if "source_arm" in new else len(new)==expected_new and new.photo_id.is_unique,
            "Unexpected new source photo cohort or acquisition identity")
    require(not (source_photo_set(new)&source_photo_set(old)),
            "Previously processed original photo ID was re-inferred in the new 93m wave")
    require(int(new.classified.sum())==receipt["classifiable"],
            "Actual new four-hue ROI labels disagree with frozen original method receipt")
    require(int(new.technical_failure.fillna(False).sum())==receipt["acquisition_failed"]
            and int((~new.classified).sum())==
                receipt["acquisition_failed"]+receipt["unclassifiable_acquired"],
            "Original photo acquisition failures and technically unknown ROIs were modified")
    require(new.photo_id.str.fullmatch(r"\d+").all()
            and new.inat_taxon_id.ne("").all(),
            "Original photo ID or source taxon was not preserved")
    old_obs=set(old.loc[old.observation_id.ne(""),"observation_id"])
    new_obs=set(new.loc[new.observation_id.ne(""),"observation_id"])
    old_any_taxa=set(old.inat_taxon_id)-{""}
    new_any_taxa=set(new.inat_taxon_id)-{""}
    old_named_hues=set(old.loc[old.classified,"inat_taxon_id"])-{""}
    new_named_hues=set(new.loc[new.classified,"inat_taxon_id"])-{""}

    raw_old=normalise_source(old)
    raw_new=normalise_source(new)
    union=pd.concat([raw_old,raw_new],ignore_index=True)
    require(len(union)==expected_old+expected_new
            and union.photo_id.is_unique
            and union.iloc[:len(raw_old)].reset_index(drop=True).equals(
                raw_old.reset_index(drop=True)),
            "Old original RGB source outcome rows mutated or image IDs duplicated")
    require(int(union.classified.sum())==int(old.classified.sum())+
            receipt["classifiable"],
            "Photo hue totals were not preserved when appending genuine original RGB")
    report={
        "schema":"mois_genuine_full93m_first_probability_photo_RGB_append_to_350399_v1",
        "photographic_scope":"all original 350399 true pixel attempts + genuinely processed original 93m-source p=.0002 observation-photo subset, NO ecological modelling",
        "full_metadata_source_candidate_original_photo_ids":93209351,
        "preexisting_actual_photo_ROIs":int(len(old)),
        "new_actual_original_photo_attempts":int(len(new)),
        "total_actual_original_photo_attempts":int(len(union)),
        "preexisting_actual_four_hue_classified":int(old.classified.sum()),
        "new_actual_four_hue_classified":int(new.classified.sum()),
        "total_actual_four_hue_classified":int(union.classified.sum()),
        "new_original_photo_attempts_unclassifiable_or_failed":int((~new.classified).sum()),
        "new_source_acquisition_failures":int(receipt["acquisition_failed"]),
        "new_source_acquired_unknown_hue":int(receipt["unclassifiable_acquired"]),
        "new_source_downloaded_photo_attempts":int(receipt["downloaded"]),
        "new_source_reused_prior_photo_ids_not_refetched":expected_reused,
        "old_photo_IDs_altered":False,
        "photo_ID_duplicate_count":0,
        "preexisting_observation_IDs_in_new_photos_for_future_confirmation":
            len(new_obs&old_obs),
        "old_distinct_observation_ids":len(old_obs),
        "new_distinct_observation_ids":len(new_obs),
        "old_nominal_taxa_all_ROIs":len(old_any_taxa),
        "new_nominal_taxa_all_ROIs":len(new_any_taxa),
        "new_nominal_taxa_not_previously_photographed":len(new_any_taxa-old_any_taxa),
        "old_nominal_taxa_four_hue":len(old_named_hues),
        "new_nominal_taxa_four_hue":len(new_named_hues),
        "new_nominal_taxa_four_hue_not_previously_classified":len(new_named_hues-old_named_hues),
        "new_actual_pixel_source_hue_counts":{
            hue:int(new.morph.eq(hue).sum()) for hue in HUES},
        "new_actual_pixel_source_original_unknown_label_counts":{
            str(k):int(v) for k,v in
            new.loc[~new.classified,"morph"].value_counts().items()},
        "source_signed_files_sha256":{
            "first_93m_real_photo_parquet":sha256(new_parquet),
            "first_93m_measurement_receipt":sha256(new_receipt),
            "full93m_source_hash_selection_manifest":sha256(source_manifest),
            "old_original_pixels_core":sha256(core),
            "old_original_pixels_broad":sha256(broad),
            "old_original_pixels_wave2":sha256(wave2),
            "old_original_pixels_cached8192":sha256(cached8192),
            "old_original_pixels_tail16384":sha256(tail16384)
        },
        "frozen_method":receipt["frozen_original_ROI_method"],
        "photo_equal_four_hue_fraction_not_census":{
            "old":float(old.white.mean()) if len(old) else None,
            "new_classified_only":float(new.loc[new.classified,"white"].mean())
                if new.classified.any() else None,
            "full_classified_only":float(union.loc[union.classified,"white"].mean())
                if union.classified.any() else None
        },
        "still_no_independent_botanical_flower_open_petal_confirmation":True,
        "no_native_island_or_environmental_effect_estimated":True,
        "historic_source_arm_probabilities_not_equivalent":True,
        "all_original_licensed_image_bytes_remain_off_release":True,
    }
    report["photo_equal_four_hue_fraction_not_census"]["old"]=(
        float(old.loc[old.classified,"white"].mean()) if old.classified.any() else None)
    return union,report


def run(core,broad,wave2,cached8192,tail16384,new_parquet,new_receipt,
        source_manifest,output):
    output=Path(output)
    require(not output.exists() or not any(output.iterdir()),
            "Do not rewrite prior source-signed authentic image results")
    union,report=join(core,broad,wave2,cached8192,tail16384,
                      new_parquet,new_receipt,source_manifest)
    output.mkdir(parents=True,exist_ok=True)
    union.to_parquet(output/OUTPUT_PHOTOS,index=False,compression="zstd")
    frames=[]
    for name,x in union.groupby("arm",observed=True,sort=False):
        classified=x.loc[x.classified]
        frames.append({
            "source_arm":name,
            "original_photo_attempts":len(x),
            "four_hue_classified":len(classified),
            "technical_unknown_ROIs":int((~x.classified).sum()),
            "acquisition_failures_with_source_record":int(x.technical_failure.sum()),
            "white":int(classified.morph.eq("white").sum()),
            "yellow_orange":int(classified.morph.eq("yellow_orange").sum()),
            "red_pink":int(classified.morph.eq("red_pink").sum()),
            "blue_purple":int(classified.morph.eq("blue_purple").sum()),
            "classifiable_fraction":float(x.classified.mean()),
            "white_fraction_among_four_hue":(
                float(classified.white.mean()) if len(classified) else None),
            "nominal_taxa_with_classifiable_ROIs":int(
                classified.inat_taxon_id.loc[
                    classified.inat_taxon_id.ne("")].nunique()),
            "distinct_source_observation_ids":int(
                x.observation_id.loc[x.observation_id.ne("")].nunique())
        })
    pd.DataFrame(frames).to_csv(output/OUTPUT_COHORT,index=False)
    report["output_sha256"]={
        OUTPUT_PHOTOS:sha256(output/OUTPUT_PHOTOS),
        OUTPUT_COHORT:sha256(output/OUTPUT_COHORT),
    }
    (output/OUTPUT_REPORT).write_text(json.dumps(
        report,ensure_ascii=False,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("MOIS_FULL93M_FIRST_REAL_RGB_ADDED_TO_350399",json.dumps({
        k:report[k] for k in (
            "total_actual_original_photo_attempts",
            "new_actual_original_photo_attempts",
            "new_actual_four_hue_classified",
            "total_actual_four_hue_classified",
            "new_original_photo_attempts_unclassifiable_or_failed",
            "new_nominal_taxa_four_hue_not_previously_classified",
            "preexisting_observation_IDs_in_new_photos_for_future_confirmation")
    },sort_keys=True),flush=True)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("core","broad","wave2","cached8192","tail16384",
                "new-parquet","new-receipt","source-manifest","output"):
        p.add_argument("--"+key,type=Path,required=True)
    a=p.parse_args()
    run(a.core,a.broad,a.wave2,a.cached8192,a.tail16384,
        a.new_parquet,a.new_receipt,a.source_manifest,a.output)


if __name__=="__main__":
    main()
