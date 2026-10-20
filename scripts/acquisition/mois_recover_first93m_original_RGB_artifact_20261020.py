#!/usr/bin/env python3
"""Restore COMPLETED p=.0002 MOIS original RGB photo artifacts after only a shell release error.

The 19 real-photo worker groups already completed and an original measured
Parquet/receipt/source manifest were uploaded in GitHub Actions run
39232152549. Do not download or remeasure any photo, change a photo hue,
treat source metadata as colour, or count FCP-overlapping photos again here.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
import pandas as pd

from scripts.analysis.mois_actual_photo_cohort_season_audit_20261019 import (
    HUES, require, sha256,
)

RUN_ID=39232152549
RUN_HEAD_SHA="95ff0248947c73b17cd2ff6fc5a5e6f47c038831"
PHOTO_COUNT=18462
CLASSIFIED_COUNT=7082
DOWNLOADED_COUNT=18436
ACQUISITION_FAILED=26
UNKNOWN_ACQUIRED=11354
FIRST_SOURCE_PHOTOS=18533
ALREADY_MEASURED=71
UPSTREAM_PHOTOS=93209351
GROUPS=19
SCHEMA="mois_global93m_observation_baseline_real_RGB_photo_results_v1"
SOURCE_SCHEMA="mois_full_source_observation_sampling_cohort_v1"
ARTIFACT_NAME="mois-93m-observation-baseline-genuine-original-RGB"


def audit(jobs_path, measured_path, receipt_path, source_path):
    jobs=json.loads(Path(jobs_path).read_text(encoding="utf-8"))
    receipt=json.loads(Path(receipt_path).read_text(encoding="utf-8"))
    source=json.loads(Path(source_path).read_text(encoding="utf-8"))
    names={}
    for j in jobs.get("jobs", []):
        m=re.fullmatch(r"measure-original-RGB \((\d+)\)", str(j.get("name","")))
        if m:
            index=int(m.group(1))
            require(index not in names,"Duplicate original ROI worker group in GitHub run")
            names[index]=j.get("conclusion")
    require(names=={i:"success" for i in range(GROUPS)},
            "Not all 19 actual original-photo ROI worker groups succeeded")
    require(any(j.get("name")=="freeze-source-queue" and
                j.get("conclusion")=="success" for j in jobs["jobs"]),
            "Original full-index observation random-source queue did not complete")
    require(any(j.get("name")=="persist-authenticated-RGB" and
                j.get("conclusion")=="failure" for j in jobs["jobs"]),
            "This pinned recovery is only for the original post-ROI shell syntax failure")
    require(receipt.get("schema")==SCHEMA and source.get("schema")==SOURCE_SCHEMA,
            "Original-photo result or source-only manifest has changed schemas")
    require(sha256(measured_path)==receipt.get("measured_photo_parquet_sha256")
            and sha256(source_path)==receipt.get("source_plan_sha256"),
            "Measured source pixels or pre-colour photo hash roster differs from original receipt")
    require(receipt.get("attempted")==PHOTO_COUNT
            and receipt.get("classifiable")==CLASSIFIED_COUNT
            and receipt.get("downloaded")==DOWNLOADED_COUNT
            and receipt.get("acquisition_failed")==ACQUISITION_FAILED
            and receipt.get("unclassifiable_acquired")==UNKNOWN_ACQUIRED
            and receipt.get("full_original_upstream_photo_candidates")==UPSTREAM_PHOTOS
            and receipt.get("matching_prior_measured_photos_reused")==ALREADY_MEASURED
            and receipt.get("observation_level_selection_probability")==0.0002,
            "Frozen original first full93m photo-image terminal totals changed")
    require(CLASSIFIED_COUNT+UNKNOWN_ACQUIRED+ACQUISITION_FAILED==PHOTO_COUNT
            and DOWNLOADED_COUNT+ACQUISITION_FAILED==PHOTO_COUNT,
            "Photo classification, acquisition and downloaded denominators disagree")
    require(source.get("selection_probability")==0.0002
            and source.get("selection_unit")=="observation_uuid"
            and source.get("sampled_photo_ids_before_prior_dedup")==FIRST_SOURCE_PHOTOS
            and source.get("previously_measured_photo_ids_skipped")==ALREADY_MEASURED
            and source.get("new_photo_attempts_planned")==PHOTO_COUNT
            and source.get("selected_new_photos_unsupported_by_frozen_loader")==0
            and source.get("full_frame_photo_ids")==UPSTREAM_PHOTOS
            and source.get("source_original_images_measured") is False
            and source.get("colour_measured") is False
            and source.get("no_per_taxon_or_observation_photo_count_cap") is True
            and source.get("all_unselected_source_photos_eligible_at_higher_probability") is True,
            "Not the frozen before-colour full93m source observation draw")
    require(receipt.get("technical_unknowns_not_assumed_nonwhite") is True
            and receipt.get("ecological_analysis_run") is False
            and receipt.get("raw_licensed_photo_pixels_redistributed") is False,
            "Source indicates fabricated hues, ecological analysis or raw images")
    method=receipt.get("frozen_original_ROI_method",{})
    require(isinstance(method, dict) and bool(method.get("source_commit"))
            and bool(method.get("detector_sha256")),
            "Original SHA-pinned image classifier method provenance missing")

    photos=pd.read_parquet(measured_path)
    mandatory={"photo_id","measurement_id","inat_taxon_id","observation_id",
               "morph","measurement_status","classified_four_hue","source_arm"}
    require(mandatory<=set(photos), "Missing frozen source photo ROI identity/status fields")
    require(len(photos)==PHOTO_COUNT and photos.photo_id.is_unique
            and photos.measurement_id.is_unique
            and photos.photo_id.astype(str).str.fullmatch(r"\d+").all(),
            "Actual original photo IDs duplicated or incomplete")
    require(photos["source_arm"].eq("global93m_probability_sample_observations").all(),
            "Photo not from original selected observation probability arm")
    require(photos.observation_id.astype(str).str.len().gt(0).all()
            and photos.inat_taxon_id.astype(str).str.len().gt(0).all(),
            "Original source observation or taxon ID missing")
    known=photos.morph.isin(HUES)
    stated=photos.classified_four_hue.astype(str).str.casefold().eq("true")
    failed=photos.measurement_status.eq("image_acquisition_failed")
    require(known.eq(stated).all() and
            known.eq(photos.measurement_status.eq("classified_four_state_morph")).all()
            and int(known.sum())==CLASSIFIED_COUNT
            and int(failed.sum())==ACQUISITION_FAILED,
            "Original ROI white/nonwhite and unknown status changed after worker completion")
    require((~known & ~failed).sum()==UNKNOWN_ACQUIRED,
            "Unclassifiable but downloaded actual original photos were not retained")
    require(photos.loc[~known,"morph"].eq("mixed_uncertain").all(),
            "Unclassified real RGB photos must NEVER be called chromatic")
    return {
        "schema":"mois_first93m_p0002_completed_original_RGB_release_recovery_v1",
        "recovered_original_github_actions_run":RUN_ID,
        "original_run_head_sha":RUN_HEAD_SHA,
        "original_worker_groups_succeeded":len(names),
        "original_terminal_measurement_artifact_name":ARTIFACT_NAME,
        "original_actual_photo_attempts":len(photos),
        "original_classified_four_hue_photos":int(known.sum()),
        "original_downloaded_photo_attempts":DOWNLOADED_COUNT,
        "original_acquisition_failed_photo_attempts":int(failed.sum()),
        "original_unknown_acquired_photo_attempts":UNKNOWN_ACQUIRED,
        "source_selected_photo_ids_prior_dedup":FIRST_SOURCE_PHOTOS,
        "source_reused_original_measured_photo_ids":ALREADY_MEASURED,
        "source_photo_candidate_count_not_rgb_attempts":UPSTREAM_PHOTOS,
        "original_measurement_parquet_sha256":sha256(measured_path),
        "original_measurement_receipt_sha256":sha256(receipt_path),
        "original_precolour_source_sha256":sha256(source_path),
        "frozen_roi_method":method,
        "new_photo_original_taxa_all_attempts":int(photos.inat_taxon_id.nunique()),
        "new_photo_distinct_observations":int(photos.observation_id.nunique()),
        "new_photo_four_hue_counts":{h:int(photos.morph.eq(h).sum()) for h in HUES},
        "all_19_real_original_image_ROI_groups_succeeded":True,
        "failure_only_last_release_shell_fi_missing":True,
        "photo_pixels_not_refetched_or_reclassified":True,
        "native_and_climate_analysis_deferred":True,
        "fcp_original_photo_overlap_with_this_new_cohort_not_evaluated_here":True
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("jobs","measured","receipt","source","report"):
        p.add_argument("--"+name,type=Path,required=True)
    x=p.parse_args()
    record=audit(x.jobs,x.measured,x.receipt,x.source)
    require(not x.report.exists(),"Do not rewrite signed recovery evidence")
    x.report.parent.mkdir(parents=True,exist_ok=True)
    x.report.write_text(json.dumps(record,ensure_ascii=False,indent=2,
                                  sort_keys=True,allow_nan=False)+"\n")
    print("MOIS_COMPLETE_REAL18462_FROM_ALREADY_FINISHED_GITHUB_ARTIFACT",
          json.dumps({k:record[k] for k in (
              "original_actual_photo_attempts","original_classified_four_hue_photos",
              "new_photo_distinct_observations","new_photo_original_taxa_all_attempts",
              "original_acquisition_failed_photo_attempts")},sort_keys=True),
          flush=True)


if __name__=="__main__":
    main()
