"""Contracts for publication recovery only; toy data contain no botanical evidence."""
import json
import pandas as pd
import pytest
from scripts.acquisition import mois_recover_first93m_original_RGB_artifact_20261020 as m
from scripts.analysis.mois_actual_photo_cohort_season_audit_20261019 import sha256


def make(tmp_path, monkeypatch):
    for key, value in {
        "PHOTO_COUNT":3,"CLASSIFIED_COUNT":2,
        "DOWNLOADED_COUNT":2,"ACQUISITION_FAILED":1,
        "UNKNOWN_ACQUIRED":0,"FIRST_SOURCE_PHOTOS":4,
        "ALREADY_MEASURED":1,"GROUPS":2
    }.items():
        monkeypatch.setattr(m,key,value)
    parquet=tmp_path/"measured_photos.parquet"
    rows=[
        dict(photo_id="100",measurement_id="one",inat_taxon_id="51",
             observation_id="oA",morph="white",
             measurement_status="classified_four_state_morph",
             classified_four_hue="true",
             source_arm="global93m_probability_sample_observations"),
        dict(photo_id="101",measurement_id="two",inat_taxon_id="52",
             observation_id="oA",morph="blue_purple",
             measurement_status="classified_four_state_morph",
             classified_four_hue="true",
             source_arm="global93m_probability_sample_observations"),
        dict(photo_id="102",measurement_id="three",inat_taxon_id="52",
             observation_id="oB",morph="mixed_uncertain",
             measurement_status="image_acquisition_failed",
             classified_four_hue="false",
             source_arm="global93m_probability_sample_observations")
    ]
    pd.DataFrame(rows).to_parquet(parquet,index=False)
    source=tmp_path/"source.json"
    source.write_text(json.dumps({
        "schema":m.SOURCE_SCHEMA,"selection_probability":.0002,
        "selection_unit":"observation_uuid",
        "sampled_photo_ids_before_prior_dedup":4,
        "previously_measured_photo_ids_skipped":1,
        "new_photo_attempts_planned":3,
        "selected_new_photos_unsupported_by_frozen_loader":0,
        "full_frame_photo_ids":m.UPSTREAM_PHOTOS,
        "source_original_images_measured":False,
        "colour_measured":False,
        "no_per_taxon_or_observation_photo_count_cap":True,
        "all_unselected_source_photos_eligible_at_higher_probability":True,
    }))
    receipt=tmp_path/"receipt.json"
    receipt.write_text(json.dumps({
        "schema":m.SCHEMA,"attempted":3,
        "classifiable":2,"downloaded":2,
        "acquisition_failed":1,"unclassifiable_acquired":0,
        "full_original_upstream_photo_candidates":m.UPSTREAM_PHOTOS,
        "matching_prior_measured_photos_reused":1,
        "observation_level_selection_probability":.0002,
        "measured_photo_parquet_sha256":sha256(parquet),
        "source_plan_sha256":sha256(source),
        "technical_unknowns_not_assumed_nonwhite":True,
        "ecological_analysis_run":False,
        "raw_licensed_photo_pixels_redistributed":False,
        "frozen_original_ROI_method":{
            "source_commit":"same_model","detector_sha256":"abc"
        }
    }))
    jobs=tmp_path/"jobs.json"
    jobs.write_text(json.dumps({"jobs":[
        {"name":"freeze-source-queue","conclusion":"success"},
        {"name":"measure-original-RGB (0)","conclusion":"success"},
        {"name":"measure-original-RGB (1)","conclusion":"success"},
        {"name":"persist-authenticated-RGB","conclusion":"failure"}
    ]}))
    return jobs,parquet,receipt,source


def verify(a):
    return m.audit(*a)


def test_signed_workers_finished_only_release_shell_failed(tmp_path, monkeypatch):
    a=make(tmp_path, monkeypatch)
    r=verify(a)
    assert r["original_actual_photo_attempts"]==3
    assert r["original_classified_four_hue_photos"]==2
    assert r["original_acquisition_failed_photo_attempts"]==1
    assert r["new_photo_distinct_observations"]==2
    assert r["new_photo_four_hue_counts"]["white"]==1
    assert r["photo_pixels_not_refetched_or_reclassified"]


def test_missing_successful_worker_group_rejects_resume(tmp_path, monkeypatch):
    a=make(tmp_path,monkeypatch)
    d=json.loads(a[0].read_text())
    d["jobs"][2]["conclusion"]="failure"
    a[0].write_text(json.dumps(d))
    with pytest.raises(ValueError,match="19 actual original-photo ROI worker"):
        verify(a)


def test_mutated_measured_original_photo_pixel_rejected_by_byte_hash(tmp_path,monkeypatch):
    a=make(tmp_path,monkeypatch)
    d=pd.read_parquet(a[1])
    d.loc[0,"morph"]="red_pink"
    d.to_parquet(a[1],index=False)
    with pytest.raises(ValueError,match="Measured source pixels"):
        verify(a)


def test_unclassifiable_photo_cannot_be_mislabeled_as_nonwhite(tmp_path,monkeypatch):
    a=make(tmp_path,monkeypatch)
    d=pd.read_parquet(a[1])
    d.loc[2,"morph"]="red_pink"
    d.to_parquet(a[1],index=False)
    r=json.loads(a[2].read_text())
    r["measured_photo_parquet_sha256"]=sha256(a[1])
    a[2].write_text(json.dumps(r))
    with pytest.raises(ValueError,match="Original ROI white/nonwhite"):
        verify(a)


def test_duplicate_original_source_photo_ids_fail_even_if_signed(tmp_path,monkeypatch):
    a=make(tmp_path,monkeypatch)
    d=pd.read_parquet(a[1])
    d.loc[2,"photo_id"]="100"
    d.to_parquet(a[1],index=False)
    r=json.loads(a[2].read_text())
    r["measured_photo_parquet_sha256"]=sha256(a[1])
    a[2].write_text(json.dumps(r))
    with pytest.raises(ValueError,match="duplicated or incomplete"):
        verify(a)
