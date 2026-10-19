"""Synthetic source-identity checks for append of newly measured 93m ROI photos.

Synthetic classifiable white labels here are tests ONLY, never field results.
"""
import json

import pandas as pd
import pytest

from scripts.analysis.mois_actual_photo_cohort_season_audit_20261019 import (
    load_pixels, sha256,
)
from scripts.analysis import mois_append_full93m_first_original_RGB_photo_wave_20261019 as mod


def row(pid, tax, obs, morph, status="classified_four_state_morph", source_arm=None):
    d={
        "photo_id":str(pid),"inat_taxon_id":str(tax),
        "species":"Campanula punctata" if tax==10 else "Cirsium brevicaule",
        "observation_id":obs,"latitude":37.1,"longitude":140.7,
        "observed_on":"2025-07-14","morph":morph,"measurement_status":status
    }
    if source_arm is not None:
        d["source_arm"]=source_arm
        d["classified_four_hue"]=str(morph in mod.HUES).lower()
    return d


def setup(tmp_path,monkeypatch):
    old_path=tmp_path/"old.parquet"
    new_path=tmp_path/"new.parquet"
    receipt_path=tmp_path/"receipt.json"
    source_path=tmp_path/"source.json"
    old=[
        row(1,10,"oldA","white"),
        row(2,11,"oldB","red_pink")
    ]
    new=[
        row(3,10,"oldA","blue_purple",
            source_arm="global93m_probability_sample_observations"),
        row(4,20,"newB","mixed_uncertain","image_acquisition_failed",
            source_arm="global93m_probability_sample_observations"),
    ]
    pd.DataFrame(old).to_parquet(old_path,index=False)
    pd.DataFrame(new).to_parquet(new_path,index=False)
    manifest={
        "schema":mod.EXPECTED_QUEUE_SCHEMA,
        "full_frame_photo_ids":93209351,
        "selection_unit":"observation_uuid",
        "selection_probability":0.0002,
        "sampled_photo_ids_before_prior_dedup":3,
        "previously_measured_photo_ids_skipped":1,
        "new_photo_attempts_planned":2,
        "selected_new_photos_unsupported_by_frozen_loader":0
    }
    source_path.write_text(json.dumps(manifest))
    receipt={
        "schema":mod.EXPECTED_RESULT_SCHEMA,
        "measured_photo_parquet_sha256":sha256(new_path),
        "source_plan_sha256":sha256(source_path),
        "observation_level_selection_probability":0.0002,
        "attempted":2,"classifiable":1,
        "acquisition_failed":1,"unclassifiable_acquired":0,
        "full_original_upstream_photo_candidates":93209351,
        "matching_prior_measured_photos_reused":1,
        "technical_unknowns_not_assumed_nonwhite":True,
        "ecological_analysis_run":False,
        "downloaded":1,
        "frozen_original_ROI_method":{"blob":"fixed-original-method"}
    }
    receipt_path.write_text(json.dumps(receipt))
    old_frame=load_pixels(old_path,"core")
    def fake_assemble(*args,**kwargs):
        assert kwargs.get("strict") is True
        return old_frame,[old_frame]
    monkeypatch.setattr(mod,"assemble",fake_assemble)
    return old_path,new_path,receipt_path,source_path


def invoke(args):
    old,new,receipt,plan=args
    return mod.join(old,old,old,old,old,new,receipt,plan,
                    expected_old=2,expected_new=2,
                    expected_candidate=3,expected_reused=1)


def test_exact_new_original_photo_join_and_old_source_immutable(tmp_path,monkeypatch):
    args=setup(tmp_path,monkeypatch)
    result,report=invoke(args)
    assert result.photo_id.tolist()==["1","2","3","4"]
    assert report["new_actual_four_hue_classified"]==1
    assert report["total_actual_four_hue_classified"]==3
    assert report["new_original_photo_attempts_unclassifiable_or_failed"]==1
    assert report["new_nominal_taxa_four_hue_not_previously_classified"]==0
    assert report["new_source_acquisition_failures"]==1
    assert report["preexisting_observation_IDs_in_new_photos_for_future_confirmation"]==1
    assert report["new_nominal_taxa_not_previously_photographed"]==1
    assert not result.duplicated("photo_id").any()
    assert result.loc[result.photo_id.eq("4"),"classified"].iloc[0] == False
    assert result.loc[result.photo_id.eq("4"),"white"].iloc[0] == False
    assert result.loc[result.photo_id.eq("4"),"technical_failure"].iloc[0] == True


def test_previously_measured_original_photo_id_is_rejected(tmp_path,monkeypatch):
    old,new,receipt,plan=setup(tmp_path,monkeypatch)
    d=pd.read_parquet(new);d.loc[0,"photo_id"]="1"
    d.to_parquet(new,index=False)
    r=json.loads(receipt.read_text());r["measured_photo_parquet_sha256"]=sha256(new)
    receipt.write_text(json.dumps(r))
    with pytest.raises(ValueError,match="Previously processed original photo ID"):
        invoke((old,new,receipt,plan))


def test_changed_parquet_checksum_cannot_be_published(tmp_path,monkeypatch):
    old,new,receipt,plan=setup(tmp_path,monkeypatch)
    d=pd.read_parquet(new);d.loc[0,"morph"]="white"
    d.to_parquet(new,index=False)
    with pytest.raises(ValueError,match="SHA-pinned selection plan changed"):
        invoke((old,new,receipt,plan))


def test_zero_technical_missingness_must_match_real_image_receipt(tmp_path,monkeypatch):
    old,new,receipt,plan=setup(tmp_path,monkeypatch)
    r=json.loads(receipt.read_text());r["acquisition_failed"]=0
    r["unclassifiable_acquired"]=1
    receipt.write_text(json.dumps(r))
    with pytest.raises(ValueError,match="technical"):
        invoke((old,new,receipt,plan))


def test_source_plan_sample_size_and_probability_are_not_flexible(tmp_path,monkeypatch):
    old,new,receipt,plan=setup(tmp_path,monkeypatch)
    d=json.loads(plan.read_text());d["selection_probability"]=0.1
    plan.write_text(json.dumps(d))
    r=json.loads(receipt.read_text());r["source_plan_sha256"]=sha256(plan)
    receipt.write_text(json.dumps(r))
    with pytest.raises(ValueError,match="first observation-level p=0.0002"):
        invoke((old,new,receipt,plan))
