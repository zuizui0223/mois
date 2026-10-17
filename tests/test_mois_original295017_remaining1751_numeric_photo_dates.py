"""Tests of strict source photo+numeric observation recovery and no date overwrite."""
from __future__ import annotations
import pandas as pd
import pytest

from scripts.acquisition import mois_original295017_remaining1751_numeric_photo_dates as m


def fake_frames(monkeypatch):
    monkeypatch.setattr(m,"ORIGINAL_PHOTO_COHORT",4)
    monkeypatch.setattr(m,"NUMERIC_SOURCE_MISSING",2)
    monkeypatch.setattr(m,"UUID_NONNUMERIC_MISSING",1)
    d=pd.DataFrame([
        {"photo_id":"100","observation_id":"1100","inat_taxon_id":"44","morph":"white",
         "observed_on":"","qualified_date":"","old_mois_or_new_fcp":"historical",
         "original_source_date_provenance":"unknown_no_observed_on"},
        {"photo_id":"101","observation_id":"1101","inat_taxon_id":"44","morph":"mixed_uncertain",
         "observed_on":"","qualified_date":"","old_mois_or_new_fcp":"fcp_additional_ROI",
         "original_source_date_provenance":"unknown_no_observed_on"},
        {"photo_id":"102","observation_id":"uuid-alpha","inat_taxon_id":"44","morph":"white",
         "observed_on":"","qualified_date":"","old_mois_or_new_fcp":"mois_new_original_ROI",
         "original_source_date_provenance":"unknown_no_observed_on"},
        {"photo_id":"103","observation_id":"1103","inat_taxon_id":"45","morph":"red_pink",
         "observed_on":"","qualified_date":"2024-02-21","old_mois_or_new_fcp":"historical",
         "original_source_date_provenance":"historic_s3"}
    ])
    selected=d.iloc[:2].copy()
    found=pd.DataFrame([
        {"original_photo_id":"100","original_observation_id":"1100","source_taxon_id":"44",
         "provenance_status":"source_photo_and_observation_verified_with_current_public_date",
         "retrieved_observed_on":"2026-07-18"},
        {"original_photo_id":"101","original_observation_id":"1101","source_taxon_id":"44",
         "provenance_status":"original_photo_verified_but_current_observation_date_unavailable_or_invalid",
         "retrieved_observed_on":None},
    ])
    return d,selected,found


def test_new_public_original_numeric_photo_date_only_and_preserve_existing(monkeypatch):
    a,s,found=fake_frames(monkeypatch)
    q,counts=m.merge(a,s,found)
    assert counts["newly_verified_source_photo_numeric_original_obs_dates"]==1
    assert counts["newly_dated_classifiable_four_hue_original_source_photos"]==1
    assert q.loc[0,"qualified_date"]=="2026-07-18"
    assert q.loc[1,"qualified_date"]==""
    assert q.loc[2,"qualified_date"]==""
    assert q.loc[3,"qualified_date"]=="2024-02-21"
    assert q[["photo_id","morph","inat_taxon_id"]].equals(
        a[["photo_id","morph","inat_taxon_id"]])


def test_source_photo_identity_change_must_fail(monkeypatch):
    a,s,found=fake_frames(monkeypatch)
    changed=found.copy()
    changed.loc[0,"original_observation_id"]="other"
    with pytest.raises(ValueError,match="numeric observation"):
        m.merge(a,s,changed)
    changed=found.copy()
    changed.loc[0,"original_photo_id"]="999999"
    with pytest.raises(ValueError,match="source ID"):
        m.merge(a,s,changed)


def test_existing_public_day_never_overwritten(monkeypatch):
    a,s,found=fake_frames(monkeypatch)
    s=s.copy()
    s.loc[1,["photo_id","observation_id","inat_taxon_id"]]=["103","1103","45"]
    found=found.copy()
    found.loc[1,"original_photo_id"]="103"
    found.loc[1,"original_observation_id"]="1103"
    found.loc[1,"source_taxon_id"]="45"
    found.loc[1,"provenance_status"]="source_photo_and_observation_verified_with_current_public_date"
    found.loc[1,"retrieved_observed_on"]="2025-02-21"
    with pytest.raises(ValueError,match="not be overwritten"):
        m.merge(a,s,found)
