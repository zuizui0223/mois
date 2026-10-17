"""The only permitted changes to original photos are 1732/4096 sourced days.

Synthetic fixtures test source identity; CI separately checks actual numeric
295k/FCP218/full3430502 source SHA receipts and measured outcomes.
"""
from __future__ import annotations
import pandas as pd
import pytest

from scripts.analysis import mois_reconcile_original302242_all_three_verified_photo_date_sources as m


def mock_actual_source(monkeypatch):
    for attr,value in {
        "PIXEL_N":5,"CLASSIFIABLE_N":4,"BASE_DATED_ORIG":2,
        "BASE_DATED_HUES":2,"OLD_N":3,"API_NEW_DATES":1,
        "API_NEW_HUES":1,"UUID_NEW_DATES":1,
        "UUID_NEW_HUES":1,"FINAL_DATED_ORIG":4,
        "FINAL_DATED_HUES":4,
    }.items():
        monkeypatch.setattr(m,attr,value)
    original=pd.DataFrame([
        dict(photo_id="100",inat_taxon_id="5",morph="white",
             observation_id="10",observed_on="",qualified_date="2020-07-01",
             source_original_measured="historical"),
        dict(photo_id="101",inat_taxon_id="5",morph="red_pink",
             observation_id="11",observed_on="",qualified_date="",
             source_original_measured="historical"),
        dict(photo_id="102",inat_taxon_id="6",morph="blue_purple",
             observation_id="uuid-1",observed_on="",qualified_date="",
             source_original_measured="mois_wave1_new_original_ROI_v4"),
        dict(photo_id="103",inat_taxon_id="7",morph="yellow_orange",
             observation_id="12",observed_on="2026-07-04",qualified_date="2026-07-04",
             source_original_measured="fcp_wave4_partial_original_ROI_v4"),
        dict(photo_id="104",inat_taxon_id="7",morph="mixed_uncertain",
             observation_id="13",observed_on="",qualified_date="",
             source_original_measured="fcp_wave4_partial_original_ROI_v4"),
    ])
    original["original_source_date_provenance"]="old"
    numeric=original.iloc[:3].copy()
    numeric.loc[1,"qualified_date"]="2024-08-21"
    numeric.loc[1,"original_source_date_provenance"]=(
        "current_iNat_exact_numeric_observation_original_photo_date_verified")
    uuid=pd.DataFrame([{
      "original_photo_id":"102",
      "original_first4096_observation_source_UUID":"uuid-1",
      "original_first4096_taxon_id":"6",
      "source_taxon_ID_matches_original_ROI_photo":"True",
      "current_FCP_source_observation_ID":"202",
      "verified_original_full_source_observed_on":"2025-06-18",
      "source_metadata_snapshot":"2026-09-27",
      "photo_source_maximum_abs_geocoordinate_difference":"0.0000",
      "original_source_photo_hue":"blue_purple",
      "date_identity_status":"photo_ID_original_taxon_and_snapshot_date_verified"
    }])
    return original,numeric,uuid


def test_exact_three_way_original_photo_date_recovery_only(monkeypatch):
    orig,numeric,uuid=mock_actual_source(monkeypatch)
    out,stats=m.apply_all_verified_days(orig,numeric,uuid)
    assert out.photo_id.tolist()==orig.photo_id.tolist()
    assert out.qualified_date.tolist()==[
        "2020-07-01","2024-08-21","2025-06-18","2026-07-04",""]
    assert stats["new_exact_numeric_observation_photo_ID_dates"]==1
    assert stats["new_original_UUID_only_observation_photo_dates"]==1
    assert stats["final_real_photo_IDs_with_four_hue_and_verified_observed_on"]==4
    assert orig.loc[1,"qualified_date"]==""
    assert out.morph.tolist()==orig.morph.tolist()


def test_deny_overwrite_existing_S3_date(monkeypatch):
    orig,numeric,uuid=mock_actual_source(monkeypatch)
    numeric.loc[1,"photo_id"]="100"
    with pytest.raises(ValueError,match="source.*photo|date"):
        m.apply_all_verified_days(orig,numeric,uuid)


def test_reject_unverified_photo_to_uuid_id_or_taxon(monkeypatch):
    orig,numeric,uuid=mock_actual_source(monkeypatch)
    uuid.loc[0,"original_first4096_observation_source_UUID"]="invented_numeric_obs"
    with pytest.raises(ValueError,match="UUID-only photograph"):
        m.apply_all_verified_days(orig,numeric,uuid)
    uuid.loc[0,"original_first4096_observation_source_UUID"]="uuid-1"
    uuid.loc[0,"original_first4096_taxon_id"]="9"
    with pytest.raises(ValueError,match="UUID-only photograph"):
        m.apply_all_verified_days(orig,numeric,uuid)


def test_original_pixels_and_colour_never_change_due_to_date_patch(monkeypatch):
    orig,numeric,uuid=mock_actual_source(monkeypatch)
    numeric.loc[1,"morph"]="white"
    with pytest.raises(ValueError,match="Original photo"):
        m.apply_all_verified_days(orig,numeric,uuid)


def test_no_double_date_source_for_one_original_photo(monkeypatch):
    orig,numeric,uuid=mock_actual_source(monkeypatch)
    uuid.loc[0,"original_photo_id"]="101"
    with pytest.raises(ValueError,match="double-count"):
        m.apply_all_verified_days(orig,numeric,uuid)
