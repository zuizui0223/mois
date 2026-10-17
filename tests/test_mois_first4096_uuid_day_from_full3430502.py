"""Small source-ledger fixture: photo_ID not observation UUID predicts source day.

These fake three rows test correctness only, never reported as real biology.
"""
import csv
import hashlib
from pathlib import Path

import pandas as pd
import pytest

from scripts.acquisition import mois_first4096_undated_uuid_photo_days_from_original3430502 as m


def original_source_fixture(tmp_path:Path,monkeypatch):
    path=tmp_path/"authoritative_photo_source.csv"
    items=[
        {"photo_id":"111","inat_taxon_id":"220","observation_id":"777",
         "observed_on":"2025-06-19","source_snapshot_date":"2026-09-27",
         "latitude":"42.0000","longitude":"140.0000"},
        {"photo_id":"222","inat_taxon_id":"221","observation_id":"778",
         "observed_on":"2024-04-09","source_snapshot_date":"2026-09-27",
         "latitude":"43.0000","longitude":"141.0000"},
        {"photo_id":"333","inat_taxon_id":"222","observation_id":"779",
         "observed_on":"invalid","source_snapshot_date":"2026-09-27",
         "latitude":"44.0000","longitude":"142.0000"},
    ]
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=items[0].keys())
        w.writeheader()
        w.writerows(items)
    monkeypatch.setattr(m,"FCP_AUTHORITATIVE_ROWS",3)
    monkeypatch.setattr(m,"FCP_AUTHORITATIVE_FULL_SOURCE_SHA",
                        hashlib.sha256(path.read_bytes()).hexdigest())
    target=pd.DataFrame([
        {"photo_id":"111","inat_taxon_id":"220","observation_id":"uuid-xx",
         "latitude":"42","longitude":"140","morph":"white"},
        {"photo_id":"333","inat_taxon_id":"222","observation_id":"uuid-yy",
         "latitude":"44","longitude":"142","morph":"red_pink"},
    ])
    return path,target


def test_real_source_photo_id_and_taxon_map_uuid_to_photo_day(monkeypatch,tmp_path):
    path,target=original_source_fixture(tmp_path,monkeypatch)
    found,proof=m.scan_full_source(path,target)
    assert proof["FCP_source_rows_fully_scanned"]==3
    assert set(found)=={"111","333"}
    assert found["111"]["verified_original_full_source_observed_on"]=="2025-06-19"
    assert found["111"]["current_FCP_source_observation_ID"]=="777"
    assert found["111"]["original_first4096_observation_source_UUID"]=="uuid-xx"
    assert found["333"]["verified_original_full_source_observed_on"] is None


def test_source_historical_taxon_mismatch_is_not_joined(monkeypatch,tmp_path):
    path,target=original_source_fixture(tmp_path,monkeypatch)
    target.loc[0,"inat_taxon_id"]="different_taxon"
    found,_=m.scan_full_source(path,target)
    assert found["111"]["verified_original_full_source_observed_on"] is None
    assert found["111"]["source_taxon_ID_matches_original_ROI_photo"] is False


def test_source_drift_fails_even_when_target_ids_all_match(monkeypatch,tmp_path):
    path,target=original_source_fixture(tmp_path,monkeypatch)
    monkeypatch.setattr(m,"FCP_AUTHORITATIVE_FULL_SOURCE_SHA","0"*64)
    with pytest.raises(ValueError,match="source SHA"):
        m.scan_full_source(path,target)


def test_source_invalid_calendar_day_is_not_imputed():
    assert m.valid_day("2025-03-01")=="2025-03-01"
    assert m.valid_day("2025-02-30") is None
    assert m.valid_day("") is None
    assert m.valid_day("1800-01-01")=="1800-01-01"
    assert m.valid_day("1700-01-01") is None


def test_reused_original_photo_duplicate_ID_must_fail(monkeypatch,tmp_path):
    path,target=original_source_fixture(tmp_path,monkeypatch)
    with path.open("a") as f:
        f.write("111,220,777,2025-06-19,2026-09-27,42.0,140.0\n")
    monkeypatch.setattr(m,"FCP_AUTHORITATIVE_ROWS",4)
    monkeypatch.setattr(m,"FCP_AUTHORITATIVE_FULL_SOURCE_SHA",
                        hashlib.sha256(path.read_bytes()).hexdigest())
    with pytest.raises(ValueError,match="repeats"):
        m.scan_full_source(path,target)
