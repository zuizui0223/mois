"""Strict photo-ID, UUID original-source day and preserved numeric-date tests."""
import csv
import hashlib
import json
from pathlib import Path
import pandas as pd
import pytest

from scripts.analysis import mois_full302242_uuid4096_numeric1732_exact_observed_on as m


def example():
    d=pd.DataFrame([
       dict(photo_id="101",observation_id="1101",inat_taxon_id="40",
            latitude="40.125",longitude="141.25",morph="white",
            qualified_date="",original_source_date_provenance="old",
            old_mois_or_new_fcp="historical"),
       dict(photo_id="102",observation_id="1102",inat_taxon_id="40",
            latitude="40.125",longitude="141.25",morph="red_pink",
            qualified_date="2025-04-17",original_source_date_provenance="historic_S3",
            old_mois_or_new_fcp="historical"),
       dict(photo_id="103",observation_id="aabbccdd-1111-2222-3333-444455556666",
            inat_taxon_id="41",latitude="34.1",longitude="135.15",
            morph="blue_purple",qualified_date="",
            original_source_date_provenance="unknown",
            old_mois_or_new_fcp="mois_new_original_ROI"),
       dict(photo_id="104",observation_id="aabbccdd-1111-2222-3333-777788889999",
            inat_taxon_id="42",latitude="34.2",longitude="135.15",
            morph="mixed_uncertain",qualified_date="",
            original_source_date_provenance="unknown",
            old_mois_or_new_fcp="mois_new_original_ROI"),
    ])
    n=d.iloc[:2].copy()
    n.loc[0,"qualified_date"]="2026-06-17"
    n.loc[0,"original_source_date_provenance"]="numeric_api_exact_id"
    return d,n


def test_exact_numeric_photo_date_only_no_old_hue_overwrite():
    old,numeric=example()
    out,r=m.numeric_patch(old,numeric,first_n=2)
    assert r=={"new_original_numeric_photo_observation_date_IDs":1,
               "new_original_numeric_four_hue_date_photo_IDs":1}
    assert out.qualified_date.tolist()==[
        "2026-06-17","2025-04-17","",""]
    assert out.photo_id.tolist()==old.photo_id.tolist()
    assert out.morph.tolist()==old.morph.tolist()
    assert old.qualified_date.iloc[0]==""


def test_numerical_source_cannot_change_preexisting_hue_or_date():
    old,numeric=example()
    numeric.loc[0,"morph"]="yellow_orange"
    with pytest.raises(ValueError,match="pixels/taxon/hue"):
        m.numeric_patch(old,numeric,first_n=2)
    _,numeric=example()
    numeric.loc[1,"qualified_date"]="2024-04-17"
    with pytest.raises(ValueError,match="replace"):
        m.numeric_patch(old,numeric,first_n=2)


def test_original_UUID_photo_date_and_unclassifiable_photo_remain_distinct(monkeypatch):
    monkeypatch.setattr(m,"UUID_MISSING",2)
    monkeypatch.setattr(m,"UUID_NEW_HUE",1)
    old,numeric=example()
    patched,_=m.numeric_patch(old,numeric,first_n=2)
    missing=m.missing_original_UUID_rows(patched)
    assert len(missing)==2
    uuids=pd.DataFrame([
      dict(original_photo_id="103",
           matched_source_original_observation_UUID=old.iloc[2].observation_id,
           matched_source_original_taxon_id="41",verified_original_observed_on="2025-05-16"),
      dict(original_photo_id="104",
           matched_source_original_observation_UUID=old.iloc[3].observation_id,
           matched_source_original_taxon_id="42",verified_original_observed_on="2026-09-27"),
    ])
    result,counts=m.apply_uuid_days(patched,uuids)
    assert counts["new_source_UUID_original_photo_day_IDs"]==2
    assert counts["new_UUID_source_classifiable_original_hues"]==1
    assert result.qualified_date.str.len().eq(10).all()
    assert result.morph.eq("mixed_uncertain").sum()==1


def test_source_UUID_mismatch_or_existing_date_rejected():
    old,numeric=example()
    p,_=m.numeric_patch(old,numeric,first_n=2)
    wrong=pd.DataFrame([
        dict(original_photo_id="103",
             matched_source_original_observation_UUID="different-UUID",
             matched_source_original_taxon_id="41",
             verified_original_observed_on="2025-05-16")
    ])
    with pytest.raises(ValueError,match="UUID/taxon"):
        m.apply_uuid_days(p,wrong)
    wrong.loc[0,"matched_source_original_observation_UUID"]=old.iloc[2].observation_id
    wrong.loc[0,"original_photo_id"]="102"
    wrong.loc[0,"matched_source_original_taxon_id"]="40"
    with pytest.raises(ValueError,match="overwrite"):
        m.apply_uuid_days(p,wrong)


def test_full_original_provider_source_metadata_scan_signed_SHA_and_identity(monkeypatch,tmp_path):
    old,_=example()
    rows=[
      dict(photo_id="101",observation_id="1101",inat_taxon_id="40",latitude="40.125",longitude="141.25",observed_on="2026-06-17",source_snapshot_date="2026-09-27"),
      dict(photo_id="103",observation_id=old.iloc[2].observation_id,inat_taxon_id="41",latitude="34.1",longitude="135.15",observed_on="2025-05-16",source_snapshot_date="2026-09-27"),
      dict(photo_id="104",observation_id=old.iloc[3].observation_id,inat_taxon_id="42",latitude="34.2",longitude="135.15",observed_on="2026-09-27",source_snapshot_date="2026-09-27"),
    ]
    source=tmp_path/"original.csv"
    with source.open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)
    hashval=hashlib.sha256(source.read_bytes()).hexdigest()
    monkeypatch.setattr(m,"SOURCE_SHA",hashval)
    monkeypatch.setattr(m,"TOTAL_ORIGINAL_SOURCE_ROWS",3)
    r=tmp_path/"signed_source.json"
    r.write_text(json.dumps({
       "schema":"fcp_rgfca_targeted_original_20260927_photo_ID_census_v1",
       "selected_all_photo_IDs_csv_sha256":hashval,
       "n_selected_original_photo_IDs":3,
       "source_snapshot_date":"2026-09-27"
    }))
    result,stats=m.scan_all_full_FCP_source(source,r,old.iloc[2:])
    assert len(result)==2
    assert stats["original_3430502_source_full_CSV_SHA256"]==hashval
    assert result.verified_original_observed_on.tolist()==[
        "2025-05-16","2026-09-27"]
    rows[2]["observed_on"]="2027-01-01"
    with source.open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)
    # Re-pinning the current mutated source never authorizes future-dated observations.
    hashval=hashlib.sha256(source.read_bytes()).hexdigest()
    monkeypatch.setattr(m,"SOURCE_SHA",hashval)
    r.write_text(json.dumps({
       "schema":"fcp_rgfca_targeted_original_20260927_photo_ID_census_v1",
       "selected_all_photo_IDs_csv_sha256":hashval,
       "n_selected_original_photo_IDs":3,
       "source_snapshot_date":"2026-09-27"
    }))
    with pytest.raises(ValueError,match="day missing or invalid"):
        m.scan_all_full_FCP_source(source,r,old.iloc[2:])
