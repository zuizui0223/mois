"""Date coverage tests for ALL immutable original RGFCA source-photo IDs.

Synthetic metadata only. This is not a flowering or white-pigment result.
"""
from __future__ import annotations

import csv
import json

import pytest

from scripts.acquisition.audit_fcp_original_rgfca_monthly_photo_opportunities import (
    execute, precise_day, source_coordinates, digest, REQUIRED,
)


def fixture(tmp_path):
    root=tmp_path/"source.csv"
    receipt=tmp_path/"source_receipt.json"
    entries=[
        (1,101,"Species A","o1","u1",35.0,135.0,150,"2026-01-10"),
        (2,101,"Species A","o1","u1",35.0,135.0,150,"2026-01-10"),
        (3,101,"Species A","o2","u2",35.01,135.01,250,"2026-01-15"),
        (4,101,"Species A","o3","u3",35.0,135.0,500,"2026-02-10"),
        (5,101,"Species A","o4","u4",35.01,135.01,500,"2026-02-11"),
        (6,101,"Species A","o4","u4",35.01,135.01,500,"2026-02-11"),
        (7,101,"Species A","o5","u5",35.0,135.0,0,"2026-03-10"),
        (8,104,"Species B","o6","u6",-34.0,20.0,300,""),
        (9,104,"Species B","o7","u7",-34.0,20.0,25000,"2026-09-27"),
        (10,104,"Species B","o8","u8",-34.0,20.0,50,"2026-10-01"),
    ]
    rows=[]
    for pid,tax,name,obs,user,lat,lon,accuracy,day in entries:
        rows.append(dict(
            inat_taxon_id=str(tax),species=name,photo_id=str(pid),
            observation_id=obs,observer_id=user,latitude=str(lat),
            longitude=str(lon),positional_accuracy_m=str(accuracy),
            observed_on=day,source_snapshot_date="2026-09-27"))
    with root.open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=REQUIRED,lineterminator="\n")
        writer.writeheader();writer.writerows(rows)
    receipt.write_text(json.dumps({
        "status":"COMPLETE_TARGETED_1000_MONTHLY_SOURCE_IDS_COUNT_MATCHED_NO_PHOTO_PIXELS",
        "source_snapshot_date":"2026-09-27",
        "selected_all_photo_IDs_csv_sha256":digest(root),
        "original_image_pixels_opened":False,
        "historical_flower_colour_CSV_imported":False,
        "all_enrolled_species_match_independent_full_source_photo_ID_counts":True,
        "n_selected_original_photo_IDs":len(rows),
        "enrolled_outcome_blind_taxa":2,
    }))
    return root,receipt


def csv_rows(path):
    with path.open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))


def test_only_full_original_observation_day_and_defensible_site_precision():
    assert precise_day("2026-01-10")=="2026-01"
    for invalid in ("2026-01","2026-02-30","2026-10-01","",None,"2026/02/10"):
        assert precise_day(invalid)==""
    assert source_coordinates(35.0,135.0,"150")[-1] is True
    assert source_coordinates(35.0,135.0,"0")[-1] is False
    assert source_coordinates(35.0,135.0,"")[-1] is False
    assert source_coordinates(-34.0,20.0,"25000")[-1] is False
    with pytest.raises(ValueError,match="coordinates"):
        source_coordinates(95.,135.,"50")


def test_real_original_photo_source_date_census_has_no_flower_outcome(tmp_path):
    inputs=fixture(tmp_path)
    out=tmp_path/"results"
    r=execute(*inputs,out,_small_test=True)
    assert r["status"]=="SYNTHETIC_ONLY_NOT_BIOLOGICAL_PHOTOS"
    assert r["n_taxa"]==2
    assert r["n_original_photo_ids"]==10
    assert r["n_photo_ids_with_full_observed_on_date"]==8
    assert r["n_photo_ids_without_valid_observed_on_date"]==2
    assert r["n_dated_source_photos_5km_location_accuracy"]==6
    assert r["n_taxa_with_original_dated_photo"]==2
    assert r["n_taxa_with_repeated_source_calendar_months"]==1
    assert r["n_taxa_with_supported_0p25_degree_multiple_month_same_year"]==1
    assert r["n_supported_geocell_years_two_months_two_obs_two_observers_each"]==1
    assert r["source_photo_opportunity_is_not_flowering_detection"] is True
    assert r["all_original_image_pixels_classified"] is False
    assert r["new_flower_colour_or_white_frequency_estimated"] is False
    monthly=csv_rows(out/"all_original_1000_source_photo_observation_months.csv")
    assert len(monthly)==4
    jan=next(x for x in monthly if x["inat_taxon_id"]=="101" and
             x["year_month"]=="2026-01")
    assert jan["all_original_date_known_photo_ids"]=="3"
    assert jan["distinct_original_observations"]=="2"
    assert jan["distinct_original_observers"]=="2"
    sites=csv_rows(out/"all_original_1000_source_month_0p25degree_5km_accuracy.csv")
    assert len(sites)==2
    feb=next(x for x in sites if x["year_month"]=="2026-02")
    assert feb["all_source_photo_ids_with_5km_accuracy"]=="3"
    assert feb["distinct_source_observations"]=="2"
    assert feb["distinct_source_observers"]=="2"
    taxa=csv_rows(out/"all_original_1000_taxa_date_and_site_sampling_coverage.csv")
    assert len(taxa)==2
    a=next(x for x in taxa if x["inat_taxon_id"]=="101")
    assert a["all_original_photo_ids"]=="7"
    assert a["photos_with_valid_original_observed_day"]=="7"
    assert a["distinct_original_months_with_dated_candidates"]=="3"
    assert a["independent_site_years_with_two_supported_observation_months"]=="1"
    b=next(x for x in taxa if x["inat_taxon_id"]=="104")
    assert b["all_original_photo_ids"]=="3"
    assert b["photos_without_valid_original_observed_day"]=="2"
    assert b["independent_site_years_with_two_supported_observation_months"]=="0"
    assert r["original_source_months_file_sha256"]==digest(
        out/"all_original_1000_source_photo_observation_months.csv")
    with pytest.raises(ValueError,match="overwrite"):
        execute(*inputs,out,_small_test=True)


def test_tampered_source_manifest_is_rejected_even_with_same_count(tmp_path):
    manifest,receipt=fixture(tmp_path)
    manifest.write_text(manifest.read_text().replace("2026-01-10","2026-01-11",1))
    with pytest.raises(ValueError,match="exact complete dated"):
        execute(manifest,receipt,tmp_path/"bad",_small_test=True)


def test_fake_duplicate_original_photo_id_is_rejected_even_with_new_sha(tmp_path):
    manifest,receipt=fixture(tmp_path)
    rows=csv_rows(manifest)
    rows[4]["photo_id"]="1"
    with manifest.open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=REQUIRED)
        writer.writeheader();writer.writerows(rows)
    report=json.loads(receipt.read_text())
    report["selected_all_photo_IDs_csv_sha256"]=digest(manifest)
    receipt.write_text(json.dumps(report))
    with pytest.raises(ValueError,match="photo ID repeated"):
        execute(manifest,receipt,tmp_path/"bad",_small_test=True)


def test_full_original_1000_source_gate_blocks_synthetic_as_real(tmp_path):
    manifest,receipt=fixture(tmp_path)
    with pytest.raises(ValueError,match="Full original 1000 species"):
        execute(manifest,receipt,tmp_path/"bad",_small_test=False)


def test_missing_licence_and_observed_source_identity_blocks_reuse(tmp_path):
    manifest,receipt=fixture(tmp_path)
    report=json.loads(receipt.read_text())
    report["all_enrolled_species_match_independent_full_source_photo_ID_counts"]=False
    receipt.write_text(json.dumps(report))
    with pytest.raises(ValueError,match="exact complete dated"):
        execute(manifest,receipt,tmp_path/"bad",_small_test=True)
