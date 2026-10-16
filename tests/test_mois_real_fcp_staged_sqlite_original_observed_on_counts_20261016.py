"""Synthetic exact FCP original permanent SQL date count; no image pixels."""
from __future__ import annotations

import csv
import importlib.util
import json
import sqlite3
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
module=importlib.util.spec_from_file_location(
    "mois_source_full_date_histogram",
    ROOT/"scripts/acquisition/mois_real_fcp_staged_sqlite_original_observed_on_counts_20261016.py")
m=importlib.util.module_from_spec(module)
module.loader.exec_module(m)


def fixture(tmp):
    stage=tmp/"original_staged"
    stage.mkdir()
    photo=stage/"_working_fcp_snapshot.sqlite"
    with sqlite3.connect(photo) as db:
        db.execute("""CREATE TABLE observations(
            uuid TEXT PRIMARY KEY,taxon INTEGER,observer TEXT,observed_on TEXT,
            lat REAL,lon REAL,accuracy TEXT)""")
        db.execute("""CREATE TABLE photos(
            photo_id TEXT PRIMARY KEY,uuid TEXT,license TEXT,extension TEXT,
            photographer TEXT,taxon INTEGER NOT NULL,
            ambiguous INTEGER NOT NULL DEFAULT 0)""")
        db.executemany("INSERT INTO observations VALUES (?,?,?,?,?,?,?)",[
            ("obsa",11,"a","2025-04-03",34.,140.,"10"),
            ("obsb",11,"b","2025-08-23",34.,140.,"10"),
            ("obsc",20,"c","2025-12-10",-34.,151.,"10"),
            ("obsd",20,"d","",0.,22.,"10"),
            ("obse",30,"e","2027-01-03",0.,22.,"10"),
            ("obsf",30,"e","2021-02-29",12.,22.,"10")
        ])
        db.executemany("INSERT INTO photos VALUES (?,?,?,?,?,?,?)",[
            ("101","obsa","cc-by","jpg","a",11,0),
            ("102","obsa","cc-by","jpg","a",11,0),
            ("103","obsb","cc-by","jpg","b",11,0),
            ("104","obsc","cc-by","jpg","c",20,0),
            ("105","obsd","cc-by","jpg","d",20,0),
            ("106","obse","cc-by","jpg","e",30,0),
            ("107","obsf","cc-by","jpg","f",30,0),
            # Quarantined duplicate photo ID source does NOT count as
            # unique licensed original source opportunity.
            ("108","obsf","cc-by","jpg","f",30,1)
        ])
        db.commit()
    counts=tmp/"original_photo_counts.csv"
    with counts.open("w",newline="") as f:
        w=csv.writer(f)
        w.writerow(("inat_taxon_id","species","n_unique_source_photo_ids","snapshot_date"))
        w.writerows([
            (11,"Alpha beta",3,"2026-09-27"),
            (20,"Gamma delta",2,"2026-09-27"),
            (30,"Omega xi",2,"2026-09-27"),
            (40,"Zero original photos",0,"2026-09-27"),
        ])
    doc={
        "schema":m.STAGE_SCHEMA,"status":m.STAGE_STATUS,
        "sqlite_basename":photo.name,
        "sqlite_bytes":photo.stat().st_size,
        "sqlite_sha256":m.sha(photo),
        "historical_flower_colour_csv_used":False,
        "source_image_pixels_opened":False,
        "original_photo_ids_deduplicated_globally":True,
        "full_source_census_complete":False,
        "full_study_original_images_measured":False,
        "source_snapshot_date":"2026-09-27",
        "source_species":4,
        "source_species_frame_sha256":"SYNTHETIC_FRAME",
        "n_target_taxon_licensed_photos":7,
        "shards":2,
        "verified_original_dated_archive":{"source_archive_sha256":"SYNTHETIC_ONLY"}
    }
    (stage/"fcp_staged_sqlite_source_receipt.json").write_text(json.dumps(doc))
    return stage,counts


def test_full_source_sqlite_original_date_counts_preserve_each_photo_even_same_observation(tmp_path):
    stage,counts=fixture(tmp_path)
    outfile=tmp_path/"all_source_months.json"
    result=m.audit_staged_original_dates(stage,counts,outfile,synthetic=True,min_public=10)
    assert result["status"]==m.STATUS
    assert result["n_original_source_named_taxa"]==4
    assert result["n_original_photo_opportunities_with_globally_unique_ids"]==7
    assert result["n_original_snapshot_observed_on_date_known"]==4
    assert result["n_original_snapshot_observed_on_date_unknown_or_invalid"]==3
    assert result["original_observed_on_status_counts"]=={
        "ORIGINAL_SOURCE_DATE_KNOWN":4,
        "SOURCE_DATE_NOT_VALID_ISO_DAY":1,
        "SOURCE_DATE_IMPOSSIBLE_CALENDAR_DAY":1,
        "SOURCE_DATE_AFTER_FROZEN_SNAPSHOT":1
    }
    assert result["n_original_named_taxa_with_valid_observed_on_photos"]==2
    assert result["n_source_taxa_with_no_valid_original_observation_dates"]==2
    assert len(result["original_latitude_domain_observed_month_photo_opportunity_counts"])==3
    assert sum(v["public_suppressed_below"]==10 for v in result[
        "original_latitude_domain_observed_month_photo_opportunity_counts"])==3
    assert result["all_original_photo_ids_inclusion_cap"] is None
    assert result["all_original_images_opened_or_flowering_roi_confirmed"] is False
    assert result["n_classifiable_flowering_photographs_estimated_here"]==0
    assert result["real_original_photo_white_fraction_estimated_here"] is False
    assert result["synthetic_test_only"] is True
    assert outfile.is_file()
    assert "photo_id" not in json.loads(outfile.read_text())
    assert "observation_uuid" not in json.loads(outfile.read_text())
    with pytest.raises(ValueError,match="cannot be overwritten"):
        m.audit_staged_original_dates(stage,counts,outfile,synthetic=True)


def test_wrong_original_complete_sqlite_byte_or_taxon_ledger_rejected(tmp_path):
    stage,counts=fixture(tmp_path)
    db=stage/"_working_fcp_snapshot.sqlite"
    with sqlite3.connect(db) as connection:
        connection.execute("UPDATE photos SET ambiguous=1 WHERE photo_id='101'")
        connection.commit()
    with pytest.raises(ValueError,match="not byte authenticated"):
        m.audit_staged_original_dates(stage,counts,tmp_path/"fail.json",synthetic=True)
    assert not (tmp_path/"fail.json").exists()


def test_true_production_42k_not_created_from_synthetic_fcp_sqlite(tmp_path):
    stage,counts=fixture(tmp_path)
    with pytest.raises(ValueError,match="Not genuine complete"):
        m.audit_staged_original_dates(stage,counts,tmp_path/"false_world.json",synthetic=False)
    assert not (tmp_path/"false_world.json").exists()


def test_date_month_conservation_does_not_assume_unannotated_photographs_have_no_flowers(tmp_path):
    stage,counts=fixture(tmp_path)
    result=m.audit_staged_original_dates(stage,counts,tmp_path/"month.json",synthetic=True)
    assert result["all_original_photo_ids_permanently_manifested_in_this_report"] is False
    assert result["historic_100543_photo_colour_outcomes_reused"] is False
    assert result["source_original_date_not_created_or_uploaded_at"] is True
    assert result["snapshot_photo_count_denominator_not_wild_blooming_individuals"] is True
    assert result["n_original_snapshot_observed_on_date_known"]==4
