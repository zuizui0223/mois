"""Synthetic FCP all-photo original observed_on pre-ROI source-admission tests.

Test one original sealed FCP per-species photo shard. Source PHOTO ID counts,
distinct known observation dates, full original photo count and outcome-free
ROI queue remain independently checked. No real photo pixels or network calls.
"""
from __future__ import annotations
import csv
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
S=importlib.util.spec_from_file_location("mois_full_source_dated",
 ROOT/"scripts/acquisition/mois_full_source_dated_flower_candidates_20261016.py")
m=importlib.util.module_from_spec(S)
S.loader.exec_module(m)

TAXA={11:("Alpha beta",5),20:("Gamma delta",4),30:("Omega xi",0)}


def example_shard(tmp:Path):
    counts=tmp/"all_42111_counts_synthetic.csv"
    with counts.open("w",newline="") as f:
        w=csv.writer(f)
        w.writerow(["inat_taxon_id","species","n_unique_source_photo_ids","snapshot_date"])
        for tid,(name,n) in TAXA.items():w.writerow([tid,name,n,m.SNAPSHOT])
    folder=tmp/"shard_000_of_128"
    folder.mkdir()
    (folder/"source_contract.json").write_text(json.dumps({
        "schema":"fcp_uncapped_shard_source_contract_v1",
        "source_species_frame_sha256":"synthetic_frame_sha",
        "photo_cap":None,"outcome_blind":True,
        "created_d2":m.SNAPSHOT,"shard":0,"nshards":128
    }))
    (folder/"shard_source_receipt.json").write_text(json.dumps({
        "schema":"fcp_frozen_source_shard_emission_v1",
        "status":m.SOURCE_STATUS,
        "source_species_frame_sha256":"synthetic_frame_sha",
        "source_species_ids":list(TAXA),
        "n_source_species":len(TAXA),
        "partition":0,"n_shards":128,
        "n_original_source_photo_ids":9,
        "historical_colour_outcomes_used":False,
        "photo_pixels_opened":False,"per_species_photo_cap":None
    }))
    dated=["2020-04-16","2020-07-09","", "2020-13-01","2027-01-01",
           "2021-03-01","2021-10-05","2019-02-12","none"]
    for tid,(species,n) in TAXA.items():
        directory=folder/"species"/str(tid)
        directory.mkdir(parents=True)
        photo=directory/"photo_manifest.ndjson"
        with photo.open("w") as f:
            for k in range(n):
                idx=(k if tid==11 else 5+k)
                row={
                    "inat_taxon_id":tid,"species":species,
                    "observation_id":f"synthetic_obs_{idx}",
                    "photo_id":str(101+idx),
                    "photo_url_large":f"https://inaturalist-open-data.s3.amazonaws.com/photos/{101+idx}/large.jpg",
                    "photo_license":"cc-by","attribution":"synthetic only",
                    "observer_id":str(10+idx),
                    "observed_on":dated[idx],
                    "latitude":(37 if tid==11 else -35),
                    "longitude":(140 if tid==11 else 138),
                    "positional_accuracy_m":5,
                    "source_snapshot_date":m.SNAPSHOT,
                    "source_query":"inaturalist_open_dataset_monthly_snapshot"
                }
                f.write(json.dumps(row,sort_keys=True)+"\n")
        (directory/"checkpoint.json").write_text(json.dumps({
            "taxon_id":tid,"status":"complete","source_photo_rows":n,
            "photo_cap":None,"colour_pixels_opened":False,
            "metadata_sha256":m.sha(photo)
        }))
    return folder,counts


def test_entire_source_shard_dates_no_flower_claim_and_no_species_cap(tmp_path):
    shard,counts=example_shard(tmp_path)
    private=tmp_path/"protected_dated.ndjson"
    public=tmp_path/"safe_public.json"
    r=m.audit_shard(shard,counts,private,public,synthetic=True,min_public_cell=10)
    assert r["status"]==m.STATUS
    assert r["n_source_taxa_in_this_independent_shard"]==3
    assert r["n_original_eligible_source_photographs_examined"]==9
    assert r["n_original_photo_observation_dates_usable"]==5
    assert r["n_original_photo_observation_dates_excluded"]==4
    assert r["date_screen_state"]=={
        "ORIGINAL_SOURCE_DATE_KNOWN":5,
        "SOURCE_DATE_AFTER_FROZEN_SNAPSHOT":1,
        "SOURCE_DATE_IMPOSSIBLE_CALENDAR_DAY":1,
        "SOURCE_DATE_NOT_VALID_ISO_DAY":2
    }
    assert r["source_photo_cap_per_taxon"] is None
    assert r["metadata_rejected_by_observer_added_Flowers_annotation"]==0
    assert r["source_annotation_term_Flowers_available_in_original_four_tables"] is False
    assert r["actual_flower_photographs_verified_from_pixels"]==0
    assert r["photo_hue_classification_attempted"] is False
    assert r["is_an_all_fcp_128_shard_complete_census"] is False
    assert r["source_original_counts_ledger_sha256"]==m.sha(counts)
    dated=[json.loads(z) for z in private.read_text().splitlines()]
    assert len(dated)==5
    assert all(z["original_photo_flowering_state"]=="NOT_YET_PHOTO_VERIFIED"
               for z in dated)
    assert all("observed_on" in z and "photo_id" in z for z in dated)
    assert all("photo_id" not in json.dumps(z) for z in [r])
    assert all("public_suppressed_below" in z for z in
               r["source_observation_month_domain_counts"])
    with pytest.raises(ValueError,match="overwrite"):
        m.audit_shard(shard,counts,private,public,synthetic=True)


def test_invalid_dates_are_not_guessed_using_photographer_upload_times():
    assert m.source_observation_day("2020-02-29")[0].month==2
    for raw in ("2021-02-29","2020-00-01","2020/03/15",
                "2020-05","2020-05-12T08:30:00","","nan",None):
        day,status=m.source_observation_day(raw)
        assert day is None
        assert status.startswith("SOURCE_DATE_")
    assert m.source_observation_day("2026-09-28")[0] is None
    assert m.climatic_season(-35,1)==("SOUTH_EXTRATROPICAL","SUMMER")
    assert m.climatic_season(0,1)==("TROPICAL","CALENDAR_Q1")
    assert m.climatic_season(35,1)==("NORTH_EXTRATROPICAL","WINTER")


def test_original_manifest_changed_or_fake_photo_row_blocked(tmp_path):
    shard,counts=example_shard(tmp_path)
    file=shard/"species"/"11"/"photo_manifest.ndjson"
    original=file.read_text()
    file.write_text(original.replace("2020-04-16","2020-06-16"))
    with pytest.raises(ValueError,match="SHA/count"):
        m.audit_shard(shard,counts,tmp_path/"a.ndjson",tmp_path/"a.json",
                      synthetic=True)
    assert not (tmp_path/"a.json").exists()
    file.write_text(original)
    counts.write_text(counts.read_text().replace("Gamma delta,4","Gamma delta,5"))
    with pytest.raises(ValueError,match="original dated metadata source shard|taxon checkpoint"):
        m.audit_shard(shard,counts,tmp_path/"a.ndjson",tmp_path/"a.json",
                      synthetic=True)
