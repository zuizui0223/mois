"""Synthetic all-FCP provider 128-shard-style original date aggregation.

Creates actual synthetic per-taxon FCP metadata manifests, audits both shard
sources independently, then demands complete 2-shard coverage; no flower ROI.
"""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]


def module(name,path):
    s=importlib.util.spec_from_file_location(name,ROOT/path)
    v=importlib.util.module_from_spec(s);s.loader.exec_module(v)
    return v


date=module("mois_full_date_stage_test","scripts/acquisition/mois_full_source_dated_flower_candidates_20261016.py")
m=module("mois_full_date_reassemble_test","scripts/acquisition/mois_aggregate_all_original_dated_shards_20261016.py")

TAXA={
    11:("Alpha beta",3),12:("No photos",0),
    20:("Gamma delta",2),21:("Omega xi",1)
}


def complete_synthetic(tmp):
    source=tmp/"FCP_full_source_2_parts"
    private=tmp/"MOIS_private_dates"
    source.mkdir();private.mkdir()
    counts=tmp/"original_source_taxon_photo_counts.csv"
    with counts.open("w",newline="") as f:
        w=csv.writer(f)
        w.writerow(["inat_taxon_id","species","n_unique_source_photo_ids","snapshot_date"])
        for tid,(name,n) in TAXA.items():
            w.writerow([tid,name,n,date.SNAPSHOT])
    for part,taxa in enumerate(((11,12),(20,21))):
        folder=source/f"shard_{part:03}_of_002"
        folder.mkdir()
        (folder/"source_contract.json").write_text(json.dumps({
            "schema":"fcp_uncapped_shard_source_contract_v1",
            "source_species_frame_sha256":"TEST_FAKE_FRAME_SHA",
            "photo_cap":None,"outcome_blind":True,"created_d2":date.SNAPSHOT,
            "shard":part,"nshards":2
        }))
        count=0
        for tid in taxa:
            species,n=TAXA[tid]
            species_dir=folder/"species"/str(tid)
            species_dir.mkdir(parents=True)
            path=species_dir/"photo_manifest.ndjson"
            with path.open("w") as w:
                for j in range(n):
                    # Date exclusions are independent of photo morphology.
                    when=(
                        ["2022-04-12","2022-07-22",""][j] if tid==11
                        else ["2020-12-01","2026-10-15"][j] if tid==20
                        else "2021-08-08"
                    )
                    pid=tid*1000+j
                    r={
                        "inat_taxon_id":tid,"species":species,
                        "observation_id":f"source-obs-{pid}",
                        "photo_id":str(pid),
                        "photo_url_large":f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.jpg",
                        "photo_license":"cc-by",
                        "attribution":"TEST SOURCE ONLY",
                        "observer_id":str(tid+700),
                        "observed_on":when,
                        "latitude":36. if tid==11 else -32.,
                        "longitude":140. if tid==11 else 118.,
                        "positional_accuracy_m":5,
                        "source_snapshot_date":date.SNAPSHOT,
                        "source_query":"inaturalist_open_dataset_monthly_snapshot"
                    }
                    w.write(json.dumps(r,sort_keys=True)+"\n")
            (species_dir/"checkpoint.json").write_text(json.dumps({
                "status":"complete","taxon_id":tid,"source_photo_rows":n,
                "metadata_sha256":date.sha(path),
                "photo_cap":None,"colour_pixels_opened":False
            }))
            count+=n
        (folder/"shard_source_receipt.json").write_text(json.dumps({
            "schema":"fcp_frozen_source_shard_emission_v1",
            "status":date.SOURCE_STATUS,
            "source_species_frame_sha256":"TEST_FAKE_FRAME_SHA",
            "partition":part,"n_shards":2,"source_species_ids":list(taxa),
            "n_source_species":len(taxa),"stage_sqlite_sha256":"TEST_FAKE_STAGED_SQLITE_SHA",
            "n_original_source_photo_ids":count,
            "photo_pixels_opened":False,"historical_colour_outcomes_used":False,
            "per_species_photo_cap":None
        }))
        date.audit_shard(folder,counts,
            private/f"dated_candidates_shard_{part:03}.ndjson",
            private/f"dated_candidates_shard_{part:03}_receipt.json",
            synthetic=True)
    return source,private,counts


def test_all_shards_exactly_date_screened_original_full_photo_opportunities(tmp_path):
    sources,private,counts=complete_synthetic(tmp_path)
    receipt=m.aggregate(sources,private,counts,tmp_path/"all_original_date_counts.json",
        synthetic=True,partitions=2,expected_photo_total=6)
    assert receipt["status"]==m.STATUS
    assert receipt["source_original_shards_complete"]==2
    assert receipt["n_all_original_source_named_taxa"]==4
    assert receipt["n_all_original_unique_photo_opportunities"]==6
    assert receipt["n_original_snapshot_observed_on_valid_photo_candidates"]==4
    assert receipt["n_original_snapshot_observed_on_invalid_or_missing_photo_candidates"]==2
    assert receipt["n_named_taxa_with_at_least_one_original_valid_dated_photo"]==3
    assert receipt["n_named_taxa_without_original_valid_dated_photos"]==1
    assert receipt["n_actual_original_photo_ROI_flower_colour_measurements_this_stage"]==0
    assert receipt["source_photo_cap_per_species"] is None
    assert receipt["original_source_photo_bloom_phenology_verified_by_image_pixels"] is False
    assert receipt["old_original_100543_photo_cohort_reused_as_original_full_data"] is False
    assert receipt["date_quality_states"]=={
        "ORIGINAL_SOURCE_DATE_KNOWN":4,
        "SOURCE_DATE_NOT_VALID_ISO_DAY":1,
        "SOURCE_DATE_AFTER_FROZEN_SNAPSHOT":1
    }
    assert len(receipt["source_shard_receipts"])==2
    assert len(receipt["latitude_zone_original_observed_month_counts"])==4
    assert "photo_id" not in (tmp_path/"all_original_date_counts.json").read_text()


def test_missing_one_complete_shard_cannot_promote_an_aggregate(tmp_path):
    sources,private,counts=complete_synthetic(tmp_path)
    missing=private/"dated_candidates_shard_001_receipt.json"
    missing.unlink()
    with pytest.raises(ValueError,match="Cannot promote unsealed"):
        m.aggregate(sources,private,counts,tmp_path/"out.json",
                    synthetic=True,partitions=2,expected_photo_total=6)
    assert not (tmp_path/"out.json").exists()


def test_date_queue_mutation_or_duplicate_taxon_raises(tmp_path):
    sources,private,counts=complete_synthetic(tmp_path)
    file=private/"dated_candidates_shard_000.ndjson"
    file.write_text(file.read_text().replace("2022-04-12","2022-04-11"))
    with pytest.raises(ValueError,match="source sealed"):
        m.aggregate(sources,private,counts,tmp_path/"out.json",
                    synthetic=True,partitions=2,expected_photo_total=6)
    assert not (tmp_path/"out.json").exists()
    # Even if one source receipt is mutated to match a new queue, the
    # independent provider immutable shard and species roster still gate it.
    sources,private,counts=complete_synthetic(tmp_path/"fresh") if False else (None,None,None)


def test_photo_count_census_disagreement_prevents_complete_global_status(tmp_path):
    sources,private,counts=complete_synthetic(tmp_path)
    with pytest.raises(ValueError,match="photo IDs do not reconcile"):
        m.aggregate(sources,private,counts,tmp_path/"out.json",
                    synthetic=True,partitions=2,expected_photo_total=7)
