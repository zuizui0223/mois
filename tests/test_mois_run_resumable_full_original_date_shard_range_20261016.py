"""Source-photo date resumability: exact original metadata, no source cap."""
from __future__ import annotations
import importlib.util
import json
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]


def module(name,path):
    s=importlib.util.spec_from_file_location(name,ROOT/path)
    x=importlib.util.module_from_spec(s);s.loader.exec_module(x)
    return x


m=module("mois_source_resumable_date_range",
         "scripts/acquisition/mois_run_resumable_full_original_date_shard_range_20261016.py")
fixtures=module("mois_all_fcp_original_photo_shard_fixture",
         "tests/test_mois_aggregate_all_original_dated_shards_20261016.py")


def test_complete_2_shard_source_receipts_reverified_without_replacing(tmp_path):
    shards,private,counts=fixtures.complete_synthetic(tmp_path)
    result=m.run(shards,counts,private,tmp_path/"partial_block_0_2.json",
                 start=0,stop=2,partitions=2,synthetic=True)
    assert result["status"]==m.STATUS
    assert result["n_shards_verified_in_this_execution"]==2
    assert result["photo_opportunities_in_this_partial_source_block"]==6
    assert result["dated_candidates_in_this_partial_source_block"]==4
    assert all(x["state"]=="PREVIOUSLY_COMPLETED_REVERIFIED"
               for x in result["per_original_shard"])
    assert result["source_photo_sampling_cap_per_taxon"] is None
    assert result["photo_image_pixels_opened"] is False
    assert result["flowering_photos_counted_or_inferred"] is False
    assert result["worldwide_complete_42111_taxon_93million_date_census"] is False
    assert "photo_id" not in (tmp_path/"partial_block_0_2.json").read_text()
    with pytest.raises(ValueError,match="Do not overwrite"):
        m.run(shards,counts,private,tmp_path/"partial_block_0_2.json",
              start=0,stop=2,partitions=2,synthetic=True)


def test_resumption_of_new_original_shard_without_erasing_prior_seal(tmp_path):
    shards,private,counts=fixtures.complete_synthetic(tmp_path)
    before=private/"dated_candidates_shard_000.ndjson"
    digest=m.sha(before)
    (private/"dated_candidates_shard_001.ndjson").unlink()
    (private/"dated_candidates_shard_001_receipt.json").unlink()
    r=m.run(shards,counts,private,tmp_path/"one_new.json",
            start=1,stop=2,partitions=2,synthetic=True)
    assert r["per_original_shard"][0]["state"]=="NEW_SOURCE_DATES_VALIDATED"
    assert r["dated_candidates_in_this_partial_source_block"]==2
    assert m.sha(before)==digest
    recheck=m.run(shards,counts,private,tmp_path/"reverified.json",
                  start=0,stop=2,partitions=2,synthetic=True)
    assert all(x["state"]=="PREVIOUSLY_COMPLETED_REVERIFIED"
               for x in recheck["per_original_shard"])


def test_missing_one_half_finished_shard_fails_closed_without_guessing_photo_dates(tmp_path):
    shards,private,counts=fixtures.complete_synthetic(tmp_path)
    (private/"dated_candidates_shard_001_receipt.json").unlink()
    with pytest.raises(ValueError,match="Partially written"):
        m.run(shards,counts,private,tmp_path/"not_completed.json",
              start=1,stop=2,partitions=2,synthetic=True)
    assert not (tmp_path/"not_completed.json").exists()


def test_bounded_original_work_shard_range_not_biological_sample(tmp_path):
    shards,private,counts=fixtures.complete_synthetic(tmp_path)
    with pytest.raises(ValueError,match="four original source photo shards"):
        m.run(shards,counts,private,tmp_path/"not_completed.json",
              start=0,stop=5,partitions=6,synthetic=True)
