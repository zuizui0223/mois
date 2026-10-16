"""Synthetic complete source 2-shard date-first image ROI worldwide audit.

The full provider photo frame is fake and SMALL, correctly marked synthetic.
No real global 93M count, source flowering pixels, or ecological coefficients
are claimed. Every fake dated photo has one original-FCP-compatible ROI
terminal, and photograph/observation/species monthly rates remain distinct.
"""
from __future__ import annotations
import importlib.util
import json
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]


def module(name,path):
    spec=importlib.util.spec_from_file_location(name,ROOT/path)
    obj=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


fixtures=module("mois_full_2shard_date_source",
 "tests/test_mois_aggregate_all_original_dated_shards_20261016.py")
packet=module("mois_full_source_date_dispatch",
 "scripts/acquisition/mois_prepare_all_dated_original_roi_packets_20261016.py")
worker=module("mois_full_original_roi_worker",
 "scripts/acquisition/mois_execute_all_dated_original_roi_workunit_20261016.py")
backend=module("mois_fake_original_photo_worker_for_end_to_end",
 "tests/test_mois_execute_all_dated_original_roi_workunit_20261016.py")
shard=module("mois_full_original_roi_shard_reassembly",
 "scripts/analysis/mois_reassemble_dated_original_roi_shard_20261016.py")
m=module("mois_full_date_first_bloom_global",
 "scripts/analysis/mois_complete_original_date_first_blooming_photo_census_20261016.py")
climate_test=module("mois_climate_full_source_date_first_fixture",
 "tests/test_mois_join_monthly_dynamic_climate_20261016.py")
climate_analysis=module("mois_global_observation_species_climate",
 "scripts/analysis/mois_full_bloom_monthly_TSPVPD_observation_species_FE_20261016.py")


def prepare_full(tmp,*,skip_last_shard=False):
    original,dates,counts=fixtures.complete_synthetic(tmp)
    audit=fixtures.m.aggregate(original,dates,counts,
        tmp/"all_original_date_receipt.json",
        synthetic=True,partitions=2,expected_photo_total=6)
    assert audit["n_original_snapshot_observed_on_valid_photo_candidates"]==4
    consolidated=tmp/"completed_roi_shards"
    consolidated.mkdir()
    sealed=tmp/"sealed_date_roi"
    all_packets=tmp/"source_roi_packets"
    sealed.mkdir();all_packets.mkdir()
    for i in range(1 if skip_last_shard else 2):
        name=f"shard_{i:03}_of_002"
        queue=dates/f"dated_candidates_shard_{i:03}.ndjson"
        receipt=dates/f"dated_candidates_shard_{i:03}_receipt.json"
        src=all_packets/name
        plan=packet.dispatch(queue,receipt,src,
            packet_size=2,synthetic=True)
        sub=sealed/name
        for k in range(plan["n_bounded_original_roi_packets"]):
            worker.execute(src/f"packet_{k:06d}",src,sub,k,0,
                backend=backend.FakeFCP(),
                source_script=tmp/"original_roi_script.py",
                original_detector=tmp/"original_yolo.pt",
                original_sam=tmp/"original_sam",synthetic=True)
        dest=consolidated/name
        dest.mkdir()
        proof=shard.reassemble(src,sub,
            dest/"positive_actual_bloom_photos.ndjson",
            dest/"roi_shard_receipt.json",
            synthetic=True)
        assert proof["n_original_date_known_photo_ROI_terminals"]==2
    return tmp/"all_original_date_receipt.json",consolidated


def test_full_original_2shard_photo_date_eligibility_bloom_classification_counts(tmp_path):
    audit,shards=prepare_full(tmp_path)
    dst=tmp_path/"protected_global_blooms.ndjson"
    receipt=tmp_path/"public_global_flower_white.json"
    r=m.aggregate_all(audit,shards,dst,receipt,
        synthetic=True,partitions=2,expected_photos=6,expected_taxa=4)
    assert r["schema"]==m.SCHEMA
    assert r["status"]==m.STATUS
    assert r["full_42111_original_photo_snapshot_source_count"]==6
    assert r["n_original_source_date_known_photo_candidates"]==4
    assert r["n_original_date_unknown_photo_opportunities_not_analysed_for_bloom"]==2
    assert r["n_original_date_known_photo_ROI_image_terminals"]==4
    assert r["n_original_date_and_site_known_ROI_bloom_hues"]==4
    assert r["n_original_date_known_white_flower_photos"]==4
    assert r["n_original_distinct_photo_observation_UUIDs_with_flower"]==4
    assert r["n_original_source_taxa_with_classifiable_date_known_flower"]==3
    assert r["n_original_source_genera_with_classifiable_date_known_flower"]==3
    assert r["source_shards_completed"]==2
    assert r["full_42111_original_source_not_selectively_date_sampled"] is True
    assert r["all_missing_date_source_photos_were_excluded_before_roi"] is True
    assert r["all_original_date_eligible_photo_pixels_attempted"] is True
    assert r["no_old_100543_photo_ROI_hues_reused"] is True
    assert r["original_photo_population_random_sample_of_all_wild_blooms"] is False
    assert r["synthetic_test_only"] is True
    positives=[json.loads(x) for x in dst.read_text().splitlines()]
    assert len(positives)==4 and all(
        row["actual_photo_contains_classifiable_flower"] for row in positives)
    assert len(set(row["original_photo_id"] for row in positives))==4
    for x in r["source_original_domain_month_classified_photo_white_and_species_equal"]:
        assert x["public_suppressed_below"]==10
    assert "original_photo_id" not in json.loads(receipt.read_text())
    assert "original_observation_uuid" not in json.loads(receipt.read_text())
    with pytest.raises(ValueError,match="overwritten"):
        m.aggregate_all(audit,shards,dst,receipt,
           synthetic=True,partitions=2,expected_photos=6,expected_taxa=4)


def test_missing_one_completed_roi_source_shard_cannot_fake_globally_complete(tmp_path):
    audit,shards=prepare_full(tmp_path,skip_last_shard=True)
    with pytest.raises(ValueError,match="not source-complete"):
        m.aggregate_all(audit,shards,tmp_path/"positive.ndjson",
            tmp_path/"false_global_receipt.json",
            synthetic=True,partitions=2,expected_photos=6,expected_taxa=4)
    assert not (tmp_path/"false_global_receipt.json").exists()


def test_genuine_completed_date_first_photo_positive_source_feeds_monthly_worldclim(tmp_path):
    audit,shards=prepare_full(tmp_path)
    positive=tmp_path/"positive_all.ndjson"
    proof=tmp_path/"global_bloom_receipt.json"
    m.aggregate_all(audit,shards,positive,proof,
        synthetic=True,partitions=2,expected_photos=6,expected_taxa=4)
    climate=climate_test.normal(tmp_path/"weather") if False else None
    # The weather-source fixture writes ZIPs in an existing test directory.
    weather=tmp_path/"fake_monthly_climate"
    weather.mkdir()
    climate=climate_test.normal(weather)
    result=climate_analysis.evaluate(
        positive,proof,climate,tmp_path/"private_species_site_month.csv",
        tmp_path/"public_monthly_flower_white_TSPVPD.json",
        synthetic=True,min_repeat_cells=2)
    assert result["status"]==climate_analysis.STATUS
    assert result["n_original_FCP_all_photo_opportunities_processed_with_ROI"] is None
    assert result["n_original_FCP_all_photo_opportunities_reconciled"]==6
    assert result["n_original_date_known_photo_ROI_terminals"]==4
    assert result["n_original_bloom_and_source_date_qualified_photo_IDs"]==4
    assert result["n_original_white_bloom_date_known_photos"]==4
    assert result["n_distinct_original_blooming_photo_observation_UUIDs"]==4
    assert result["source_photo_sampling_cap_per_species"] is None
    assert result["n_bloom_photos_no_monthly_worldclim"]==0
    assert result["monthly_weather_is_1970_2000_climatological_not_photo_year"] is True


def test_all_date_eligible_roi_count_incomplete_cannot_be_relabelled_complete(tmp_path):
    audit,shards=prepare_full(tmp_path)
    q=shards/"shard_001_of_002/roi_shard_receipt.json"
    old=json.loads(q.read_text())
    old["n_original_date_known_photo_ROI_terminals"]=1
    q.write_text(json.dumps(old))
    with pytest.raises(ValueError,match="double-positive counts unsealed"):
        m.aggregate_all(audit,shards,tmp_path/"positive.ndjson",
            tmp_path/"public.json",synthetic=True,
            partitions=2,expected_photos=6,expected_taxa=4)
