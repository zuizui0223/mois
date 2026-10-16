"""Validate completed original-date-known PHOTO ROI-v4 source shard reassembly.

No real plant images, usernames or external URLs requested; verified
synthetic original source emits valid-day plants, a subset of photo-pixel
colour-classifiable blooming evidence and complete unit source receipts.
"""
from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]


def module(name,path):
    spec=importlib.util.spec_from_file_location(name,ROOT/path)
    v=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(v)
    return v


fake=module("mois_source_date_first_fcp_worker_fixture",
    "tests/test_mois_execute_all_dated_original_roi_workunit_20261016.py")
m=module("mois_shard_actual_date_bloom_terminal_reassemble",
    "scripts/analysis/mois_reassemble_dated_original_roi_shard_20261016.py")


def prepared(tmp_path,*,complete=True,permanent_failure=False):
    packet_root,original=fake.setup(tmp_path)
    seals=tmp_path/"roi_sealed"
    backend=fake.FakeFCP(fail_one=permanent_failure)
    stop=len(original["per_original_packet"]) if complete else 1
    for i in range(stop):
        fake.m.execute(
            packet_root/f"packet_{i:06d}",packet_root,
            seals,i,0,backend=backend,
            source_script=tmp_path/"original_method.py",
            original_detector=tmp_path/"model.pt",
            original_sam=tmp_path/"sam",synthetic=True)
    return packet_root,seals,original


def test_all_source_original_dated_photos_must_be_measured_before_complete_shard(tmp_path):
    source,seal,original=prepared(tmp_path,complete=True)
    positive=tmp_path/"private_original_bloom.ndjson"
    report=tmp_path/"public_shard.json"
    r=m.reassemble(source,seal,positive,report,synthetic=True)
    assert r["status"]==m.STATUS
    assert r["original_photo_opportunities_in_source_shard"]==9
    assert r["n_original_date_known_photos_dispatched"]==5
    assert r["n_original_date_known_photo_ROI_terminals"]==5
    assert r["n_original_bloom_four_hue_classified_photo_pixels"]==5
    assert r["n_original_bloom_white_four_hue_classified_photo_pixels"]==5
    assert r["n_original_distinct_observation_uuids_with_classifiable_flower"]==5
    assert r["n_not_date_known_source_photo_opportunities_not_image_processed"]==4
    assert r["n_finalized_original_photo_roi_packets"]==2
    assert r["old_100543_photo_outcomes_ever_reused"] is False
    assert r["observation_level_iNaturalist_Flowers_annotations_used_as_photo_positives"] is False
    assert r["source_photo_per_species_sample_limit"] is None
    assert r["this_shard_alone_is_NOT_complete_global_42111_source"] is True
    assert r["actual_all_128_original_photo_shards_complete"] is False
    assert r["synthetic_test_only"] is True
    blooms=[json.loads(line) for line in positive.read_text().splitlines()]
    assert len(blooms)==5
    assert len({r["original_photo_id"] for r in blooms})==5
    assert all(r["original_ROI_v4_hue"]=="white" for r in blooms)
    assert all(r["actual_photo_contains_classifiable_flower"] is True for r in blooms)
    assert all(r["original_snapshot_20260927_date_proven"] is True for r in blooms)
    month_path=tmp_path/"original_roi_month_effort_private.json"
    assert month_path.is_file()
    month_data=json.loads(month_path.read_text(encoding="utf8"))
    assert month_data["schema"]=="mois_private_original_source_photo_month_roi_effort_v1"
    assert month_data["n_original_dated_photo_roi_terminals"]==5
    assert month_data["n_original_roi_four_hue_classifiable"]==5
    assert month_data["n_original_roi_white"]==5
    assert sum(x["n_original_date_known_photo_roi_terminals"]
               for x in month_data["latitude_domain_month_photo_effort"])==5
    assert r["private_all_month_roi_effort_sha256"]==m.sha(month_path)
    p=json.loads(report.read_text())
    assert "original_photo_id" not in p
    assert "original_observation_uuid" not in p
    assert "latitude" not in p and "longitude" not in p
    with pytest.raises(ValueError,match="overwrite"):
        m.reassemble(source,seal,positive,report,synthetic=True)


def test_missing_one_original_packet_must_not_count_any_complete_shard(tmp_path):
    packet_root,terminals,_=prepared(tmp_path,complete=False)
    public=tmp_path/"shard_should_not_exist.json"
    private=tmp_path/"positive_not_published.ndjson"
    with pytest.raises(ValueError,match="unmeasured ROI"):
        m.reassemble(packet_root,terminals,private,public,synthetic=True)
    assert not public.exists()
    assert not private.exists()
    assert private.with_name(private.name+".INCOMPLETE").exists()


def test_permanent_unavailable_photo_counts_as_unclassifiable_not_fake_nonwhite(tmp_path):
    packets,seals,_=prepared(tmp_path,permanent_failure=True)
    result=m.reassemble(packets,seals,tmp_path/"positive.ndjson",
                        tmp_path/"summary.json",synthetic=True)
    # The explicitly synthetic backend fails only its FIRST URL globally,
    # even though the source contains two independent compute packets.
    assert result["n_original_date_known_photo_ROI_terminals"]==5
    assert result["n_original_permanent_photo_acquisition_failures"]==1
    assert result["n_original_bloom_four_hue_classified_photo_pixels"]==4
    assert result["n_original_bloom_white_four_hue_classified_photo_pixels"]==4
    assert result["n_original_image_roi_other_unclassifiable"]==0
    assert len((tmp_path/"positive.ndjson").read_text().splitlines())==4


def test_altered_photo_terminal_result_cannot_seal_valid_original_blooms(tmp_path):
    packet_root,terminals,_=prepared(tmp_path)
    photo=terminals/"photo_roi_packet_000000_batch_000/original_fcp_date_valid_photo_ROI_v4_terminal.csv"
    photo.write_text(photo.read_text().replace("white","red_pink"))
    with pytest.raises(ValueError,match="not sealed"):
        m.reassemble(packet_root,terminals,tmp_path/"bloom.ndjson",
                     tmp_path/"report.json",synthetic=True)
    assert not (tmp_path/"report.json").exists()
