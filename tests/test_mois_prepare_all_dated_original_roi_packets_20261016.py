"""Uncapped full-source date-FIRST original-photo ROI packet builder tests.

Uses original FCP-style 2026-09-27 photo-ID, date, species, licence and
observation UUID tables, synthetic only. Tests photo source identity and
worker location-blinding; never downloads a photograph or infers blooming.
"""
from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]


def load(name,path):
    s=importlib.util.spec_from_file_location(name,ROOT/path)
    x=importlib.util.module_from_spec(s);s.loader.exec_module(x)
    return x


date_test=load("mois_all_source_date_candidate_fixture_for_packets",
   "tests/test_mois_full_source_dated_flower_candidates_20261016.py")
m=load("mois_full_date_first_packets",
   "scripts/acquisition/mois_prepare_all_dated_original_roi_packets_20261016.py")


def fixture(tmp):
    shard,counts=date_test.example_shard(tmp)
    private=tmp/"original_private_dated.ndjson"
    receipt=tmp/"original_dated_receipt.json"
    date_test.m.audit_shard(shard,counts,private,receipt,synthetic=True)
    return private,receipt


def test_every_original_date_known_photo_dispatched_without_flower_annotation(tmp_path):
    private,receipt=fixture(tmp_path)
    packets=tmp_path/"blind_original_packets"
    r=m.dispatch(private,receipt,packets,packet_size=2,synthetic=True)
    assert r["status"]==m.STATUS
    assert r["n_original_photo_opportunities_in_sealed_shard"]==9
    assert r["n_date_known_original_photo_ids_dispatched_to_roi"]==5
    assert r["n_original_independent_observation_uuids_dispatched"]==5
    assert r["n_bounded_original_roi_packets"]==3
    assert [x["n_source_dated_photos"] for x in r["per_original_packet"]]==[2,2,1]
    assert r["n_original_photo_ids_omitted_despite_valid_source_date"]==0
    assert r["n_original_photos_with_undetermined_flowering_state"]==5
    assert r["n_original_photo_pixels_measured_at_dispatch"]==0
    assert r["all_photo_ids_with_source_date_included_without_sampling_cap"] is True
    assert r["photos_from_source_without_date_not_sent_to_image_evaluator"] is True
    assert r["protected_source_photo_ids_and_coordinates_in_public_receipt"] is False
    mids=set();photoids=set();observationids=set()
    for k in range(3):
        folder=packets/f"packet_{k:06}"
        proof=json.loads((folder/"packet_receipt.json").read_text())
        assert proof["status"]==m.PACKET_STATUS
        assert proof["worker_has_observation_location_species_or_hue"] is False
        assert proof["original_flower_roi_verified"] is False
        assert proof["original_photo_selection_cap"] is None
        assert proof["n_compute_batches_of_max_64_photos"]==1
        with (folder/"worker.csv").open(newline="") as f:
            worker=list(csv.DictReader(f))
        with (folder/"acquisition.csv").open(newline="") as f:
            acquisition=list(csv.DictReader(f))
        with (folder/"private_source.csv").open(newline="") as f:
            source=list(csv.DictReader(f))
        assert set(worker[0])==set(m.WORKER)
        assert "photo_id" not in worker[0]
        assert "species" not in worker[0]
        assert "latitude" not in worker[0]
        assert "observed_on" not in worker[0]
        assert "morph" not in worker[0]
        assert len(worker)==len(acquisition)==len(source)
        for w,a,original in zip(worker,acquisition,source):
            assert w["measurement_id"]==a["measurement_id"]==original["measurement_id"]
            assert original["bloom_state_before_roi"]=="NOT_YET_PHOTO_VERIFIED"
            assert original["observed_on"]
            assert w["measurement_id"]==m.blind_id(original["original_photo_id"])
            assert a["photo_url_large"].startswith(
                "https://inaturalist-open-data.s3.amazonaws.com/photos/")
            mids.add(w["measurement_id"])
            photoids.add(original["original_photo_id"])
            observationids.add(original["original_observation_uuid"])
    assert len(mids)==len(photoids)==len(observationids)==5
    public=json.loads((packets/"full_dated_roi_dispatch_receipt.json").read_text())
    assert "original_photo_id" not in public
    assert "original_observation_uuid" not in public
    assert "latitude" not in public and "longitude" not in public
    assert "photo_url_large" not in public
    assert all(set(z)=={"packet_index","packet_receipt_sha256","n_source_dated_photos"}
               for z in public["per_original_packet"])
    with pytest.raises(ValueError,match="overwrite"):
        m.dispatch(private,receipt,packets,packet_size=2,synthetic=True)


def test_multiple_photos_from_one_original_observation_are_all_dispatched(tmp_path):
    private,receipt=fixture(tmp_path)
    lines=[json.loads(s) for s in private.read_text().splitlines()]
    assert len(lines)==5
    # An original observation can have more than one photo, but the
    # observer, source time, taxon and original site must be identical.
    second=lines[1]
    first=lines[0]
    second["observation_id"]=first["observation_id"]
    second["inat_taxon_id"]=first["inat_taxon_id"]
    second["species"]=first["species"]
    second["observed_on"]=first["observed_on"]
    second["observed_month"]=first["observed_month"]
    second["observed_year"]=first["observed_year"]
    second["latitude"]=first["latitude"]
    second["longitude"]=first["longitude"]
    second["observer_id"]=first["observer_id"]
    second["latitude_domain"]=first["latitude_domain"]
    second["local_season"]=first["local_season"]
    private.write_text("\n".join(json.dumps(row) for row in lines)+"\n")
    d=json.loads(receipt.read_text())
    d["private_full_dated_source_photo_queue_sha256"]=m.sha(private)
    receipt.write_text(json.dumps(d))
    outcome=m.dispatch(private,receipt,tmp_path/"out_packets",packet_size=4,synthetic=True)
    assert outcome["n_date_known_original_photo_ids_dispatched_to_roi"]==5
    assert outcome["n_original_independent_observation_uuids_dispatched"]==4


def test_shard_missing_date_or_spoofed_flower_labels_do_not_enter_worker(tmp_path):
    private,receipt=fixture(tmp_path)
    rows=[json.loads(s) for s in private.read_text().splitlines()]
    rows[0]["original_photo_flowering_state"]="FLOWERING_FROM_FAKE_ANNOTATION"
    private.write_text("\n".join(json.dumps(z) for z in rows)+"\n")
    d=json.loads(receipt.read_text())
    d["private_full_dated_source_photo_queue_sha256"]=m.sha(private)
    receipt.write_text(json.dumps(d))
    with pytest.raises(ValueError,match="flower-unverified"):
        m.dispatch(private,receipt,tmp_path/"bad_packets",packet_size=2,synthetic=True)


def test_reject_tampered_original_snapshot_sha_source_queuing(tmp_path):
    private,receipt=fixture(tmp_path)
    private.write_text(private.read_text().replace("2020-04-16","2020-08-16"))
    with pytest.raises(ValueError,match="source packet"):
        m.dispatch(private,receipt,tmp_path/"bad_packets",packet_size=2,synthetic=True)
    with pytest.raises(ValueError,match="4096"):
        m.dispatch(private,receipt,tmp_path/"bad_packets",packet_size=4097,synthetic=True)



def test_more_than_4096_original_photos_of_one_species_are_not_biologically_capped(tmp_path):
    """4096 is a compute packet boundary; ALL 4097 IDs must be dispatched."""
    private=tmp_path/"source_full_dated.ndjson"
    record={
        "inat_taxon_id":101,"species":"Campanula synthetica",
        "observation_id":"the_same_original_observation_uuid",
        "photo_license":"cc-by",
        "observer_id":"501",
        "observed_on":"2020-05-16",
        "latitude":35.,
        "longitude":139.,
        "observed_month":5,
        "observed_year":2020,
        "source_snapshot_date":"2026-09-27",
        "latitude_domain":"NORTH_EXTRATROPICAL",
        "local_season":"SPRING",
        "original_photo_flowering_state":"NOT_YET_PHOTO_VERIFIED",
        "origin":"ORIGINAL_DATED_FCP_SHARD_NO_PHOTO_PIXELS",
    }
    with private.open("w") as f:
        for i in range(4097):
            pid=str(100000+i)
            one={**record,
                 "photo_id":pid,
                 "photo_url_large":
                   "https://inaturalist-open-data.s3.amazonaws.com/photos/"+pid+"/large.jpg"}
            f.write(json.dumps(one,sort_keys=True)+"\n")
    report={
        "schema":"mois_original_42111_all_dated_photo_candidate_shard_v1",
        "status":"ORIGINAL_DATED_PHOTO_OPPORTUNITIES_BLOOM_UNVERIFIED_ROI_NOT_MEASURED",
        "private_full_dated_source_photo_queue_sha256":m.sha(private),
        "source_photo_cap_per_taxon":None,
        "photo_hue_classification_attempted":False,
        "actual_flower_photographs_verified_from_pixels":0,
        "metadata_date_known_is_not_blooming_photo_confirmation":True,
        "synthetic_test_only":True,
        "shard_number":0,
        "n_original_eligible_source_photographs_examined":4097,
        "n_original_photo_observation_dates_usable":4097,
    }
    receipt=tmp_path/"source_dated_receipt.json"
    receipt.write_text(json.dumps(report))
    result=m.dispatch(private,receipt,tmp_path/"photo_packets",
                      packet_size=4096,synthetic=True)
    assert result["n_date_known_original_photo_ids_dispatched_to_roi"]==4097
    assert result["n_original_independent_observation_uuids_dispatched"]==1
    assert result["n_bounded_original_roi_packets"]==2
    assert [x["n_source_dated_photos"] for x in result["per_original_packet"]]==[4096,1]
    assert result["n_original_photo_ids_omitted_despite_valid_source_date"]==0
    assert result["all_photo_ids_with_source_date_included_without_sampling_cap"] is True
