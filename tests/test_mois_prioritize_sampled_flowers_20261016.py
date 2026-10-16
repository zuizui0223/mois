"""Synthetic-only annotation-first ROI scheduling; never changes photo sample membership."""
import csv
import importlib.util
import json
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]

def load(name,path):
    spec=importlib.util.spec_from_file_location(name,ROOT/path)
    obj=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj

packet_test=load("mois_prioritize_fixture","tests/test_mois_sample_roi_handoff_20261016.py")
m=load("mois_prioritize","scripts/acquisition/mois_prioritize_sampled_flowers_20261016.py")


def fixture(tmp_path):
    sample, receipt, _ = packet_test.setup_sample(tmp_path)
    origin=tmp_path/"original_packet"
    packet_test.m.assemble(sample,receipt,origin,synthetic=True)
    overlay=tmp_path/"annotation.csv"
    with overlay.open("w",newline="") as f:
        w=csv.writer(f)
        w.writerow(["observation_id","annotation_term_id","annotation_value_ids",
                    "lookup_status","annotation_source","lookup_snapshot_date"])
        w.writerows([
            ["o1","12","13","LOOKED_UP","inat_api_v1","2026-10-16"],
            ["o2","12","","LOOKED_UP","inat_api_v1","2026-10-16"],
            ["o9","12","21","LOOKED_UP","inat_api_v1","2026-10-16"],
        ])
    return origin,overlay


def photo_ids(path):
    out=[]
    for p in range(128):
        with (path/f"worker_partition_{p:03}.csv").open(newline="") as f:
            out.extend(row["measurement_id"] for row in csv.DictReader(f))
    return sorted(out)


def test_flowers_first_within_partition_and_global_batch_schedule(tmp_path):
    source,overlay=fixture(tmp_path)
    output=tmp_path/"prioritized"
    z=m.reorder(source,overlay,m.sha256_file(overlay),output)
    assert z["selected_photo_ids"]==3
    assert z["phenology_priority_counts"]["FLOWERS"]==1
    assert z["phenology_priority_counts"]["NO_ANNOTATION"]==1
    assert z["phenology_priority_counts"]["NO_FLOWERS_FRUITS"]==1
    assert z["photographs_dropped_or_added"]==0
    assert z["sampling_inclusion_probabilities_changed"] is False
    assert z["original_pixels_measured"] is False
    assert z["file_sha256"]["private_metadata_join.csv"] == m.sha256_file(
        source/"private_metadata_join.csv")
    assert photo_ids(source)==photo_ids(output)
    batches=list(csv.DictReader((output/"private_phenology_roi_work_order.csv").open()))
    assert sum(int(r["n_selected"]) for r in batches)==3
    assert batches[0]["first_priority"]=="0"
    assert z["phenology_first_batch_schedule_preserves_all_sample_units"] is True
    # Existing location-blind subset executor can consume these resealed
    # packet hashes without opening any pixels.
    original=load("mois_check_original_packet","scripts/acquisition/mois_execute_sample_roi_v4_20261016.py")
    for p,n in enumerate(z["partition_photo_counts"]):
        if n:
            rc,workers,acq=original.read_packet(output,p,0,64)
            assert len(workers)==n and len(acq)==n
            assert set(workers[0])==set(m.WORKER)
    with pytest.raises(ValueError,match="Refuse"):
        m.reorder(source,overlay,m.sha256_file(overlay),output)


def test_partial_annotation_crosswalk_not_promoted_to_no_flower(tmp_path):
    source,overlay=fixture(tmp_path)
    rows=list(csv.DictReader(overlay.open()))
    with overlay.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows[:-1])
    with pytest.raises(ValueError,match="unverified original observation"):
        m.reorder(source,overlay,m.sha256_file(overlay),tmp_path/"prioritized")


def test_unverified_annotations_still_remain_in_complete_work_queue(tmp_path):
    source,overlay=fixture(tmp_path)
    rows=list(csv.DictReader(overlay.open()))
    rows[2]["annotation_value_ids"]=""
    rows[2]["lookup_status"]="NOT_LOOKED_UP"
    with overlay.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    report=m.reorder(source,overlay,m.sha256_file(overlay),tmp_path/"prioritized")
    assert report["phenology_priority_counts"]["NOT_VERIFIED"]==1
    assert photo_ids(source)==photo_ids(tmp_path/"prioritized")
