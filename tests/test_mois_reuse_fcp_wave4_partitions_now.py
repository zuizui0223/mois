"""Audit that ONLY actually finished, source-matched FCP original RGB terminals count."""
import csv
import json
from pathlib import Path
import pytest

from scripts.analysis.mois_reuse_fcp_wave4_completed_roi_partitions import (
    PATTERN, HUES, MIN_COMPLETE_PARTITIONS, EXPECTED_MODEL_PARTITIONS,
    COLS, sealed_terminals, digest
)


def make_partitions(tmp_path:Path,wrong_id=False):
    src,assignments={},{}
    for p in range(MIN_COMPLETE_PARTITIONS):
        mid=f"measurement_{p:02d}"
        src[mid]={
            "measurement_id":mid,"photo_id":str(10000+p),
            "inat_taxon_id":"42","species":"Testus sp.",
            "latitude":"42.0","longitude":"-70","photo_license":"cc-by",
            "observation_id":str(p),"observed_on":"2026-04-01",
        }
        assignments[mid]=(0,p,0)
        row={
            "measurement_id":("unknown" if wrong_id and p==0 else mid),
            "morph":("white" if p%2==0 else "mixed_uncertain"),
            "measurement_status":(
                "classified_four_state_morph" if p%2==0
                else "not_evaluable_roi_or_flip_gate"
            ),
            "image_sha256":"deadbeef" if p%2==0 else "",
        }
        stem=tmp_path/f"b0_partition_s{p:02d}_p00"
        with (stem.with_suffix(".csv")).open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=list(row))
            w.writeheader();w.writerow(row)
        (stem.with_suffix(".json")).write_text(json.dumps({
            "status":"complete_random_photo_first_terminal_partition",
            "semantic_shard":p,"compute_partition":0,
            "terminal_rows":1,"selected_rows":1,
            "classified_rows":int(p%2==0),
            "mixed_uncertain_rows":int(p%2!=0),
            "source_urls_present":False,
            "species_present":False,
            "coordinates_present":False,
        }))
    return src,assignments


def test_exact_original_partition_names_and_scientific_response_boundary():
    assert PATTERN.fullmatch("b0_partition_s13_p03.csv")
    assert PATTERN.fullmatch("b1_partition_s31_p00.json")
    assert not PATTERN.fullmatch("b2_partition_s01_p01.csv")
    assert EXPECTED_MODEL_PARTITIONS==256
    assert "morph" in COLS and "photo_id" in COLS
    assert {"D","W","Ppoly","Fst"}.isdisjoint(COLS)


def test_only_completed_original_rgb_rows_allowed(tmp_path:Path):
    src,assign=make_partitions(tmp_path)
    rows,proof=sealed_terminals(tmp_path,src,assign)
    assert len(rows)==MIN_COMPLETE_PARTITIONS
    assert len(proof)==MIN_COMPLETE_PARTITIONS
    assert sum(r["morph"] in HUES for r in rows)==8
    assert all(r["source_original_measured"]=="fcp_wave4_partial_original_ROI_v4" for r in rows)
    assert all(len(r)==len(COLS) for r in rows)
    assert all(p["terminal_csv_sha256"]==digest(
        tmp_path/f'b0_partition_s{p["semantic_shard"]:02d}_p00.csv')
        for p in proof)


def test_completed_worker_cannot_change_original_photo_identity(tmp_path:Path):
    src,assign=make_partitions(tmp_path,wrong_id=True)
    with pytest.raises(ValueError,match="not selected"):
        sealed_terminals(tmp_path,src,assign)


def test_missing_completed_receipt_cannot_be_mistaken_for_negative_white(tmp_path:Path):
    src,assign=make_partitions(tmp_path)
    (tmp_path/"b0_partition_s00_p00.json").unlink()
    with pytest.raises(ValueError,match="incomplete"):
        sealed_terminals(tmp_path,src,assign)
