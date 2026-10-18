"""No sample-depth drift when admitting global breadth source real original RGB."""
import csv
import hashlib
import json
import pandas as pd
import pytest
from scripts.analysis import mois_measure_broad_angiosperm_original_photo_roi_v4 as worker
from scripts.analysis import mois_integrate_broad_angiosperm_real4096 as agg
from scripts.acquisition import mois_global_angiosperm_breadth_source as broad
from scripts.analysis.mois_reassemble_original_new4096_measurements import (
    HUES,METHOD,MEASURED_SCHEMA,
)


def make_plan(n=256):
    rows=[]
    for i in range(n):
        mid="MOIS-"+str(i).zfill(24)
        rows.append({
            "photo_id":str(10000+i),"inat_taxon_id":str(400+i),
            "species":"Novelus herb"+str(i),
            "measurement_id":mid,
            "image_filename":mid+".jpg",
            "observer_id":str(20000+i),
            "observation_id":"uuid"+str(i),
            "source_snapshot_date":"2026-10-18",
            "photo_license":"cc-by",
            "photo_url_large":
                 f"https://inaturalist-open-data.s3.amazonaws.com/photos/{10000+i}/large.jpg",
            "shard":"0",
            "latitude":"35.1","longitude":"132.9",
        })
    return rows


def test_broad_photo_worker_accepts_pinned_dynamic_snapshot_not_old_FCP_20260927(tmp_path):
    source=tmp_path/"source.csv"
    rows=make_plan()
    with source.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)
    digest=hashlib.sha256(source.read_bytes()).hexdigest()
    plan=tmp_path/"source.json"
    plan.write_text(json.dumps({
        "schema":"mois_open_global_angiosperm_breadth_first4096_v1",
        "status":"SOURCE_ONLY_4096_RESERVED_NO_NEW_PIXEL_HUES",
        "source_snapshot_date":"2026-10-18",
        "first_wave_original_source_photo_ids":4096,
        "no_photo_colour_or_fitness_or_bloom_label_used":True,
        "first_wave_original_4096_csv_sha256":"broad-test-first4096",
        "per_shard_sha256":{"0":digest},
    }))
    selection,r=worker.validate(source,plan,0)
    assert len(selection)==256
    assert r["source_snapshot_date"]=="2026-10-18"
    r["source_snapshot_date"]="2026-09-27"
    plan.write_text(json.dumps(r))
    with pytest.raises(ValueError,match="Original scientific source"):
        worker.validate(source,plan,0)


def test_all_classifier_four_hues_and_source_model_method_remain_frozen():
    assert HUES==("white","yellow_orange","red_pink","blue_purple")
    assert MEASURED_SCHEMA=="mois_256_source_original_roi_v4_four_hue_terminal_v1"
    assert METHOD["location_blind"] is True
    assert "model" not in broad.COLUMNS
    assert "morph" not in broad.COLUMNS
    assert agg.NEW==4096
    assert agg.EXPECTED==313535


def test_breadth_source_coverage_not_deep_FCP_photo_volume():
    assert broad.FIRST_WAVE==4096
    assert broad.N_SHARDS==16
    assert broad.PER_SHARD==256
    assert "1000taxon" not in worker.__doc__.lower()
    assert "within-species" not in agg.__doc__.lower()


def test_old_cohort_cannot_be_overwritten_and_source_only_photo_not_labelled(tmp_path):
    r={
        "schema":"mois_open_global_angiosperm_breadth_first4096_v1",
        "status":"SOURCE_ONLY_4096_RESERVED_NO_NEW_PIXEL_HUES",
        "first_wave_original_source_photo_ids":4096,
        "first_wave_distinct_photo_ids":4096,
        "new_photo_pixels_actually_measured":False,
        "first_wave_original_4096_csv_sha256":"mismatch",
    }
    (tmp_path/broad.RECEIPT).write_text(json.dumps(r))
    (tmp_path/broad.FIRST).write_text("photo_id\n123\n")
    with pytest.raises(ValueError,match="Source-only global"):
        agg.verified_blind_plan(tmp_path)
