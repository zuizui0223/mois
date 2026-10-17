"""Synthetic contract tests. Real MOIS photo counts are checked by source replay CI."""
import csv
from pathlib import Path
import pandas as pd
import pytest
from scripts.analysis.mois_scale_original_precolour_4096 import (
    measurement, SOURCE_SHA, ORIGINAL_SELECTED, ELIGIBLE_NEW,
    WAVE_PHOTOS, SHARDS, PHOTOS_PER_SHARD
)
from scripts.analysis.mois_measure_original_new_photos_roi_v4 import (
    PhotoSource,validate
)


def test_original_sample_is_full_photo_source_opportunity_not_colours():
    assert ORIGINAL_SELECTED == 106376
    assert ELIGIBLE_NEW == 104904
    assert WAVE_PHOTOS==SHARDS*PHOTOS_PER_SHARD==4096
    assert len(SOURCE_SHA)==64
    assert measurement("123") == measurement("123")
    assert measurement("123") != measurement("456")
    assert measurement("123").startswith("MOIS-")


def test_fixed_64_photo_blind_worker_uses_no_geography_or_morph():
    rows=[]
    for k in range(256):
        rows.append({
          "measurement_id":measurement(str(k)),"image_filename":measurement(str(k))+".jpg",
          "photo_license":"cc-by","photo_url_large":
            f"https://inaturalist-open-data.s3.amazonaws.com/photos/{k}/large.jpg",
          "latitude":"30","longitude":"100","inat_taxon_id":"123","observer_id":str(k),
          "species":"Testus flower","photo_id":str(k),"observation_id":str(k)
        })
    p=PhotoSource(rows,5)
    for batch in range(4):
        worker,acq,num=p.get_batch(5,batch,64)
        assert num==4 and len(worker)==len(acq)==64
        assert tuple(worker)==("measurement_id","image_filename","photo_license")
        assert "latitude" not in worker and "species" not in worker
        assert "photo_url_large" not in worker and "morph" not in worker
        assert "photo_url_large" in acq
    with pytest.raises(ValueError):
        p.get_batch(5,4,64)


def test_mutated_original_source_not_accepted(tmp_path):
    meta=tmp_path/"source.csv"
    meta.write_text("photo_id,morph,shard\n1,white,0\n")
    receipt=tmp_path/"source.json"
    receipt.write_text('{"status":"FAKE_SOURCE"}')
    with pytest.raises((KeyError,ValueError)):
        validate(meta,receipt,0)
