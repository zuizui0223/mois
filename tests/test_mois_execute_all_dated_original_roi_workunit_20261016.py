"""Synthetic end-to-end, no-network ROI-v4 source metadata blinding tests.

The 42,111-source taxon/photo eligibility is date-first/no-cap; actual image
classification runs only through original FCP method-compatible interface.
"""
from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path

import pandas as pd
import pytest

ROOT=Path(__file__).resolve().parents[1]


def load(name,path):
    s=importlib.util.spec_from_file_location(name,ROOT/path)
    o=importlib.util.module_from_spec(s);s.loader.exec_module(o)
    return o


test_source=load("mois_original_date_packets_source_test",
    "tests/test_mois_prepare_all_dated_original_roi_packets_20261016.py")
m=load("mois_date_full_fcp_roi_batch",
    "scripts/acquisition/mois_execute_all_dated_original_roi_workunit_20261016.py")
original_date_builder=test_source.m


class FakeFCP:
    BIO=m.BIO
    class PermanentOriginalPhotoError(Exception):pass

    def __init__(self, fail_one=False):
        self.fail_one=fail_one
        self.downloaded=[]
        self.classifier_worker_rows=[]

    def verify_frozen_method(self,*args):
        return {
            "synthetic_test_only":True,
            "location_blind":True,
            "no_old_measurement_csv_inputs":True
        }

    def image_loader(self,url):
        self.downloaded.append(url)
        if self.fail_one and url==self.downloaded[0]:
            raise self.PermanentOriginalPhotoError("source image genuinely missing")
        return b"synthetic bytes not real original photograph"

    def evaluate_with_frozen_source(self,worker,images,source,detector,sam,out):
        with worker.open(newline="") as f:
            r=csv.DictReader(f)
            assert tuple(r.fieldnames)==m.WORKER
            entries=list(r)
        assert all("latitude" not in e and "observed_on" not in e
                   and "inat_taxon_id" not in e and "photo_url_large" not in e
                   and "morph" not in e for e in entries)
        self.classifier_worker_rows.extend(entries)
        measurements=[]
        for item in entries:
            assert (images/item["image_filename"]).read_bytes()==(
                b"synthetic bytes not real original photograph")
            hue="white"
            measurements.append({
                "measurement_id":item["measurement_id"],
                "morph":hue,"measurement_status":"classified_four_state_morph",
                **{"flower_fraction_"+c:(1. if c=="white" else 0.)
                   for c in m.BIO}
            })
        return pd.DataFrame(measurements)

    def validate_terminal(self,rows,expected):
        assert list(rows["measurement_id"]) or expected==[]
        assert len(rows)==len(expected)
        assert set(rows["measurement_id"])==set(expected)
        assert all(
            (r["measurement_status"]=="classified_four_state_morph")
            ==(r["morph"] in {"white","yellow_orange","red_pink","blue_purple"})
            for _,r in rows.iterrows()
        )


def setup(tmp):
    private,receipt=test_source.fixture(tmp)
    packets=tmp/"blind_packets"
    source=original_date_builder.dispatch(
        private,receipt,packets,packet_size=3,synthetic=True)
    assert source["n_date_known_original_photo_ids_dispatched_to_roi"]==5
    return packets,source


def test_date_first_photo_roi_actually_routes_only_hash_IDs_into_fcp_worker(tmp_path):
    packets,source=setup(tmp_path)
    backend=FakeFCP()
    result=m.execute(packets/"packet_000000",packets,
          tmp_path/"sealed_original_roi",0,0,backend=backend,
          source_script=tmp_path/"fake_fcp_method.py",
          original_detector=tmp_path/"fake_detector.pt",
          original_sam=tmp_path/"fake_sam",synthetic=True)
    assert result["status"]==m.STATUS
    assert result["n_original_date_known_unverified_bloom_photos_submitted"]==3
    assert result["n_source_roi_four_colour_classifiable_flower_photos"]==3
    assert result["n_source_roi_white_classifiable_photos"]==3
    assert result["n_source_image_acquisition_failed"]==0
    assert result["original_fcp_image_pixels_saved_after_evaluation"] is False
    assert result["missing_photo_flower_hues_imputed"] is False
    assert len(backend.classifier_worker_rows)==3
    assert len(backend.downloaded)==3
    expected=(tmp_path/"sealed_original_roi/photo_roi_packet_000000_batch_000")
    assert expected.is_dir()
    assert not (expected/"images").exists()
    terminal=list(csv.DictReader(
        (expected/"original_fcp_date_valid_photo_ROI_v4_terminal.csv").open()))
    assert len(terminal)==3 and all(r["morph"]=="white" for r in terminal)
    assert "latitude" not in json.dumps(result)
    with pytest.raises(ValueError,match="overwrite"):
        m.execute(packets/"packet_000000",packets,
          tmp_path/"sealed_original_roi",0,0,backend=backend,
          source_script=tmp_path/"fake_fcp_method.py",
          original_detector=tmp_path/"fake_detector.pt",
          original_sam=tmp_path/"fake_sam",synthetic=True)


def test_source_permanent_failure_does_not_become_false_nonwhite_bloom(tmp_path):
    packets,_=setup(tmp_path)
    backend=FakeFCP(fail_one=True)
    result=m.execute(packets/"packet_000000",packets,
          tmp_path/"sealed_roi",0,0,backend=backend,
          source_script=tmp_path/"source.py",
          original_detector=tmp_path/"detector.pt",
          original_sam=tmp_path/"sam",synthetic=True)
    assert result["n_source_image_acquisition_failed"]==1
    assert result["n_source_roi_four_colour_classifiable_flower_photos"]==2
    assert result["n_source_roi_white_classifiable_photos"]==2
    data=list(csv.DictReader(
        (tmp_path/"sealed_roi/photo_roi_packet_000000_batch_000/original_fcp_date_valid_photo_ROI_v4_terminal.csv").open()))
    nonclass=[x for x in data if x["measurement_status"]=="image_acquisition_failed"]
    assert len(nonclass)==1 and nonclass[0]["morph"]=="mixed_uncertain"
    assert all(nonclass[0]["flower_fraction_"+c] in ("0.0","0.0") for c in m.BIO)


def test_original_bloom_packet_source_identity_and_packet_number_cannot_shift(tmp_path):
    packets,_=setup(tmp_path)
    with pytest.raises(ValueError,match="beyond"):
        m.read_batch(packets/"packet_000000",packets,0,1)
    with pytest.raises(ValueError,match="receipt"):
        m.read_batch(packets/"packet_000001",packets,0,0)
    target=packets/"packet_000000/acquisition.csv"
    target.write_text(target.read_text().replace("cc-by","cc0"))
    with pytest.raises(ValueError,match="bytes drifted"):
        m.read_batch(packets/"packet_000000",packets,0,0)


def test_real_image_worker_requires_source_access_coordination(tmp_path,monkeypatch):
    packets,_=setup(tmp_path)
    monkeypatch.delenv("FCP_PROVIDER_GLOBAL_RATE_FILE",raising=False)
    with pytest.raises(ValueError,match="coordinated global pacing"):
        m.execute(packets/"packet_000000",packets,
           tmp_path/"would_download_real",0,0,backend=FakeFCP(),
           source_script=tmp_path/"original.py",
           original_detector=tmp_path/"original.pt",
           original_sam=tmp_path/"original_sam",synthetic=False)
    assert not (tmp_path/"would_download_real").exists()
