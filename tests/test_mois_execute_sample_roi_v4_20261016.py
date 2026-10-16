"""Synthetic-only unit tests for Mois subset original ROI-v4 executor.

Uses injected fake backend, no downloads, real FCP model assets or old colours.
"""
import csv
import importlib.util
import json
from pathlib import Path

import pandas as pd
import pytest

ROOT=Path(__file__).resolve().parents[1]


def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

samples=module("mois_test_packets",ROOT/"tests/test_mois_sample_roi_handoff_20261016.py")
m=module("mois_subset_executor",ROOT/"scripts/acquisition/mois_execute_sample_roi_v4_20261016.py")


class FakeBackend:
    class PermanentOriginalPhotoError(ValueError):
        pass
    BIO=("white","yellow","orange","red","pink","magenta","purple","blue","bronze")
    MORPHS=("white","yellow_orange","red_pink","blue_purple")
    @staticmethod
    def verify_frozen_method(source,detector,sam):
        return {"synthetic_test_only":True,"location_blind":True,"no_old_measurement_csv_inputs":True}
    @staticmethod
    def image_loader(url):
        if "/12/" in url:
            raise FakeBackend.PermanentOriginalPhotoError("synthetic 404")
        return b"fake-image-bytes"
    @staticmethod
    def evaluate_with_frozen_source(worker,images,source,detector,sam,out):
        rows=[]
        with worker.open(newline="") as f:
            reader=csv.DictReader(f)
            assert reader.fieldnames == m.WORKER
            for r in reader:
                assert (images/r["image_filename"]).read_bytes()==b"fake-image-bytes"
                rows.append({"measurement_id":r["measurement_id"],
                             "morph":"white","measurement_status":"classified_four_state_morph",
                             **{"flower_fraction_"+x:1.0 if x=="white" else 0.0 for x in FakeBackend.BIO}})
        d=pd.DataFrame(rows)
        d.to_csv(out,index=False)
        return d
    @staticmethod
    def validate_terminal(records, expected):
        assert len(records)==len(expected)
        assert set(records["measurement_id"])==set(expected)
        assert "species" not in records and "latitude" not in records
        for r in records.to_dict("records"):
            assert r["measurement_status"] in ("classified_four_state_morph","image_acquisition_failed")


def create(root):
    sample,selection_receipt,_=samples.setup_sample(root)
    private=root/"packets"
    packets=samples.m.assemble(sample,selection_receipt,private,synthetic=True)
    return private,packets


def run_all(private,out,meta):
    completed=[]
    for part,n in enumerate(meta["partition_photo_counts"]):
        if n:
            done=m.execute(private,out,part,0,64,backend=FakeBackend,
               original_measure=Path("/nonexistent"),original_detector=Path("/nonexistent"),
               efficient_sam=Path("/nonexistent"),synthetic=True)
            completed.append(done)
    return completed


def test_sampled_pixels_are_sealed_without_identity_leak(tmp_path):
    private,packets=create(tmp_path)
    out=tmp_path/"terminals"
    completed=run_all(private,out,packets)
    assert sum(x["batch_photo_count"] for x in completed)==3
    assert sum(x["photo_ids_classified"] for x in completed)==2
    assert sum(x["photo_ids_acquisition_failed"] for x in completed)==1
    for sealed in out.glob("terminal_p*_b*"):
        data=list(csv.DictReader((sealed/"original_roi_sample_terminal.csv").open()))
        assert len(data)>=1
        assert not {"inat_taxon_id","photo_id","observation_id","latitude","longitude","bio12_mm"} & set(data[0])
        assert m.sha256_file(sealed/"original_roi_sample_terminal.csv")==json.loads(
            (sealed/"work_unit_receipt.json").read_text())["sealed_terminal_csv_sha256"]
    for x in completed:
        part=x["partition"]
        with pytest.raises(ValueError,match="Existing sealed"):
            m.execute(private,out,part,0,64,backend=FakeBackend,
               original_measure=Path("/none"),original_detector=Path("/none"),
               efficient_sam=Path("/none"),synthetic=True)


def test_location_blind_packets_reject_unsealed_change(tmp_path):
    private,meta=create(tmp_path)
    selected=next(i for i,n in enumerate(meta["partition_photo_counts"]) if n)
    file=private/f"worker_partition_{selected:03}.csv"
    file.write_text(file.read_text()+"broken line\n")
    with pytest.raises(ValueError,match="source packet changed"):
        m.read_packet(private,selected,0,64)


def test_excessive_budget_missing_provider_gate_and_out_of_range(tmp_path,monkeypatch):
    private,meta=create(tmp_path)
    selected=next(i for i,n in enumerate(meta["partition_photo_counts"]) if n)
    with pytest.raises(ValueError,match="1-64 photo"):
        m.read_packet(private,selected,0,65)
    with pytest.raises(ValueError,match="outside"):
        m.read_packet(private,selected,999,64)
    monkeypatch.delenv("FCP_PROVIDER_GLOBAL_RATE_FILE",raising=False)
    with pytest.raises(ValueError,match="RATE_FILE"):
        m.execute(private,tmp_path/"output",selected,0,64,backend=FakeBackend,
                  original_measure=Path("/none"),original_detector=Path("/none"),
                  efficient_sam=Path("/none"))


def test_provider_transient_aborts_and_preserves_incomplete(tmp_path):
    private,meta=create(tmp_path)
    class NetworkDown(FakeBackend):
        @staticmethod
        def image_loader(url):
            raise ConnectionError("simulated provider outage")
    selected=next(i for i,n in enumerate(meta["partition_photo_counts"]) if n)
    with pytest.raises(ConnectionError,match="outage"):
        m.execute(private,tmp_path/"output",selected,0,64,backend=NetworkDown,
            original_measure=Path("/none"),original_detector=Path("/none"),
            efficient_sam=Path("/none"),synthetic=True)
    assert list((tmp_path/"output").glob("*.INCOMPLETE"))
    assert not list((tmp_path/"output").glob("terminal_p*_b*"))
