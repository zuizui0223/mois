"""Synthetic-only checks for offline import of iNaturalist observation annotations."""
import csv
import importlib.util
import json
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location(
    "mois_annotation_import",ROOT/"scripts/acquisition/mois_annotation_overlay_from_inat_export_20261016.py")
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def make(tmp):
    photo=tmp/"selected.csv"
    with photo.open("w",newline="") as f:
        w=csv.writer(f);w.writerow(["observation_id","inat_taxon_id","photo_id"])
        w.writerows([("u1","1","111"),("u2","2","222"),("u3","3","333")])
    source=tmp/"inat.jsonl"
    source.write_text(
        json.dumps({"uuid":"u1","annotations":[
            {"controlled_attribute_id":12,"controlled_value_id":13,"vote_score":1},
            {"controlled_attribute_id":12,"controlled_value_id":15,"vote_score":2}]})+"\n"+
        json.dumps({"uuid":"u2","annotations":[]})+"\n",encoding="utf8")
    return photo,source


def test_original_annotation_vs_unlooked_up_are_distinct(tmp_path):
    sample,data=make(tmp_path)
    out=tmp_path/"private.csv";receipt=tmp_path/"receipt.json"
    doc=m.build(sample,data,m.sha(data),"2026-10-16",out,receipt,synthetic=True)
    assert doc["n_observations_annotation_LOOKED_UP"]==2
    assert doc["n_observations_NOT_LOOKED_UP"]==1
    assert doc["n_checked_with_accepted_Flowers_value_13"]==1
    rows={r["observation_id"]:r for r in csv.DictReader(out.open())}
    assert rows["u1"]["annotation_value_ids"]=="13|15"
    assert rows["u1"]["lookup_status"]=="LOOKED_UP"
    assert rows["u2"]["lookup_status"]=="LOOKED_UP"
    assert rows["u2"]["annotation_value_ids"]==""
    assert rows["u3"]["lookup_status"]=="NOT_LOOKED_UP"
    assert rows["u3"]["annotation_value_ids"]==""
    with pytest.raises(ValueError,match="No overwrite"):
        m.build(sample,data,m.sha(data),"2026-10-16",out,receipt,synthetic=True)


def test_disputed_annotation_is_not_positive(tmp_path):
    sample,data=make(tmp_path)
    rows=[
        {"uuid":"u1","annotations":[{"controlled_attribute_id":12,
                 "controlled_value_id":13,"vote_score":0}]},
        {"uuid":"u2","annotations":[]},
    ]
    data.write_text("".join(json.dumps(v)+"\n" for v in rows))
    doc=m.build(sample,data,m.sha(data),"2026-10-16",
                tmp_path/"out.csv",tmp_path/"out.json",synthetic=True)
    assert doc["n_checked_with_accepted_Flowers_value_13"]==0


def test_unknown_annotation_and_wrong_source_hash_fail(tmp_path):
    sample,data=make(tmp_path)
    with pytest.raises(ValueError,match="SHA256"):
        m.build(sample,data,"0"*64,"2026-10-16",
                tmp_path/"out.csv",tmp_path/"out.json",synthetic=True)
    data.write_text(json.dumps({"uuid":"u1","annotations":[
        {"controlled_attribute_id":12,"controlled_value_id":999,"vote_score":1}]})+"\n")
    with pytest.raises(ValueError,match="Unrecognized"):
        m.build(sample,data,m.sha(data),"2026-10-16",
                tmp_path/"out.csv",tmp_path/"out.json",synthetic=True)


def test_missing_annotations_field_cannot_be_called_unannotated(tmp_path):
    sample,data=make(tmp_path)
    data.write_text(json.dumps({"uuid":"u1"})+"\n")
    with pytest.raises(ValueError,match="complete annotations"):
        m.build(sample,data,m.sha(data),"2026-10-16",
                tmp_path/"out.csv",tmp_path/"out.json",synthetic=True)
