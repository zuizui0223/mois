"""Synthetic anti-drift tests for the actual 295017 original photo maxT re-fit.

Fixtures are a deliberately small source and DO NOT represent ecological
observations; full authentic source and water numbers are computed in Actions.
"""
import gzip
import hashlib
import io
import json
from pathlib import Path

import pandas as pd
import pytest

from scripts.analysis import mois_refit_original_water_295017_actual_3wave_FCP180 as m


def gzip_source(path:Path,data:pd.DataFrame)->str:
    raw=data.to_csv(index=False,lineterminator="\n").encode()
    with path.open("wb") as f:
        with gzip.GzipFile(fileobj=f,mode="wb",filename="",mtime=0) as gz:
            gz.write(raw)
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fixtures(monkeypatch,tmp_path):
    old=pd.DataFrame([
        {"photo_id":"1","inat_taxon_id":"101","species":"Alpha beta",
         "morph":"white","original_source":"historical"},
        {"photo_id":"2","inat_taxon_id":"102","species":"Beta gamma",
         "morph":"red_pink","original_source":"new_fcp"},
    ])
    current=pd.DataFrame([
        {"photo_id":"1","inat_taxon_id":"101","species":"Alpha beta",
         "morph":"white","source_original_measured":"historical",
         "observed_on":"","latitude":35,"longitude":140},
        {"photo_id":"2","inat_taxon_id":"102","species":"Beta gamma",
         "morph":"red_pink","source_original_measured":"fcp_wave4",
         "observed_on":"2025-05-20","latitude":35,"longitude":140},
        {"photo_id":"3","inat_taxon_id":"103","species":"Gamma delta",
         "morph":"yellow_orange","source_original_measured":"mois_wave3_new_original_ROI_v4",
         "observed_on":"2025-07-20","latitude":36,"longitude":141},
    ])
    o=tmp_path/"old.csv.gz"
    p=tmp_path/"new.csv.gz"
    so=gzip_source(o,old)
    sp=gzip_source(p,current)
    receipt=tmp_path/"full_receipt.json"
    receipt.write_text(json.dumps({
        "schema":"mois_295017_genuine_all_original_photos_three4096_waves_FCP180_date_enriched_v1",
        "original_measured_photographs":3,
        "original_measured_photo_four_hue_labels":3,
        "files_SHA256":{m.PHOTO_NAME:sp},
        "source_shared_conflicts":{"inat_taxon_id_source_conflicts":0,
                                   "morph_source_conflicts":0,
                                   "observed_on_source_conflicts":0}}))
    old_receipt=tmp_path/"old_receipt.json"
    old_receipt.write_text(json.dumps({
        "n_original_photo_ids":2,"n_four_hue_classified":2,
        "full_photo_static_monthly_normal_csv_sha256":so}))
    for key,v in {"OLD_N":2,"OLD_HUES":2,"NEW_N":1,"NEW_HUES":1,
                  "FULL_N":3,"FULL_HUES":3,"FULL_SHA":sp,
                  "OLD_ENV_SHA":so}.items():
        monkeypatch.setattr(m,key,v)
    return o,p,old_receipt,receipt


def test_all_true_original_photo_ids_retained_with_two_finished_hues(monkeypatch,tmp_path):
    old,photo,oldproof,photoproof=fixtures(monkeypatch,tmp_path)
    r=m.verify(photo,photoproof,old,oldproof)
    assert r["all_original_unique_photos"]==3
    assert r["all_four_hue_classified_photos"]==3
    assert r["new_original_source_photos"]==1
    assert r["new_four_hue_measured_photo_ids"]==1
    assert r["original_model_SHA_verified_for_every_photo"]


def test_any_prior_source_photo_colour_reclassification_fails(monkeypatch,tmp_path):
    old,photo,oldproof,photoproof=fixtures(monkeypatch,tmp_path)
    d=pd.read_csv(photo,compression="gzip",dtype={"photo_id":str})
    d.loc[d.photo_id=="1","morph"]="blue_purple"
    digest=gzip_source(photo,d)
    monkeypatch.setattr(m,"FULL_SHA",digest)
    raw=json.loads(photoproof.read_text())
    raw["files_SHA256"][m.PHOTO_NAME]=digest
    photoproof.write_text(json.dumps(raw))
    with pytest.raises(ValueError,match="morph"):
        m.verify(photo,photoproof,old,oldproof)


def test_unknown_new_original_photo_not_counted_as_a_classified_hue(monkeypatch,tmp_path):
    old,photo,oldproof,photoproof=fixtures(monkeypatch,tmp_path)
    d=pd.read_csv(photo,compression="gzip",dtype={"photo_id":str})
    d.loc[d.photo_id=="3","morph"]="mixed_uncertain"
    digest=gzip_source(photo,d)
    monkeypatch.setattr(m,"FULL_SHA",digest)
    photoproof.write_text(json.dumps({
        **json.loads(photoproof.read_text()),
        "files_SHA256":{m.PHOTO_NAME:digest}}))
    with pytest.raises(ValueError,match="denominators|pixel"):
        m.verify(photo,photoproof,old,oldproof)


def test_a_foreign_new_photo_source_is_not_pretended_measured_mois_wave(monkeypatch,tmp_path):
    old,photo,oldproof,photoproof=fixtures(monkeypatch,tmp_path)
    d=pd.read_csv(photo,compression="gzip",dtype={"photo_id":str})
    d.loc[d.photo_id=="3","source_original_measured"]="unknown_future_imputed"
    digest=gzip_source(photo,d)
    monkeypatch.setattr(m,"FULL_SHA",digest)
    photoproof.write_text(json.dumps({
        **json.loads(photoproof.read_text()),
        "files_SHA256":{m.PHOTO_NAME:digest}}))
    with pytest.raises(ValueError,match="8191|photo"):
        m.verify(photo,photoproof,old,oldproof)
