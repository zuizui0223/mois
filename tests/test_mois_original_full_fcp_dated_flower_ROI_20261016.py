"""Synthetic end-to-end FCP dated original 42k all-photo flowering ROI gate.

No real botanical pixels or 93M-run claimed. These tests confirm exact date
and per-image four-colour status must BOTH be positive before a record
contributes to the seasonal observed-white denominator.
"""
from __future__ import annotations
import csv
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location(
    "mois_full_dated_42k_blooms",
    ROOT/"scripts/analysis/mois_original_full_fcp_dated_flower_ROI_20261016.py")
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def fixture(tmp):
    photo=tmp/"synthetic_newly_remeasured_photos.csv"
    model=tmp/"synthetic_method_attestation.json"
    model.write_text(json.dumps({
        "location_blind":True,"synthetic_only":True,
        "source_commit":"NOT_FCP_ORIGINAL"
    }))
    rows=[]
    for i in range(40):
        if i<14:date_="2020-04-10";lat=35;lon=130
        elif i<28:date_="2020-07-20";lat=35;lon=130
        else:date_="2020-12-12";lat=-32;lon=150
        if i in (1,15,29):date_=""  # must fail month classification
        if i==4:date_="2026-12-19"  # after frozen original snapshot
        if i==5:date_="2020-02-30"  # impossible
        if i%3==0:
            hue="mixed_uncertain"
            status="not_evaluable_roi_or_flip_gate"
        else:
            hue="white" if i%5==0 else "blue_purple"
            status="classified_four_state_morph"
        value=dict.fromkeys([f"flower_fraction_{c}" for c in m.BIO],0.0)
        if status=="classified_four_state_morph":
            value["flower_fraction_white" if hue=="white" else
                  "flower_fraction_purple"]=1.0
        rows.append({
            "inat_taxon_id":str(101+i%4),"species":f"S_{i%4}",
            "photo_id":str(i+1000),"observation_id":f"synthetic_{i}",
            "observer_id":str(i%6+200),"observed_on":date_,
            "latitude":str(lat),"longitude":str(lon),
            "morph":hue,"measurement_status":status,
            "global_classifiable":str(status=="classified_four_state_morph").lower(),
            **value,
        })
    with photo.open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)
    receipt=tmp/"complete_fcp_synthetic.json"
    nclass=sum(r["measurement_status"]=="classified_four_state_morph" for r in rows)
    receipt.write_text(json.dumps({
        "schema":"fcp_fresh_complete_remeasurement_v1",
        "status":m.SOURCE_STATUS,
        "study_id":"fcp_fresh_uncapped_20261014_v1",
        "n_frozen_licensed_photo_ids":len(rows),
        "n_terminal_measurements":len(rows),
        "n_classifiable":nclass,
        "n_unclassifiable":len(rows)-nclass,
        "n_source_taxa":4,
        "historical_measurements_reused":False,
        "source_firewall_sha256":"f"*64,
        "new_photo_table_sha256":m.sha(photo),
        "method_attestation_sha256":m.sha(model),
    }))
    return photo,receipt,model,rows


def test_final_flowers_require_real_source_date_and_classifiable_ROI(tmp_path):
    photo,receipt,method,rows=fixture(tmp_path)
    out=tmp_path/"private_positive_original_photo_blooms.ndjson"
    public=tmp_path/"public_full_month.json"
    d=m.analyze(photo,receipt,method,out,public,synthetic=True)
    assert d["schema"]==m.SCHEMA
    assert d["n_original_FCP_all_source_photo_ROI_results"]==40
    assert d["n_original_ROI_classifiable_flower_hues_all_dates"]==sum(
        r["measurement_status"]=="classified_four_state_morph" for r in rows)
    qualified=[r for r in rows if r["measurement_status"]=="classified_four_state_morph"
               and m.source_observation_day(r["observed_on"])[0] is not None]
    assert d["n_original_date_and_site_known_ROI_bloom_hues"]==len(qualified)
    assert d["n_original_date_known_white_flower_photos"]==sum(
        r["morph"]=="white" for r in qualified)
    actual=[json.loads(r) for r in out.read_text().splitlines()]
    assert len(actual)==len(qualified)
    assert all(r["actual_photo_contains_classifiable_flower"] is True for r in actual)
    assert all(r["original_snapshot_20260927_date_proven"] is True for r in actual)
    assert not any(r["original_ROI_v4_hue"]=="mixed_uncertain" for r in actual)
    assert all("API_" not in r for r in actual)
    assert d["full_42111_original_source_not_selectively_date_sampled"] is True
    assert d["per_species_source_photo_cap"] is None
    assert d["unannotated_original_photo_rejected_as_nonflower_only"] is False
    assert d["unclassifiable_photographs_assigned_fake_nonwhite"] is False
    assert d["WorldClim_monthly_T_S_P_VPD_source_site_link_not_yet_run"] is True
    assert d["synthetic_test_only"] is True
    assert "original_photo_id" not in public.read_text()
    assert "original_observation_uuid" not in public.read_text()
    for group in d["species_equal_and_photo_weighted_calendar_month_source_white"]:
        if "public_suppressed_below" in group:continue
        assert group["n_original_roi_white"]<=group["n_original_photo_ROI_classifiable_bloom"]
        assert group["n_original_photo_ROI_classifiable_bloom"]<=group[
            "n_original_photos_with_original_date_and_location"]
        assert 0<=group["white_fraction_given_actual_photo_bloom_ROI"]<=1
        if group["n_observed_named_taxa_with_classifiable_flower"]:
            assert 0<=group["species_equal_white_fraction_given_actual_photo_bloom_ROI"]<=1
    with pytest.raises(ValueError,match="overwrite"):
        m.analyze(photo,receipt,method,out,public,synthetic=True)


def test_non_complete_original_source_and_source_method_mutation_are_refused(tmp_path):
    photo,receipt,method,rows=fixture(tmp_path)
    with pytest.raises(ValueError,match="uncapped"):
        m.verify_complete(photo,receipt,method,synthetic=False)
    bad=json.loads(receipt.read_text())
    bad["status"]="FAKE_COMPUTED_93M_HUE_STATUS"
    receipt.write_text(json.dumps(bad))
    with pytest.raises(ValueError,match="fully sealed"):
        m.analyze(photo,receipt,method,tmp_path/"a.ndjson",
                  tmp_path/"b.json",synthetic=True)


def test_ambiguous_four_hue_and_incorrect_nine_palette_bloom_fails(tmp_path):
    photo,receipt,method,rows=fixture(tmp_path)
    rows[0]["morph"]="white"  # not classifiable: never silently convert!
    with photo.open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)
    doc=json.loads(receipt.read_text());doc["new_photo_table_sha256"]=m.sha(photo)
    receipt.write_text(json.dumps(doc))
    with pytest.raises(ValueError,match="classification outcome"):
        m.analyze(photo,receipt,method,tmp_path/"a.ndjson",
                  tmp_path/"b.json",synthetic=True)
