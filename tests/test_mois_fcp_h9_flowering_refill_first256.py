"""FCP H9 flowering-photo metadata rescue after ROI failure, MOIS taxon-first."""
import csv
import json
from pathlib import Path

import pandas as pd
import pytest

from scripts.acquisition import mois_fcp_h9_flowering_refill_first256 as m


def failed_source(tid="4567",photo="991"):
    return {
        "inat_taxon_id":tid,
        "species":"Campanula punctata",
        "photo_id_source":photo,
        "observation_id":"original-pilot-uuid",
        "measurement_status":"not_evaluable_roi_or_flip_gate",
    }


def observation(pid="55201",tid=4567,lic="cc-by",annotation=True):
    return {
        "id":10025,
        "uuid":"new-h9-observation-uuid",
        "taxon":{"id":tid,"rank":"species","name":"Campanula punctata"},
        "quality_grade":"research",
        "obscured":False,"geoprivacy":None,
        "geojson":{"coordinates":[136.40,37.18]},
        "positional_accuracy":250,
        "observed_on":"2025-07-20",
        "user":{"id":6583},
        "photos":[{"id":int(pid),"url":f"https://static.inaturalist.org/photos/{pid}/square.jpg",
                   "license_code":lic}],
        "annotations":[{"controlled_attribute_id":12,"controlled_value_id":13}]
           if annotation else [],
    }


def test_FCP_H9_params_exactly_flowers_not_fruits_or_buds():
    q=m.fcp_h9_params("4567")
    assert (q["term_id"],q["term_value_id"])==(12,13)
    assert q["per_page"]==200 and q["page"]==1
    assert q["order_by"]=="random" and q["quality_grade"]=="research"
    assert q["photos"]=="true" and q["geo"]=="true"
    assert q["acc_below"]==5000 and q["obscuration"]=="none"
    assert "cc-by" in q["photo_license"]
    with pytest.raises(ValueError,match="numeric"):
        m.fcp_h9_params("not-a-taxon")


def test_H9_licensed_flower_observation_has_original_photo_but_NO_hue():
    item=m.candidate_from_response(observation(),failed_source(),{"991"})
    assert item["photo_id"]=="55201"
    assert item["source_annotation"].endswith("(12:13)")
    assert item["photo_url_large"].endswith("/55201/large.jpg")
    assert item["original_failed_photo_id"]=="991"
    assert "morph" not in item and "white" not in item
    assert item["source_sampling_arm"].startswith("fcp_h9_flowering")


def test_H9_api_source_validation_rejects_bad_geo_species_and_licence():
    o=observation()
    assert m.candidate_from_response(o,failed_source(),{"55201"}) is None
    assert m.candidate_from_response(observation(lic="all-rights-reserved"),
                                     failed_source(),set()) is None
    bad=observation(tid=9)
    assert m.candidate_from_response(bad,failed_source(),set()) is None
    bad=observation()
    bad["positional_accuracy"]=15000
    assert m.candidate_from_response(bad,failed_source(),set()) is None
    bad=observation()
    bad["quality_grade"]="casual"
    assert m.candidate_from_response(bad,failed_source(),set()) is None
    bad=observation()
    bad["obscured"]=True
    assert m.candidate_from_response(bad,failed_source(),set()) is None


def test_outcome_blind_choice_does_not_depend_on_image_hue():
    first=observation(pid="775")
    second=observation(pid="9917")
    rows=[first,second]
    a=m.choose_one_h9_photo(rows,failed_source(),set())
    first["morph"]="red_pink"
    second["flower_colour"]="white"
    b=m.choose_one_h9_photo(rows,failed_source(),set())
    assert a["photo_id"]==b["photo_id"]
    assert a["photo_id"] not in {"991",""}


class DummyClient:
    def __init__(self):
        self.queries=[]
    def observations(self,params):
        self.queries.append(params)
        tid=params["taxon_id"]
        return {"total_results":1,"results":[observation(
             str(100000+tid),tid=tid)]}


def test_source_only_uses_failed_photo_species_once_and_keeps_signature(
        tmp_path,monkeypatch):
    q=pd.DataFrame([
        {**failed_source(tid=str(i),photo=str(i+10)),
         "classified_four_hue":"false","morph":"mixed_uncertain"}
        for i in range(500,503)])
    core=tmp_path/"core";core.mkdir()
    source=tmp_path/"source";source.mkdir()
    pilot=tmp_path/"pilot";pilot.mkdir()
    fifth=tmp_path/"fifth";fifth.mkdir()
    season=tmp_path/"season";season.mkdir()
    pd.DataFrame({"photo_id":["91234","95432"]}).to_csv(
        source/m.broad.POOL,index=False)
    rr={"source_snapshot_date":"2026-10-18",
        "original_broad_all_future_source_pool_SHA256":
        m.sha(source/m.broad.POOL),
        "first_wave_original_4096_csv_sha256":"sha-global"}
    monkeypatch.setattr(m,"load_pilot",lambda *_:(q,{"99"},rr,{
        "source_outcomes_csv_SHA256":"pilot-hue-sha"}))
    monkeypatch.setattr(m.broad,"reservation",lambda *a,**k:(
        set(),{"first_wave_original_4096_csv_sha256":"test-reservation-sha"}))
    client=DummyClient()
    result=m.run(core,pilot,source,fifth,season,tmp_path/"out",
                 client=client,maximum_species=3)
    assert result["rescued_species_API_attempts"]==3
    assert result["rescued_flowering_annotated_new_photo_SOURCE_ONLY_IDs"]==3
    assert len(client.queries)==3
    assert not result["new_photo_pixels_actually_measured"]
    data=pd.read_csv(tmp_path/"out"/m.FROZEN_SOURCE)
    assert len(data)==3
    assert data.inat_taxon_id.nunique()==3
    assert data.photo_id.nunique()==3
    assert "morph" not in data.columns
    assert result["source_new_photo_candidates_SHA256"]==m.sha(
        tmp_path/"out"/m.FROZEN_SOURCE)
