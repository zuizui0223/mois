"""Pure tests: real iNat API pilot count contract, never network or photo bytes."""
from __future__ import annotations
import importlib.util
import json
from pathlib import Path
from urllib.parse import parse_qs, urlsplit

import pytest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location("mois_real_inat_count_test",
  ROOT/"scripts/analysis/mois_real_inat_phenology_count_bias_pilot_20261016.py")
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def case_count(url):
    q=parse_qs(urlsplit(url).query)
    assert q["taxon_id"]==["47126"]
    assert q["photos"]==["true"]
    assert q["per_page"]==["1"]
    assert q["d1"] and q["d2"]
    if "without_term_id" in q:return {"n_observations":70,"raw_response_sha256":"f"*64}
    if q.get("term_value_id")==["13"]:
        return {"n_observations":12,"raw_response_sha256":"a"*64}
    if q.get("term_value_id")==["21"]:
        return {"n_observations":2,"raw_response_sha256":"b"*64}
    if q.get("term_id")==["12"]:
        return {"n_observations":30,"raw_response_sha256":"c"*64}
    return {"n_observations":100,"raw_response_sha256":"d"*64}


def test_observed_plant_photo_source_and_annotation_not_synonym_for_no_flower(tmp_path):
    slept=[]
    cohorts=m.COHORTS[:2]
    out=tmp_path/"counts.json"
    r=m.run(out,cohorts=cohorts,get_count=case_count,sleep=slept.append)
    assert r["status"]==m.STATUS
    assert r["queries_total"]==10
    assert slept==[m.QUERY_MIN_DELAY_SEC]*9
    assert r["no_iNaturalist_image_pixels_opened"] is True
    assert r["actual_original_source_FCP_42111_taxon_frame_sampled"] is False
    assert r["original_93209351_photo_ID_manifest_used"] is False
    assert r["annotations_assumed_frozen_as_of_2026_09_27"] is False
    a=r["cohorts"]["world_2025_observed"]
    assert a["total_observed_plant_photo_records"]==100
    assert a["with_any_phenology_annotation"]==30
    assert a["explicit_flowers_observations"]==12
    assert a["explicit_no_flowers_fruits_observations"]==2
    assert a["without_phenology_annotation_filter_records"]==70
    assert a["any_plus_without_partitions_exactly"] is True
    assert a["unannotated_is_not_nonflowering"] is True
    raw=json.loads(out.read_text())
    assert raw["schema"]==m.SCHEMA
    assert "raw_observation" not in out.read_text()
    assert "taxon_species_id" not in out.read_text()
    with pytest.raises(ValueError,match="No overwrite"):
        m.run(out,cohorts=cohorts,get_count=case_count,sleep=lambda _:None)


def test_source_without_term_mismatch_is_diagnostic_not_nonflowering():
    z={"all_plant_photo_observations":100,
       "any_flowers_fruits_annotation":30,
       "explicit_flowers_annotation":15,
       "explicit_no_flowers_fruits_annotation":1,
       "without_flowers_fruits_annotation":65}
    result=m.summarize(z)
    assert result["any_plus_without_partitions_exactly"] is False
    assert result["any_plus_without_minus_total"]==-5
    assert result["explicit_flowers_fraction_observed"]==.15
    assert result["unannotated_is_not_nonflowering"] is True
    z["explicit_flowers_annotation"]=31
    with pytest.raises(ValueError,match="larger"):
        m.summarize(z)


def test_api_response_parser_never_retains_individual_observations():
    class FakeResponse:
        status=200
        def __init__(self):self.body=json.dumps({
           "total_results":1004,"results":[{"id":102,"observer":"sensitive","photos":[{"url":"PRIVATE"}]}]}).encode()
        def __enter__(self):return self
        def __exit__(self,*args):pass
        def read(self,n):return self.body
    data=m.query_count(m.uri_query(m.COHORTS[0],{}),
                       opener=lambda *_args,**_kwargs:FakeResponse())
    assert data["n_observations"]==1004
    assert "sensitive" not in json.dumps(data)
    assert set(data)=={"n_observations","raw_response_sha256"}


def test_no_high_volume_or_unsafe_external_query_in_count_pilot():
    assert len(m.COHORTS)*len(m.FILTERS)==35
    assert m.QUERY_MIN_DELAY_SEC>=1
    m.validate_cohorts(m.COHORTS)
    with pytest.raises(ValueError,match="Exceeded"):
        m.validate_cohorts(m.COHORTS+[dict(m.COHORTS[0],name="one_more")])
    with pytest.raises(ValueError,match="Partial geographic"):
        m.validate_cohorts([{"name":"bad","d1":"2025-01-01","d2":"2025-03-01","swlat":30}])
