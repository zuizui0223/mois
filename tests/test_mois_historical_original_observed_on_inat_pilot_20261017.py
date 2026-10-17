"""Source blind original-photo selection and current API verification tests."""
import pandas as pd
import pytest
import requests

from scripts.acquisition.mois_historical_original_observed_on_inat_pilot_20261017 import (
    deterministic_blind_sample, normalize_public_observation,
    public_observations_batch, PILOT_N,
)


def original():
    return {"photo_id":"253367174","observation_id":"147261418",
            "inat_taxon_id":"400312"}


def observation(photos=True,flowering=True):
    return {"id":147261418,"observed_on":"2022-10-18",
            "taxon":{"id":400312},
            "photos":[{"id":253367174}, {"id":253367186}] if photos else [{"id":999}],
            "annotations":[{"controlled_value":{"label":"Flowers"}}] if flowering else []}


def test_recover_only_exact_original_photo_in_original_observation():
    r=normalize_public_observation(original(),observation())
    assert r["retrieved_observed_on"]=="2022-10-18"
    assert r["provenance_status"]=="source_photo_and_observation_verified_with_current_public_date"
    assert r["taxon_id_now_matches_source"]
    assert r["has_current_observation_level_flowering_annotation"] is True
    assert r["source_photo_reappeared_in_api"]


def test_observation_date_and_flowering_tag_blocked_without_matching_original_photo():
    r=normalize_public_observation(original(),observation(photos=False))
    assert r["retrieved_observed_on"] is None
    assert r["provenance_status"]=="observation_found_but_original_photo_not_in_current_photo_set"
    assert not r["has_current_observation_level_flowering_annotation"]
    assert normalize_public_observation(original(),None)["retrieved_observed_on"] is None
    with pytest.raises(ValueError,match="wrong observation"):
        normalize_public_observation(original(),{"id":999})


def test_source_taxon_reidentification_is_recorded_not_silently_reassigned():
    obs=observation()
    obs["taxon"]["id"]=777
    r=normalize_public_observation(original(),obs)
    assert r["retrieved_observed_on"]=="2022-10-18"
    assert r["source_taxon_id"]=="400312"
    assert r["current_public_taxon_id"]=="777"
    assert r["taxon_id_now_matches_source"] is False


def test_blind_sampling_selection_uses_only_original_ids_not_colour():
    a=pd.DataFrame([{"photo_id":str(i),"observation_id":str(i*2),
                     "inat_taxon_id":str(i//4),"morph":"white" if i%2 else "red_pink",
                     "observed_on":None} for i in range(250)])
    x=deterministic_blind_sample(a,200)
    b=a.copy()
    b["morph"]=b.morph.map({"white":"red_pink","red_pink":"white"})
    assert deterministic_blind_sample(b,200).photo_id.tolist()==x.photo_id.tolist()
    assert len(x)==PILOT_N
    with pytest.raises(ValueError,match="bounded"):
        deterministic_blind_sample(a,250)


class FakeResponse:
    status_code=200
    def raise_for_status(self): pass
    def json(self): return {"total_results":2,"results":[
        {"id":42, "photos":[{"id":4242}]},
        {"id":43, "photos":[{"id":4343}]}]}


class FakeSession:
    def get(self,url,params,headers,timeout):
        assert params["id"]=="42,43"
        assert params["per_page"]==2
        assert timeout==30
        return FakeResponse()


def test_batch_public_source_filters_observation_ids_and_hashes_original_payload():
    got,meta=public_observations_batch(["42","43"],session=FakeSession())
    assert set(got)=={"42","43"} and meta["returned_observations"]==2
    assert len(meta["original_public_API_response_SHA256"])==64
    with pytest.raises(ValueError,match="over-limit"):
        public_observations_batch([str(i) for i in range(51)],session=FakeSession())
