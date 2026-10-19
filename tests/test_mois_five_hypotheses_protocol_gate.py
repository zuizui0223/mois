"""Evidence-state and source-disjointness safeguards for MOIS future photo tests.

These synthetic observations are NOT ecological or floral colour data.
"""
import copy
import csv
import hashlib
import json
from pathlib import Path

import pytest

from scripts.analysis.mois_five_hypotheses_protocol_gate import (
    FRAME_PHOTOS, REQUIRED_HYPOTHESES, gate, sha256, validate_protocol,
)

ROOT=Path(__file__).resolve().parents[1]
SPEC=ROOT/"analysis_specs"/"mois_five_global_photo_hypotheses_future_holdout_v1.json"


def test_real_five_hypotheses_protocol_explicitly_exploratory_previously():
    p=validate_protocol(SPEC)
    assert p["photo_universe"]["photo_candidates"]==FRAME_PHOTOS
    assert len(p["directional_hypotheses"])==5
    assert tuple(x["id"] for x in p["directional_hypotheses"])==REQUIRED_HYPOTHESES
    assert all(x["previously_explored"] for x in p["directional_hypotheses"])
    assert not p["external_source_gates"]["WCVP"]["same_island_native_assertion_permitted"]
    assert not p["external_source_gates"]["bee_richness"]["visitation_causality_claim_permitted"]


def test_cannot_change_historical_discovery_to_preregistered():
    p=json.loads(SPEC.read_text())
    p["status"]="already_preregistered"
    with pytest.raises(ValueError,match="retrospective"):
        validate_protocol(p)
    p=json.loads(SPEC.read_text())
    p["directional_hypotheses"][2]["previously_explored"]=False
    with pytest.raises(ValueError,match="pre-existing"):
        validate_protocol(p)


def write_json(path, d):
    path.write_text(json.dumps(d,indent=2,sort_keys=True)+"\n")


def write_csv(path, records, fields=("photo_id","observation_id")):
    with path.open("w",newline="",encoding="utf8") as f:
        writer=csv.DictWriter(f,fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)


def fixtures(tmp_path):
    # Small finished mock source index to test behavior without 93 million IDs.
    index_path=tmp_path/"index.json"
    queue_path=tmp_path/"queue.json"
    prior_path=tmp_path/"prior.csv"
    future_path=tmp_path/"future.csv"
    prior_receipt=tmp_path/"prior.json"
    index={
        "schema":"mois_full_42111_93209351_uncapped_original_photo_observation_index_v1",
        "photos_emitted":5,"taxa_in_frame":2,
        "photo_caps":None,"per_observation_photo_cap":None,
        "photo_colours_measured":False,
        "shards":[{"photo_count":2,"observation_count":2,
                   "photos_sha256":"a"*64,"observations_sha256":"b"*64},
                  {"photo_count":3,"observation_count":2,
                   "photos_sha256":"c"*64,"observations_sha256":"d"*64}]
    }
    write_json(index_path,index)
    queue={
        "schema":"mois_full_source_observation_sampling_cohort_v1",
        "full_frame_photo_ids":5,
        "full_source_manifest_sha256":sha256(index_path),
        "selection_unit":"observation_uuid",
        "selection_probability":.0002,
        "sampled_before_any_hue_measurement":True,
        "colour_measured":False,
        "source_original_images_measured":False,
        "all_unselected_source_photos_eligible_at_higher_probability":True,
        "no_per_taxon_or_observation_photo_count_cap":True,
        "new_photo_attempts_planned":2,
        "selected_observations_with_new_photos":1,
        "ecological_analysis_run":False
    }
    write_json(queue_path,queue)
    write_csv(prior_path,[{"photo_id":"1","observation_id":"OLD"},
                          {"photo_id":"2","observation_id":"OLD"},
                          {"photo_id":"3","observation_id":"OLD2"}])
    write_csv(future_path,[{"photo_id":"4","observation_id":"NEW"},
                           {"photo_id":"5","observation_id":"NEW"}])
    receipt={
        "schema":"mois_prior_exposed_photo_and_observation_roster_v1",
        "includes_all_mois_original_attempts":True,
        "includes_other_repositories_and_failed_original_attempts":True,
        "photo_id_observation_id_csv_sha256":sha256(prior_path)
    }
    write_json(prior_receipt,receipt)
    return index_path,queue_path,future_path,prior_path,prior_receipt


def mock_gate(tmp_path, **changes):
    index,queue,future,prior,receipt=fixtures(tmp_path)
    options={"expected_photos":5,"expected_taxa":2,"expected_shards":2,
             "min_prior_attempts":3}
    options.update(changes)
    return gate(SPEC,index,queue,future,prior,receipt,**options)


def test_mock_future_two_photos_from_same_unseen_observation_admitted(tmp_path):
    r=mock_gate(tmp_path)
    assert r["status"]=="OUTCOME_BLIND_PHOTO_ROSTER_ELIGIBLE_NOT_ECOLOGICAL_CONFIRMATION"
    assert r["future_photos"]==2
    assert r["future_distinct_observations"]==1
    assert r["previously_exposed_observation_ids"]==2
    assert r["independent_source_observation_and_photo_IDs_verified"]
    assert r["hypothesis_confirmatory_inference_permitted"] is False
    assert r["native_UVB_bee_and_botanical_review_external_data_ready"] is False


def test_different_photo_of_previously_seen_observation_is_leakage(tmp_path):
    index,queue,future,prior,receipt=fixtures(tmp_path)
    write_csv(future,[{"photo_id":"4","observation_id":"OLD"},
                      {"photo_id":"5","observation_id":"OLD"}])
    with pytest.raises(ValueError,match="OBSERVATION"):
        gate(SPEC,index,queue,future,prior,receipt,
             expected_photos=5,expected_taxa=2,expected_shards=2,
             min_prior_attempts=3)


def test_known_photo_leakage_rejected(tmp_path):
    index,queue,future,prior,receipt=fixtures(tmp_path)
    write_csv(future,[{"photo_id":"2","observation_id":"NEW"},
                      {"photo_id":"5","observation_id":"NEW"}])
    with pytest.raises(ValueError,match="photo IDs"):
        gate(SPEC,index,queue,future,prior,receipt,
             expected_photos=5,expected_taxa=2,expected_shards=2,
             min_prior_attempts=3)


def test_source_only_photo_ids_do_not_become_known_colour_hues(tmp_path):
    index,queue,future,prior,receipt=fixtures(tmp_path)
    write_csv(future,[{"photo_id":"4","observation_id":"NEW","morph":"white"},
                      {"photo_id":"5","observation_id":"NEW","morph":"red_pink"}],
              fields=("photo_id","observation_id","morph"))
    with pytest.raises(ValueError,match="colour or ROI outcomes"):
        gate(SPEC,index,queue,future,prior,receipt,
             expected_photos=5,expected_taxa=2,expected_shards=2,
             min_prior_attempts=3)


def test_full_index_must_be_signed_and_complete_not_counts_only(tmp_path):
    index,queue,future,prior,receipt=fixtures(tmp_path)
    d=json.loads(index.read_text())
    d["photo_colours_measured"]=True
    write_json(index,d)
    with pytest.raises(ValueError,match="has not finished"):
        gate(SPEC,index,queue,future,prior,receipt,
             expected_photos=5,expected_taxa=2,expected_shards=2,
             min_prior_attempts=3)


def test_frozen_queue_must_be_source_sha_matched_precolour(tmp_path):
    index,queue,future,prior,receipt=fixtures(tmp_path)
    d=json.loads(queue.read_text())
    d["sampled_before_any_hue_measurement"]=False
    write_json(queue,d)
    with pytest.raises(ValueError,match="pre-colour"):
        gate(SPEC,index,queue,future,prior,receipt,
             expected_photos=5,expected_taxa=2,expected_shards=2,
             min_prior_attempts=3)


def test_prior_known_observations_must_include_other_repositories_and_unknowns(tmp_path):
    index,queue,future,prior,receipt=fixtures(tmp_path)
    d=json.loads(receipt.read_text())
    d["includes_other_repositories_and_failed_original_attempts"]=False
    write_json(receipt,d)
    with pytest.raises(ValueError,match="all repositories"):
        gate(SPEC,index,queue,future,prior,receipt,
             expected_photos=5,expected_taxa=2,expected_shards=2,
             min_prior_attempts=3)
