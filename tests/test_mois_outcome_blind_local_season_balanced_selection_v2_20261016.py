"""Synthetic tests of prospective season-balanced V2 Mois photo selection.

No original photo pixels, retrospective hue labels, annotation-selection
feedback or unsealed global source photos are used to choose any photo IDs.
"""
import csv
import importlib.util
import io
import json
from collections import Counter
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]


def module(name,relpath):
    s=importlib.util.spec_from_file_location(name,ROOT/relpath)
    v=importlib.util.module_from_spec(s)
    s.loader.exec_module(v)
    return v


v2=module("mois_v2_season_balanced","scripts/acquisition/mois_outcome_blind_local_season_balanced_selection_v2_20261016.py")
v1=module("mois_v1_sampling","scripts/acquisition/mois_stratified_photo_selection_20261016.py")
direct_test=module("mois_direct_test_source","tests/test_mois_direct_fcp_source_sampling_20261016.py")
direct=module("mois_direct_source","scripts/acquisition/mois_sample_direct_from_fcp_source_20261016.py")


def make_rows():
    rows=[]
    for month,count in ((4,184),(7,8),(10,8)):
        for k in range(count):
            photo=len(rows)+11
            rows.append({
                "inat_taxon_id":"100","genus":"Alpha",
                "photo_id":str(photo),"observation_id":f"obs-{photo}",
                "latitude":"40","longitude":"130","bio12_mm":"1000",
                "photo_license":"cc-by",
                "observed_on":f"2024-{month:02}-12",
                "photo_url_large":f"https://inaturalist-open-data.s3.amazonaws.com/photos/{photo}/large.jpg",
                "observer_id":f"observer-{k % 8}","attribution":"synthetic source",
                "positional_accuracy_m":"30",
            })
    return rows


def sampled(func,rows):
    stream=io.StringIO()
    writer=csv.DictWriter(stream,fieldnames=v1.PRIVATE_FIELDS)
    writer.writeheader()
    result=func(100,"Alpha",rows,writer)
    stream.seek(0)
    return result,list(csv.DictReader(stream))


def seasons(rows):
    return Counter(
        direct_test.s.local_season(x["latitude"],x["observed_on"])[1]
        if hasattr(direct_test.s,"local_season") else
        v2.local_season(x["latitude"],x["observed_on"])[1]
        for x in rows
    )


def test_seasonal_sample_improves_temporal_support_without_changing_total_opportunities():
    rows=make_rows()
    old,old_rows=sampled(v1.sample_species,rows)
    new,new_rows=sampled(v2.sample_species_season_v2,rows)
    assert old[:3]==new[:3]==(200,200,40)
    assert len({x["observation_id"] for x in new_rows})==40
    assert len({x["photo_id"] for x in new_rows})==40
    assert all(x["latitude"]=="40" and x["bio12_mm"]=="1000" for x in new_rows)
    early=seasons(old_rows)
    seasonal=seasons(new_rows)
    assert seasonal["SPRING"]>=3
    assert seasonal["SUMMER"]>=3
    assert seasonal["AUTUMN"]>=3
    assert sum(seasonal.values())==sum(early.values())==40
    assert sum(seasonal[k]>=3 for k in ("SPRING","SUMMER","AUTUMN")) >= (
        sum(early[k]>=3 for k in ("SPRING","SUMMER","AUTUMN"))
    )
    # They are intentionally different frames, NOT independent repetitions.
    assert {x["photo_id"] for x in new_rows}!={x["photo_id"] for x in old_rows}
    for x in new_rows:
        pi=float(x["observation_inclusion_probability"])
        assert 0<pi<=1
        assert abs(pi*float(x["photo_design_weight"])-1)<1e-9
    assert new[4]["__nested_observations__"]==40
    assert new[4].get("__fallback_observations__",0)==0


def test_nested_phase_quotas_never_force_unobserved_months_or_zero_probabilities():
    x=v2.nested_season_quotas({"SPRING":184,"SUMMER":8,"AUTUMN":8},40)
    assert sum(x.values())==40
    assert all(x[s]>=3 for s in x)
    assert v2.nested_season_quotas({"SPRING":184,"SUMMER":8,"AUTUMN":8},2) is None
    assert v2.nested_season_quotas({"SPRING":1,"SUMMER":1},2)=={"SPRING":1,"SUMMER":1}
    with pytest.raises(ValueError,match="Invalid"):
        v2.nested_season_quotas({"SPRING":1},5)


def test_proper_two_stage_photo_probability_for_multiple_photos_one_observation():
    rows=make_rows()
    # Add a second original photo attached to the SAME observation; source
    # photo-level probability must divide by that observation's multiplicity.
    extra=dict(rows[0],photo_id="100000",
               photo_url_large="https://inaturalist-open-data.s3.amazonaws.com/photos/100000/large.jpg")
    rows.append(extra)
    result,selected=sampled(v2.sample_species_season_v2,rows)
    assert result[:3]==(201,200,40)
    for x in selected:
        obs_pi=float(x["observation_inclusion_probability"])
        photo_pi=float(x["photo_inclusion_probability"])
        m=int(x["source_photos_in_observation"])
        assert abs(photo_pi-obs_pi/m)<1e-10
        assert abs(float(x["photo_design_weight"])-1/photo_pi)<1e-7
        if x["observation_id"]==rows[0]["observation_id"]:
            assert m==2


def test_unverified_observation_date_and_positive_season_fallback():
    assert v2.nested_season_quotas({"UNKNOWN_SEASON":3,"SPRING":5},1) is None
    rows=make_rows()
    rows[0]["observed_on"]=""
    result,chosen=sampled(v2.sample_species_season_v2,rows)
    assert result[2]==40
    assert len(chosen)==40
    assert all(float(x["observation_inclusion_probability"])>0 for x in chosen)


def test_same_observation_with_changed_date_is_rejected_not_blindly_duplicated():
    rows=make_rows()
    extra=dict(rows[0],photo_id="100000",observed_on="2024-08-15")
    with pytest.raises(ValueError,match="inconsistent site/season/observer"):
        sampled(v2.sample_species_season_v2,[*rows,extra])


def test_sealed_direct_FCP_128_shard_replay_has_explicit_v2_protocol(tmp_path):
    source,counts,admission=direct_test.fixture(tmp_path)
    selection=tmp_path/"private_season_balanced.csv"
    receipt=tmp_path/"receipt_v2.json"
    result=direct.run(source,counts,admission,direct_test.FakeBio12(),
                      selection,receipt,calibration={"synthetic":True},
                      synthetic=True,species_sampler=v2.sample_species_season_v2,
                      selection_protocol=v2.STRATEGY)
    assert result["schema"]==v2.SCHEMA
    assert result["sampling_protocol"]==v2.STRATEGY
    assert result["source_species"]==3
    assert result["source_photo_ids_authenticated"]==9
    assert result["source_observations"]==8
    assert result["selected_original_photo_ids"]==8
    assert result["selected_local_season_nested_observations"]==8
    assert result["selected_local_season_fallback_observations"]==0
    assert result["pixels_opened_or_colours_measured"] is False
    assert direct.sha256_file(selection)==result["selected_private_photo_csv_sha256"]
    assert len(list(csv.DictReader(selection.open())))==8


def test_v2_protocol_not_accepted_with_accidentally_missing_sampler(tmp_path):
    source,counts,admission=direct_test.fixture(tmp_path)
    with pytest.raises(ValueError,match="Alternate sampling protocol must be named"):
        direct.run(source,counts,admission,direct_test.FakeBio12(),
                   tmp_path/"out.csv",tmp_path/"out.json",
                   calibration={"synthetic":True},synthetic=True,
                   selection_protocol=v2.STRATEGY)



def test_v2_sealed_photo_selection_keeps_original_design_weights_through_roi_and_season_preflight(tmp_path):
    handoff=module("mois_v2_roi_handoff",
                   "scripts/acquisition/mois_build_sample_roi_handoff_20261016.py")
    season_support=module("mois_v2_season_support",
                          "scripts/analysis/mois_seasonal_sample_support_preflight_20261016.py")
    source,counts,admission=direct_test.fixture(tmp_path)
    selected=tmp_path/"v2_private.csv"
    source_receipt=tmp_path/"v2_receipt.json"
    original=direct.run(
        source,counts,admission,direct_test.FakeBio12(),
        selected,source_receipt,calibration={"synthetic":True},
        synthetic=True,species_sampler=v2.sample_species_season_v2,
        selection_protocol=v2.STRATEGY
    )
    packets=tmp_path/"private_packets"
    verified=handoff.assemble(selected,source_receipt,packets,synthetic=True)
    assert verified["selected_photo_ids"]==original["selected_original_photo_ids"]
    assert verified["original_source_selection_schema"]==v2.SCHEMA
    assert verified["original_source_selection_protocol"]==v2.STRATEGY
    assert verified["original_pixels_measured"] is False
    preflight=season_support.audit(
        selected,source_receipt,tmp_path/"private_season.json",
        tmp_path/"public_season.json",synthetic=True)
    assert preflight["source_sampling_protocol"]==v2.STRATEGY
    assert preflight["preclassification_selected_unique_observations"]==8
    assert preflight["real_classifiable_photos_measured"]==0
    assert preflight["synthetic_test_only"] is True


def test_v2_roi_handoff_refuses_forged_or_missing_seasonal_strategy(tmp_path):
    handoff=module("mois_v2_roi_handoff_guard",
                   "scripts/acquisition/mois_build_sample_roi_handoff_20261016.py")
    source,counts,admission=direct_test.fixture(tmp_path)
    selected=tmp_path/"v2_private.csv"
    receipt=tmp_path/"v2_receipt.json"
    direct.run(source,counts,admission,direct_test.FakeBio12(),
               selected,receipt,calibration={"synthetic":True},
               synthetic=True,species_sampler=v2.sample_species_season_v2,
               selection_protocol=v2.STRATEGY)
    forged=json.loads(receipt.read_text())
    forged["sampling_protocol"]="UNKNOWN_FLOWER_OUTCOME_SELECTED_MODEL"
    receipt.write_text(json.dumps(forged))
    with pytest.raises(ValueError,match="season-balanced sampler provenance"):
        handoff.assemble(selected,receipt,tmp_path/"private_packets",synthetic=True)
