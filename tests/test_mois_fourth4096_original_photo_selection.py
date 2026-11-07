"""Real fourth-wave sampling identity safeguards, synthetic source-only fixtures.

Tests do not qualify any fake photos as measured original pixels or petals.
"""
import copy
import json
from collections import Counter

import pytest
from scripts.analysis import mois_original_fourth4096_after_309439_source_freeze as m


def fake_source(n_tax=3,n_each=12):
    rows=[]
    for tax in range(1,n_tax+1):
        for i in range(n_each):
            pid=str(1000000+tax*1000+i)
            rows.append({
              "photo_id":pid,"inat_taxon_id":str(tax),"species":f"Species {tax}",
              "observer_id":f"observer_{tax}_{i}",
              "observation_id":str(150000000+tax*1000+i),
              "observed_on":"2025-06-18","latitude":str(35+tax*0.5),
              "longitude":"138.95","photo_license":"cc-by",
              "photo_url_large":f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.jpg",
              "source_snapshot_date":"2026-09-27",
              "positional_accuracy_m":"40",
            })
    return rows


def test_true_original_photo_and_observer_id_disjoint_selection(monkeypatch):
    monkeypatch.setattr(m,"SELECTION_N",12)
    observed=fake_source()
    old={observed[0]["photo_id"],observed[15]["photo_id"]}
    chosen,stats=m.choose(iter(observed),old,n=12,expected_taxa=3)
    assert len(chosen)==12 and len({x["photo_id"] for x in chosen})==12
    assert len({x["observer_id"] for x in chosen})==12
    assert len({x["observation_id"] for x in chosen})==12
    assert len({x["inat_taxon_id"] for x in chosen})==3
    assert not set(x["photo_id"] for x in chosen)&old
    assert stats["source_original_photo_IDs_seen"]==len(observed)
    assert all(x["measurement_id"].startswith("MOIS-") for x in chosen)
    assert all(x["wave_sequence"].isdigit() for x in chosen)


def test_rank_and_selection_photo_source_only_not_pixel_hue_or_climate(monkeypatch):
    monkeypatch.setattr(m,"SELECTION_N",12)
    src=fake_source()
    first,_=m.choose(src,set(),n=12,expected_taxa=3)
    second=copy.deepcopy(src)
    second.reverse()
    for x in second:
        x["flower_hue"]="white"
        x["petal_quality"]="fake_bloom"
        x["monthly_precipitation"]=9999
        x["local_pollinator_density"]=0
    again,_=m.choose(second,set(),n=12,expected_taxa=3)
    assert [x["photo_id"] for x in first]==[x["photo_id"] for x in again]
    assert all("flower_hue" not in x for x in again)
    assert all("monthly_precipitation" not in x for x in again)


def test_reject_duplicate_observers_and_invalid_photo_url(monkeypatch):
    monkeypatch.setattr(m,"SELECTION_N",12)
    src=fake_source()
    for row in src:
        row["observer_id"]="same_observer_for_every_photo"
    with pytest.raises(ValueError,match="unique photo/observation/observer"):
        m.choose(src,set(),n=12,expected_taxa=3)
    bad=fake_source()
    bad[0]["photo_url_large"]="https://other-provider/photos/1/large.jpg"
    _,count=m.choose(bad,set(),n=12,expected_taxa=3)
    assert count["invalid_original_licensed_photo_URL_no_pixel_inference"]==1


def test_invalid_date_original_not_selected_even_if_valid_photo_ID(monkeypatch):
    monkeypatch.setattr(m,"SELECTION_N",12)
    src=fake_source()
    src[0]["observed_on"]="2025-02-30"
    src[1]["source_snapshot_date"]="2025-01-01"
    _,count=m.choose(src,set(),n=12,expected_taxa=3)
    assert count["invalid_source_metadata_no_pixel_inference"]==2


def test_never_represent_source_only_4096_as_classified_measurement():
    assert m.MEASURED_N==309439
    assert m.MEASURED_HUES==149662
    assert m.MEASURED_DATED_HUES==149653
    assert len(m.FULL_SHA)==len(m.MEASURED_SHA)==64
    assert m.SELECTION_N==16*256==4096
    assert m.SOURCE_TAXA==1000
    assert "WAVE4" in m.SALT or "WAVE4" in m.SALT.upper()
