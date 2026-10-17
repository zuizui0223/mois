"""Reproducibility without inventing botanical open-petal ecological results."""
from __future__ import annotations
from collections import Counter
from datetime import date
import json
import hashlib
import pandas as pd
import pytest

from scripts.analysis import mois_target_season_cell_original_photo_source_only_20261017 as s
from scripts.analysis import mois_admit_targeted_256_real_pixels_and_siteyear_support_20261017 as m
from scripts.analysis.mois_reassemble_original_new4096_measurements import BASE_COLS,HUES


def photo(pid,tax="100",observer="observer-0",obs="obs",day="2025-01-15",lat=40.1,
          lon=139.1):
    p={k:"" for k in s.COLS}
    p.update(photo_id=pid,inat_taxon_id=tax,species="Example species",
             observer_id=observer,observation_id=obs,photo_license="cc-by",
             photo_url_large=f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.jpg",
             latitude=str(lat),longitude=str(lon),observed_on=day,
             source_snapshot_date="2026-09-27")
    return p


def old_photo(pid,day,season,obs):
    d={k:"" for k in BASE_COLS}
    d.update(photo_id=pid,inat_taxon_id="100",species="Example species",
             observation_id=obs,morph="white" if season=="summer" else "red_pink",
             observed_on=day,latitude="40.1",longitude="139.1")
    d["qualified_date"]=day
    d["original_source_date_provenance"]="original"
    return d


def test_old_seasonyear_cell_pivot_and_december_rollover():
    jan=s.locus("2025-01-18","40.1","139.1","100")
    dec=s.locus("2024-12-20","40.1","139.1","100")
    assert jan and jan[0]==dec[0] and jan[1]==dec[1]=="winter"
    assert s.locus("2025-07-08",40.1,139.1,"100")[1]=="summer"
    assert s.locus("2025-07-08",-40.1,139.1,"100")[1]=="winter"
    assert s.locus("2026-11-12",40.1,139.1,"100") is None
    assert s.locus("2025-01-18",0,139.1,"100") is None


def test_candidate_priority_preexisting_other_season_with_no_cue_colour(monkeypatch,tmp_path):
    monkeypatch.setattr(s,"SOURCE_N",5)
    monkeypatch.setattr(s,"FULL_SHA","fake_digest")
    monkeypatch.setattr(s,"sha",lambda *args:"fake_digest")
    oldloc,season=s.locus("2025-01-12",40.1,139.1,"100")
    assert season=="winter"
    pair={oldloc:{"summer":2,"winter":0}}
    rows=[
        photo("1001",observer="u1",obs="o1",day="2025-01-12"),
        photo("1002",observer="u2",obs="o2",day="2025-07-12"),
        photo("1003",observer="u3",obs="o3",day="2025-02-12",tax="101"),
        photo("1004",observer="u4",obs="o4",day="2025-01-12",lat=-1),
        photo("1005",observer="u5",obs="o5",day="2025-08-12")]
    ]
    src=tmp_path/"test.csv"
    with src.open("w",newline="") as f:
        import csv
        writer=csv.DictWriter(f,fieldnames=s.COLS)
        writer.writeheader();writer.writerows(rows)
    tops,commons,counts=s.scan(src,set(),pair)
    assert counts["source_photo_ID_candidates_scanned"]==5
    assert counts["strict_source_candidate_potential_repeat_photo_IDs"]==3
    assert any(x[2]["photo_id"]=="1001" and x[2]["target_priority"]=="1"
               for key,heap in tops.items() for x in heap)
    assert all("morph" not in x[2] for heap in tops.values() for x in heap)


def test_select_never_dup_photo_or_observer_or_observation(monkeypatch):
    monkeypatch.setattr(s,"RESERVE_N",3)
    monkeypatch.setattr(s,"SOURCE_TAXA",2,raising=False)
    top={}
    globalpool={}
    for tax in ("100","101"):
        heap=[]
        for i in range(4):
            p=photo(f"10{tax}{i}",tax=tax,
                    observer=f"u-{tax}-{i}",obs=f"o-{tax}-{i}")
            s.push_smallest(heap,s.identity(p["photo_id"]),p,4)
        globalpool[tax]=heap
    out,st=s.pick(top,globalpool,reserved=set(),n=3)
    assert len(out)==3
    assert len({p["photo_id"] for p in out})==3
    assert len({p["observer_id"] for p in out})==3
    assert len({p["observation_id"] for p in out})==3
    assert [int(p["shard"]) for p in out]==[0,0,0]


def test_targeting_replication_is_photographic_source_only():
    rows=[
        old_photo("1","2025-07-12","summer","o1"),
        old_photo("2","2025-01-12","winter","o2"),
    ]
    d=pd.DataFrame(rows)
    b,_=m.support(d)
    assert b["north_extratropical"][
        "same_original_taxon_025degree_cell_same_seasonyear_both_summer_winter_min1_obs"]==1
    assert b["north_extratropical"][
        "same_original_taxon_025degree_cell_same_seasonyear_both_summer_winter_min2_obs"]==0
    fresh=old_photo("3","2025-01-16","winter","o3")
    merged=pd.concat([d,pd.DataFrame([fresh])],ignore_index=True)
    a,_=m.support(merged)
    assert a["north_extratropical"][
        "same_original_taxon_025degree_cell_same_seasonyear_both_summer_winter_min2_obs"]==0
    # Need both summer and winter 2 observations, never treat extra winter alone as >=2.
    merged=pd.concat([merged,pd.DataFrame([old_photo(
        "4","2025-07-17","summer","o4")])],ignore_index=True)
    a,_=m.support(merged)
    assert a["north_extratropical"][
        "same_original_taxon_025degree_cell_same_seasonyear_both_summer_winter_min2_obs"]==1


def test_model_attestation_not_faked_by_source_identifiers():
    assert m.METHOD["location_blind"] is True
    assert s.N_SHARDS==16 and s.SHARD_SIZE==256
    assert s.RESERVE_N==4096
    assert m.SHARD_SIZE==256
    assert not hasattr(s,"morph")
