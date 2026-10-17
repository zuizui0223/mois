"""Source-independent synthetic locality/season-year photograph audit.

These fixtures test code identities, NOT real flower colour ecology.
"""
import numpy as np
import pandas as pd
import pytest

from scripts.analysis.mois_seasonal_267k_climate_normals import labelled_calendar
from scripts.analysis.mois_original_photo_community_season_geoscale import (
    matched_community, cell_cluster_interval,
)


def example(count=5, *, winter_year=2025, summer_lat=40.01, winter_lat=40.01):
    rows=[]
    ix=0
    for season, month, year, latitude, morph in (
        ("winter",1,winter_year,winter_lat,"red_pink"),
        ("summer",7,2025,summer_lat,"white")):
        for k in range(count):
            ix+=1
            rows.append({"photo_id":str(ix),"observation_id":str(ix),
              "inat_taxon_id":str(k%2+100),"morph":morph,
              "observed_on":f"{year}-{month:02d}-04",
              "latitude":latitude,"longitude":138.01})
    return labelled_calendar(pd.DataFrame(rows))


def test_local_community_can_compare_different_season_taxa_not_individual_pigments():
    d=example()
    x=matched_community(d,width=.25,min_obs=5,min_taxa=2)
    assert len(x)==1
    assert x.summer_photos.iloc[0]==5 and x.winter_photos.iloc[0]==5
    assert x.summer_minus_winter_white_pp.iloc[0]==pytest.approx(100)
    assert matched_community(d,width=.25,min_obs=6).empty
    # Site-year sample does not require a named taxon photographed in both seasons.


def test_distinct_observations_not_multiple_photos_from_one_observation():
    d=example()
    d.loc[d.season.eq("summer"),"observation_id"]="same_observation"
    assert matched_community(d,width=.25,min_obs=3,min_taxa=2).empty


def test_community_repeats_do_not_span_incompatible_years():
    assert matched_community(example(winter_year=2024),width=.25,min_obs=5).empty


def test_increasing_geocell_width_creates_apparent_local_comparison():
    d=example(summer_lat=40.01,winter_lat=40.30)
    assert matched_community(d,width=.25,min_obs=5).empty
    assert len(matched_community(d,width=.5,min_obs=5))==1


def test_bootstrap_returns_no_interval_from_one_cell_and_reproducible_intervals():
    one=matched_community(example(),width=.25,min_obs=5)
    assert cell_cluster_interval(one)["bootstrap_95pp"] is None
    multi=pd.concat([one.assign(latitude_cell=int(one.latitude_cell.iloc[0])+j)
                     for j in range(10)],ignore_index=True)
    r1=cell_cluster_interval(multi,draws=30,seed=111)
    r2=cell_cluster_interval(multi,draws=30,seed=111)
    assert r1==r2
    assert r1["n_unique_geocells"]==10
    assert r1["bootstrap_95pp"][0]==pytest.approx(100)
    assert r1["bootstrap_95pp"][1]==pytest.approx(100)


def test_unknown_photo_hues_do_not_count_as_white_nor_winter_photos():
    d=example(count=5)
    d.loc[(d.season=="summer") & (d.photo_id=="10"),"morph"]="mixed_uncertain"
    d["classified"]=d.morph.isin(("white","yellow_orange","red_pink","blue_purple"))
    d["white"]=d.morph.eq("white")
    assert matched_community(d,width=.25,min_obs=5).empty
    q=matched_community(d,width=.25,min_obs=3)
    assert len(q)==1
    assert int(q.summer_photos.iloc[0])==4
