"""Synthetic support/denominator regression tests, not biological observations."""
import pandas as pd
import pytest

from scripts.analysis import mois_historical9600_taxonomic_breadth_vs_local_season_support as m


def example():
    records=[]
    for pid,oid,tax,when,morph,lat,source in [
        ("1","11","A",None,"white",40,"historical"),
        ("2","12","B",None,"mixed_uncertain",-40,"historical"),
        ("3","13","C",None,"red_pink",40,"historical"),
        ("4","14","A","2025-07-03","white",40,"new_fcp"),
        ("5","15","A","2025-01-03","red_pink",40,"new_fcp"),
    ]:
        records.append(dict(photo_id=pid,observation_id=oid,
          inat_taxon_id=tax,morph=morph,observed_on=when,
          latitude=lat,longitude=135.01,old_mois_or_new_fcp=source))
    return pd.DataFrame(records)


def recovered():
    return pd.DataFrame([
       dict(original_photo_id="1",retrieved_observed_on="2024-07-12",
         provenance_status="source_photo_and_observation_verified_with_current_public_date"),
       dict(original_photo_id="2",retrieved_observed_on=None,
         provenance_status="not_found_in_current_public_api"),
       dict(original_photo_id="3",retrieved_observed_on="2026-02-12",
         provenance_status="source_photo_and_observation_verified_with_current_public_date"),
    ])


def test_exact_source_photo_dates_added_but_unclassified_remain_unknown(monkeypatch):
    monkeypatch.setattr(m,"MOIS_ROWS",5)
    source=example()
    old=m.make_panel(source,recovered(),False)
    new=m.make_panel(source,recovered(),True)
    assert int(old.date.notna().sum())==2
    assert int(new.date.notna().sum())==4
    assert int((new.classified&new.date.notna()).sum())==4
    assert new.loc[new.photo_id=="2","date"].isna().all()
    assert new.loc[new.photo_id=="3","morph"].iloc[0]=="red_pink"


def test_no_invented_same_population_or_year_summer_winter_pair(monkeypatch):
    monkeypatch.setattr(m,"MOIS_ROWS",5)
    data=m.make_panel(example(),recovered(),True)
    broad=m.spatial_support(data,width=.25,same_year=False,minimum=1,
                            zone="north_extratropical")
    strict=m.spatial_support(data,width=.25,same_year=True,minimum=1,
                             zone="north_extratropical")
    assert broad["matched_groups"]==1
    assert strict["matched_groups"]==1  # original new_fcp 2025 summer and winter
    assert m.spatial_support(data,width=.25,same_year=True,minimum=2,
                             zone="north_extratropical")["matched_groups"]==0


def test_source_photo_date_collision_is_rejected(monkeypatch):
    monkeypatch.setattr(m,"MOIS_ROWS",5)
    changed=recovered().copy()
    changed.loc[0,"original_photo_id"]="4"
    with pytest.raises(ValueError,match="overwrite"):
        m.make_panel(example(),changed,True)


def test_southern_empty_original_source_both_seasons_returns_no_match(monkeypatch):
    monkeypatch.setattr(m,"MOIS_ROWS",5)
    panel=m.make_panel(example(),recovered(),True)
    assert m.spatial_support(panel,width=.25,same_year=True,minimum=1,
                             zone="south_extratropical")["matched_groups"]==0
