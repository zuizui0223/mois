"""Synthetic verification of prospective 9600 authentic original photo-date rejoin."""
import hashlib
import numpy as np
import pandas as pd
import pytest
from scripts.acquisition import mois_historical9600_source_date_expansion as m


def tiny_photo_frame():
    return pd.DataFrame([
        {"photo_id":"101","observation_id":"1","inat_taxon_id":"20","morph":"white",
         "old_mois_or_new_fcp":"historical","observed_on":None,"latitude":40.1,"longitude":138.1},
        {"photo_id":"102","observation_id":"2","inat_taxon_id":"20","morph":"mixed_uncertain",
         "old_mois_or_new_fcp":"historical","observed_on":None,"latitude":40.1,"longitude":138.1},
        {"photo_id":"103","observation_id":"3","inat_taxon_id":"21","morph":"red_pink",
         "old_mois_or_new_fcp":"historical","observed_on":None,"latitude":-41.1,"longitude":139.1},
        {"photo_id":"200","observation_id":"4","inat_taxon_id":"11","morph":"white",
         "old_mois_or_new_fcp":"new_fcp","observed_on":"2025-06-18","latitude":40.1,"longitude":138.1},
    ])


def test_original_photo_hash_selection_is_color_and_environment_blind():
    frame=pd.DataFrame([{"photo_id":str(i),"observation_id":str(i*10),
                        "morph":"white" if i%2 else "red_pink",
                        "latitude":30 if i%2 else -30} for i in range(100)])
    chosen=m.fixed_prefix(frame,40)
    another=frame.sample(frac=1,random_state=22).copy()
    another["morph"]="mixed_uncertain"
    another["latitude"]=0
    assert chosen.photo_id.tolist()==m.fixed_prefix(another,40).photo_id.tolist()
    assert chosen.observation_id.is_unique
    with pytest.raises(ValueError,match="bounded"):
        m.fixed_prefix(frame,9601)


def test_original_photo_hash_selection_guard_duplicate_observation():
    a=pd.DataFrame([{"photo_id":"1","observation_id":"2"},
                    {"photo_id":"3","observation_id":"2"}])
    with pytest.raises(ValueError,match="identifier uniqueness"):
        m.fixed_prefix(a,2)


def test_photo_dates_rejoin_only_validated_exact_original_phenotype_and_other_source(monkeypatch):
    monkeypatch.setattr(m,"MOIS_ROWS",4)
    monkeypatch.setattr(m,"HISTORICAL_ROWS",3)
    source=tiny_photo_frame()
    selected=source.iloc[:3].copy()
    found=pd.DataFrame([
        {"original_photo_id":"101","retrieved_observed_on":"2026-01-21",
         "provenance_status":"source_photo_and_observation_verified_with_current_public_date"},
        {"original_photo_id":"102","retrieved_observed_on":"2025-05-11",
         "provenance_status":"source_photo_and_observation_verified_with_current_public_date"},
        {"original_photo_id":"103","retrieved_observed_on":None,
         "provenance_status":"observation_found_but_original_photo_not_in_current_photo_set"}
    ])
    out, receipt=m.normal_season_table(source,selected,found)
    assert receipt["original_classified_with_original_observed_on"]==1
    assert receipt["new_classified_historical_source_dates_from_verified_original_photo_ID"]==1
    assert receipt["expanded_classified_original_photo_day_records"]==2
    assert receipt["historic_original_photo_dates_still_unverified"]==1
    hist=out[out.photo_source_group=="historical_sample_with_verified_current_date"]
    # Unknown hue contributes an attempt, not a fake white or chromatic photo.
    assert hist.photo_attempts.sum()==2
    assert hist.classified_hues.sum()==1


def test_refuse_overwrite_existing_original_photo_day(monkeypatch):
    monkeypatch.setattr(m,"MOIS_ROWS",4)
    monkeypatch.setattr(m,"HISTORICAL_ROWS",3)
    source=tiny_photo_frame()
    source.loc[0,"observed_on"]="2020-12-10"
    selected=source.iloc[:3].copy()
    found=pd.DataFrame([
        {"original_photo_id":"101","retrieved_observed_on":"2021-12-10",
         "provenance_status":"source_photo_and_observation_verified_with_current_public_date"},
        {"original_photo_id":"102","retrieved_observed_on":None,
         "provenance_status":"not_found_in_current_public_api"},
        {"original_photo_id":"103","retrieved_observed_on":None,
         "provenance_status":"not_found_in_current_public_api"}
    ])
    with pytest.raises(ValueError,match="collide"):
        m.normal_season_table(source,selected,found)


def test_refuse_an_out_of_order_or_replaced_photograph(monkeypatch):
    monkeypatch.setattr(m,"MOIS_ROWS",4)
    monkeypatch.setattr(m,"HISTORICAL_ROWS",3)
    source=tiny_photo_frame()
    selected=source.iloc[:3].copy()
    wrong=pd.DataFrame([
        {"original_photo_id":"102","retrieved_observed_on":"2026-01-21",
         "provenance_status":"source_photo_and_observation_verified_with_current_public_date"},
        {"original_photo_id":"101","retrieved_observed_on":None,
         "provenance_status":"not_found_in_current_public_api"},
        {"original_photo_id":"103","retrieved_observed_on":None,
         "provenance_status":"not_found_in_current_public_api"}
    ])
    with pytest.raises(ValueError,match="photo-ID order"):
        m.normal_season_table(source,selected,wrong)


def test_source_prefix_hash_matches_old_algorithm():
    source=pd.DataFrame([{"photo_id":str(i),"observation_id":str(i+100)} for i in range(100)])
    picked=m.fixed_prefix(source,15)
    keys=[(hashlib.sha256(("mois_historical_dated_source_pilot_v1:"+str(i)).encode()).hexdigest(),
           str(i)) for i in range(100)]
    expect=[p for _,p in sorted(keys)[:15]]
    assert picked.photo_id.tolist()==expect
