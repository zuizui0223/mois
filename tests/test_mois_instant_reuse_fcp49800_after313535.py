"""Reuse completed FCP photos without repeating pixel inference or inflating MOIS source."""
import numpy as np
import pandas as pd
import pytest
from scripts.analysis import mois_instant_reuse_fcp49800_after313535 as m


def original():
    return pd.DataFrame([
        {"photo_id":"1","inat_taxon_id":"7","species":"Cirsium a",
         "morph":"white","observation_id":"10","qualified_date":"2025-03-18",
         "old_mois_or_new_fcp":"historical","latitude":33., "longitude":135.,
         "classified":True,"white":True},
        {"photo_id":"2","inat_taxon_id":"9","species":"Campanula b",
         "morph":"mixed_uncertain","observation_id":"11","qualified_date":"2025-08-19",
         "old_mois_or_new_fcp":"historical","latitude":38.,"longitude":137.,
         "classified":False,"white":False},
    ])


def fcp_row(pid,morph,*,tax="8",date="2025-07-18"):
    return {"photo_id":pid,"inat_taxon_id":tax,"species":"Cirsium c",
            "morph":morph,"observation_id":"obs"+pid,
            "latitude":"32.1","longitude":"128.8",
            "measurement_status":("classified_four_state_morph"
                                  if morph in m.HUES else "mixed_uncertain"),
            "observed_on":date}


def prepare(monkeypatch):
    old=original()
    monkeypatch.setattr(m.prior,"N_PHOTOS",2)
    monkeypatch.setattr(m.prior,"N_HUES",1)
    monkeypatch.setattr(m.prior,"WHITE_N",1)
    monkeypatch.setattr(m,"EXPECTED_ORIGINAL",2)
    monkeypatch.setattr(m,"EXPECTED_FCP_PHOTOS",5)
    monkeypatch.setattr(m.prior,"load_source",lambda root:(old,{"source":"fake"}))
    monkeypatch.setattr(m.pd,"read_parquet",lambda path:old)
    rows=[fcp_row("1","red_pink",tax="7"),
          fcp_row("2","mixed_uncertain",tax="9"),
          fcp_row("3","white"),
          fcp_row("4","red_pink"),
          fcp_row("5","yellow_orange")]
    monkeypatch.setattr(m,"read_fcp_wave4",lambda root:(rows,{
        "measurement_original_fixed_colour_worker":True,
        "measured_sha256":"5"*64,
        "source_metadata_sha256":"6"*64,
    }))
    monkeypatch.setattr(m,"reservation_ids",
        lambda a,b,c:({"4"},{"5"},{"first_wave_original_4096_csv_sha256":"x"},
                            {"first_wave_original_4096_csv_sha256":"y"}))
    return old


def test_real_source_reuse_prevents_prior_and_reserved_photo_overlap(monkeypatch,tmp_path):
    old=prepare(monkeypatch)
    before,added,whole,counts,receipt,manifest,fifth,season=m.main_source(
        tmp_path,tmp_path,tmp_path,tmp_path)
    assert len(before)==2
    assert len(added)==1
    assert added.photo_id.tolist()==["3"]
    assert len(whole)==3
    assert whole.iloc[:2].reset_index(drop=True).equals(old)
    assert whole.photo_id.is_unique
    assert counts["already_in_primary313535"]==2
    assert counts["deferred_to_running_fifth4096_wave"]==1
    assert counts["excluded_SI_season_targeted4096"]==1
    assert counts["overlap_hue_disagreement"]==1
    assert counts["actually_novel_photo_four_hue_ROIs"]==1


def test_classifier_technical_unknown_not_assigned_coloured(monkeypatch,tmp_path):
    old=prepare(monkeypatch)
    allrows=[fcp_row("1","white",tax="7"),fcp_row("2","mixed_uncertain",tax="9"),
             fcp_row("3","mixed_uncertain"),fcp_row("4","red_pink"),
             fcp_row("5","blue_purple")]
    monkeypatch.setattr(m,"read_fcp_wave4",lambda root:(allrows,{
        "measurement_original_fixed_colour_worker":True,
        "measured_sha256":"5"*64,
        "source_metadata_sha256":"6"*64}))
    _,added,combined,counts,*_=m.main_source(tmp_path,tmp_path,tmp_path,tmp_path)
    assert len(added)==1
    assert added.morph.tolist()==["mixed_uncertain"]
    assert not bool(added.classified.iloc[0])
    assert not bool(added.white.iloc[0])
    assert counts["actually_novel_ROI_processed_photo_IDs"]==1
    assert counts["actually_novel_photo_four_hue_ROIs"]==0
    assert combined.photo_id.is_unique


def test_source_failure_prevents_imputed_new_pixel_count(monkeypatch,tmp_path):
    old=prepare(monkeypatch)
    monkeypatch.setattr(m,"read_fcp_wave4",lambda root:([fcp_row("3","white")],{
        "measurement_original_fixed_colour_worker":True}))
    with pytest.raises(ValueError,match="49800 real source"):
        m.main_source(tmp_path,tmp_path,tmp_path,tmp_path)


def test_source_unknowns_date_and_principle():
    assert "anthocyanin" not in m.SOURCE_FCP
    assert "wave4_49800" in m.SOURCE_FCP
    assert m.DATE.fullmatch("2025-07-18")
    assert not m.DATE.fullmatch("2025-07")
    assert m.SOURCE_FCP!="historical"
