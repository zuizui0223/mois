"""Synthetic checks of new genuine ROI date/pixel merge, no ecology inferred."""
from __future__ import annotations

import pandas as pd
import pytest

from scripts.analysis.mois_reassemble_original_new4096_measurements import BASE_COLS,HUES
from scripts.analysis import mois_replay_fourth4096_actual_ROI_pixels_and_source_dates as m


def fake_photo(pid,kind,day,morph):
    r={key:"" for key in BASE_COLS}
    r.update(photo_id=pid,inat_taxon_id="456",species="Campanula sample",
             morph=morph,source_original_measured=kind,
             old_mois_or_new_fcp="mois_new_original_ROI",
             measurement_status="classified_four_state_morph" if morph in HUES
                     else "not_evaluable_no_biological_palette_mass",
             latitude="40.12",longitude="139.71",
             photo_license="cc-by",
             observation_id="obs-"+pid,observed_on=day)
    return r


def test_append_actual_new_original_ROI_and_observed_on_no_reclassification(monkeypatch):
    monkeypatch.setattr(m,"MEASURED_N",2)
    monkeypatch.setattr(m,"MEASURED_HUES",1)
    monkeypatch.setattr(m,"MEASURED_DATED_HUES",1)
    monkeypatch.setattr(m,"NEW_PIXEL_N",3)
    old=pd.DataFrame([fake_photo("10","old","2025-06-18","white"),
                      fake_photo("11","old","","mixed_uncertain")])
    old["qualified_date"]=["2025-06-18",""]
    old["original_source_date_provenance"]=["old","missing"]
    fresh=[
        fake_photo("12","mois_wave4_new_original_ROI_v4","2025-07-11","red_pink"),
        fake_photo("13","mois_wave4_new_original_ROI_v4","2024-07-11","mixed_uncertain"),
        fake_photo("14","mois_wave4_new_original_ROI_v4","2025-12-30","blue_purple")
    ]
    before=old.copy(deep=True)
    result,summary=m.append(old,fresh)
    assert len(result)==5
    assert result.photo_id.is_unique
    assert result.iloc[:2].reset_index(drop=True).equals(before)
    assert summary["new_photo_four_hue_original_ROI_classified"]==2
    assert summary["new_photo_unknown_pixel_OR_ROI_color"]==1
    assert summary["expanded_original_four_hue_ROI_photo_ids"]==3
    assert summary["expanded_actual_date_and_four_hue_ROI_photo_ids"]==3
    assert result.loc[result.photo_id.eq("13"),"qualified_date"].iloc[0]=="2024-07-11"
    assert result.loc[result.photo_id.eq("11"),"qualified_date"].iloc[0]==""
    assert result.loc[result.photo_id.eq("10"),"morph"].iloc[0]=="white"


def test_reject_original_photo_duplicate(monkeypatch):
    monkeypatch.setattr(m,"MEASURED_N",2)
    monkeypatch.setattr(m,"NEW_PIXEL_N",1)
    old=pd.DataFrame([fake_photo("10","old","2025-06-18","white"),
                      fake_photo("11","old","","mixed_uncertain")])
    old["qualified_date"]=["2025-06-18",""]
    old["original_source_date_provenance"]=["old","missing"]
    with pytest.raises(ValueError,match="overlap"):
        m.append(old,[fake_photo("10","wave4","2025-01-01","white")])


def test_reject_undated_unmeasured_or_unclassified_source(monkeypatch):
    monkeypatch.setattr(m,"MEASURED_N",2)
    monkeypatch.setattr(m,"NEW_PIXEL_N",1)
    old=pd.DataFrame([fake_photo("10","old","2025-06-18","white"),
                      fake_photo("11","old","","mixed_uncertain")])
    old["qualified_date"]=["2025-06-18",""]
    old["original_source_date_provenance"]=["old","missing"]
    with pytest.raises(ValueError,match="source calendar"):
        m.append(old,[fake_photo("20","wave4","","red_pink")])


def test_no_false_completed_photo_labels_from_blind_source_only():
    assert m.NEW_PIXEL_N==4096
    assert m.MEASURED_N==309439
    assert m.MEASURED_HUES==149662
    assert m.MEASURED_DATED_HUES==149653
