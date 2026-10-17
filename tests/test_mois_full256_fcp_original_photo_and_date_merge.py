"""Synthetic source integrity tests for finished original 49,800-photo import.

Synthetic fixtures are strictly validation inputs; no simulated photo count or
white fraction is ever published as an actual ecological measurement.
"""
from __future__ import annotations
import pandas as pd
import pytest
from scripts.analysis.mois_full256_fcp_original_photo_and_date_merge import (
    append_finished_pixels, _date_real,
)
from scripts.analysis.mois_import_full_FCP_wave4_completed_original_photos import (
    BASE_COLS,
)

HUES=("white","yellow_orange","red_pink","blue_purple")


def source(pid, hue="white", *, taxon="501", day="2026-08-10", raw_blank=False):
    r={col:"" for col in BASE_COLS}
    r.update(photo_id=str(pid),inat_taxon_id=taxon,species="source specimen",
       morph=hue,source_original_measured="fcp_wave4_partial_original_ROI_v4",
       old_mois_or_new_fcp="fcp_additional_ROI",
       observation_id=str(int(pid)+150),observed_on="" if raw_blank else day,
       latitude="40.1",longitude="136.4",photo_license="cc-by-nc",
       measurement_status="classified_four_state_morph" if hue in HUES else
                          "not_evaluable_roi_or_flip_gate")
    return r


def originals():
    old=pd.DataFrame([source("10","red_pink",raw_blank=True),
        source("11","white"),source("12","mixed_uncertain")])
    old["qualified_date"]=["2026-08-10","2026-08-10","2026-08-10"]
    old["original_source_date_provenance"]=["historic_s3","old","old"]
    whole=pd.DataFrame([
        source("10","red_pink"),
        source("11","white"),
        source("12","mixed_uncertain"),
        source("13","white"),
        source("14","mixed_uncertain")])
    whole["global_classifiable"]=["True","True","False","True","False"]
    return old,whole


def test_completed_photo_import_preserves_all_original_pixel_date_rows():
    old,whole=originals()
    orig=old.copy(deep=True)
    merged,stats=append_finished_pixels(old,whole,strict=False)
    assert len(merged)==5 and merged.photo_id.is_unique
    assert merged.iloc[:3].reset_index(drop=True).equals(orig)
    assert stats["new_original_genuinely_measured_photo_pixels"]==2
    assert stats["new_original_four_hue_ROI_photo_labels"]==1
    assert stats["new_date_and_four_hue_original_photo_ROIs"]==1
    assert stats["original_FCP_pixels_overlapping_previous_source_photo_IDs"]==3
    assert stats["source_original_observed_on_blank_but_independently_qualified_days_agree"]==1
    assert merged.iloc[-1].morph=="mixed_uncertain"
    assert merged.iloc[-1].qualified_date=="2026-08-10"


@pytest.mark.parametrize("col,replacement",[
    ("morph","white"),("inat_taxon_id","777"),
    ("observation_id","999991"),("observed_on","2026-08-11"),
    ("latitude","39.5")])
def test_refuse_same_photo_different_original_identity(col,replacement):
    old,whole=originals()
    whole.loc[whole.photo_id.eq("10"),col]=replacement
    with pytest.raises(ValueError):
        append_finished_pixels(old,whole,strict=False)


def test_provider_date_may_be_missing_in_old_raw_but_must_match_verified_day():
    old,whole=originals()
    old.loc[0,"qualified_date"]="2026-09-01"
    with pytest.raises(ValueError,match="independently"):
        append_finished_pixels(old,whole,strict=False)


def test_source_pixel_photo_must_have_valid_actual_day_and_status():
    old,whole=originals()
    whole.loc[whole.photo_id.eq("13"),"observed_on"]="2026-02-30"
    with pytest.raises(ValueError,match="no original source date"):
        append_finished_pixels(old,whole,strict=False)
    old,whole=originals()
    whole.loc[whole.photo_id.eq("13"),"global_classifiable"]="False"
    with pytest.raises(ValueError,match="classification status"):
        append_finished_pixels(old,whole,strict=False)


def test_duplicate_photo_ids_and_no_new_pixels_are_not_data_growth():
    old,whole=originals()
    repeated=pd.concat([whole,whole.iloc[:1]],ignore_index=True)
    with pytest.raises(ValueError,match="duplicated"):
        append_finished_pixels(old,repeated,strict=False)
    with pytest.raises(ValueError,match="do not add"):
        append_finished_pixels(old,whole.iloc[:3],strict=False)


def test_photo_observation_date_is_not_a_weather_or_flowering_annotation():
    assert _date_real("2026-09-20")
    assert not _date_real("2028-12-05")
    assert not _date_real("2026-02-30")
    assert not _date_real("arbitrary month")
