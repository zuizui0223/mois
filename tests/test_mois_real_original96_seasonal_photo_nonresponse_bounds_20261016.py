"""Original 96 photo seasonal classifier success and unclassified photo bounds.

Real original pixel classes, not fictitious new flower labels. Hypothetical
unclassified-photo binary assignment MUST NOT be called a wild-flower result.
"""
from __future__ import annotations
import importlib.util
import json
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location(
    "mois_actual_96_mnar_photo_seasons",
    ROOT/"scripts/analysis/mois_real_original96_seasonal_photo_nonresponse_bounds_20261016.py")
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
ORIGINAL=ROOT/"results/mois_actual_original96_photo_monthly_tavg_srad_prec_vpd_technical_pilot_20261016.json"
PUBLIC=ROOT/"results/mois_real_original96_seasonal_photo_nonresponse_bounds_20261016.json"


def test_real_source_replays_exact_photo_classifiability_and_hypothetical_missingness():
    first=json.loads(ORIGINAL.read_text(encoding="utf8"))
    stored=json.loads(PUBLIC.read_text(encoding="utf8"))
    recompute=m.derive(first)
    assert recompute==stored
    assert stored["all_original_photographs"]==96
    assert stored["original_roi_classifiable"]==20
    assert stored["original_roi_white"]==3
    assert stored["photographs_without_classifiable_4hue"]==76
    assert stored["missing_photo_binary_colour_assignments_not_valid_plant_flowering_traits"] is True
    assert stored["overall_image_pilot_is_not_random_42111_taxon_photographic_census"] is True
    assert "photo_id" not in json.dumps(stored)
    assert "latitude" not in json.dumps(stored)


def test_spring_vs_summer_observer_and_ROI_selection_ambiguous_colour_contrast():
    stored=json.loads(PUBLIC.read_text(encoding="utf8"))
    group={(z["domain"],z["local_season"]):z for z in stored[
        "source_season_counts_and_logical_missing_photo_bounds"]}
    spring=group["NORTH_EXTRATROPICAL","SPRING"]
    summer=group["NORTH_EXTRATROPICAL","SUMMER"]
    assert (spring["source_photo_opportunities"],
            spring["four_hue_classifiable_photographs"],
            spring["classified_white_photographs"])==(45,7,1)
    assert (summer["source_photo_opportunities"],
            summer["four_hue_classifiable_photographs"],
            summer["classified_white_photographs"])==(25,9,1)
    assert spring["roi_four_colour_classifiability_fraction"]==pytest.approx(7/45)
    assert summer["roi_four_colour_classifiability_fraction"]==pytest.approx(9/25)
    old=stored["north_spring_minus_summer"]
    assert old["source_white_given_classifiable_difference_percentage_points"]==pytest.approx(
        100*(1/7-1/9))
    assert old["worst_case_all_photo_binary_assignment_difference_percentage_points"]==pytest.approx([
        100*(1/45-17/25),
        100*(39/45-1/25),
    ])
    assert old["both_signs_of_raw_all_photo_binary_assignment_difference_possible"] is True
    assert old["is_comparison_of_same_taxa_and_same_geographic_cell"] is False
    assert old["no_temperature_solar_rainfall_VPD_causal_effect_estimated"] is True


def test_original_image_status_or_count_forgery_fails():
    original=json.loads(ORIGINAL.read_text())
    original["n_roi_original_classified_white"]=4
    with pytest.raises(ValueError,match="not admitted"):
        m.derive(original)
    original=json.loads(ORIGINAL.read_text())
    original["observed_season_rows"][0]["classified_white"]=1
    with pytest.raises(ValueError,match="denominator"):
        m.derive(original)
