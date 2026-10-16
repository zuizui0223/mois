"""Verify original photographic classifiability planning; no image-model inference."""
import importlib.util
import json
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
s=importlib.util.spec_from_file_location(
    "mois_two_season_roi_yield_design",
    ROOT/"scripts/analysis/mois_seasonal_photo_classifiability_power_20261016.py")
m=importlib.util.module_from_spec(s)
s.loader.exec_module(m)

SOURCE=ROOT/"results/mois_seasonal_photo_classifiability_power_20261016.json"


def test_existing_source_count_receipt_is_exact_and_calculations_replay():
    data=json.loads(SOURCE.read_text())
    calculated=m.report()
    assert data["status"]==calculated["status"]==m.STATUS
    assert data["source_historical_photographs_unique"]==100543
    assert data["historical_four_hue_classifiable_photographs"]==45953
    assert data["historical_photo_records_unclassifiable"]==54590
    assert data["target_classifiable_in_each_season"]==3
    assert data["proposal_min_original_selected_observations_per_occupied_season_when_quota_allows"]==10
    assert data["no_new_photo_or_phenology_annotation_labels_measured"] is True
    assert data["no_seasonal_white_climate_effect_computed"] is True
    for key in calculated["same_historical_classifiability_rate_scenarios"]:
        current=data["same_historical_classifiability_rate_scenarios"][key]
        expectation=calculated["same_historical_classifiability_rate_scenarios"][key]
        assert current["previously_classifiable_photo_fraction"]==pytest.approx(
            expectation["previously_classifiable_photo_fraction"],abs=1e-12)
        assert current["minimum_equal_original_observations_per_season_for_pair_success"]==expectation[
            "minimum_equal_original_observations_per_season_for_pair_success"]
        for a,b in zip(current["scenarios"],expectation["scenarios"]):
            assert a["independent_original_observations_per_season"]==b[
                "independent_original_observations_per_season"]
            assert a["probability_both_seasons_have_at_least_three_classifiable"]==pytest.approx(
                b["probability_both_seasons_have_at_least_three_classifiable"],abs=1e-12)
    current=data["unequal_original_dry_wet_classifiability_scenarios"]
    expected=calculated["unequal_original_dry_wet_classifiability_scenarios"]
    assert current["minimum_equal_source_observations_per_dry_and_wet_phase"]==expected[
        "minimum_equal_source_observations_per_dry_and_wet_phase"]
    for a,b in zip(current["scenarios"],expected["scenarios"]):
        assert a["probability_both_have_at_least_three_classifiable"]==pytest.approx(
            b["probability_both_have_at_least_three_classifiable"],abs=1e-12)


def test_historical_source_ratio_scenarios_justify_a_ten_per_phase_floor():
    p=45953/100543
    assert m.probability_two_seasons(3,3,p,p)==pytest.approx(.009,abs=.001)
    assert m.probability_two_seasons(10,10,p,p)>0.8
    assert m.probability_two_seasons(9,9,p,p)<0.8
    assert m.minimum_per_phase(p,p,.8)==10
    assert m.minimum_per_phase(p,p,.9)==12
    assert m.minimum_per_phase(3904/9990,3904/9990,.9)==14
    assert m.minimum_per_phase(4981/10084,4981/10084,.8)==9
    # Wet and dry REGION deciles are not necessarily the same as the local
    # rainy/dry season at an individual original plant-photo location.
    p_dry=4981/10084
    p_wet=3904/9990
    assert m.minimum_per_phase(p_dry,p_wet,.8)==11
    assert m.minimum_per_phase(p_dry,p_wet,.9)==13
    assert m.probability_two_seasons(10,10,p_dry,p_wet)<.8
    assert m.probability_two_seasons(13,13,p_dry,p_wet)>.9


def test_photo_classifiability_probability_is_monotone_and_not_a_biological_pvalue():
    for p in (.3,.457,.7):
        ps=[m.probability_two_seasons(n,n,p,p) for n in range(3,21)]
        assert all(0<=x<=1 for x in ps)
        assert ps==sorted(ps)
    with pytest.raises(ValueError):
        m.classified_probability(3,-0.01)
    with pytest.raises(ValueError):
        m.minimum_per_phase(.4,.5,1)
