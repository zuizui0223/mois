"""Validate the committed REAL 127504-photo result receipt and figure source.

This test uses already source-frozen numeric evidence, never invented
synthetic data or a new downloaded iNaturalist source observation.
"""
import csv
import json
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/"results/mois_real127504_photo_white_classifiability_season_selection_20261017.json"
FIG=ROOT/"results/mois_real127504_photo_white_classifiability_season_figure_data_20261017.csv"


def sources():
    receipt=json.loads(RECEIPT.read_text(encoding="utf-8"))
    with FIG.open(newline="",encoding="utf-8") as f:
        data=list(csv.DictReader(f))
    return receipt,data


def test_frozen_real_source_photo_denominators_not_missing_equals_nonwhite():
    j,_=sources()
    assert j["original_source_photo_csv_sha256"]=="fbb54b58dc91e40e2e166081212588acc18adb19af4bdbcd1c7fe1778bee28d9"
    assert j["total_original_attempted_photo_ids"]==267425
    assert j["total_four_hue_original_photo_ids"]==129795
    assert j["date_original_photo_attempts"]==261578
    assert j["date_and_classifiable_photo_hues"]==127504
    assert j["date_and_technical_unknown_photo_hues"]==134074
    assert j["date_and_classifiable_photo_hues"]+j["date_and_technical_unknown_photo_hues"]==j["date_original_photo_attempts"]


@pytest.mark.parametrize("hemisphere",["north_extratropical","south_extratropical"])
def test_white_detected_per_photo_attempt_symmetric_identity_exact(hemisphere):
    j,_=sources()
    g=j["global_hemisphere"][hemisphere]
    s_attempt=g["summer_original_photo_attempts"]
    w_attempt=g["winter_original_photo_attempts"]
    s_label=g["summer_four_hue_classifiable"]
    w_label=g["winter_four_hue_classifiable"]
    s_white=g["summer_white_labeled"]
    w_white=g["winter_white_labeled"]
    c_s=s_label/s_attempt
    c_w=w_label/w_attempt
    p_s=s_white/s_label
    p_w=w_white/w_label
    total=100*(c_s*p_s-c_w*p_w)
    hue=100*(c_s+c_w)/2*(p_s-p_w)
    classifiability=100*(p_s+p_w)/2*(c_s-c_w)
    assert total == pytest.approx(g["summer_minus_winter_recorded_white_labels_per_attempt_pp"],abs=1e-5)
    assert hue == pytest.approx(g["recorded_white_label_yield_symmetric_photo_hue_component_pp"],abs=1e-5)
    assert classifiability == pytest.approx(g["recorded_white_label_yield_symmetric_classifiability_component_pp"],abs=1e-5)
    assert total == pytest.approx(hue+classifiability)


def test_eight_figure_rows_from_signed_photo_counts_and_independent_matched_localities():
    j,table=sources()
    assert len(table)==8
    assert set(r["photographic_estimand"] for r in table)=={
        "technical_four_hue_classifiability",
        "white_fraction_conditional_on_classified"}
    for zone in ("north_extratropical","south_extratropical"):
        rows=[r for r in table if r["hemisphere"]==zone]
        assert len(rows)==4
        g=j["global_hemisphere"][zone]
        l=j["matched_0p25_degree_cell_same_meteorological_season_year_min5_classified_observations_min2_taxa_each_season"][zone]
        for r in rows:
            if r["source_scope"]=="unpaired_global_source":
                key=("summer_minus_winter_classifiability_pp" if
                    r["photographic_estimand"]=="technical_four_hue_classifiability"
                    else "summer_minus_winter_white_share_conditional_classified_pp")
                assert float(r["summer_minus_winter_pp"])==pytest.approx(g[key])
                assert int(r["dated_original_photo_attempts"])==(g["summer_original_photo_attempts"]+g["winter_original_photo_attempts"])
                assert not r["bootstrap95_low_pp"] and not r["bootstrap95_high_pp"]
            else:
                assert int(r["site_season_year_comparisons"])==l["eligible_site_years"]
                assert int(r["dated_original_photo_attempts"])==l["all_attempted_photos_in_paired_seasons"]
                assert int(r["classifiable_original_photos"])==l["classifiable_hue_photo_ids_in_paired_seasons"]
                if r["photographic_estimand"]=="technical_four_hue_classifiability":
                    assert float(r["summer_minus_winter_pp"])==pytest.approx(l["classifiability_summer_minus_winter_pp"])
                    expected=l["classifiability_geocell_bootstrap999_percentile95_pp"]
                else:
                    assert float(r["summer_minus_winter_pp"])==pytest.approx(l["white_given_classified_summer_minus_winter_pp"])
                    expected=l["white_given_classified_geocell_bootstrap999_percentile95_pp"]
                assert float(r["bootstrap95_low_pp"])==pytest.approx(expected[0],abs=1e-5)
                assert float(r["bootstrap95_high_pp"])==pytest.approx(expected[1],abs=1e-5)
