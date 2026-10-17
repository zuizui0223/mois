"""Scientific source checksum and inference guard for actual 309439 water fit.

The complete 309439 original-photo/photo-site 999-genus run is archived in the
successful original MOIS GitHub Actions; this test verifies the compact
source-locked, authentic numeric receipt and manuscript all share one source.
"""
from __future__ import annotations
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/"results/mois_real309439_fullFCP256_six_climate_maxT999_result_20261017.json"
PAPER=ROOT/"docs/MOIS_GLOBAL_FLORAL_COLOUR_MOISTURE_AND_SEASONAL_SAMPLING_MANUSCRIPT_20261017.md"


def source():
    return json.loads(RECEIPT.read_text(encoding="utf-8"))


def test_complete_photo_and_worldclim_original_source_population_sha():
    r=source()
    assert r["genuinely_measured_original_photo_ids"]==309439
    assert r["four_hue_original_ROI_photo_outcomes"]==149662
    assert r["actual_original_dated_four_hue_photo_id_count"]==149653
    assert r["new_real_photo_ids_point_sampled_on_old_295017_worldclim"]==14422
    assert r["new_classifiable_four_hue_original_photo_ids"]==7124
    assert r["genus_resample_draws"]==999
    assert r["simultaneous_hue_environment_comparisons"]==24
    assert r["source_original_309439_measured_photo_hues_GZ_SHA256"]=="546266e2bbaa340a2523864c28c5a877842eb4b0904a7e0dff90870c6519a47e"
    assert r["source_old_295017_photo_original_climate_GZ_SHA256"]=="84c651aeaacf0e317445593c249a85940d230a19ec8b7f7223ec0e46553013f7"
    assert r["actual_full_original_309439_photo_climate_GZ_SHA256"]=="e8cf9ab1c935b906731b6c783a45193c9c13cca8cb226fe46c4335d0e1c75cee"
    assert r["full_actual_309439_model_decision_JSON_SHA256"]=="c925fd8d318e21169989f7667e9f87993714adefbabfde293b61c68758555da9"
    assert len(r["prespecified_environment_anchors"])==6
    assert len(r["photographic_hue_outcomes"])==4


def test_within_one_maxT_family_positive_globally_but_not_geographically():
    r=source()
    global_anchors=(
        "white_annual_precipitation_global_genus_equal",
        "white_atmospheric_vapour_annual_global_genus_equal",
    )
    conditioned=(
        "white_annual_precipitation_coarse_geography_conditioned",
        "white_annual_vapour_coarse_geography_conditioned",
    )
    for key in global_anchors:
        v=r[key]
        assert v["point_pp_per_original_env_SD"]>0
        assert 0 < v["simultaneous_maxT_95CI_pp"][0] < v["point_pp_per_original_env_SD"] < v["simultaneous_maxT_95CI_pp"][1]
        assert v["maxT_excludes_zero"] is True
    for key in conditioned:
        v=r[key]
        low,high=v["simultaneous_maxT_95CI_pp"]
        assert low < 0 < high
        assert v["maxT_excludes_zero"] is False
    assert round(r[global_anchors[0]]["point_pp_per_original_env_SD"],3)==1.988
    assert round(r[global_anchors[1]]["point_pp_per_original_env_SD"],3)==3.161


def test_actual_full309439_manuscript_includes_both_source_and_geography_effects():
    s=PAPER.read_text(encoding="utf-8")
    assert "309,439" in s and "149,662" in s and "149,653" in s
    assert "### 3.1." in s
    assert "1.988" in s and "3.161" in s
    assert "−2.635" in s and "−2.640" in s
    assert "999" in s and "1970–2000" in s
    assert "not" in s.lower()
    assert "FCP" in s
