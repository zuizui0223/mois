"""Source-result guard: real 96 image ROI×48 monthly climate, strictly pilot."""
from __future__ import annotations
import json
import math
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
F=ROOT/"results/mois_actual_original96_photo_monthly_tavg_srad_prec_vpd_technical_pilot_20261016.json"
DYN=ROOT/"results/mois_worldclim21_dynamic_48_real_source_receipt_20261016.json"

def test_actual_original_image_roi_method_and_true_source_monthly_climate_ids():
    r=json.loads(F.read_text(encoding="utf8"))
    source=json.loads(DYN.read_text(encoding="utf8"))
    assert r["status"]=="ACTUAL_96_IMAGE_ROI_TECHNICAL_PILOT_CLIMATE_JOIN_NOT_GLOBAL_SEASON_EFFECT"
    assert r["source_original_FCP_96_pilot_csv_sha256"]=="911ec9edd10a7126bee28cfa67738cef0ff4d275b1f0ca1c06e3f26e7e290039"
    assert r["source_48_worldclim_monthly_manifest_sha256"]==source["original_48_raster_manifest_sha256"]
    assert r["verified_original_FCP_photo_climate_execution_url"]=="https://github.com/zuizui0223/mois/actions/runs/38849479497"
    assert r["n_original_source_photo_IDs"]==96
    assert r["n_original_distinct_observations"]==96
    assert r["n_original_taxa"]==4
    assert r["n_roi_original_four_colour_classifiable"]==20
    assert r["n_roi_original_classified_white"]==3
    assert r["n_monthly_source_site_climate_joined"]==96
    assert r["n_monthly_source_site_climate_missing"]==0
    assert r["roi_status_counts"]["not_evaluable_roi_or_flip_gate"]==67
    assert r["roi_status_counts"]["not_evaluable_ambiguous_palette_composition"]==9
    assert r["original_photo_selection_is_nonrepresentative_four_taxon_technical_pilot"] is True
    assert r["original_source_population_42111_taxa_photo_IDs_joined"] is False
    assert r["seasonal_climate_white_association_or_fitness_estimated"] is False

def test_original_pilot_16_classified_climate_cells_but_only_three_month_repeats():
    r=json.loads(F.read_text(encoding="utf8"))
    assert r["real_10m_source_climate_cells_with_any_original_photo"]==60
    assert r["real_10m_source_climate_cells_with_original_photos_in_two_or_more_months"]==11
    assert r["real_10m_source_climate_cells_with_classifiable_photos"]==16
    assert r["real_10m_source_climate_cells_with_classifiable_photos_in_two_or_more_months"]==3
    assert r["cell_FE_primary_four_dynamic_months_estimable_in_this_pilot"] is False
    assert r["prespecified_site_cell_season_FE_within_month_estimated"] is False
    assert r["synthetic_test_only"] is False
    assert "photo_id" not in r and "coordinates" not in r

def test_real_photo_site_seasonal_white_tavg_srad_prec_vpd_are_descriptive_only():
    r=json.loads(F.read_text(encoding="utf8"))
    rows=r["observed_season_rows"]
    assert sum(v["sampled_original_photo_opportunities"] for v in rows)==96
    assert sum(v["four_hue_classifiable"] for v in rows)==20
    assert sum(v["classified_white"] for v in rows)==3
    north={v["local_season"]:v for v in rows if v["domain"]=="NORTH_EXTRATROPICAL"}
    assert {k:(v["four_hue_classifiable"],v["classified_white"]) for k,v in north.items()}=={
        "AUTUMN":(0,0),"SPRING":(7,1),"SUMMER":(9,1),"WINTER":(3,1)
    }
    assert north["SPRING"]["white_fraction_among_classifiable"]==pytest.approx(1/7)
    assert north["SUMMER"]["white_fraction_among_classifiable"]==pytest.approx(1/9)
    assert north["WINTER"]["white_fraction_among_classifiable"]==pytest.approx(1/3)
    s=north["SPRING"]["source_classified_monthly_climate_means"]
    u=north["SUMMER"]["source_classified_monthly_climate_means"]
    w=north["WINTER"]["source_classified_monthly_climate_means"]
    assert s["monthly_tavg_c"]==pytest.approx(13.538253852299281)
    assert u["monthly_tavg_c"]==pytest.approx(17.57780541314019)
    assert w["monthly_tavg_c"]==pytest.approx(10.171749750773111)
    assert s["monthly_srad_kj_m2_day"]==pytest.approx(22395.14285714286)
    assert u["monthly_srad_kj_m2_day"]==pytest.approx(26025.88888888889)
    assert w["monthly_srad_kj_m2_day"]==pytest.approx(10831.333333333334)
    assert s["monthly_prec_mm"]==pytest.approx(28.428571428571427)
    assert u["monthly_prec_mm"]==pytest.approx(9.)
    assert w["monthly_prec_mm"]==pytest.approx(87.)
    assert s["monthly_vpd_kpa"]==pytest.approx(.7313202705294432)
    assert u["monthly_vpd_kpa"]==pytest.approx(1.146363719791248)
    assert w["monthly_vpd_kpa"]==pytest.approx(.5222332167229968)
    # Means from different seasons' photos do NOT identify an adjusted
    # within-cell, causal, or interspecific global climate-colour slope.
    assert r["seasonal_climate_white_association_or_fitness_estimated"] is False
