"""Actual original WorldClim 48-band geographic cell/month covariance receipt."""
import json
import math
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
F=ROOT/"results/mois_real_worldclim48_within_cell_monthly_climate_collinearity_20261016.json"


def test_actual_original_global_48_raster_monthly_climate_collinearity_not_botany():
    r=json.loads(F.read_text("utf8"))
    assert r["status"]=="REAL_MONTHLY_CLIMATE_ONLY_FIXED_CELL_SEASON_COLLINEARITY_NO_WHITE_PHOTOS"
    assert r["schema"]=="mois_real_worldclim48_cell_centered_monthly_climate_correlations_v1"
    assert r["verified_original_48_worldclim_climate_covariance_run"]=="https://github.com/zuizui0223/mois/actions/runs/38849962223"
    assert r["source_original_48_month_raster_manifest_sha256"]=="9d0cf37a64534542d7788debef28c13f106a15c3bf7bf91be0653f3ce106baae"
    assert r["source_complete_all_48_monthly_climate_cells"]==808053
    assert r["all_climate_axes"]==["tavg_c","srad_kj_m2_day","prec_mm","derived_vpd_kpa"]
    assert r["original_fcp_sampled_flower_photographs_joined"] is False
    assert r["real_13_env_genus_spatial_colour_model_reestimated"] is False
    assert r["observed_flowering_plant_frequency_measured"] is False
    assert r["synthetic_test_only"] is False


def test_true_temperature_solar_monthly_correlation_depends_on_time_and_cell_controls():
    r=json.loads(F.read_text("utf8"))
    b=r["latitude_domains"]
    assert sum(v["n_complete_48layer_original_source_cells"] for v in b.values())==808053
    expected={
        "NORTH_EXTRATROPICAL":(.8015069526430434,.2970699482472602),
        "SOUTH_EXTRATROPICAL":(.9288024713382301,.8412613991588237),
        "TROPICAL":(.676804667423698,.6669089206513864),
    }
    for name,(cell_only,cell_month) in expected.items():
        p=b[name]["model_A_after_exact_cell_fixed_effects"]
        q=b[name]["model_B_after_cell_and_latitude_zone_calendar_month_fixed_effects"]
        assert p["n_approx_independent_dynamic_axes"]==4
        assert q["n_approx_independent_dynamic_axes"]==4
        assert p["correlation_matrix"][0][1]==pytest.approx(cell_only,abs=1e-10)
        assert q["correlation_matrix"][0][1]==pytest.approx(cell_month,abs=1e-10)
        assert p["condition_number_correlation"]>1
        assert q["condition_number_correlation"]>1
        for m in (p["correlation_matrix"],q["correlation_matrix"]):
            for i in range(4):
                assert m[i][i]==pytest.approx(1,abs=1e-10)
                for j in range(4):
                    assert m[i][j]==pytest.approx(m[j][i],abs=1e-10)
                    assert -1.00000001<=m[i][j]<=1.00000001
    # These different climate-only slopes are NOT photo-hue partial slopes.
    assert r["original_fcp_sampled_flower_photographs_joined"] is False
