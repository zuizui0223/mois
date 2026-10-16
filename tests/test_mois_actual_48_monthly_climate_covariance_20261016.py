"""Test geography-centred seasonal climate collinearity, no iNaturalist pixels."""
import importlib.util
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,ROOT/path)
    o=importlib.util.module_from_spec(spec);spec.loader.exec_module(o);return o
fixture=load("mois_cov_fake48_source","tests/test_mois_join_monthly_dynamic_climate_20261016.py")
m=load("mois_worldclim48_actual_climate_cov","scripts/analysis/mois_actual_48_monthly_climate_covariance_20261016.py")


def test_all_four_climate_variables_area_weighted_within_grid_and_zone_month(tmp_path):
    normals=fixture.normal(tmp_path)
    report=m.audit(normals,synthetic=True)
    assert report["status"]==m.STATUS
    assert report["source_complete_all_48_monthly_climate_cells"]==23
    assert report["original_fcp_sampled_flower_photographs_joined"] is False
    assert report["real_13_env_genus_spatial_colour_model_reestimated"] is False
    assert report["month_climate_normal_not_observed_weather"] is True
    assert set(report["latitude_domains"])==set(m.DOMAINS)
    assert report["latitude_domains"]["NORTH_EXTRATROPICAL"][
        "n_complete_48layer_original_source_cells"]==6
    assert report["latitude_domains"]["TROPICAL"][
        "n_complete_48layer_original_source_cells"]==12
    assert report["latitude_domains"]["SOUTH_EXTRATROPICAL"][
        "n_complete_48layer_original_source_cells"]==5
    for region in m.DOMAINS:
        z=report["latitude_domains"][region]
        within=z["model_A_after_exact_cell_fixed_effects"]
        adjusted=z["model_B_after_cell_and_latitude_zone_calendar_month_fixed_effects"]
        assert within["n_approx_independent_dynamic_axes"]<=4
        assert within["status"]=="RANK_DEFICIENT_DYNAMIC_CLIMATE_PROFILES"
        assert within["correlation_matrix"][0][1]==pytest.approx(1,abs=1e-9)
        assert adjusted["status"]=="RANK_DEFICIENT_AFTER_FIXED_EFFECTS"
        assert adjusted["n_approx_independent_dynamic_axes"]==0
        assert z["all_months_weighted_equally_not_photo_observer_calendar"] is True
    assert report["synthetic_test_only"] is True


def test_actual_48_band_climatology_provenance_must_precede_global_climate_statistics(tmp_path):
    normals=fixture.normal(tmp_path)
    with pytest.raises(ValueError,match="byte-pinned"):
        m.audit(normals,synthetic=False)
