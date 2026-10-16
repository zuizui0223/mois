"""Synthetic latitudinal climate-seasonality area audit; no plant-photo results."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]


def source(name,relpath):
    spec=importlib.util.spec_from_file_location(name,ROOT/relpath)
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


m=source("mois_real_monthly_atlas","scripts/analysis/mois_real_worldclim_monthly_climate_atlas_20261016.py")
grid=source("mois_climatology_grid_fixture","tests/test_mois_monthly_rainfall_flower_colour_20261016.py")


def test_original_geo_grid_area_and_hemispheric_seasonality(tmp_path):
    normals,_,_=grid.raster_fixture(tmp_path)
    result=m.climate_atlas(normals,synthetic=True)
    assert result["status"]==m.STATUS
    assert result["source_grid_n_valid_complete_12month_cells"]==23
    assert result["source_grid_n_columns"]==6
    assert result["source_grid_n_rows"]==4
    zones=result["latitude_domains"]
    assert zones["NORTH_EXTRATROPICAL"]["n_valid_grid_cells"]==6
    assert zones["TROPICAL"]["n_valid_grid_cells"]==12
    assert zones["SOUTH_EXTRATROPICAL"]["n_valid_grid_cells"]==5
    assert zones["NORTH_EXTRATROPICAL"]["n_seasonal_grid_cells"]==5
    assert zones["TROPICAL"]["n_seasonal_grid_cells"]==10
    assert zones["SOUTH_EXTRATROPICAL"]["n_seasonal_grid_cells"]==4
    assert sum(x["n_valid_grid_cells"] for x in zones.values())==23
    assert sum(x["valid_cell_area_million_km2"] for x in zones.values())==pytest.approx(
        result["source_complete_valid_grid_cell_area_million_km2"],rel=1e-13)
    for x in zones.values():
        assert x["mean_annual_precipitation_mm_area_weighted"]==pytest.approx(660)
        assert x["seasonal_area_fraction"]+x["aseasonal_area_fraction"]==pytest.approx(1)
        assert sum(x["area_fraction_of_seasonal_cells_peak_rainfall_in_calendar_quarter"].values())==pytest.approx(1)
        assert x["area_fraction_of_seasonal_cells_peak_rainfall_in_calendar_quarter"]["Q3"]==pytest.approx(.5)
        assert x["area_fraction_of_seasonal_cells_peak_rainfall_in_calendar_quarter"]["Q4"]==pytest.approx(.5)
    north=zones["NORTH_EXTRATROPICAL"]
    south=zones["SOUTH_EXTRATROPICAL"]
    assert north["seasonal_area_fraction_with_wettest_month_in_local_summer"]==pytest.approx(1/3)
    assert north["seasonal_area_fraction_with_wettest_month_in_local_winter"]==pytest.approx(1/6)
    assert south["seasonal_area_fraction_with_wettest_month_in_local_summer"]==pytest.approx(1/6)
    assert south["seasonal_area_fraction_with_wettest_month_in_local_winter"]==pytest.approx(1/3)
    assert north["seasonal_area_fraction_summer_rainfall_over_winter_by_annual_10pp"]==pytest.approx(1)
    assert south["seasonal_area_fraction_winter_rainfall_over_summer_by_annual_10pp"]==pytest.approx(1)
    assert result["normal_precipitation_not_observed_weather_at_photo_date"] is True
    assert result["flora_phenotype_occurrences_or_fitness_observed"] is False
    assert result["iNaturalist_photo_colour_or_annotation_counts_observed"] is False
    assert result["synthetic_test_only"] is True


def test_real_climate_source_identity_cannot_be_forged_by_a_synthetic_grid(tmp_path):
    normals,_,_=grid.raster_fixture(tmp_path)
    with pytest.raises(ValueError,match="climate manifest drift"):
        m.climate_atlas(normals)
    normals.monthly=np.full((12,4,6),np.nan)
    with pytest.raises(ValueError,match="No valid source WorldClim"):
        m.climate_atlas(normals,synthetic=True)


def test_monthly_ties_are_split_evenly_not_allocated_to_earliest_month(tmp_path):
    normals,_,_=grid.raster_fixture(tmp_path)
    out=m.climate_atlas(normals,synthetic=True)
    tropics=out["latitude_domains"]["TROPICAL"]
    # Original synthetic has a six-month rainy plateau (July--December).
    # Arbitrary argmax would count all peak land area in Q3.
    assert tropics["area_fraction_of_seasonal_cells_peak_rainfall_in_calendar_quarter"]["Q3"]==pytest.approx(.5)
    assert tropics["area_fraction_of_seasonal_cells_peak_rainfall_in_calendar_quarter"]["Q4"]==pytest.approx(.5)
    assert out["wettest_month_ties_split_equally_across_tied_months_for_aggregate_area"] is True
