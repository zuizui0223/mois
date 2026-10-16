"""Regression of real, original-source WorldClim 2.1 monthly climate atlas.

The data were actually calculated from 12 original GeoTIFFs in
https://github.com/zuizui0223/mois/actions/runs/38840741791.
It is climate CELL-AREA inference, NOT iNaturalist photo frequencies or
a newly measured flower-white phenomenon.
"""
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
ATLAS=ROOT/"results/mois_real_worldclim_monthly_climate_atlas_20261016.json"


def test_real_monthly_climate_receipt_has_exact_original_source_and_no_flowers():
    r=json.loads(ATLAS.read_text(encoding="utf-8"))
    assert r["schema"]=="mois_global_historical_monthly_precipitation_land_area_seasonality_v1"
    assert r["status"]=="REAL_WORLDCLIM_ONLY_NO_PLANT_PHOTO_COLOUR_OR_OBSERVER_INFERENCE"
    assert r["source_original_monthly_manifest_sha256"]=="07f45ac2973589d8951ac3eb0fa1c81b700018c5ab7322e2346dcf7d6723ee23"
    assert r["original_source_execution_url"]=="https://github.com/zuizui0223/mois/actions/runs/38840741791"
    assert r["source_grid_n_rows"]==1080 and r["source_grid_n_columns"]==2160
    assert r["source_grid_n_valid_complete_12month_cells"]==808053
    assert len(r["original_12_monthly_raster_sha256"])==12
    assert set(r["original_12_monthly_raster_sha256"])=={f"{m:02}" for m in range(1,13)}
    assert r["normal_precipitation_not_observed_weather_at_photo_date"] is True
    assert r["flora_phenotype_occurrences_or_fitness_observed"] is False
    assert r["iNaturalist_photo_colour_or_annotation_counts_observed"] is False
    assert r["source_42111_taxa_13predictor_refitted"] is False
    assert r["synthetic_test_only"] is False


def test_original_tropical_and_extratropical_seasonality_is_area_weighted():
    result=json.loads(ATLAS.read_text(encoding="utf-8"))
    dom=result["latitude_domains"]
    assert set(dom)=={"NORTH_EXTRATROPICAL","SOUTH_EXTRATROPICAL","TROPICAL"}
    assert sum(x["n_valid_grid_cells"] for x in dom.values())==808053
    assert math.isclose(
        sum(x["valid_cell_area_million_km2"] for x in dom.values()),
        result["source_complete_valid_grid_cell_area_million_km2"],rel_tol=1e-12
    )
    assert 150<result["source_complete_valid_grid_cell_area_million_km2"]<160
    expected={
        "TROPICAL":(.8873479301693107,158740),
        "NORTH_EXTRATROPICAL":(.756361611709869,386365),
        "SOUTH_EXTRATROPICAL":(.5509914470224445,262948),
    }
    for key,(fraction,n) in expected.items():
        row=dom[key]
        assert row["n_valid_grid_cells"]==n
        assert math.isclose(row["seasonal_area_fraction"],fraction,abs_tol=1e-11)
        assert math.isclose(row["seasonal_area_fraction"]+row["aseasonal_area_fraction"],1.,abs_tol=1e-11)
        assert math.isclose(sum(row["area_fraction_of_seasonal_cells_peak_rainfall_in_calendar_quarter"].values()),1.,abs_tol=1e-11)
    trop=dom["TROPICAL"]["area_fraction_of_seasonal_cells_peak_rainfall_in_calendar_quarter"]
    assert trop["Q1"]>.3 and trop["Q3"]>.3
    assert dom["NORTH_EXTRATROPICAL"]["seasonal_area_fraction_with_wettest_month_in_local_summer"]>.7
    assert dom["SOUTH_EXTRATROPICAL"]["seasonal_area_fraction_with_wettest_month_in_local_summer"]<.4
    assert result["wettest_month_ties_split_equally_across_tied_months_for_aggregate_area"] is True
