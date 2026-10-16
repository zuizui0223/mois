#!/usr/bin/env python3
"""REAL global latitude-zone climate seasonality on byte-pinned WorldClim normals.

This is an independent climate reference frame, NOT sampled iNaturalist flower
colour or plant occurrence inference. Do not convert area-weighted geographic
cell fractions into frequency of plants, flowers, species or photograph uploads.
Only original historical monthly normals for 1970–2000 are analysed.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))
from scripts.analysis.mois_monthly_rainfall_flower_colour_20261016 import (
    MonthlyNormals, sha256_file, NORMALS
)

EXPECTED_MANIFEST_SHA256="07f45ac2973589d8951ac3eb0fa1c81b700018c5ab7322e2346dcf7d6723ee23"
SCHEMA="mois_global_historical_monthly_precipitation_land_area_seasonality_v1"
STATUS="REAL_WORLDCLIM_ONLY_NO_PLANT_PHOTO_COLOUR_OR_OBSERVER_INFERENCE"
RADIUS_KM=6371.0088
CV_LIMIT=0.15
MONTHLY_SPAN_LIMIT_MM=30
DOMAINS=("TROPICAL","NORTH_EXTRATROPICAL","SOUTH_EXTRATROPICAL")


def _finite_fraction(numerator:float,denominator:float):
    return numerator/denominator if denominator else None


def climate_atlas(normals, *, synthetic=False):
    """Return only coarse area-weighted summaries of true monthly climate.

    Valid mask follows the same original land/no-data pixels in all 12
    provider layers. Each cell is weighted by spherical area, not by raw
    pixel count or iNaturalist photographic observation density.
    """
    data=np.asarray(normals.monthly,dtype=np.float64)
    if data.ndim!=3 or data.shape[0]!=12 or data.shape[1:]!=(normals.height,normals.width):
        raise ValueError("Expected 12 aligned source monthly rain GeoTIFF arrays")
    if not synthetic and normals.manifest_sha256!=EXPECTED_MANIFEST_SHA256:
        raise ValueError("Original verified 2026-10-16 12-raster climate manifest drift")
    transform=normals.transform
    if abs(transform.b)>1e-12 or abs(transform.d)>1e-12 or transform.a<=0 or transform.e>=0:
        raise ValueError("Only north-up original WGS84 source rainfall grids are admitted")
    valid=np.all(np.isfinite(data)&(data>=0)&(data<=10000),axis=0)
    if not np.any(valid):
        raise ValueError("No valid source WorldClim land pixels with 12 months")
    nrows,ncols=valid.shape
    lat_north=transform.f+np.arange(nrows)*transform.e
    lat_south=lat_north+transform.e
    lat_mid=(lat_north+lat_south)/2
    if np.any((lat_north>90.00001)|(lat_south< -90.00001)):
        raise ValueError("Source global climate grid extends outside WGS84 latitude")
    dlon=math.radians(transform.a)
    row_area=(RADIUS_KM**2)*dlon*np.abs(
        np.sin(np.deg2rad(lat_north))-np.sin(np.deg2rad(lat_south)))
    # Broadcast row area to all lon cells without over-weighting high latitude.
    km2=np.broadcast_to(row_area[:,None],valid.shape)
    annual=np.sum(data,axis=0)
    mean=annual/12
    maximum=np.max(data,axis=0)
    minimum=np.min(data,axis=0)
    cv=np.std(data,axis=0)/np.where(mean>0,mean,1)
    aseasonal=(cv<CV_LIMIT)|((maximum-minimum)<MONTHLY_SPAN_LIMIT_MM)
    seasonal=valid&(~aseasonal)
    # Equal highest-precipitation normal months are common where the source
    # rounds monthly data. Share a cell's area equally among peak months;
    # do NOT spuriously mark the first tied month as the rainy season.
    ties=(data==maximum[None,:,:])&valid[None,:,:]
    n_tied=np.maximum(1,np.sum(ties,axis=0))
    weighted_tied_area=km2/n_tied
    tropic=np.abs(lat_mid)<23.5
    north=lat_mid>=23.5
    south=lat_mid<=-23.5
    domain_masks={
        "TROPICAL":tropic,
        "NORTH_EXTRATROPICAL":north,
        "SOUTH_EXTRATROPICAL":south
    }
    all_valid_area=float(np.sum(km2[valid]))
    outputs={}
    for domain in DOMAINS:
        mask=valid&domain_masks[domain][:,None]
        mask_seasonal=seasonal&mask
        n=int(np.count_nonzero(mask))
        area=float(np.sum(km2[mask]))
        seasonal_area=float(np.sum(km2[mask_seasonal]))
        count_seasonal=int(np.count_nonzero(mask_seasonal))
        result={
            "n_valid_grid_cells":n,
            "valid_cell_area_million_km2":area/1e6,
            "share_of_all_valid_source_cell_area":_finite_fraction(area,all_valid_area),
            "n_seasonal_grid_cells":count_seasonal,
            "n_aseasonal_grid_cells":n-count_seasonal,
            "seasonal_area_fraction":_finite_fraction(seasonal_area,area),
            "aseasonal_area_fraction":_finite_fraction(area-seasonal_area,area),
            "mean_annual_precipitation_mm_area_weighted":_finite_fraction(
                float(np.sum(annual[mask]*km2[mask])),area),
            "area_fraction_of_seasonal_cells_peak_rainfall_in_calendar_quarter":{},
        }
        for quarter in range(1,5):
            peak_share=np.sum(ties[(quarter-1)*3:quarter*3],axis=0)
            result["area_fraction_of_seasonal_cells_peak_rainfall_in_calendar_quarter"][
                f"Q{quarter}"]=_finite_fraction(
                    float(np.sum((weighted_tied_area*peak_share)[mask_seasonal])),
                    seasonal_area)
        if domain!="TROPICAL":
            if domain=="NORTH_EXTRATROPICAL":
                summer=(6,7,8)
                winter=(12,1,2)
            else:
                summer=(12,1,2)
                winter=(6,7,8)
            summer_share=np.sum(ties[np.array(summer)-1],axis=0)
            winter_share=np.sum(ties[np.array(winter)-1],axis=0)
            summer_area=float(np.sum(
                (weighted_tied_area*summer_share)[mask_seasonal]))
            winter_area=float(np.sum(
                (weighted_tied_area*winter_share)[mask_seasonal]))
            result["seasonal_area_fraction_with_wettest_month_in_local_summer"]=_finite_fraction(
                summer_area,seasonal_area)
            result["seasonal_area_fraction_with_wettest_month_in_local_winter"]=_finite_fraction(
                winter_area,seasonal_area)
            rain_summer=np.sum(data[np.array(summer)-1],axis=0)
            rain_winter=np.sum(data[np.array(winter)-1],axis=0)
            # Require more than ten percentage points of normal annual rain
            # above the other 3-month period; classification of climate, not
            # season-specific observation or flowering.
            summer_dominant=mask_seasonal&(
                (rain_summer-rain_winter)>0.1*annual)
            winter_dominant=mask_seasonal&(
                (rain_winter-rain_summer)>0.1*annual)
            result["seasonal_area_fraction_summer_rainfall_over_winter_by_annual_10pp"]=_finite_fraction(
                float(np.sum(km2[summer_dominant])),seasonal_area)
            result["seasonal_area_fraction_winter_rainfall_over_summer_by_annual_10pp"]=_finite_fraction(
                float(np.sum(km2[winter_dominant])),seasonal_area)
        outputs[domain]=result
    if sum(x["n_valid_grid_cells"] for x in outputs.values())!=int(np.count_nonzero(valid)):
        raise AssertionError("Latitude-domain area partition failed")
    if not math.isclose(sum(x["valid_cell_area_million_km2"] for x in outputs.values()),
                        all_valid_area/1e6,rel_tol=1e-12):
        raise AssertionError("Latitude-domain source geospatial area did not sum")
    return {
        "schema":SCHEMA,
        "status":STATUS,
        "source_original_monthly_manifest_sha256":normals.manifest_sha256,
        "original_12_monthly_raster_sha256":normals.raster_hashes,
        "climate_source":NORMALS,
        "climatology_years":"1970-2000",
        "grid_resolution_degrees":abs(float(transform.a)),
        "area_weighting":"spherical geographic grid area; excludes missing source cells",
        "source_grid_n_columns":ncols,
        "source_grid_n_rows":nrows,
        "source_grid_n_valid_complete_12month_cells":int(np.count_nonzero(valid)),
        "source_complete_valid_grid_cell_area_million_km2":all_valid_area/1e6,
        "aseasonal_cv_threshold":CV_LIMIT,
        "aseasonal_max_minus_min_monthly_mm_threshold":MONTHLY_SPAN_LIMIT_MM,
        "region_definitions":"TROPICAL if abs(lat)<23.5; northern extratropical if lat>=23.5; southern extratropical if lat<=-23.5",
        "wettest_month_ties_split_equally_across_tied_months_for_aggregate_area":True,
        "latitude_domains":outputs,
        "normal_precipitation_not_observed_weather_at_photo_date":True,
        "flora_phenotype_occurrences_or_fitness_observed":False,
        "iNaturalist_photo_colour_or_annotation_counts_observed":False,
        "source_42111_taxa_13predictor_refitted":False,
        "synthetic_test_only":bool(synthetic),
        "note":"Independent real rainfall-only spatial availability evidence for seasonal study design. Climate cell areas do NOT approximate the frequency of species or actual blooming flowers."
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--worldclim-monthly-manifest",type=Path,required=True)
    p.add_argument("--manifest-sha256",required=True)
    p.add_argument("--out-aggregate",type=Path,required=True)
    a=p.parse_args()
    if a.out_aggregate.exists():
        raise ValueError("Do not overwrite a source-audited global seasonal climatology aggregate")
    if a.manifest_sha256!=EXPECTED_MANIFEST_SHA256:
        raise ValueError("Original monthly climatology manifest must match admitted 12-raster source")
    normals=MonthlyNormals(a.worldclim_monthly_manifest,a.manifest_sha256)
    result=climate_atlas(normals)
    a.out_aggregate.parent.mkdir(parents=True,exist_ok=True)
    a.out_aggregate.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("MOIS_REAL_WORLDCLIM_GLOBAL_SEASONALITY_ATLAS",
          json.dumps(result,sort_keys=True))


if __name__=="__main__":
    main()
