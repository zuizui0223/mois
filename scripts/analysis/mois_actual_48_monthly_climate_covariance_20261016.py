#!/usr/bin/env python3
"""REAL historical monthly climate collinearity AFTER exact geographic-cell FE.

WorldClim v2.1 ORIGINAL 48 monthly tavg/srad/prec/vapr rasters (1970--2000).
Area-weight climate CELLS, not observed flowers or iNaturalist uploads.
Never use a photo-derived calendar distribution to create climate normals.

Compute two distinct four-axis correlation matrices by latitude zone:
 A: residualize each variable on exact grid cell (12-month site mean);
 B: additionally residualize each variable on region × calendar month.
VPD is an approximation from monthly mean temperature and mean vapour pressure.
A correlation/condition-number check does NOT identify climate's effect on
white-flower frequencies; it supports the predeclared future model design.
"""
from __future__ import annotations
import argparse
import json
import math
import sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from scripts.analysis.mois_join_monthly_dynamic_climate_20261016 import (
    DynamicNormals,saturation_pressure_kpa,REAL_MANIFEST_SHA
)

SCHEMA="mois_real_worldclim48_cell_centered_monthly_climate_correlations_v1"
STATUS="REAL_MONTHLY_CLIMATE_ONLY_FIXED_CELL_SEASON_COLLINEARITY_NO_WHITE_PHOTOS"
DOMAINS=("TROPICAL","NORTH_EXTRATROPICAL","SOUTH_EXTRATROPICAL")
AXES=("tavg_c","srad_kj_m2_day","prec_mm","derived_vpd_kpa")
VALID_LIMITS={
    "tavg":(-95.,65.),"srad":(0.,50000.),
    "prec":(0.,10000.),"vapr":(0.,15.)
}
RADIUS_KM=6371.0088


def weighted_covariance(delta:np.ndarray,area:np.ndarray):
    """delta=[4 climate axes,12 months,N valid cells]; weights=area per cell."""
    if delta.ndim!=3 or delta.shape[0]!=4 or delta.shape[1]!=12 \
       or delta.shape[2]!=len(area) or not len(area) or np.any(area<=0):
        raise ValueError("Expected four climatological axes over twelve measured source months")
    total=12*float(np.sum(area))
    cov=np.einsum("atn,btn,n->ab",delta,delta,area,optimize=True)/total
    spread=np.sqrt(np.maximum(np.diag(cov),0))
    if np.any(spread<=1e-12):
        return {"status":"RANK_DEFICIENT_AFTER_FIXED_EFFECTS",
                "covariance_matrix_original_units":cov.tolist(),
                "correlation_matrix":None,
                "condition_number_correlation":None,
                "n_approx_independent_dynamic_axes":int(np.sum(spread>1e-12))}
    correlation=cov/spread[:,None]/spread[None,:]
    rank=int(np.linalg.matrix_rank(correlation,tol=1e-8))
    condition=float(np.linalg.cond(correlation))
    return {
        "status":("MOMENTS_COMPUTED_CORRELATIONS_NOT_HUE_RESPONSE" if
                  rank==4 else "RANK_DEFICIENT_DYNAMIC_CLIMATE_PROFILES"),
        "covariance_matrix_original_units":cov.tolist(),
        "correlation_matrix":correlation.tolist(),
        "condition_number_correlation":(
            condition if rank==4 and math.isfinite(condition) else None),
        "n_approx_independent_dynamic_axes":rank,
        "source_weighted_within_grid_monthly_SD":spread.tolist(),
    }


def audit(source:DynamicNormals, *, synthetic=False):
    if source.synthetic is not synthetic or not synthetic and source.manifest_sha!=REAL_MANIFEST_SHA:
        raise ValueError("Not the original byte-pinned 48-monthly-source climate archive")
    t=source.arrays
    if set(t)!=set(VALID_LIMITS):
        raise ValueError("Missing required actual T/S/P/vapour climate axes")
    h,w=source.height,source.width
    valid=np.ones((h,w),dtype=bool)
    for var,(low,high) in VALID_LIMITS.items():
        arr=t[var]
        if arr.shape!=(12,h,w):
            raise ValueError("Original monthly raster shape drift")
        valid &= np.all(np.isfinite(arr)&(arr>=low)&(arr<=high),axis=0)
    latnorth=source.transform.f+np.arange(h)*source.transform.e
    latsouth=latnorth+source.transform.e
    latmid=(latnorth+latsouth)/2
    area_by_row=(RADIUS_KM**2)*math.radians(source.transform.a)*np.abs(
        np.sin(np.deg2rad(latnorth))-np.sin(np.deg2rad(latsouth)))
    area=np.broadcast_to(area_by_row[:,None],valid.shape)
    groups={
        "TROPICAL":np.abs(latmid)<23.5,
        "NORTH_EXTRATROPICAL":latmid>=23.5,
        "SOUTH_EXTRATROPICAL":latmid<=-23.5,
    }
    result={}
    for domain in DOMAINS:
        mask=valid&groups[domain][:,None]
        n=int(np.count_nonzero(mask))
        if n<1:
            result[domain]={"status":"NO_COMPLETE_FOUR_VARIABLE_SOURCE_LAND_CELLS"}
            continue
        weight=np.asarray(area[mask],dtype=np.float64)
        # Avoid building a 48*global-grid float64 source copy. Gather ONLY
        # the current latitude area's complete 12-month land cells.
        temp=np.asarray(t["tavg"][:,mask],dtype=np.float64)
        solar=np.asarray(t["srad"][:,mask],dtype=np.float64)
        rain=np.asarray(t["prec"][:,mask],dtype=np.float64)
        vap=np.asarray(t["vapr"][:,mask],dtype=np.float64)
        raw_vpd=saturation_pressure_kpa(temp)-vap
        clipped=np.maximum(raw_vpd,0.)
        climate=np.stack((temp,solar,rain,clipped),axis=0)
        delta=climate-climate.mean(axis=1,keepdims=True)
        stageA=weighted_covariance(delta,weight)
        month_mean=np.einsum("atn,n->at",delta,weight,optimize=True)/float(np.sum(weight))
        stageB=weighted_covariance(delta-month_mean[:,:,None],weight)
        result[domain]={
            "n_complete_48layer_original_source_cells":n,
            "approx_complete_source_grid_area_million_km2":float(weight.sum()/1e6),
            "percent_valid_cells_with_any_negative_untruncated_sat_minus_actual_vapour":
                float(np.mean(np.any(raw_vpd<0,axis=0))*100),
            "model_A_after_exact_cell_fixed_effects":stageA,
            "model_B_after_cell_and_latitude_zone_calendar_month_fixed_effects":stageB,
            "source_climate_month_effects_weighted_annual_global_profile":month_mean.tolist(),
            "all_months_weighted_equally_not_photo_observer_calendar":True,
        }
    n_all=sum(result[z].get("n_complete_48layer_original_source_cells",0) for z in DOMAINS)
    if n_all!=int(np.count_nonzero(valid)):
        raise ValueError("Source 48-month raster valid geography partition changed")
    return {
        "schema":SCHEMA,"status":STATUS,
        "climatology":"WorldClim 2.1 1970–2000, 10 arc-minute",
        "source_original_48_month_raster_manifest_sha256":source.manifest_sha,
        "source_complete_all_48_monthly_climate_cells":n_all,
        "all_climate_axes":list(AXES),
        "original_climate_units":["°C","kJ/m2/day","mm/month","kPa, derived VPD"],
        "cell_month_climate_only":True,
        "original_fcp_sampled_flower_photographs_joined":False,
        "real_13_env_genus_spatial_colour_model_reestimated":False,
        "observed_flowering_plant_frequency_measured":False,
        "month_climate_normal_not_observed_weather":True,
        "synthetic_test_only":bool(synthetic),
        "latitude_domains":result,
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--verified-worldclim48-root",type=Path,required=True)
    p.add_argument("--verified-worldclim48-sha256",required=True)
    p.add_argument("--out-public-climate-collinearity",type=Path,required=True)
    a=p.parse_args()
    if a.out_public_climate_collinearity.exists():
        raise ValueError("Existing real climate-only audit cannot be overwritten")
    a.out_public_climate_collinearity.parent.mkdir(parents=True,exist_ok=True)
    climate=DynamicNormals(a.verified_worldclim48_root,a.verified_worldclim48_sha256)
    report=audit(climate)
    a.out_public_climate_collinearity.write_text(
        json.dumps(report,indent=2,sort_keys=True)+"\n",encoding="utf8")
    print("MOIS_ACTUAL_48_WORLDCLIM_ORIGINAL_CELL_CENTERED_MONTHLY_CLIMATE_ONLY",json.dumps(report,sort_keys=True))


if __name__=="__main__":main()
