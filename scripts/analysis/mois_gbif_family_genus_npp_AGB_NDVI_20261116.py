#!/usr/bin/env python3
"""Independent NPP increment for GBIF/WCVP taxonomic assemblage segregation.

Outcome: observed genus Jaccard minus exact native/introduction label shuffle
expectation conditional on WCVP family and local species/genus counts.
Response does NOT use floral photos or hue proxies.
Predictors: 2001-2023 authentic MOD17A3HGF.061-derived annual NPP (0.5deg
open ICDC product DOI 10.25592/uhhfdm.14635), woody ESA AGB 2020, MODIS
NDVI 2001-21, WorldClim present, LGM temperature/precip change, latitude,
species counts, recorded GBIF effort. Same actual 100km GBIF spatial cells.
NPP is kg C m^-2 yr^-1 flux; AGB is dry woody mass stock; NDVI greenness.
No causality, local native individual identity, or LGM refugia inferred.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

from scripts.analysis.mois_gbif_family_excess_biomass_ndvi_20261116 import (
    signed_parquet,build_comparison,normalized_model_matrix,
    FAMILY_SHA,SPATIAL_SHA,NATIVE_SHA,NDVI_SHA,
    FAMILY_MEMBER,SPATIAL_MEMBER,NATIVE_MEMBER,NDVI_MEMBER
)
from scripts.analysis.mois_lgm_climate_history_species_20261116 import (
    require,sha
)

NPP_SHA="68b2cfa1aeac9dc5d02cd310c9b2f45553a5be4db3319a1e40897612ad9b1bef"
NPP_MEMBER="MOD17_C061_NPP_2001_2023_18577_true_GBIF100km_cells.parquet"
NPP_FIELD="NPP_2001_2023_mean_kgC_m2_yr"
MIN_YEARS=18

def read_sources(family,spatial,native,ndvi,npp):
    f=signed_parquet(family,FAMILY_SHA,FAMILY_MEMBER,"output_SHA256")
    s=signed_parquet(spatial,SPATIAL_SHA,SPATIAL_MEMBER,"output_sha256")
    n=signed_parquet(native,NATIVE_SHA,NATIVE_MEMBER,"output_sha256")
    v=signed_parquet(ndvi,NDVI_SHA,NDVI_MEMBER,"output_SHA256")
    p=signed_parquet(npp,NPP_SHA,NPP_MEMBER,"output_SHA256")
    require(len(p)==18577 and not p[["grid_x","grid_y"]].duplicated().any(),
            "true NPP series global GBIF cell coverage altered")
    cols=["grid_x","grid_y",NPP_FIELD,"NPP_valid_years_2001_2023",
          "NPP_2001_2023_sd_kgC_m2_yr",
          "NPP_2001_2023_cv","NPP_2020_kgC_m2_yr",
          "NPP_2001_kgC_m2_yr","NPP_2023_kgC_m2_yr"]
    require(set(cols).issubset(p.columns),"signed 23-annual NPP fields missing")
    base=build_comparison(f,s,n,v)
    d=base.merge(p[cols],on=["grid_x","grid_y"],
                 how="left",validate="many_to_one")
    require(len(d)==len(base),"NPP join duplicated WCVP GBIF cells")
    d["valid_NPP_18of23"]=(d.NPP_valid_years_2001_2023>=MIN_YEARS)
    d["valid_NPP_all23"]=d.NPP_valid_years_2001_2023.eq(23)
    d["complete_3vegetation"]=(d.complete_woodyNDVI
              &d.valid_NPP_18of23
              &np.isfinite(d[NPP_FIELD])
              &d.NPP_2001_2023_sd_kgC_m2_yr.notna())
    return d

def matrix_with_NPP(d,yearly_npp=False):
    x,sd=normalized_model_matrix(d)
    col="NPP_2020_kgC_m2_yr" if yearly_npp else NPP_FIELD
    values=d[col].astype(float)
    deviation=values.std()
    require(deviation>1e-9 and np.isfinite(deviation),
            "annual NPP species-grid cells have no usable variance")
    x["NPP_annual_carbon_flux"]=(values-values.mean())/deviation
    return x, {**sd,"NPP_annual_carbon_flux":float(deviation)}

def fit(d,response="residual_after_family_conditioning",
        block="500",npp=True,agb=True,ndvi=True,support_weighted=False,
        year2020_only=False):
    require(len(d)>300 and d.block_500km.nunique()>40,
            "spatial replicates insufficient")
    x,rawsd=matrix_with_NPP(d,yearly_npp=year2020_only)
    if not npp:x=x.drop(columns=["NPP_annual_carbon_flux"])
    if not agb:x=x.drop(columns=["log_woody_AGB"])
    if not ndvi:x=x.drop(columns=["NDVI_mean"])
    y=d[response].astype(float)
    b=d.block_500km.astype(str) if block=="500" else d.block_1000km.astype(str)
    x=x-x.groupby(b).transform("mean")
    y=y-y.groupby(b).transform("mean")
    useful=x.abs().sum(axis=1)>0
    x=x.loc[useful]
    y=y.loc[useful]
    groups=d.loc[useful,"block_500km"].astype(str)
    require(groups.nunique()>30 and len(y)>300,
            "too little within-block independent exposure support")
    if support_weighted:
        w=d.loc[useful,["native_species_n","introduced_species_n"]].min(axis=1)
        model=sm.WLS(y,x,weights=w/w.median()).fit(
            cov_type="cluster",cov_kwds={"groups":groups},use_t=False)
    else:
        model=sm.OLS(y,x).fit(
            cov_type="cluster",cov_kwds={"groups":groups},use_t=False)
    effects={}
    for name in ("NPP_annual_carbon_flux","log_woody_AGB","NDVI_mean",
                 "LGM_warming"):
        if name in model.params:
            low,high=model.conf_int().loc[name].tolist()
            effects[name]={"beta_per_global_SD":float(model.params[name]),
                "95CI":[float(low),float(high)],
                "p":float(model.pvalues[name])}
    return {"n_cells":len(y),"n_independent_500km_blocks":int(groups.nunique()),
            "R2_in_sample":float(model.rsquared),
            "response":response,"fixed_block_km":int(block),
            "npp_years":"2020" if year2020_only else "2001-2023",
            "AGB_included":bool(agb),"NDVI_included":bool(ndvi),
            "NPP_included":bool(npp),
            "weights":"min_WCVP_status_group_species_count" if support_weighted else "equal_cells",
            "effects":effects,"raw_global_SD":{
                k:float(rawsd[k]) for k in (
                    "NPP_annual_carbon_flux","log_woody_AGB",
                    "NDVI_mean","LGM_warming")},
            "NPP_is_annual_carbon_flux_not_biomass":True,
            "R2_is_in_sample_not_external_predictive_validity":True}

def evaluate(d,mode):
    z=d[d.region_mode.eq(mode)].copy()
    require(len(z)==({"strict":3106,"moderate":3602}[mode]),
            "same WCVP native-introduced GBIF region support changed")
    usable=z[z.complete_3vegetation].copy()
    require(len(usable)>=2000,
            "not enough matched 23year-NPP, AGB and NDVI support")
    results={
       "n_total_original_WCVP_cells":len(z),
       "n_cells_18plus_NPP_same_GBIF_full_environment":len(usable),
       "n_500km_spatial_clusters":int(usable.block_500km.nunique()),
       "NPP_annual_median_kgC_m2_year":float(usable[NPP_FIELD].median()),
       "NPP_annual_mean_kgC_m2_year":float(usable[NPP_FIELD].mean()),
       "AGB2020_median_Mg_ha":float(usable.ESA_AGB_woody_Mg_ha_2020.median()),
       "NDVI_median":float(usable.MODIS_NDVI_2001_2021_mean.median()),
       "corr_spearman":{
         "NPP_AGB":float(usable[NPP_FIELD].corr(
             usable.ESA_AGB_woody_Mg_ha_2020,method="spearman")),
         "NPP_NDVI":float(usable[NPP_FIELD].corr(
             usable.MODIS_NDVI_2001_2021_mean,method="spearman")),
         "NDVI_AGB":float(usable.MODIS_NDVI_2001_2021_mean.corr(
             usable.ESA_AGB_woody_Mg_ha_2020,method="spearman"))
       },
       "models":{}
    }
    for outcome in ("residual_after_family_conditioning","observed_jaccard"):
        fit_same={
            "500km_full_NPP_AGB_NDVI":fit(usable,outcome),
            "500km_no_NPP":fit(usable,outcome,npp=False),
            "500km_no_AGB":fit(usable,outcome,agb=False),
            "500km_no_NDVI":fit(usable,outcome,ndvi=False),
            "500km_NPP_only_vegetation":fit(usable,outcome,agb=False,ndvi=False),
            "1000km_full_NPP_AGB_NDVI":fit(usable,outcome,block="1000"),
            "500km_full_support_weighted":fit(usable,outcome,support_weighted=True)
        }
        if usable.valid_NPP_all23.sum()>=1500:
            y23=usable[usable.valid_NPP_all23].copy()
            fit_same["500km_full_all23years"]=fit(y23,outcome)
        if usable.NPP_2020_kgC_m2_yr.notna().sum()>=1500:
            y2020=usable[usable.NPP_2020_kgC_m2_yr.notna()].copy()
            fit_same["500km_2020_NPP_single_year"]=fit(y2020,outcome,year2020_only=True)
        fit_same["R2_increment_from_NPP_on_identical_cells"]=(
            fit_same["500km_full_NPP_AGB_NDVI"]["R2_in_sample"]
            -fit_same["500km_no_NPP"]["R2_in_sample"])
        results["models"][outcome]=fit_same
    return results,usable

def run(family,spatial,native,ndvi,npp,output):
    out=Path(output)
    require(not out.exists(),"refuse overwrite signed original-NPP regression")
    d=read_sources(family,spatial,native,ndvi,npp)
    results={};tables=[]
    for mode in ("strict","moderate"):
        results[mode],usable=evaluate(d,mode)
        t=usable[["grid_x","grid_y","region_mode",
             "residual_after_family_conditioning",
             "observed_jaccard",
             "ESA_AGB_woody_Mg_ha_2020",
             "MODIS_NDVI_2001_2021_mean",
             "NPP_valid_years_2001_2023",
             NPP_FIELD,
             "NPP_2001_2023_sd_kgC_m2_yr",
             "threeGCM_delta_temperature_C_median",
             "block_500km"]].copy()
        tables.append(t)
        fit0=results[mode]["models"]["residual_after_family_conditioning"]
        effects=fit0["500km_full_NPP_AGB_NDVI"]["effects"]
        print("MOIS_GBIF_23YR_REAL_NPP_VS_FAMILY_FIXED_GENUS_NULL",
              json.dumps({"mode":mode,"cells":len(usable),
                "NPP_beta":effects["NPP_annual_carbon_flux"]["beta_per_global_SD"],
                "NPP_p":effects["NPP_annual_carbon_flux"]["p"],
                "NPP_CI":effects["NPP_annual_carbon_flux"]["95CI"],
                "LGM_p":effects["LGM_warming"]["p"],
                "AGB_p":effects["log_woody_AGB"]["p"],
                "NDVI_p":effects["NDVI_mean"]["p"],
                "R2_increment_in_sample":fit0["R2_increment_from_NPP_on_identical_cells"],
                "correlations":results[mode]["corr_spearman"]},
                sort_keys=True),flush=True)
    out.mkdir(parents=True)
    table=out/"GBIF_WCVP_family_null_NPP_2001_2023_AGB_NDVI_same_cells.parquet"
    pd.concat(tables,ignore_index=True).to_parquet(table,index=False)
    models=out/"models.json"
    models.write_text(json.dumps(results,indent=2,sort_keys=True)+"\n")
    audit={"schema":"mois_GBIF_WCVP_family_null_23year_NPP_AGB_NDVI_vs_LGM_v1",
          "source_GBIF_accepted_species":37104,
          "source_GBIF_species_x_100km_pairs":5493931,
          "native_is_broad_WCVP_TDWG3_not_wild_individual":True,
          "LGM_is_paleoclimate_at_current_GBIF_occupancy_not_refugia":True,
          "response_uses_no_photo_hue":True,
          "NPP_is_annual_carbon_flux_not_woody_biomass":True,
          "NPP_source_0point5deg_is_not_original_500m":True,
          "NPP_is_not_flower_allocation":True,
          "sources_zip_sha256":{"FAMILY":FAMILY_SHA,"SPATIAL":SPATIAL_SHA,
                 "NATIVE":NATIVE_SHA,"NDVI":NDVI_SHA,"NPP":NPP_SHA},
          "results":results,"output_SHA256":{
                 table.name:sha(table),models.name:sha(models)}}
    (out/"receipt.json").write_text(json.dumps(audit,indent=2,sort_keys=True)+"\n")
    return audit

if __name__=="__main__":
    p=argparse.ArgumentParser()
    for item in ("family","spatial","native","ndvi","npp","output"):
        p.add_argument("--"+item,required=True)
    a=p.parse_args()
    run(a.family,a.spatial,a.native,a.ndvi,a.npp,a.output)
