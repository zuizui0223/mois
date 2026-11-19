#!/usr/bin/env python3
"""MOIS: verify whether global photo9–palaeorainfall associations survive
conservative WCVP region-native filtering, and compare native/introduced
assemblages in exactly matched current GBIF 100km cells.

Signed original sources (all source image classes, WCVP and GBIF occurrences
unchanged); source photo ROI modal colour is not independently verified petal
colour. WCVP TDWG3 status is regional, not individual native provenance.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from sklearn.model_selection import GroupKFold
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.preprocessing import StandardScaler

from scripts.analysis.mois_hue9_assembly_environment_spatial_20261118 import (
    require, sha, spatial_fit, original_GBIF_record_effort,
)

GEO_SHA = "cbe1eef381e3e1a0062dbf5134b047f947338b7c499d2ecc0cd22ee8855909fd"
ENV_SHA = "4f5d8fc3a308482fc3645d796a779d1a4725c51ed2afbf6a06e61fd6213f3ee6"
GBIF_SHA = "422f6d6489e61c6e57c5d51db614e47f613039781a901e3fd3481c13b2b085c2"
GEO_CSV = "global_WCVP_TDWG3_native_only_photo9_island_alpine_polar_south_cells.csv.gz"
ENV_CSV = "all_current_GBIF_100km_cells_genus_breadth_entropy_ecology.csv.gz"
NPP = "NPP_2001_2023_mean_kgC_m2_yr"
LGM_TEMP = "threeGCM_delta_temperature_C_median"
LGM_PRECIP = "threeGCM_delta_log_precip_median"
ECOLOGY = ("temperature","log_precip","LGM_temp_delta","LGM_precip_delta",
           "NPP","NDVI","log_woody_AGB")
NO_LGM = ("temperature","log_precip","NPP","NDVI","log_woody_AGB")
GEO_BASE = ("log_status_richness","expected_status_null_Shannon",
            "log_all_GBIF_records","abs_latitude","longitude_sine","longitude_cosine")
PAIRED_BASE = ("log_native_richness","log_intro_richness",
               "difference_expected_null_Shannon","log_all_GBIF_records",
               "abs_latitude","longitude_sine","longitude_cosine")
SEED = 20261119

def source_csv(archive,expected,member,digest_key):
    require(sha(archive)==expected,"Source ZIP digest changed: "+member)
    with ZipFile(archive) as z:
        receipt=json.loads(z.read("receipt.json"))
        for filename,digest in receipt[digest_key].items():
            require(sha(z.read(filename))==digest,
                    "Signed prior result member digest drift: "+filename)
        d=pd.read_csv(io.BytesIO(z.read(member)),compression="gzip")
    require(not d[["grid_x","grid_y"]].duplicated().any(),
            "Prior WCVP/GBIF environment output has repeated geographic cell")
    return d,receipt

def load_sources(args):
    geo,rg=source_csv(args.geography,GEO_SHA,GEO_CSV,"output_sha256")
    env,re=source_csv(args.environment,ENV_SHA,ENV_CSV,"numeric_outputs_sha256")
    require(len(geo)==18577 and len(env)==12295,
            "Global WCVP grid or environment cohort mismatch")
    require(rg["source_37104_accepted_photo_hue_taxa"]==37104
            and rg["source_full_GBIF_species_cell_pairs"]==5493931
            and rg["source_null_199_same_global_genus_by_width_colour_label_permutations"]
            and rg["source_WCVP_exact_match_and_TDWG3_status"][
               "n_WCVP_accepted_unique_exact_matched_species"]==36762,
            "Previous WCVP source not the validated unique-exact region-native cohort")
    require(rg["results"]["strict_native"]["geographic_groups"][
               "all_WCVP_status_species_richness_ge30"]["n_real_100km_cells"]==7231,
            "Previous WCVP strict-native count drift")
    effort=original_GBIF_record_effort(args.GBIF)
    original=env.merge(effort,on=["grid_x","grid_y"],validate="one_to_one")
    require(len(original)==12295 and original.GBIF_occurrence_count.gt(0).all(),
            "Original 396M full GBIF observation effort not recovered")
    # The identical full-flora environmental model must be exactly reproducible.
    previous=re["models"]["primary_500km_richness30plus"]
    reproduced=spatial_fit(original)
    require(reproduced["n_fit_cells"]==previous["n_fit_cells"]==12239
            and reproduced["n_spatial_clusters"]==previous["n_spatial_clusters"]==804,
            "Historic full-flora FE denominator drift")
    for key,val in previous["coefficients"].items():
        require(abs(reproduced["coefficients"][key]["beta_Shannon_per_global_SD"]-
                    val["beta_Shannon_per_global_SD"])<1e-9,
                "Historic full-flora coefficient not reproduced: "+key)
    d=original.merge(geo,on=["grid_x","grid_y"],validate="one_to_one")
    require(len(d)==12295 and d[["grid_x","grid_y"]].is_unique,
            "WCVP native table and original climate cells did not join exactly")
    require(d.strict_native_n_photo9_species.ge(30).sum()==7201
            and d.moderate_native_n_photo9_species.ge(30).sum()==8119,
            "Source strict/moderate native30 climate-complete sample drift")
    require(((d.strict_native_n_photo9_species>=30)&
             (d.strict_introduced_n_photo9_species>=30)).sum()==3106,
            "Source 3106 strict paired WCVP cells drift")
    return d,re,rg

def standardized(raw):
    require(np.isfinite(raw.to_numpy(float)).all(),"Incomplete source predictors")
    mean=raw.mean(axis=0);sd=raw.std(axis=0,ddof=1)
    require(sd.gt(1e-10).all(),"Unidentified covariate within selected WCVP cohort")
    return (raw-mean)/sd,{
        k:{"mean":float(mean[k]),"SD":float(sd[k])} for k in raw.columns
    }

def climate_columns(d):
    return {
      "temperature":d.present_temperature_C.astype(float),
      "log_precip":np.log1p(d.present_precip_mm.astype(float)),
      "LGM_temp_delta":d[LGM_TEMP].astype(float),
      "LGM_precip_delta":d[LGM_PRECIP].astype(float),
      "NPP":d[NPP].astype(float),
      "NDVI":d.MODIS_NDVI_2001_2021_mean.astype(float),
      "log_woody_AGB":np.log1p(d.ESA_AGB_woody_Mg_ha_2020.astype(float)),
      "log_all_GBIF_records":np.log1p(d.GBIF_occurrence_count.astype(float)),
      "abs_latitude":d.lat.abs().astype(float),
      "longitude_sine":np.sin(np.deg2rad(d.lon.astype(float))),
      "longitude_cosine":np.cos(np.deg2rad(d.lon.astype(float))),
    }

def make_design(d,status=None,paired_mode=None):
    raw=climate_columns(d)
    if status is not None:
        raw["log_status_richness"]=np.log1p(d[f"{status}_n_photo9_species"].astype(float))
        raw["expected_status_null_Shannon"]=d[f"{status}_Shannon_null"].astype(float)
        cols=list(ECOLOGY)+list(GEO_BASE)
        y=d[f"{status}_Shannon_anomaly"].astype(float)
    else:
        require(paired_mode in ("strict","moderate"),"Unknown WCVP paired mode")
        raw["log_native_richness"]=np.log1p(d[f"{paired_mode}_native_n_photo9_species"].astype(float))
        raw["log_intro_richness"]=np.log1p(d[f"{paired_mode}_introduced_n_photo9_species"].astype(float))
        raw["difference_expected_null_Shannon"]=(
           d[f"{paired_mode}_native_Shannon_null"]-
           d[f"{paired_mode}_introduced_Shannon_null"]).astype(float)
        cols=list(ECOLOGY)+list(PAIRED_BASE)
        y=(d[f"{paired_mode}_native_Shannon_anomaly"]-
           d[f"{paired_mode}_introduced_Shannon_anomaly"]).astype(float)
    raw=pd.DataFrame(raw,index=d.index)
    raw=raw[cols]
    x,scale=standardized(raw)
    require(np.isfinite(y).all(),"Source response missing after WCVP species minimum")
    return x,y,scale

def fit_spatial(d,block_km,status=None,paired_mode=None):
    x,y,scale=make_design(d,status,paired_mode)
    b=d[f"block_{block_km}km"].astype(str)
    xc=x-x.groupby(b).transform("mean")
    yc=y-y.groupby(b).transform("mean")
    keep=xc.abs().sum(axis=1).gt(1e-10)
    xx=xc.loc[keep]; yy=yc.loc[keep]; groups=b.loc[keep]
    require(len(xx)>650 and groups.nunique()>90,
            "Insufficient WCVP native/introduced spatially informative groups")
    require(np.linalg.matrix_rank(xx.to_numpy())==len(xx.columns),
            "Within-block WCVP climate model singular")
    fit=sm.OLS(yy,xx).fit(cov_type="cluster",
                          cov_kwds={"groups":groups},use_t=False)
    p=[float(fit.pvalues[k]) for k in ECOLOGY]
    q=multipletests(p,method="fdr_bh")[1]
    coefs={}
    for k in x.columns:
        ci=fit.conf_int().loc[k]
        coefs[k]={"beta_Shannon_per_cohort_SD":float(fit.params[k]),
                  "CI95":[float(ci.iloc[0]),float(ci.iloc[1])],
                  "p_cluster":float(fit.pvalues[k])}
    for k,v in zip(ECOLOGY,q):
        coefs[k]["BH_q7_ecological_exposures"]=float(v)
    return {"status":status,"paired_mode":paired_mode,
        "n_original_cells":int(len(d)),"n_fit_cells":int(len(yy)),
        "n_spatial_blocks":int(groups.nunique()),"block_km":block_km,
        "R2_within_block_in_sample":float(fit.rsquared),
        "term_scaling":scale,"coefficients":coefs,
        "full_GBIF_record_count_adjusts_global_effort_not_native_only_effort":True,
        "WCVP_TDWG3_native_not_proven_true_wild_local_individual":True}

def spatial_holdout(d,status):
    x,y,_=make_design(d,status=status)
    no_lgm=list(GEO_BASE)+list(NO_LGM)
    with_lgm=list(GEO_BASE)+list(ECOLOGY)
    block=d.block_1000km.astype(str).to_numpy()
    y=y.to_numpy(float)
    before=np.full(len(d),np.nan); after=np.full(len(d),np.nan)
    fold=[]
    for i,(tr,te) in enumerate(GroupKFold(n_splits=5).split(x,y,block)):
        require(set(block[tr]).isdisjoint(set(block[te])),
                "WCVP native 1000km validation has block overlap")
        for cols,pred in ((no_lgm,before),(with_lgm,after)):
            scaler=StandardScaler().fit(x.iloc[tr][cols])
            fitted=Ridge(alpha=100,solver="cholesky").fit(
                scaler.transform(x.iloc[tr][cols]),y[tr])
            pred[te]=fitted.predict(scaler.transform(x.iloc[te][cols]))
        fold.append({"fold":i,"n_test_cells":int(len(te)),
             "n_1000km_blocks":int(len(set(block[te]))),
             "no_LGM_MSE":float(mean_squared_error(y[te],before[te])),
             "with_LGM_MSE":float(mean_squared_error(y[te],after[te]))})
    require(np.isfinite(before).all() and np.isfinite(after).all(),
            "Missing native out-of-region predicted values")
    diff=(y-after)**2-(y-before)**2
    unique,idx=np.unique(block,return_inverse=True)
    sums=np.bincount(idx,weights=diff,minlength=len(unique))
    counts=np.bincount(idx,minlength=len(unique))
    rng=np.random.default_rng(SEED)
    boot=[]
    for _ in range(2000):
        draw=np.bincount(rng.integers(len(unique),size=len(unique)),
                         minlength=len(unique))
        boot.append(float(np.dot(draw,sums)/np.dot(draw,counts)))
    return {"n_cells":int(len(y)),"n_1000km_block_groups":len(unique),
       "fivefold_1000km_blocks_disjoint":True,
       "historical_climate_added_to_modern_weather_and_productivity":True,
       "with_LGM":{"MSE":float(mean_squared_error(y,after)),
                   "R2_OOF":float(r2_score(y,after))},
       "no_LGM":{"MSE":float(mean_squared_error(y,before)),
                 "R2_OOF":float(r2_score(y,before))},
       "delta_MSE_with_minus_without_LGM":float(diff.mean()),
       "delta_MSE_95CI_1000km_block_bootstrap":[float(v) for v in
                                                   np.percentile(boot,[2.5,97.5])],
       "folds":fold,"Ridge_solver":"cholesky","LGM_3GCM_median_is_not_species_palaeorange":True}

def run(args):
    root=Path(args.out)
    require(not root.exists(),"Do not overwrite signed WCVP climate replication")
    d,old,rg=load_sources(args)
    models={}
    specs=[("strict_native",30,500),("strict_native",100,500),
           ("strict_native",30,1000),("moderate_native",30,500),
           ("moderate_native",30,1000),("strict_introduced",30,500),
           ("strict_introduced",30,1000)]
    for status,n,km in specs:
        subset=d.loc[d[f"{status}_n_photo9_species"].ge(n)].reset_index(drop=True)
        models[f"{status}_species{n}_{km}km"]=fit_spatial(subset,km,status=status)
    for mode,threshold,km in [("strict",30,500),("strict",30,1000),
                              ("moderate",30,500),("strict",100,500)]:
        mask=(d[f"{mode}_native_n_photo9_species"].ge(threshold)&
              d[f"{mode}_introduced_n_photo9_species"].ge(threshold))
        subset=d.loc[mask].reset_index(drop=True)
        models[f"{mode}_native_minus_introduced_paired{threshold}_{km}km"]=fit_spatial(
             subset,km,paired_mode=mode)
    native=d.loc[d.strict_native_n_photo9_species.ge(30)].reset_index(drop=True)
    cv=spatial_holdout(native,"strict_native")
    repeat=spatial_holdout(native,"strict_native")
    require(cv["with_LGM"]["MSE"]==repeat["with_LGM"]["MSE"]
            and cv["no_LGM"]["MSE"]==repeat["no_LGM"]["MSE"],
            "Ridge solver not exactly deterministic")
    root.mkdir(parents=True)
    rows=[]
    for name,model in models.items():
        for term,val in model["coefficients"].items():
            rows.append({"model":name,"n_cells":model["n_fit_cells"],
               "n_blocks":model["n_spatial_blocks"],"term":term,
               "beta_per_cohort_SD":val["beta_Shannon_per_cohort_SD"],
               "CI_low":val["CI95"][0],"CI_high":val["CI95"][1],
               "p_cluster":val["p_cluster"],
               "q7_BH":val.get("BH_q7_ecological_exposures")})
    csv=root/"wcvp_photo9_native_vs_intro_LGM_NPP_spatial_effects.csv"
    pd.DataFrame(rows).to_csv(csv,index=False)
    meta={"source": {
        "original_photo9_COLXR_accepted_species":37104,
        "original_species_grid_pairs":5493931,
        "original_GBIF_records":396142825,
        "original_env_cells":12295,
        "original_WCVP_native_199_null_permutations_reused":True,
        "strict_native30_complete_env_cells":int(len(native)),
        "strict_WCVP_TDWG3_native_not_same_as_wild_specimen":True,
        "original_full_flora_climate_FE_coefficients_reproduced":True,
        "original_source_ZIP_SHA256":{"geography":GEO_SHA,
                                      "environment":ENV_SHA,"GBIF":GBIF_SHA},
    },
    "models":models,"strict_native_1000km_block_holdout":cv,
    "model_response":"original WCVP status-specific observed minus genus x occupancy-breadth 199-permutation expected 9-hue Shannon",
    "native_vs_introduced_paired_samecells":True,
    "not_petal_colour_or_LGM_refugia_or_pigment_fitness_causality":True,
    "no_new_original_photograph_collection":True}
    js=root/"wcvp_native_intro_paleorain_photo9_regression_and_holdout.json"
    js.write_text(json.dumps(meta,indent=2,allow_nan=False)+"\n")
    meta["numeric_outputs_sha256"]={str(f.name):sha(f) for f in (csv,js)}
    (root/"receipt.json").write_text(json.dumps(meta,indent=2,allow_nan=False)+"\n")
    print("MOIS_AUTHENTIC_NATIVE_LGM_COMPARE_SUCCESS",json.dumps({
       "n_native":len(native),
       "strict_native_500km":models["strict_native_species30_500km"]["coefficients"]["LGM_precip_delta"],
       "strict_native_1000km":models["strict_native_species30_1000km"]["coefficients"]["LGM_precip_delta"],
       "moderate_native_500km":models["moderate_native_species30_500km"]["coefficients"]["LGM_precip_delta"],
       "strict_paired500":models["strict_native_minus_introduced_paired30_500km"]["coefficients"]["LGM_precip_delta"],
       "strict_paired1000":models["strict_native_minus_introduced_paired30_1000km"]["coefficients"]["LGM_precip_delta"],
       "heldout":cv},allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--geography",required=True)
    p.add_argument("--environment",required=True)
    p.add_argument("--GBIF",required=True)
    p.add_argument("--out",required=True)
    run(p.parse_args())
