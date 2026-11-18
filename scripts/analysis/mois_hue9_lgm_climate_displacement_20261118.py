#!/usr/bin/env python3
"""MOIS true-source LGM climate-displacement analysis; not a refugia map.

Use identical 12,295 source cells and effort counts as existing 9-hue climate FE.
Primary new test: does absolute joint LGM->present displacement explain
additional colour-entropy anomaly after the *signed* climate shifts and
vegetation are already included? Source images measure photo ROI hues only.
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
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler

from scripts.analysis.mois_hue9_assembly_environment_spatial_20261118 import (
    matrix, original_GBIF_record_effort, require, sha, spatial_fit,
    LGM_TEMP, LGM_RAIN, BASE, ECOLOGY,
)

ENV_SHA = "4f5d8fc3a308482fc3645d796a779d1a4725c51ed2afbf6a06e61fd6213f3ee6"
ENV_CSV = "all_current_GBIF_100km_cells_genus_breadth_entropy_ecology.csv.gz"
SEED = 20261118


def load_environment(path):
    require(sha(path) == ENV_SHA, "Historic environment output ZIP SHA changed")
    with ZipFile(path) as z:
        receipt = json.loads(z.read("receipt.json"))
        for filename, digest in receipt["numeric_outputs_sha256"].items():
            require(sha(z.read(filename)) == digest,
                    "Historic climate output member SHA changed: " + filename)
        frame = pd.read_csv(io.BytesIO(z.read(ENV_CSV)), compression="gzip")
    require(len(frame) == 12295 and not frame[["grid_x", "grid_y"]].duplicated().any(),
            "Original complete-case 100km cell count changed")
    require(receipt["source_cohort"]["strict_complete_30plus_cells"] == 12295,
            "Source climate receipt and grid rows disagree")
    require(frame.n_accepted_species.ge(30).all(), "Source richness filter changed")
    require(np.isclose(frame.entropy_anomaly.mean(),
                       receipt["source_cohort"]["global_mean_entropy_anomaly_in_complete_cells"],
                       atol=1e-12), "Source Shannon response changed")
    return frame, receipt


def assemble(env, gbif):
    effort = original_GBIF_record_effort(gbif)  # SHA, 5.49m pairs, 396m actual records checked
    out = env.merge(effort, on=["grid_x", "grid_y"], how="left", validate="one_to_one")
    require(len(out) == len(env) and out.GBIF_occurrence_count.notna().all()
            and out.GBIF_occurrence_count.gt(0).all(),
            "Original GBIF observation effort did not join one-to-one")
    return out.reset_index(drop=True)


def displacement(d):
    """Climate-space distance, NOT kilometres/year velocity or refugial survival.

    Zero delta denotes identical local modern/LGM climates. Divide each delta
    by its between-cell source standard deviation (do not mean-centre changes).
    """
    dt = d[LGM_TEMP].to_numpy(float)
    dp = d[LGM_RAIN].to_numpy(float)
    st, sp = np.std(dt, ddof=1), np.std(dp, ddof=1)
    require(np.isfinite(dt).all() and np.isfinite(dp).all() and st > 0 and sp > 0,
            "Cannot quantify LGM-to-present climate-space shift")
    score = np.hypot(dt / st, dp / sp)
    require(np.isfinite(score).all() and np.std(score, ddof=1) > 1e-6,
            "Climate-displacement index has no independent variance")
    return score, {"delta_temperature_SD_C":float(st),
                   "delta_log_precipitation_SD":float(sp),
                   "distance_units":"Euclidean distance of signed deltas / original global across-cell SD; zero means no change",
                   "not_true_climate_velocity":True,
                   "not_species_specific_palaeosuitability_or_refugia":True}


def fit_extension(d, predictor, spatial_km=500, include_interaction=False):
    X, scales = matrix(d)
    xp = (pd.Series(predictor, index=d.index, dtype=float) - np.mean(predictor)) / np.std(predictor, ddof=1)
    X["LGM_climate_displacement"] = xp
    if include_interaction:
        # Mean-centered SD-scale predictor interaction, conditional on main terms.
        X["displacement_x_NPP"] = xp * X["NPP"]
    group = d[f"block_{spatial_km}km"].astype(str)
    residual_X = X - X.groupby(group).transform("mean")
    residual_y = d.entropy_anomaly - d.entropy_anomaly.groupby(group).transform("mean")
    keep = residual_X.abs().sum(axis=1).gt(1e-11)
    xx, yy, bb = residual_X.loc[keep], residual_y.loc[keep], group.loc[keep]
    require(len(yy) >= 9000 and bb.nunique() >= 200, "Too few independent spatial observations")
    require(np.linalg.matrix_rank(xx.to_numpy()) == len(xx.columns), "Non-identifiable climate-stability term")
    fitted = sm.OLS(yy, xx).fit(cov_type="cluster", cov_kwds={"groups":bb}, use_t=False)
    result = {"n_cells":int(len(yy)), "n_blocks":int(bb.nunique()),
              "spatial_FE_km":spatial_km, "R2_within":float(fitted.rsquared),
              "coefficients":{}}
    for name in ["LGM_climate_displacement", "LGM_precip_delta", "LGM_temp_delta",
                 "NPP", "NDVI", "log_woody_AGB"] + (["displacement_x_NPP"] if include_interaction else []):
        ci=fitted.conf_int().loc[name]
        result["coefficients"][name] = {"beta_SD":float(fitted.params[name]),
            "CI95":[float(ci.iloc[0]),float(ci.iloc[1])],
            "p_cluster":float(fitted.pvalues[name])}
    return result


def spatial_oof_increment(d):
    """Train-only climate-distance scaling; 1000 km block leakage forbidden."""
    X,_=matrix(d)
    cols=list(BASE+ECOLOGY)
    y=d.entropy_anomaly.to_numpy(float)
    blocks=d.block_1000km.astype(str).to_numpy()
    t=d[LGM_TEMP].to_numpy(float)
    p=d[LGM_RAIN].to_numpy(float)
    a=np.full(len(d),np.nan)
    b=np.full(len(d),np.nan)
    c=np.full(len(d),np.nan)
    folds=[]
    for i,(train,test) in enumerate(GroupKFold(n_splits=5).split(X,y,blocks)):
        require(set(blocks[train]).isdisjoint(set(blocks[test])),
                "1000km spatial validation blocks overlap")
        norm=StandardScaler().fit(X.iloc[train][cols])
        xx_train=norm.transform(X.iloc[train][cols]);xx_test=norm.transform(X.iloc[test][cols])
        base=Ridge(alpha=100).fit(xx_train,y[train])
        a[test]=base.predict(xx_test)
        st=np.std(t[train],ddof=1);sp=np.std(p[train],ddof=1)
        require(st>0 and sp>0,"Train-only LGM climate scaling degenerate")
        d_train=np.hypot(t[train]/st,p[train]/sp)
        d_test=np.hypot(t[test]/st,p[test]/sp)
        extra_scaler=StandardScaler().fit(d_train.reshape(-1,1))
        v_train=np.column_stack([xx_train,extra_scaler.transform(d_train.reshape(-1,1))])
        v_test=np.column_stack([xx_test,extra_scaler.transform(d_test.reshape(-1,1))])
        aug=Ridge(alpha=100).fit(v_train,y[train])
        b[test]=aug.predict(v_test)
        # Product is formed using train-only standardized NPP and displacement;
        # standardize the product on training folds again before regression.
        npp_ix=cols.index("NPP")
        product_train=(v_train[:,-1]*v_train[:,npp_ix]).reshape(-1,1)
        product_test=(v_test[:,-1]*v_test[:,npp_ix]).reshape(-1,1)
        product_scaler=StandardScaler().fit(product_train)
        interaction_train=np.column_stack([v_train,product_scaler.transform(product_train)])
        interaction_test=np.column_stack([v_test,product_scaler.transform(product_test)])
        reg_interaction=Ridge(alpha=100).fit(interaction_train,y[train])
        c[test]=reg_interaction.predict(interaction_test)
        folds.append({"fold":i,"n_test":int(len(test)),
                      "n_1000km_blocks":int(len(set(blocks[test]))),
                      "baseline_MSE":float(mean_squared_error(y[test],a[test])),
                      "plus_distance_MSE":float(mean_squared_error(y[test],b[test])),
                      "plus_NPP_interaction_MSE":float(mean_squared_error(y[test],c[test]))})
    require(np.isfinite(a).all() and np.isfinite(b).all() and np.isfinite(c).all(),
            "Incomplete spatial heldout predictions")
    ids,inv=np.unique(blocks,return_inverse=True)
    counts=np.bincount(inv,minlength=len(ids))
    rng=np.random.default_rng(SEED+909)
    changes={}
    for label,alt,ref in [("distance_vs_baseline",b,a),
                           ("NPP_interaction_vs_baseline",c,a),
                           ("NPP_interaction_vs_distance",c,b)]:
        delta=(y-alt)**2-(y-ref)**2
        sums=np.bincount(inv,weights=delta,minlength=len(ids))
        boot=[]
        for _ in range(2000):
            draws=np.bincount(rng.integers(len(ids),size=len(ids)),minlength=len(ids))
            boot.append(float(np.dot(draws,sums)/np.dot(draws,counts)))
        changes[label]={"delta_MSE":float(np.mean(delta)),
                        "CI95_1000km_block_bootstrap":[float(x) for x in np.quantile(boot,[.025,.975])]}
    return {"n_cells":int(len(d)),"n_heldout_blocks":int(len(ids)),
        "n_folds":5,"baseline_full_environment_MSE":float(mean_squared_error(y,a)),
        "plus_distance_MSE":float(mean_squared_error(y,b)),
        "plus_distance_NPP_interaction_MSE":float(mean_squared_error(y,c)),
        "baseline_R2":float(r2_score(y,a)),"plus_distance_R2":float(r2_score(y,b)),
        "plus_distance_NPP_interaction_R2":float(r2_score(y,c)),
        "paired_MSE_delta":changes,
        "folds":folds,"training_only_distance_scaling":True,
        "exploratory_interaction_train_only_scaled":True}, a,b,c


def run(args):
    root=Path(args.out)
    require(not root.exists(),"Refuse overwriting historical or new source-supported results")
    env,receipt=load_environment(args.environment)
    d=assemble(env,args.GBIF)
    baseline=spatial_fit(d)  # must reproduce the prior fully adjusted 500km model
    prior=receipt["models"]["primary_500km_richness30plus"]
    require(baseline["n_fit_cells"]==prior["n_fit_cells"]==12239
            and baseline["n_spatial_clusters"]==prior["n_spatial_clusters"]==804,
            "Original spatial FE denominator not reproduced")
    for term,v in prior["coefficients"].items():
        require(abs(baseline["coefficients"][term]["beta_Shannon_per_global_SD"]-
                    v["beta_Shannon_per_global_SD"])<1e-9,
                "Original signed climate/vegetation FE coefficient changed: "+term)
    score,metric=displacement(d)
    d["LGM_climate_displacement"] = score
    models={"primary_500km":fit_extension(d,score,500),
            "richness100_500km":fit_extension(d.loc[d.n_accepted_species.ge(100)].reset_index(drop=True),
                                            score[d.n_accepted_species.ge(100)],500),
            "spatial_1000km":fit_extension(d,score,1000),
            "NPP_interaction_exploratory":fit_extension(d,score,500,True),
            "NPP_interaction_richness100_exploratory":fit_extension(
                d.loc[d.n_accepted_species.ge(100)].reset_index(drop=True),
                score[d.n_accepted_species.ge(100)],500,True),
            "NPP_interaction_1000km_exploratory":fit_extension(d,score,1000,True)}
    cv,pred_old,pred_new,pred_interaction=spatial_oof_increment(d)
    root.mkdir(parents=True)
    outcsv=root/"100km_photo9_LGM_climate_displacement_diagnostic.csv.gz"
    keep=["grid_x","grid_y","lat","lon","n_accepted_species","block_500km",
          "block_1000km","entropy_anomaly","white_species_fraction",LGM_TEMP,LGM_RAIN,
          "LGM_climate_displacement","NPP_2001_2023_mean_kgC_m2_yr","ESA_AGB_woody_Mg_ha_2020"]
    grid=d[keep].copy();grid["CV_pred_baseline_all_env"]=pred_old
    grid["CV_pred_add_LGM_displacement"]=pred_new
    grid["CV_pred_add_NPP_interaction"]=pred_interaction
    grid.to_csv(outcsv,index=False,compression="gzip")
    summary={"cohort":{"source_photo9_species":37104,"source_species_100km_pairs":5493931,
          "same_complete_case_cells":len(d),"LGM_threeGCM_median_only":True,
          "signed_previous_coefficients_exactly_reproduced":True},
          "source_ZIP_SHA256":{"environment":ENV_SHA,
             "GBIF":"422f6d6489e61c6e57c5d51db614e47f613039781a901e3fd3481c13b2b085c2"},
          "climate_metric":metric,"models":models,"spatial_holdout":cv,
          "inference":"Photo-derived species 9-hue composition at current GBIF 100km cells; no LGM range or glacial refugia observation",
          "biomass":"MOD17 NPP is productivity and ESA AGB is woody vegetation, not flower biomass"}
    js=root/"LGM_climate_stability_photo9_ecological_models.json"
    js.write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    summary["numeric_output_sha256"]={f.name:sha(f) for f in (outcsv,js)}
    (root/"receipt.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    print("MOIS_TRUE_SOURCE_LGM_STABILITY_TEST_SUCCESS",json.dumps({
        "n_cells":len(d),"coefficients":models["primary_500km"]["coefficients"],
        "sensitivity_1000km":models["spatial_1000km"]["coefficients"]["LGM_climate_displacement"],
        "NPP_interaction_sensitivities":{
            k:v["coefficients"]["displacement_x_NPP"]
            for k,v in models.items() if "interaction" in k},
        "CV":cv,"source_reproduction":True},allow_nan=False),flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--environment",required=True)
    parser.add_argument("--GBIF",required=True)
    parser.add_argument("--out",required=True)
    run(parser.parse_args())
