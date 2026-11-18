#!/usr/bin/env python3
"""Follow-up ecological stress test of the exploratory red-photo-hue signal.

Frozen source and cohort from the original WCVP introduction-filtering analysis.
Tests red vs 8 hues, red vs white within same genera, source class imbalance,
alternative introduced threshold >=2 WCVP TDWG3 regions, native-range support,
large-genus influence and genus-disjoint out-of-fold predictive increments.
No re-selection for positive significance; all sensitivities reported.

This is NOT evidence introduced individuals evolved red flowers, nor that the
photographic ROI-red label has been botanically verified as petals.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, log_loss
from scripts.analysis.mois_wcvp_introduction_colour_filter_20261118 import (
    signed_env,signed_wcvp,prepare,design,within_genus,require,sha,SEED
)

N_BOOT=2000

def red_effect(d,label,genus_weight=False,status=True):
    require(d.source_hue9.eq("red").any(),"Red photographic hue missing")
    dd=d.copy()
    dd["is_red"]=dd.source_hue9.eq("red").astype(int)
    try:
        r=within_genus(dd,white_col="is_red",
                       weighting="genus" if genus_weight else "species",
                       photo_status=status,cohort=label)
        return r
    except ValueError as exc:
        return {"cohort":label,"status":"not_estimable",
                "reason":str(exc)[:150],"source_taxa":len(d),
                "source_red":int(dd.is_red.sum()),
                "source_introduced":int(dd.recorded_introduced.sum())}

def red_vs_white_direct(d):
    z=d[d.source_hue9.isin(("red","white"))].copy().reset_index(drop=True)
    z["is_red"]=z.source_hue9.eq("red").astype(int)
    g=z.genus.astype(str)
    # Require BOTH photo hues AND BOTH introduced states within genus:
    # Only such genera inform a conditional cross-colour comparison.
    mixed=(z.groupby("genus").is_red.transform("nunique").eq(2)
           &z.groupby("genus").recorded_introduced.transform("nunique").eq(2))
    z=z.loc[mixed].copy().reset_index(drop=True)
    require(len(z)>400 and z.genus.nunique()>30,
            "Too few genera for strict direct red-white comparison")
    x=design(z)
    x["colour_red_vs_white"]=z.is_red.astype(float)
    group=z.genus.astype(str)
    xc=x-x.groupby(group).transform("mean")
    yc=z.recorded_introduced.astype(float)
    yc=yc-yc.groupby(group).transform("mean")
    xc=xc.loc[:,xc.std()>1e-10]
    m=sm.WLS(yc,xc,weights=np.ones(len(z))).fit(
        cov_type="cluster",cov_kwds={"groups":group},use_t=False)
    ci=m.conf_int().loc["colour_red_vs_white"]
    return {"n_species":len(z),"n_genera":int(z.genus.nunique()),
            "red_n":int(z.is_red.sum()),
            "native_region_adjusted":True,
            "red_minus_white_recorded_introduction_probability_pp":
                float(100*m.params["colour_red_vs_white"]),
            "CI95_pp":[float(100*ci.iloc[0]),float(100*ci.iloc[1])],
            "p_genus_cluster":float(m.pvalues["colour_red_vs_white"])}

def cv_red_white(d):
    y=d.recorded_introduced.to_numpy(int)
    genera=d.genus.to_numpy(str)
    x=design(d).to_numpy(float)
    w=d.source_hue9.eq("white").to_numpy(float)
    r=d.source_hue9.eq("red").to_numpy(float)
    models={
      "ecology":x,
      "ecology_plus_white":np.column_stack((x,w)),
      "ecology_plus_red":np.column_stack((x,r)),
      "ecology_plus_white_and_red":np.column_stack((x,w,r))
    }
    folds=list(StratifiedGroupKFold(n_splits=3,shuffle=True,
          random_state=SEED).split(x,y,genera))
    preds={k:np.full(len(y),np.nan) for k in models}
    for i,(train,test) in enumerate(folds):
        require(set(genera[train]).isdisjoint(set(genera[test])),
                "Genus not entirely held out")
        for name,X in models.items():
            im=SimpleImputer(strategy="median").fit(X[train])
            xt=im.transform(X[train]);xv=im.transform(X[test])
            s=StandardScaler().fit(xt)
            m=LogisticRegression(C=.1,max_iter=750,solver="lbfgs")
            m.fit(s.transform(xt),y[train])
            preds[name][test]=m.predict_proba(s.transform(xv))[:,1]
    require(all(np.isfinite(a).all() for a in preds.values()),
            "Holdout probabilities are not finite")
    metrics={name:{"AUC":float(roc_auc_score(y,pr)),
                   "logloss":float(log_loss(y,pr))}
             for name,pr in preds.items()}
    baseline=preds["ecology"]
    gnames,g=np.unique(genera,return_inverse=True)
    rng=np.random.default_rng(SEED+77)
    increments=[]
    for name,p in preds.items():
        if name=="ecology":continue
        # Exact paired out-of-fold individual log-loss increments.
        p=np.clip(p,1e-10,1-1e-10)
        b=np.clip(baseline,1e-10,1-1e-10)
        v=-(y*np.log(p)+(1-y)*np.log1p(-p))+(y*np.log(b)+(1-y)*np.log1p(-b))
        sums=np.bincount(g,weights=v,minlength=len(gnames))
        n=np.bincount(g,minlength=len(gnames))
        bs=[]
        for _ in range(N_BOOT):
            c=np.bincount(rng.integers(0,len(gnames),size=len(gnames)),
                          minlength=len(gnames))
            bs.append(float(np.dot(c,sums)/np.dot(c,n)))
        increments.append({
          "new_vs_ecology":name,
          "delta_logloss":float(v.mean()),
          "CI95_genus_bootstrap":[float(z) for z in np.percentile(bs,[2.5,97.5])],
          "AUC_delta":float(metrics[name]["AUC"]-metrics["ecology"]["AUC"]),
          "N_boot":N_BOOT
        })
    return {"n_species":len(y),"n_genera":len(gnames),"model_perf":metrics,
            "paired_generalization_test":increments}

def run(env,wcvp,out):
    dst=Path(out)
    require(not dst.exists(),"Do not overwrite a verified ecology test")
    d,source=prepare(signed_env(env),*signed_wcvp(wcvp)[:2])
    d["is_red"]=d.source_hue9.eq("red").astype(int)
    grouping=d.groupby("source_hue9").agg(
       n_species=("key","size"),
       n_introduced=("recorded_introduced","sum"),
       native_regions_median=("native_region_count","median"),
       introduced_regions_mean=("introduced_region_count","mean")
    )
    grouping["intro_fraction"]=grouping.n_introduced/grouping.n_species
    grouping=grouping.reset_index()
    status_counts=pd.crosstab(d.source_hue9,d.hue_status)
    analyses=[
      red_effect(d,"primary_9hue_red_vs_other_photo_labels"),
      red_effect(d,"equal_genus_weight",genus_weight=True),
      red_effect(d,"without_photo_status",status=False),
      red_effect(d[d.native_region_count.ge(2)].reset_index(drop=True),
                 "native_WCVP_TDWG3_regions_at_least_two"),
    ]
    dd=d.copy()
    dd["recorded_introduced"]=dd.introduced_region_count.ge(2).astype(int)
    analyses.append(red_effect(dd,"introduced_in_at_least_two_TDWG3_regions"))
    # Fully documented sensitivity excluding the ten genera with most red-photo
    # taxa. This is a robustness filter not a basis for selecting a nicer result.
    red_count=d.groupby("genus").is_red.sum().sort_values(ascending=False)
    exclude=list(red_count.head(10).index)
    dd=d.loc[~d.genus.isin(exclude)].reset_index(drop=True)
    analyses.append(red_effect(dd,"excluding_top_ten_red_photo_species_genera"))
    high=d[d.hue_status.eq("sampling_consensus_single")].reset_index(drop=True)
    analyses.append(red_effect(high,"three_photo_consensus_selected"))
    contrast=red_vs_white_direct(d)
    holdout=cv_red_white(d)
    dst.mkdir(parents=True)
    grouping.to_csv(dst/"hue_introduction_frequencies.csv",index=False)
    status_counts.to_csv(dst/"hue_photo_evidence_counts.csv")
    report={
      "source":source,
      "source_env_sha256":"d8b97528188d9f8d63cb802d5b3f6e8f7f251937fd30af42aeae47e43ebad1fc",
      "red_species_distribution":grouping.to_dict("records"),
      "red_photo_evidence_class_counts":status_counts.to_dict("index"),
      "red_sensitivity":analyses,
      "red_vs_white_direct_within_genus":contrast,
      "out_of_genus_prediction":holdout,
      "ten_excluded_red_rich_genera":exclude,
      "exploratory_red_hue_found_in_first_nine_hue_screen":True,
      "primary_red_from_first_run_q_BH_nine_hues":0.0023619162305646975,
      "no_true_petal_colour_or_trait_evolution_identified":True,
      "WCVP_introduced_status_not_transplant_success":True,
      "source_accepted_species_once":True
    }
    (dst/"red_introduction_sensitivity.json").write_text(
        json.dumps(report,indent=2,allow_nan=False)+"\n")
    report["files_sha256"]={p.name:sha(p) for p in dst.iterdir()}
    (dst/"receipt.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print("MOIS_RED_INTRO_ECOLOGICAL_SENSITIVITY_DONE",json.dumps({
      "n_red":int(d.is_red.sum()),
      "red_intro_fraction":float(grouping.loc[grouping.source_hue9.eq("red"),"intro_fraction"].iloc[0]),
      "effects":analyses,
      "red_white":contrast,
      "out_of_genus":holdout
    },allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--env",required=True)
    p.add_argument("--wcvp",required=True)
    p.add_argument("--out",required=True)
    a=p.parse_args()
    run(a.env,a.wcvp,a.out)
