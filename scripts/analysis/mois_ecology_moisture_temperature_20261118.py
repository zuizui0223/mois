#!/usr/bin/env python3
"""MOIS ecological context test (18 Nov 2026): temperature × precipitation
against photo-derived species hue, not verified petal colour or refugial ranges.
Same 37104 accepted taxa in two signed sources; same species rows for all models.
Within-genus linear probability, genus-clustered CI; all 9 hues BH; genus-disjoint
OOF logistic predictive sensitivity. Ecological associations, NOT selection.
"""
from __future__ import annotations
import argparse, hashlib, io, json
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, log_loss

MAIN_SHA="1da364ff63c506f9e2c06d1d71c9609a2c48fec7a428a21ba2c3c4c6ff8f4800"
ENV_SHA="d8b97528188d9f8d63cb802d5b3f6e8f7f251937fd30af42aeae47e43ebad1fc"
MAIN_MEMBER="GBIF37104_photo9_species_NPP23years_LGM_AGB_NDVI.parquet"
ENV_MEMBER="BIO19_weather_NDVI_ESACCIv7_AGB_37104_accepted_species.parquet"
NPP="NPP_mean_GBIF_2001_2023_kgC_m2_yr"
NPP_COVER="NPP_fraction_occupied_cells_18plus_years"
LGM_T="threeGCM_delta_temperature_C_median"
LGM_P="threeGCM_delta_log_precip_median"
HUES=("white","yellow","orange","red","pink","magenta","purple","blue","bronze")
SEED=20261118

def require(ok,msg):
    if not bool(ok):raise ValueError("MOIS_ECOLOGY: "+str(msg))

def sha(inp):
    if isinstance(inp,(str,Path)):inp=Path(inp).read_bytes()
    return hashlib.sha256(inp).hexdigest()

def source(path,expected,member):
    import pyarrow.parquet as pq
    require(sha(path)==expected,"Signed source ZIP changed")
    with ZipFile(path) as z:
        buf=z.read(member)
        receipt=json.loads(z.read("receipt.json"))
    require(sha(buf)==receipt["output_SHA256"][member],"Source Parquet SHA changed")
    d=pq.read_table(io.BytesIO(buf)).to_pandas()
    d["key"]=d.key.astype(str)
    require(len(d)==37104 and d.key.is_unique,"Accepted species denominator wrong")
    return d

def assemble(main,env):
    require(set(main.key)==set(env.key),"Accepted species sources do not match")
    use=["key","source_hue9","genus","hue_status",NPP,NPP_COVER,
         LGM_T,LGM_P,"lgm_climate_coverage"]
    d=env.merge(main[use],on="key",how="inner",validate="one_to_one",
                suffixes=("","_historic"))
    require(len(d)==37104 and d.key.is_unique,"Source join damaged")
    require((d.source_hue9.astype(str)==d.source_hue9_historic.astype(str)).all(),
            "Source photo labels mismatch")
    require((d.genus.astype(str)==d.genus_historic.astype(str)).all(),
            "Source genus names mismatch")
    keys=["BIO01","BIO12","vpd_kPa_approx","srad_monthly_mean",
          "wind_monthly_mean","MODIS_NDVI_2001_2021_mean",
          "ESA_AGB_woody_Mg_ha_2020","mean_abs_latitude_over_occupied_cells",
          "records_passed_spatial_QA","occupied_grid_cells",NPP,NPP_COVER,
          LGM_T,LGM_P,"lgm_climate_coverage","hue_status"]
    require(set(keys).issubset(d.columns),"Source ecological features absent")
    keep=d[keys].notna().all(axis=1)
    keep &= d[NPP_COVER].ge(.75)&d.lgm_climate_coverage.ge(.75)
    keep &= d.BIO12.ge(0)&d.ESA_AGB_woody_Mg_ha_2020.ge(0)
    keep &= d.occupied_grid_cells.gt(0)&d.records_passed_spatial_QA.ge(0)
    d=d.loc[keep].copy().reset_index(drop=True)
    require(len(d)>35000 and d.genus.nunique()>6000,
            "Unexpected global valid species universe")
    require(set(d.source_hue9)==set(HUES),"Source nine photo-hues changed")
    return d

def design(d):
    a=pd.DataFrame({
        "temperature":d.BIO01.astype(float),
        "log_precip":np.log1p(d.BIO12.astype(float)),
        "vpd":d.vpd_kPa_approx.astype(float),
        "solar":d.srad_monthly_mean.astype(float),
        "wind":d.wind_monthly_mean.astype(float),
        "NPP":d[NPP].astype(float),
        "NDVI":d.MODIS_NDVI_2001_2021_mean.astype(float),
        "woody_AGB":np.log1p(d.ESA_AGB_woody_Mg_ha_2020.astype(float)),
        "abs_latitude":d.mean_abs_latitude_over_occupied_cells.astype(float),
        "GBIF_records":np.log1p(d.records_passed_spatial_QA.astype(float)),
        "LGM_warming":d[LGM_T].astype(float),
        "LGM_rain":d[LGM_P].astype(float)
    },index=d.index)
    require(np.isfinite(a.to_numpy()).all(),"Nonfinite controls")
    m=a.mean();sd=a.std(ddof=1)
    require(sd.gt(1e-10).all(),"No variance")
    x=(a-m)/sd
    x["temperature_X_log_precip"]=x.temperature*x.log_precip
    x["temperature_X_vpd"]=x.temperature*x.vpd
    dummy=pd.get_dummies(d.hue_status.astype(str),prefix="photo_status",
                         drop_first=True,dtype=float)
    x=pd.concat([x,dummy],axis=1)
    require(x.columns.is_unique and np.isfinite(x.to_numpy()).all(),"Model matrix bad")
    return x,{k:{"mean":float(m[k]),"SD":float(sd[k])} for k in a}

def genus_fit(d,x,hue="white",var="temperature_X_log_precip",
              weighting="species_equal",adjust_status=True):
    y=d.source_hue9.eq(hue).astype(float)
    groups=d.genus.astype(str)
    keep=y.groupby(groups).transform("nunique").eq(2)
    ii=np.flatnonzero(keep.to_numpy())
    g=groups.iloc[ii]; yy=y.iloc[ii]
    cols=[c for c in x.columns if c not in
          ("temperature_X_log_precip","temperature_X_vpd")]
    cols.append(var)
    if not adjust_status:cols=[c for c in cols if not c.startswith("photo_status_")]
    xx=x.iloc[ii][cols]
    xc=xx-xx.groupby(g).transform("mean")
    yc=yy-yy.groupby(g).transform("mean")
    xc=xc.loc[:,xc.std()>1e-9]
    require(len(ii)>100 and g.nunique()>20 and var in xc.columns,
            "Not enough mixed-colour genera")
    if weighting=="genus_equal":
        w=1/g.map(g.value_counts()).to_numpy(float)
    else:
        w=np.ones(len(ii),float)
    fit=sm.WLS(yc,xc,weights=w).fit(cov_type="cluster",
                                   cov_kwds={"groups":g},use_t=False)
    result={
      "hue":hue,"term":var,"weighting":weighting,"status_control":adjust_status,
      "n_species_mixed_genus":int(len(ii)),"n_genera_mixed_hue":int(g.nunique()),
      "beta_fraction":float(fit.params[var]),
      "CI95_fraction":[float(v) for v in fit.conf_int().loc[var]],
      "p_genus_cluster":float(fit.pvalues[var]),
      "beta_log_precip":float(fit.params["log_precip"]),
      "n_parameters":len(xc.columns),
      "ecological_main_effects":{
        k:{"beta_pp_per_SD":float(fit.params[k]*100),
           "CI95_pp":[float(v*100) for v in fit.conf_int().loc[k]],
           "p_genus_cluster":float(fit.pvalues[k])}
        for k in ("log_precip","vpd","temperature","solar","wind",
                  "NPP","NDVI","woody_AGB","LGM_warming")
        if k in fit.params
      }
    }
    if var=="temperature_X_log_precip":
        v=fit.cov_params().loc[["log_precip",var],["log_precip",var]].to_numpy()
        for q in (.2,.5,.8):
            t=float(x.temperature.quantile(q))
            c=np.array([1,t])
            est=float(c@fit.params[["log_precip",var]].to_numpy())
            err=float(np.sqrt(c@v@c))
            result[f"precip_effect_temp_quantile_{round(q*100)}"]={
              "temperature_global_z":t,
              "white_probability_pp_per_1SD_precip":100*est,
              "CI95_pp":[100*(est-1.96*err),100*(est+1.96*err)]}
    return result

def estimate_all(d,x):
    out=[]
    for weighting in ("species_equal","genus_equal"):
        row=[]
        for hue in HUES:
            r=genus_fit(d,x,hue,weighting=weighting)
            r["cohort"]="all"
            row.append(r)
        q=multipletests([r["p_genus_cluster"] for r in row],method="fdr_bh")[1]
        for r,p in zip(row,q):r["BH_q_nine_hues_this_weighting"]=float(p)
        # Prespecified ecological colour-palette check: water covariates
        # precipitation and atmospheric dryness VPD across all nine hues.
        # Correct across all 18 comparisons, not one favourable hue.
        tests=[(r,k) for r in row for k in ("log_precip","vpd")]
        qwater=multipletests([r["ecological_main_effects"][k]["p_genus_cluster"]
                             for r,k in tests],method="fdr_bh")[1]
        for (r,k),p in zip(tests,qwater):
            r["ecological_main_effects"][k]["BH_q_18_water_tests_this_weighting"]=float(p)
        out.extend(row)
    r=genus_fit(d,x,var="temperature_X_vpd");r["cohort"]="all";out.append(r)
    r=genus_fit(d,x,adjust_status=False);r["cohort"]="all";out.append(r)
    high=d.hue_status.eq("sampling_consensus_single")
    if high.sum()>1000:
        d2=d.loc[high].reset_index(drop=True)
        xx,_=design(d2)
        r=genus_fit(d2,xx);r["cohort"]="consensus_photo_selection";out.append(r)
    return out

def genus_disjoint_cv(d,x):
    y=d.source_hue9.eq("white").to_numpy(int)
    hue=d.source_hue9.to_numpy(str)
    groups=d.genus.to_numpy(str)
    A=x.drop(columns=["temperature_X_log_precip","temperature_X_vpd"]).to_numpy(float)
    B=x.drop(columns=["temperature_X_vpd"]).to_numpy(float)
    pa=np.full(len(y),np.nan);pb=np.full(len(y),np.nan)
    fold_metrics=[]
    splits=StratifiedGroupKFold(n_splits=3,shuffle=True,random_state=SEED)
    for fi,(tr,te) in enumerate(splits.split(A,hue,groups)):
        require(set(groups[tr]).isdisjoint(set(groups[te])),"Genus leak")
        for xx,p in ((A,pa),(B,pb)):
            z=StandardScaler().fit(xx[tr])
            model=LogisticRegression(C=.1,max_iter=600,solver="lbfgs")
            model.fit(z.transform(xx[tr]),y[tr])
            p[te]=model.predict_proba(z.transform(xx[te]))[:,1]
        fold_metrics.append({
            "fold":fi,"n_test_species":len(te),
            "AUC_no_interaction":float(roc_auc_score(y[te],pa[te])),
            "AUC_with_interaction":float(roc_auc_score(y[te],pb[te])),
            "logloss_no_interaction":float(log_loss(y[te],pa[te])),
            "logloss_with_interaction":float(log_loss(y[te],pb[te]))
        })
    require(np.isfinite(pa).all() and np.isfinite(pb).all(),"Holdout incomplete")
    a=-(y*np.log(np.clip(pa,1e-12,1))+(1-y)*np.log(np.clip(1-pa,1e-12,1)))
    b=-(y*np.log(np.clip(pb,1e-12,1))+(1-y)*np.log(np.clip(1-pb,1e-12,1)))
    delta=b-a
    gid,g=np.unique(groups,return_inverse=True)
    sums=np.bincount(g,weights=delta,minlength=len(gid))
    sizes=np.bincount(g,minlength=len(gid))
    rng=np.random.default_rng(SEED)
    boot=[]
    for _ in range(2000):
        counts=np.bincount(rng.integers(len(gid),size=len(gid)),minlength=len(gid))
        boot.append(float(np.dot(counts,sums)/np.dot(counts,sizes)))
    return {
        "N_species":len(y),"N_genera":len(gid),
        "AUC_no_interaction":float(roc_auc_score(y,pa)),
        "AUC_with_interaction":float(roc_auc_score(y,pb)),
        "logloss_no_interaction":float(a.mean()),
        "logloss_with_interaction":float(b.mean()),
        "delta_logloss_with_minus_without":float(delta.mean()),
        "delta_logloss_genus_bootstrap_CI95":[float(v) for v in np.percentile(boot,[2.5,97.5])],
        "folds":fold_metrics
    }

def run(mainzip,envzip,dest):
    path=Path(dest)
    require(not path.exists(),"Cannot overwrite existing signed analysis")
    main=source(mainzip,MAIN_SHA,MAIN_MEMBER)
    env=source(envzip,ENV_SHA,ENV_MEMBER)
    d=assemble(main,env)
    x,scales=design(d)
    estimates=estimate_all(d,x)
    holdout=genus_disjoint_cv(d,x)
    path.mkdir(parents=True)
    tab=pd.DataFrame([{k:v for k,v in r.items() if not k.startswith("precip_effect_temp_")}
                      for r in estimates])
    tab.to_csv(path/"hue_temperature_moisture_interactions.csv",index=False)
    (path/"all_effects.json").write_text(json.dumps(estimates,indent=2,allow_nan=False)+"\n")
    (path/"genus_disjoint_white_predictive_test.json").write_text(json.dumps(holdout,indent=2,allow_nan=False)+"\n")
    primary=estimates[0]
    report={
      "original_accepted_species":37104,"n_analyzed_species":len(d),
      "n_analyzed_genera":int(d.genus.nunique()),
      "n_white_photo_labels":int(d.source_hue9.eq("white").sum()),
      "primary_hypothesis":"Within-genus precipitation association changes with temperature; not ecological causation",
      "primary_white_result":primary,
      "genus_disjoint_holdout":holdout,
      "scales":scales,
      "input_sha256":{"NPP_LGM":MAIN_SHA,"BIO19_NDVI_AGB":ENV_SHA},
      "source_white_is_photo_proxy_not_independently_verified_petals":True,
      "species_environment_is_mean_of_current_GBIF_range":True,
      "GBIF_native_status_not_validated_here":True,
      "LGM_current_range_climate_not_fossil_refugium":True,
      "9_hue_BH_adjusted_in_separate_genus_weight_variants":True
    }
    (path/"summary.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    report["sha256_output_files"]={p.name:sha(p) for p in path.iterdir()}
    (path/"receipt.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print("MOIS_ECOLOGY_THERMAL_MOISTURE_SUCCESS",json.dumps({
       "n_species":len(d),"n_genera":int(d.genus.nunique()),"white":primary,
       "prediction":{k:holdout[k] for k in ("AUC_no_interaction","AUC_with_interaction",
                  "delta_logloss_with_minus_without","delta_logloss_genus_bootstrap_CI95")},
       "nine_hues":[{**{k:r.get(k) for k in
                    ("hue","beta_fraction","p_genus_cluster","BH_q_nine_hues_this_weighting")},
                     "precip":r["ecological_main_effects"]["log_precip"],
                     "vpd":r["ecological_main_effects"]["vpd"]}
                    for r in estimates[:9]]},allow_nan=False),flush=True)

def smoke():
    rng=np.random.default_rng(42);n=2700
    d=pd.DataFrame({
        "genus":[f"g{i//5:04d}" for i in range(n)],
        "source_hue9":[HUES[i%9] for i in range(n)],
        "hue_status":["consensus" if i%3==0 else "provisional" for i in range(n)],
        "BIO01":rng.normal(14,9,n),"BIO12":rng.uniform(100,3000,n),
        "vpd_kPa_approx":rng.uniform(.2,2.6,n),
        "srad_monthly_mean":rng.uniform(1000,15000,n),
        "wind_monthly_mean":rng.uniform(1,10,n),
        "MODIS_NDVI_2001_2021_mean":rng.uniform(.1,.8,n),
        "ESA_AGB_woody_Mg_ha_2020":rng.uniform(0,200,n),
        "mean_abs_latitude_over_occupied_cells":rng.uniform(0,67,n),
        "records_passed_spatial_QA":rng.integers(1,20000,n),
        NPP:rng.uniform(.1,2.4,n),LGM_T:rng.uniform(2,15,n),
        LGM_P:rng.normal(0,.5,n)
    })
    x,_=design(d)
    m=genus_fit(d,x);o=genus_disjoint_cv(d,x)
    require(m["n_species_mixed_genus"]>1000 and o["N_species"]==n,
            "Synthetic genus test failed")
    print("MOIS_ECOLOGY_SYNTHETIC_PASS",flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--main");p.add_argument("--env");p.add_argument("--out")
    p.add_argument("--selftest",action="store_true")
    a=p.parse_args()
    if a.selftest:smoke()
    else:
        require(a.main and a.env and a.out,"main/env/out required")
        run(a.main,a.env,a.out)
