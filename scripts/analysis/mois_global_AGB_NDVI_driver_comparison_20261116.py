#!/usr/bin/env python3
"""Which factors explain photo-colour labels? Three independent model families.

Full 36,753 GBIF-accepted plant species in genus-disjoint 3-fold validation.
Predictor groups: temperature, moisture/VPD, solar, wind, NDVI, ESA-CCI woody
aboveground biomass. Complementary within-genus linear probability model with
genus-clustered SE. Interpret all results as associations, NOT causal factors.

Training species do NOT repeat across GBIF occurrence cells. No local petal
colour, native status or validated per-species floral masks are available.
"""
from __future__ import annotations
import argparse,hashlib,io,json
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier,HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import roc_auc_score,balanced_accuracy_score,log_loss
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from scripts.analysis.mois_global_photo9_random_forest_20261116 import (
    PREDICTORS,CLIMATE,COLORS,columns_for
)

SOURCE_SHA="d8b97528188d9f8d63cb802d5b3f6e8f7f251937fd30af42aeae47e43ebad1fc"
SOURCE_MEMBER="BIO19_weather_NDVI_ESACCIv7_AGB_37104_accepted_species.parquet"
K=3
SEED=20261116
MOISTURE="precipitation_humidity"
BIOMASS="woody_biomass"
NDVI="NDVI_vegetation_greenness"
ENV_GROUPS=list(CLIMATE)+[NDVI,BIOMASS]
GROUPS={**PREDICTORS,BIOMASS:[
  "log_ESACCI_AGB_woody_Mg_ha_2020",
  "ESA_AGB_woody_positive_area_fraction_2020",
  "ESA_AGB_2020_minus_2010_Mg_ha"]}

def require(ok,msg):
    if not bool(ok):raise ValueError("MOIS_MULTIMODEL_BIOMASS_DRIVERS: "+str(msg))
def sha(b):
    if isinstance(b,(Path,str)):b=Path(b).read_bytes()
    return hashlib.sha256(b).hexdigest()

def columns(group_names):
    return list(dict.fromkeys(c for g in group_names for c in GROUPS[g]))

def read_source(path):
    import pyarrow.parquet as pq
    b=Path(path).read_bytes()
    require(sha(b)==SOURCE_SHA,"signed ESA-CCI biomass+WorldClim+NDVI ZIP changed")
    with ZipFile(io.BytesIO(b)) as z:
        raw=z.read(SOURCE_MEMBER)
        receipt=json.loads(z.read("receipt.json"))
    require(sha(raw)==receipt["output_SHA256"][SOURCE_MEMBER],
        "ESA biomass 37104 species file SHA changed")
    require(receipt["source_global_accepted_species"]==37104
        and receipt["source_GBIF_valid_occurrences"]==396142825
        and receipt["species_with_both_AGB_NDVI_BIO19"]==36753,
        "AGB global species source denominator mismatch")
    d=pq.read_table(io.BytesIO(raw)).to_pandas()
    require(len(d)==37104 and d.key.is_unique
        and d.records_passed_spatial_QA.sum()==396142825,
        "full 37104 source species and GBIF records not intact")
    require(set(d.source_hue9)==set(COLORS),"photo nine-colour labels drift")
    d["log_GBIF_occurrences"]=np.log1p(d.records_passed_spatial_QA.astype(float))
    d["log_ESACCI_AGB_woody_Mg_ha_2020"]=np.log1p(
        d.ESA_AGB_woody_Mg_ha_2020.clip(lower=0))
    d=d[(d.occupied_grid_cells>0)&d.BIO01.notna()
        &d.BIO19.notna()
        &d.MODIS_NDVI_2001_2021_mean.notna()
        &d.ESA_AGB_woody_Mg_ha_2020.notna()
        &d.mean_abs_latitude_over_occupied_cells.notna()].copy()
    require(len(d)==36753 and d.genus.nunique()>6000,
        "matched complete AGB+NDVI+WorldClim species population changed")
    for grp in ENV_GROUPS:
        require(set(GROUPS[grp]).issubset(d.columns),"model group missing actual source "+grp)
    require(d.ESA_AGB_woody_Mg_ha_2020.ge(0).all()
      and d.ESA_AGB_woody_positive_area_fraction_2020.between(0,1).all(),
      "ESA CCI woody AGB source outside legitimate values")
    d["white_photo"]=d.source_hue9.eq("white").astype(int)
    return d,receipt

def genus_folds(d):
    y=d.source_hue9.to_numpy(str)
    groups=d.genus.to_numpy(str)
    splitter=StratifiedGroupKFold(n_splits=K,shuffle=True,random_state=SEED)
    folds=list(splitter.split(np.zeros(len(d)),y,groups))
    require(len(folds)==K,"three independent heldout genus folds missing")
    for tr,te in folds:
        require(set(groups[tr]).isdisjoint(set(groups[te])),"genus leakage in split")
        require(set(y[tr])==set(COLORS) and
            len(np.unique(d.white_photo.to_numpy()[te]))==2,"photo colour lacks training/validation classes")
    require(sum(len(te) for _,te in folds)==len(d),"some source species not heldout")
    return folds

def fit_predict(kind,X,y,tr,te,fold,multiclass=False):
    imputer=SimpleImputer(strategy="median",keep_empty_features=True)
    trX=imputer.fit_transform(X[tr])
    teX=imputer.transform(X[te])
    require(np.isfinite(trX).all() and np.isfinite(teX).all(),"imputation introduced nonfinite")
    if kind=="RF":
        m=RandomForestClassifier(n_estimators=125,min_samples_leaf=30,
            class_weight="balanced_subsample",bootstrap=True,max_samples=.82,
            n_jobs=2,random_state=SEED+fold)
    elif kind=="HGB":
        m=HistGradientBoostingClassifier(max_iter=150,
            max_leaf_nodes=15,min_samples_leaf=90,learning_rate=.05,
            l2_regularization=3,random_state=SEED+fold)
    elif kind=="logistic":
        scaler=StandardScaler().fit(trX)
        trX=scaler.transform(trX);teX=scaler.transform(teX)
        m=LogisticRegression(C=.10,solver="lbfgs",max_iter=750,
            class_weight="balanced",random_state=SEED+fold)
    else:raise ValueError("Unknown model")
    m.fit(trX,y[tr])
    p=m.predict_proba(teX)
    return m,imputer,p,m.classes_

def eval_white(target,prob):
    return {"AUC":float(roc_auc_score(target,prob)),
      "balanced_accuracy":float(balanced_accuracy_score(target,prob>=.5)),
      "logloss":float(log_loss(target,np.c_[1-prob,prob],labels=[0,1]))}

def permutation_importance_full(forest,imputer,rawX,yt,feature_names,
                                groups,seed,rep=2):
    # This captures correlated variables within a group jointly; however
    # between-group correlations mean causal attribution is not identifiable.
    rng=np.random.default_rng(seed)
    classidx=list(forest.classes_).index(1)
    baseline=roc_auc_score(yt,forest.predict_proba(imputer.transform(rawX))[:,classidx])
    out=[]
    for group in groups:
        pos=[feature_names.index(c) for c in GROUPS[group]]
        samples=[]
        for i in range(rep):
            ix=rng.permutation(len(rawX))
            scrambled=rawX.copy()
            scrambled[:,pos]=rawX[ix][:,pos]
            p=forest.predict_proba(imputer.transform(scrambled))[:,classidx]
            samples.append(float(baseline-roc_auc_score(yt,p)))
        out.append({"group":group,"baseline_AUC":float(baseline),
                    "mean_drop_AUC":float(np.mean(samples)),
                    "min_drop_AUC":float(min(samples)),
                    "max_drop_AUC":float(max(samples)),
                    "repetitions":rep})
    return out

def joint_within_genus(d):
    # This specifically tests within-genus conditional associations with
    # WHITE PHOTO CLASS, not species-level RF predictive importance.
    # One-to-one species rows and genus-clustered SE.
    mixed=d.groupby("genus").white_photo.nunique().eq(2)
    s=d.loc[d.genus.map(mixed)].copy()
    require(len(s)>26000 and s.genus.nunique()>2000,"within-genus contrast unavailable")
    matrix={
      "temperature_BIO1":s.BIO01,
      "temperature_seasonality_BIO4":s.BIO04,
      "moisture_log_BIO12":np.log1p(s.BIO12.clip(lower=0)),
      "moisture_BIO15":s.BIO15,
      "moisture_vpd":s.vpd_kPa_approx,
      "moisture_vapour_pressure":s.vapr_monthly_mean,
      "solar_radiation":s.srad_monthly_mean,
      "wind":s.wind_monthly_mean,
      "ndvi_mean":s.MODIS_NDVI_2001_2021_mean,
      "ndvi_variability":s.MODIS_NDVI_2001_2021_sd,
      "woody_AGB_2020_log":s.log_ESACCI_AGB_woody_Mg_ha_2020,
      "woody_AGB_forest_fraction":s.ESA_AGB_woody_positive_area_fraction_2020,
      "woody_AGB_2010to2020_change":s.ESA_AGB_2020_minus_2010_Mg_ha,
      "log_GBIF_occurrences":s.log_GBIF_occurrences,
      "absolute_latitude":s.mean_abs_latitude_over_occupied_cells
    }
    X=pd.DataFrame(matrix,index=s.index)
    for status in sorted(s.hue_status.unique())[1:]:
        X["photo_status_"+status]=s.hue_status.eq(status).astype(int)
    # Complete cases guard the conditional inference; never replace missing
    # remotely-sensed biomass with zero woodland biomass.
    valid=X.notna().all(axis=1)
    X=X.loc[valid];s=s.loc[valid]
    Y=s.white_photo.astype(float)
    require(len(X)>25000 and s.genus.nunique()>2000,"too few complete within-genus species")
    base=X.copy()
    for col in X.columns:
        if col.startswith("photo_status_"):continue
        mean=X[col].mean();sd=X[col].std()
        require(np.isfinite(sd) and sd>0,f"non-identifiable predictor: {col}")
        X[col]=(X[col]-mean)/sd
    X=X-X.groupby(s.genus).transform("mean")
    Y=Y-Y.groupby(s.genus).transform("mean")
    cols=X.columns[X.std()>1e-8]
    X=X[cols]
    # No intercept: demeaned within genus.
    model=sm.OLS(Y,X).fit(cov_type="cluster",
        cov_kwds={"groups":s.genus},use_t=False)
    groups={
      "temperature":["temperature_BIO1","temperature_seasonality_BIO4"],
      "moisture":["moisture_log_BIO12","moisture_BIO15","moisture_vpd","moisture_vapour_pressure"],
      "solar":["solar_radiation"],"wind":["wind"],
      "NDVI":["ndvi_mean","ndvi_variability"],
      "woody_AGB":["woody_AGB_2020_log","woody_AGB_forest_fraction",
                   "woody_AGB_2010to2020_change"]
    }
    stats=[]
    for group,names in groups.items():
        require(set(names).issubset(cols),"one within-genus model column lost")
        r=np.zeros((len(names),len(cols)))
        for i,n in enumerate(names):r[i,cols.get_loc(n)]=1
        t=model.wald_test(r,scalar=True)
        stats.append({"group":group,"within_genus_joint_chi2":float(t.statistic),
                      "df":len(names),"p":float(t.pvalue)})
    q=multipletests([x["p"] for x in stats],method="fdr_bh")[1]
    for z,v in zip(stats,q):z["q_BH_six_groups"]=float(v)
    estimates=[]
    for name in ("moisture_log_BIO12","moisture_vpd","ndvi_mean","woody_AGB_2020_log",
                 "woody_AGB_forest_fraction","woody_AGB_2010to2020_change"):
        interval=model.conf_int().loc[name]
        estimates.append({"predictor":name,"linear_probability_change_per_global_1sd":
                          float(model.params[name]),"CI95_low":float(interval.iloc[0]),
                          "CI95_high":float(interval.iloc[1]),
                          "p_genus_clustered":float(model.pvalues[name])})
    return {"within_genus_mixed_white_nonwhite_species":len(s),
            "within_genus_mixed_genera":int(s.genus.nunique()),
            "joint_group_tests":stats,
            "selected_coefficients":estimates,
            "caution":"These are within-genus associations, not biological causal driver estimates."}

def evaluate(d,out):
    out=Path(out);require(not out.exists(),"refuse prior output overwrite")
    folds=genus_folds(d)
    y=d.white_photo.to_numpy(int)
    y9=d.source_hue9.to_numpy(str)
    records=[];importance=[]
    # Three model families, same genus-separated folds and exactly same species.
    jobs=[]
    for kind in ("RF","HGB","logistic"):
        jobs.append((kind,"all_environment",ENV_GROUPS))
        if kind in ("RF","HGB"):
            for group in (MOISTURE,BIOMASS,NDVI):
                jobs.append((kind,"without_"+group,[x for x in ENV_GROUPS if x!=group]))
        if kind=="RF":
            for group in ("temperature","solar_radiation","wind"):
                jobs.append((kind,"without_"+group,[x for x in ENV_GROUPS if x!=group]))
            jobs.append((kind,"plus_geography_effort",
                ENV_GROUPS+["geographic_effort"]))
    # Full 9-hue RF: if the feature group is important for all colours not
    # just whites, nine-class macro-balanced accuracy should also change.
    ninejobs=[("RF","all_environment",ENV_GROUPS),
              ("RF","without_"+MOISTURE,[x for x in ENV_GROUPS if x!=MOISTURE]),
              ("RF","without_"+BIOMASS,[x for x in ENV_GROUPS if x!=BIOMASS])]
    for kind,variant,groups in jobs:
        feature=columns(groups);X=d[feature].to_numpy(float)
        for fold,(tr,te) in enumerate(folds):
            forest,imputer,prob,classes=fit_predict(kind,X,y,tr,te,fold)
            class_idx=list(classes).index(1)
            test=eval_white(y[te],prob[:,class_idx])
            records.append({"model_type":kind,"outcome":"white_vs_eight_other_photo_hues",
               "predictor_variant":variant,"fold":fold,"heldout_species_n":len(te),
               "heldout_genus_n":len(set(d.genus.to_numpy()[te])),"features":len(feature),**test})
            if kind=="RF" and variant=="all_environment":
                perf=permutation_importance_full(forest,imputer,X[te],y[te],
                    feature,groups,seed=SEED+fold,rep=3)
                importance.extend({**row,"fold":fold} for row in perf)
        print("MOIS_AGB_RF_MODEL_VARIANT_COMPLETE",kind,variant,flush=True)
    for kind,variant,groups in ninejobs:
        feature=columns(groups);X=d[feature].to_numpy(float)
        for fold,(tr,te) in enumerate(folds):
            rf,imputer,prob,classes=fit_predict(kind,X,y9,tr,te,fold+40,multiclass=True)
            predicted=classes[prob.argmax(axis=1)]
            acc=balanced_accuracy_score(y9[te],predicted)
            records.append({"model_type":kind,"outcome":"nine_photo_hue_classes",
                "predictor_variant":variant,"fold":fold,"heldout_species_n":len(te),
                "heldout_genus_n":len(set(d.genus.to_numpy()[te])),
                "features":len(feature),"AUC":float("nan"),
                "balanced_accuracy":float(acc),
                "logloss":float(log_loss(y9[te],prob,labels=list(classes)))})
        print("MOIS_AGB_RF_NINE_COLOUR_VARIANT_COMPLETE",variant,flush=True)

    alt=joint_within_genus(d)
    out.mkdir(parents=True)
    scores=pd.DataFrame(records);groups=pd.DataFrame(importance)
    scores.to_csv(out/"AGB_NDVI_moisture_RF_HGB_logistic_3fold_genus_disjoint.csv",index=False)
    groups.to_csv(out/"RF_moisture_NDVI_AGB_grouped_permutation_all_three_folds.csv",index=False)
    sums=scores.groupby(["model_type","outcome","predictor_variant"],as_index=False).agg(
        folds=("fold","size"),mean_AUC=("AUC","mean"),
        minimum_fold_AUC=("AUC","min"),maximum_fold_AUC=("AUC","max"),
        mean_balanced_accuracy=("balanced_accuracy","mean"),
        mean_logloss=("logloss","mean"))
    sums.to_csv(out/"three_model_families_and_9hue_summary.csv",index=False)
    (out/"genus_fixed_effects_conditional_group_Wald.json").write_text(
       json.dumps(alt,indent=2,sort_keys=True)+"\n")
    for f in sums.itertuples(index=False):
        print("MOIS_BIOMASS_GROUP_DRIVER_MODEL_PERFORMANCE",
            f.model_type,f.outcome,f.predictor_variant,f.mean_AUC,f.mean_balanced_accuracy,flush=True)
    rf=sums[(sums.model_type=="RF")
          &(sums.outcome=="white_vs_eight_other_photo_hues")]
    base=float(rf.loc[rf.predictor_variant.eq("all_environment"),"mean_AUC"].iloc[0])
    rf_drops={}
    for group in (MOISTURE,BIOMASS,NDVI,"temperature","solar_radiation","wind"):
        value=float(rf.loc[rf.predictor_variant.eq("without_"+group),"mean_AUC"].iloc[0])
        rf_drops[group]=base-value
    p=groups.groupby("group").mean_drop_AUC.agg(["mean","min","max"]).to_dict("index")
    receipt={
      "schema":"MOIS_37104_9photo_BIO19_real_NDVI_real_ESA_AGB_multiple_methods_drivers_v1",
      "source_accepted_plant_taxa":37104,
      "source_eligible_complete_species":len(d),
      "distinct_source_genera":int(d.genus.nunique()),
      "source_ESA_feature_ZIP_SHA256":SOURCE_SHA,
      "model_unit":"accepted_species_one_row",
      "three_fold_genus_disjoint":True,
      "RF_white_full_environment_mean_AUC":base,
      "RF_loss_when_drop_environment_group_AUC":rf_drops,
      "RF_group_permutation_AUC_loss_three_folds":p,
      "alternative_non_RF_methods":["HistogramGradientBoosting","L2 LogisticRegression",
                                    "within-genus linear probability model, cluster-robust Wald"],
      "genus_FE_group_joint_tests":alt,
      "nine_photo_hue_RF_also_tested":True,
      "all_results_descriptive_not_causal":True,
      "NDVI_measures_green_vegetation_not_biomass":True,
      "ESACCI_AGB_measures_woody_aboveground_biomass_not_all_plants":True,
      "GBIF_species_range_not_native_only":True,
      "photo_hue_not_locally_measured_at_GBIF_coordinate":True,
      "not_validated_petal_colour_for_all_taxa":True,
      "not_ecologically_decisive_causal_predictor":True,
      "output_SHA256":{x.name:sha(x) for x in out.iterdir()}
    }
    (out/"receipt.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print("MOIS_MULTIMETHOD_TRUE_AGB_NDVI_MOISTURE_DRIVER_COMPARISON",
        json.dumps({k:receipt[k] for k in ("source_eligible_complete_species",
         "RF_white_full_environment_mean_AUC",
         "RF_loss_when_drop_environment_group_AUC",
         "RF_group_permutation_AUC_loss_three_folds",
         "alternative_non_RF_methods")},sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--AGB-features",required=True)
    parser.add_argument("--output",required=True)
    args=parser.parse_args()
    df,_=read_source(args.AGB_features)
    evaluate(df,args.output)
