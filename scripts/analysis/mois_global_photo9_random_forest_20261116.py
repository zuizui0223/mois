#!/usr/bin/env python3
"""Species-level random forests: GLOBAL photo9 vs actual BIO19, monthly meteorology, NDVI.

Models never treat the same species in multiple GBIF grid cells as independent
training observations. Genus-heldout folds avoid same-genus memorization.
One field of interest: whether real MODIS NDVI adds out-of-genus prediction
beyond WorldClim, and which *groups* of correlated features carry signal.

Model predictive association != causal flower-colour determinant. NDVI is
greenness, not standing biomass or NPP. Source hue is not verified petal truth.
"""
from __future__ import annotations
import argparse,hashlib,io,json
from collections import defaultdict
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import (
  roc_auc_score,balanced_accuracy_score,average_precision_score,log_loss,
  accuracy_score,confusion_matrix
)

FEATURE_ZIP_SHA="e08e81d1c2afd1a2ee5312f3ebc4aeb5ab7678988c2453a2a9c63d28c82eb17a"
INPUT_MEMBER="WorldClim19_weather_NDVI_37104_COLXR_photo9_species.parquet"
COLORS=("white","yellow","orange","red","pink","magenta","purple","blue","bronze")
RF_TREES=140
K=3
SEED=20261116

PREDICTORS={
 "temperature":[f"BIO{i:02d}" for i in range(1,12)]+[
    "tavg_monthly_mean","tavg_monthly_range","tavg_monthly_sd","tavg_warmest_month_C"],
 "precipitation_humidity":[f"BIO{i:02d}" for i in range(12,20)]+[
    "prec_monthly_mean","prec_monthly_range","prec_monthly_sd",
    "prec_annual_sum_of_monthly_normal_mm","vapr_monthly_mean",
    "vapr_monthly_range","vapr_monthly_sd","vpd_kPa_approx"],
 "solar_radiation":["srad_monthly_mean","srad_monthly_range","srad_monthly_sd"],
 "wind":["wind_monthly_mean","wind_monthly_range","wind_monthly_sd"],
 "NDVI_vegetation_greenness":[
    "MODIS_NDVI_2001_2021_mean","MODIS_NDVI_2001_2021_sd"],
 "geographic_effort":["log_GBIF_occurrences","mean_abs_latitude_over_occupied_cells"]
}
CLIMATE=["temperature","precipitation_humidity","solar_radiation","wind"]
FULL=["temperature","precipitation_humidity","solar_radiation","wind",
      "NDVI_vegetation_greenness","geographic_effort"]

def require(ok,msg):
    if not bool(ok): raise ValueError("MOIS_GLOBAL_RANDOM_FOREST: "+msg)
def sha(x):
    if isinstance(x,(str,Path)):x=Path(x).read_bytes()
    return hashlib.sha256(x).hexdigest()

def read_signed(path):
    import pyarrow.parquet as pq
    raw=Path(path).read_bytes()
    require(sha(raw)==FEATURE_ZIP_SHA,"full 37k real env features ZIP SHA mismatch")
    with ZipFile(io.BytesIO(raw)) as z:
        receipt=json.loads(z.read("receipt.json"))
        data=z.read(INPUT_MEMBER)
    require(sha(data)==receipt["output_SHA256"][INPUT_MEMBER],
            "full 37k features Parquet SHA mismatch")
    require(receipt["global_accepted_species"]==37104
            and receipt["unique_100km_cells"]==18577
            and receipt["WorldClim_all_19_BIO"]
            and receipt["NDVI_is_biomass"] is False,
            "source real WorldClim NDVI validation not complete")
    d=pq.read_table(io.BytesIO(data)).to_pandas()
    require(len(d)==37104 and d.key.is_unique and d.inat_taxon_id.is_unique,
            "one accepted species per row not conserved")
    require(d.records_passed_spatial_QA.sum()==396142825
            and set(d.source_hue9)==set(COLORS),"source original RGB labels or GBIF denominator changed")
    d["log_GBIF_occurrences"]=np.log1p(d.records_passed_spatial_QA.astype(float))
    for group,columns in PREDICTORS.items():
        require(set(columns).issubset(d),"missing real predictor group "+group)
    d=d[(d.occupied_grid_cells>0)
        &d.BIO01.notna() &d.BIO19.notna()
        &d.MODIS_NDVI_2001_2021_mean.notna()
        &d.mean_abs_latitude_over_occupied_cells.notna()].copy()
    require(len(d)>30000 and d.genus.nunique()>2500,"insufficient globally sampled species")
    d["white_binary"]=d.source_hue9.eq("white").astype(int)
    return d,receipt

def columns_for(groups):
    return list(dict.fromkeys(f for g in groups for f in PREDICTORS[g]))

def make_model(seed=SEED,white=True,genus_weight=False):
    return RandomForestClassifier(
        n_estimators=RF_TREES,
        max_features="sqrt",
        min_samples_leaf=20,
        class_weight="balanced_subsample",
        max_depth=None,
        n_jobs=2,
        random_state=seed,
        bootstrap=True,
        max_samples=0.82
    )

def fit_model(X,y,tr,te,white,fold,genus=None):
    imput=SimpleImputer(strategy="median",add_indicator=False,keep_empty_features=True)
    train=imput.fit_transform(X[tr])
    valid=imput.transform(X[te])
    require(np.isfinite(train).all() and np.isfinite(valid).all(),
            "non-finite predictors after training-only imputation")
    forest=make_model(seed=SEED+fold,white=white)
    # Optional inverse sqrt genus sample weights give large genera less influence
    # while retaining all species. In binary/main analysis all taxa equal weight.
    forest.fit(train,y[tr])
    pred=forest.predict_proba(valid)
    labels=forest.classes_
    return forest,imput,valid,pred,labels

def white_scores(target,prob):
    require(len(target)==len(prob) and 0<=prob.min() and prob.max()<=1,
            "out of bounds predicted white risk")
    return {"auc_ROC":float(roc_auc_score(target,prob)),
       "average_precision":float(average_precision_score(target,prob)),
       "balanced_accuracy":float(balanced_accuracy_score(target,prob>=.5)),
       "log_loss":float(log_loss(target,np.c_[1-prob,prob],labels=[0,1]))}

def multiclass_scores(target,pred,classes):
    labels=np.asarray(COLORS)
    ix={v:i for i,v in enumerate(classes)}
    padded=np.zeros((len(target),len(COLORS)),float)
    for i,label in enumerate(labels):
        if label in ix:padded[:,i]=pred[:,ix[label]]
    require((padded.sum(axis=1)>0.999).all() and (padded.sum(axis=1)<1.001).all(),
            "multiclass predicted probabilities invalid")
    guessed=labels[padded.argmax(axis=1)]
    return {
      "balanced_accuracy":float(balanced_accuracy_score(target,guessed)),
      "top1_accuracy":float(accuracy_score(target,guessed)),
      "log_loss":float(log_loss(target,padded,labels=list(COLORS)))}

def group_importance(forest,imputer,raw_holdout,y_binary,cols,groups,repeats=4,seed=SEED):
    base=forest.predict_proba(imputer.transform(raw_holdout))
    classidx=list(forest.classes_).index(1)
    base_auc=roc_auc_score(y_binary,base[:,classidx])
    rng=np.random.default_rng(seed)
    results=[]
    for label in groups:
        names=PREDICTORS[label]
        ix=[cols.index(s) for s in names if s in cols]
        if not ix:continue
        drops=[]
        for i in range(repeats):
            perm=rng.permutation(len(raw_holdout))
            changed=raw_holdout.copy()
            # One common permutation preserves within-group correlation.
            changed[:,ix]=raw_holdout[perm][:,ix]
            altered=forest.predict_proba(imputer.transform(changed))
            val=roc_auc_score(y_binary,altered[:,classidx])
            drops.append(base_auc-val)
        results.append({"group":label,"num_correlated_variables":len(ix),
          "heldout_grouped_permutation_auc_drop_mean":float(np.mean(drops)),
          "heldout_grouped_permutation_auc_drop_min":float(np.min(drops)),
          "heldout_grouped_permutation_auc_drop_max":float(np.max(drops)),
          "baseline_valid_AUC":float(base_auc),
          "repeats":repeats,
          "negative_importance_allowed":True})
    return results

def evaluate(d,output):
    output=Path(output)
    require(not output.exists(),"refuse overwrite model run")
    ybin=d.white_binary.to_numpy(int)
    y9=d.source_hue9.to_numpy(str)
    genus=d.genus.to_numpy(str)
    cv=StratifiedGroupKFold(n_splits=K,shuffle=True,random_state=SEED)
    fold_plan=list(cv.split(np.zeros(len(d)),y9,genus))
    require(all(set(genus[tr]).isdisjoint(set(genus[te]))
                for tr,te in fold_plan),"genus overlap across test and train")
    datasets={
      "BIO19_weather_only":CLIMATE,
      "BIO19_weather_plus_real_NDVI":CLIMATE+["NDVI_vegetation_greenness"],
      "full_with_geography_effort":FULL
    }
    results=[];importance=[];predictions=[]
    for dataset,groups in datasets.items():
        cols=columns_for(groups)
        X=d[cols].to_numpy(float)
        for fold,(tr,te) in enumerate(fold_plan):
            require(len(np.unique(ybin[te]))==2 and
                len(np.unique(y9[tr]))==9,"missing recorded class in split")
            forest,imp,valid,pred,classes=fit_model(X,ybin,tr,te,True,fold)
            classidx=list(classes).index(1)
            p=pred[:,classidx]
            rec={"model":"white_vs_other","predictor_set":dataset,
                 "fold":fold,"num_train_species":len(tr),"num_test_species":len(te),
                 "train_genera":len(set(genus[tr])),"test_genera":len(set(genus[te])),
                 "num_predictors":len(cols),**white_scores(ybin[te],p)}
            results.append(rec)
            predictions.extend({"key":str(key),"genus":str(g),"photo_modal9":str(h),
               "source_colour_confidence":str(status),"fold":fold,
               "predictor_set":dataset,"predicted_white_probability":float(prob),
               "observed_photo_white":int(actual)}
                for key,g,h,status,prob,actual in zip(
                   d.key.to_numpy()[te],genus[te],y9[te],d.hue_status.to_numpy()[te],p,ybin[te]))
            if dataset=="full_with_geography_effort" and fold==0:
                importance+=group_importance(forest,imp,X[te],ybin[te],cols,FULL,repeats=5)
    # Nine-color task has very rare orange class; report macro-accuracy,
    # not only ordinary accuracy that is dominated by white classes.
    for dataset in ("BIO19_weather_only","BIO19_weather_plus_real_NDVI"):
        cols=columns_for(datasets[dataset])
        X=d[cols].to_numpy(float)
        for fold,(tr,te) in enumerate(fold_plan):
            forest,imp,valid,pred,classes=fit_model(X,y9,tr,te,False,fold+42)
            vals=multiclass_scores(y9[te],pred,classes)
            results.append({"model":"nine_photo_hue_categories","predictor_set":dataset,
              "fold":fold,"num_train_species":len(tr),"num_test_species":len(te),
              "train_genera":len(set(genus[tr])),"test_genera":len(set(genus[te])),
              "num_predictors":len(cols),**vals})
    # High-confidence sensitivity: species modal hue matches mean9; do not
    # label as expert-verified flower petal truth.
    subset=d.loc[~d.modal_mean9_disagree].copy()
    require(len(subset)>27000,"quality-sensitivity subset too small")
    X=subset[columns_for(CLIMATE+["NDVI_vegetation_greenness"])].to_numpy(float)
    y=subset.white_binary.to_numpy(int);gen=subset.genus.to_numpy(str)
    sub_scores=[]
    cv2=StratifiedGroupKFold(n_splits=3,shuffle=True,random_state=SEED+1)
    for fold,(tr,te) in enumerate(cv2.split(X,subset.source_hue9.to_numpy(str),gen)):
        forest,im,valid,prob,labels=fit_model(X,y,tr,te,True,fold+88)
        sub_scores.append({"fold":fold,"species_in_quality_sample":len(subset),
          **white_scores(y[te],prob[:,list(labels).index(1)])})

    output.mkdir(parents=True)
    scores=pd.DataFrame(results)
    scores.to_csv(output/"RF_real_worldclim_MODIS_NDVI_genus_holdout_cv.csv",index=False)
    pd.DataFrame(importance).to_csv(output/"RF_grouped_feature_importance_heldout_genus.csv",index=False)
    pd.DataFrame(predictions).to_csv(output/"RF_all_species_out_of_genus_predictions.csv.gz",
       index=False,compression={"method":"gzip","mtime":0})
    pd.DataFrame(sub_scores).to_csv(output/"RF_photo_hue_modal_mean_agree_sensitivity.csv",index=False)
    agg=scores.groupby(["model","predictor_set"]).agg(
       n_folds=("fold","size"),
       auc_mean=("auc_ROC","mean"),
       auc_min=("auc_ROC","min"),
       auc_max=("auc_ROC","max"),
       balanced_acc_mean=("balanced_accuracy","mean"),
       logloss_mean=("log_loss","mean")).reset_index()
    agg.to_csv(output/"RF_model_comparison_summary.csv",index=False)
    binary=scores[scores.model.eq("white_vs_other")].pivot(index="fold",
           columns="predictor_set",values="auc_ROC")
    diffs=(binary["BIO19_weather_plus_real_NDVI"]-binary["BIO19_weather_only"])
    receipt={
      "schema":"mois_RF_37104_real_BIO19_monthly_meteorology_MODIS_NDVI_genus_holdout_v1",
      "training_species":len(d),"training_distinct_genera":d.genus.nunique(),
      "training_white_photo_species":int(ybin.sum()),
      "training_nine_hue_counts":{str(k):int(v) for k,v in d.source_hue9.value_counts().items()},
      "CV":"3_fold_stratified_genus_disjoint_species",
      "no_species_by_grid_pseudoreplication":True,
      "RF_estimators":RF_TREES,
      "RF_min_samples_leaf":20,
      "number_of_feature_groups":len(PREDICTORS),
      "worldclim19_weather_predictor_columns":columns_for(CLIMATE),
      "real_MODIS_NDVI_predictor_columns":PREDICTORS["NDVI_vegetation_greenness"],
      "NDVI_additional_out_of_genus_auc_by_fold":[float(x) for x in diffs],
      "NDVI_additional_out_of_genus_auc_mean":float(diffs.mean()),
      "strict_modal_mean_hue_agreement_subset_n":len(subset),
      "strict_modal_mean_colour_agree_is_not_botanical_validation":True,
      "GBIF_source_species_colour_measured_locally":False,
      "RANDOM_FOREST_IMPORTANCE_IS_CAUSAL_EFFECT":False,
      "actual_biomass_or_NPP_raster_integrated":False,
      "worldclim_period":"1970-2000",
      "MODIS_NDVI_period":"2001-2021",
      "feature_source_ZIP_SHA256":FEATURE_ZIP_SHA,
      "result_file_SHA256":{p.name:sha(p) for p in output.iterdir()}}
    (output/"receipt.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print("MOIS_RANDOM_FOREST_REAL_ALL_SPECIES_BIO19_WEATHER_NDVI",
          json.dumps({
           "training_species":receipt["training_species"],
           "training_distinct_genera":receipt["training_distinct_genera"],
           "NDVI_auc_gain":receipt["NDVI_additional_out_of_genus_auc_mean"],
           "models":agg.to_dict("records"),
           "group_importance":importance},sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--features",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    d,_=read_signed(a.features)
    evaluate(d,a.output)
