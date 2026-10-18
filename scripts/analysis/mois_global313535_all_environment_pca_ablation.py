#!/usr/bin/env python3
"""MOIS: all available source-measured abiotic predictors with PCA and blocked CV.

Global flower-like photographed WHITE vs chromatic (not maintained 50:50).
PCA is UNsupervised: PC1 is never declared most important for white.
Source-population matched soil, coast-distance unavailable, no imputation.
Predictive ablation refits after dropping groups/variables on unchanged
source taxon cohort; train-only PCA/scaling; held-out genus AND geography.
"""
from __future__ import annotations
import argparse,hashlib,json,re
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge

from scripts.analysis import mois_global313535_water_light_heat_habitat_partition as prev
from scripts.analysis import mois_original313535_full19_bio_vapr_terrain_lake as env

BLOCKS={
    "temperature":env.BIO[:11],
    "moisture":(*env.BIO[11:],env.VAPR),
    "solar":(prev.SOLAR,),
    "terrain":env.TERRAIN,
    "soil":prev.SOIL,
}
ABIOTIC=tuple(name for group in BLOCKS.values() for name in group)
GEO=("latitude","abs_lat","longitude_sin","longitude_cos","source_new")
OUT_MODELS="mois_313535_all_environment_group_ablation_blocked_cv.csv"
OUT_VARS="mois_313535_all_environment_individual_variable_drop_cv.csv"
OUT_PCA="mois_313535_all_environment_pca_loadings.csv"
OUT_AXES="mois_313535_all_environment_pca_axis_summary.csv"
OUT_CORR="mois_313535_all_environment_high_correlation_pairs.csv"
OUT_SUPPORT="mois_313535_all_environment_completecase_support.csv"
OUT_DECISION="mois_313535_all_environment_pca_confounded_ranking_receipt.json"
FOLDS=5
ALPHA_PER_TAXON=.01
VARIANCE=.90


def require(value,msg):
    if not value:raise ValueError(msg)


def sha(path:Path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):h.update(chunk)
    return h.hexdigest()


def read_source(source:Path,extended:Path):
    p,source_receipt=prev.load_source(source)
    doc=json.loads((extended/env.OUT_REC).read_text())
    require(doc["schema"]=="mois_actual313535_original_photo_worldclim_full19_vapr_dem_v1"
            and doc["real_original_processed_photo_ids"]==prev.N_PHOTOS
            and doc["original_four_hue_photo_labels"]==prev.N_HUES
            and doc["original_white_photo_ROIs"]==prev.WHITE_N
            and doc["previous_313535_photo_IDs_and_hues_unchanged"] is True
            and doc["new_environment_parquet_sha256"]==sha(extended/env.OUT_ENV),
            "Extended official original photo-site climate is not signed")
    a=pd.read_parquet(extended/env.OUT_ENV)
    require(len(a)==len(p) and p.photo_id.tolist()==a.photo_id.tolist()
            and len(ABIOTIC)==33 and len(set(ABIOTIC))==33,
            "Original 313535 pixel-source/19 BIO/terrain/soil misaligned")
    for name in (*env.BIO,env.VAPR,*env.TERRAIN):
        p[name]=pd.to_numeric(a[name],errors="coerce")
    p["abs_lat"]=p.latitude.abs()
    p["longitude_sin"]=np.sin(np.deg2rad(p.longitude))
    p["longitude_cos"]=np.cos(np.deg2rad(p.longitude))
    p["source_new"]=p.old_mois_or_new_fcp.astype(str).ne("historical").astype(float)
    p["genus"]=p.species.map(prev.genus)
    return p,source_receipt,doc


def source_support(p:pd.DataFrame):
    rows=[]
    for name,fields in [
        ("all_climate",(*BLOCKS["temperature"],*BLOCKS["moisture"],*BLOCKS["solar"])),
        ("climate_plus_terrain",(*BLOCKS["temperature"],*BLOCKS["moisture"],
                                *BLOCKS["solar"],*BLOCKS["terrain"])),
        ("climate_plus_soil",(*BLOCKS["temperature"],*BLOCKS["moisture"],
                             *BLOCKS["solar"],*BLOCKS["soil"])),
        ("ALL_33_ABIOTIC_SAME_PHOTO_SAMPLE",ABIOTIC),
    ]:
        good=p[p.classified & p.genus.notna()].dropna(subset=list((*GEO,*fields)))
        rows.append({
            "environment_population":name,
            "real_original_classified_photos_with_covariates":int(len(good)),
            "named_taxa":int(good.inat_taxon_id.nunique()),
            "nominal_genera":int(good.genus.nunique()),
            "source_phenotype_hue_missingness_not_imputed":True,
            "static_soil_values_not_filled_for_new_photos":True,
        })
    return pd.DataFrame(rows)


def species_panel(p:pd.DataFrame):
    # Crucial: photo-level complete cases FIRST. Never aggregate 5 variables
    # over mutually different photographed subsets within the same taxon.
    q=p[p.classified&p.genus.notna()].dropna(subset=list((*GEO,*ABIOTIC))).copy()
    q=q[q.inat_taxon_id.astype(str).str.len().gt(0)]
    numeric=list(dict.fromkeys((*GEO,*ABIOTIC)))
    ag={k:"mean" for k in numeric}
    ag.update({"white":"mean","genus":"first","photo_id":"size"})
    s=q.groupby("inat_taxon_id",sort=True).agg(ag).rename(
        columns={"white":"photographed_white_rate",
                 "photo_id":"classified_photo_count"})
    genera_counts=s.groupby("genus").size()
    s=s[s.genus.isin(genera_counts[genera_counts>=2].index)].copy()
    needcols=["latitude","longitude_sin","longitude_cos",*ABIOTIC]
    require(len(s)>1000 and s.genus.nunique()>100 and
            s[needcols].notna().all().all(),
            "No sufficient source-identical all-environment taxa/genera")
    s["latblock"]=np.floor((s.latitude+90)/10).astype(int)
    # Approximate source photo taxon centroid; not all observations
    # necessarily occupy this geographic cell.
    # longitude sin/cos means avoid naive averaging across date line.
    angles=np.rad2deg(np.arctan2(
        s.longitude_sin.to_numpy(float),s.longitude_cos.to_numpy(float)))
    s["lonblock"]=np.floor((angles+180)/10).astype(int)
    s["spatial_holdout_group"]=(s.latblock.astype(str)+"_"+
                                s.lonblock.astype(str))
    require(s.spatial_holdout_group.nunique()>=FOLDS,
            "Too few observed geographic holdout regions")
    return s,{
        "source_classified_original_photo_rows_exact_complete_case":int(len(q)),
        "source_taxa_in_regressions_after_congeneric_gate":int(len(s)),
        "source_nominal_genera_after_gate":int(s.genus.nunique()),
        "source_approximate_10degree_cell_groups":int(s.spatial_holdout_group.nunique()),
        "complete_case_photo_subset_identical_for_all_predictor_ablation":True,
        "species_means_not_species_trait_ground_truth":True,
    }


def weights_for_genus(groups):
    count=pd.Series(groups).value_counts()
    w=np.array([1/count[g] for g in groups],dtype=float)
    return w/w.mean()


def weighted_scores(y,pred,w,ref):
    pred=np.clip(pred,0.,1.)
    sse=float(np.sum(w*(y-pred)**2))
    tss=float(np.sum(w*(y-ref)**2))
    return sse,tss


def fit_block_transform(s_train,s_test,group=None):
    """Single shared transform per TRAIN fold. No global PCA leakage into CV."""
    train=np.asarray(s_train,dtype=float)
    test=np.asarray(s_test,dtype=float)
    scale=StandardScaler()
    x=scale.fit_transform(train)
    z=scale.transform(test)
    if group is None or len(group)==1:
        return x,z,1
    pca=PCA(n_components=VARIANCE,svd_solver="full")
    a=pca.fit_transform(x)
    b=pca.transform(z)
    # Limit overall number of axes to 6 per abiotic group (precommitted).
    k=min(6,a.shape[1])
    return a[:,:k],b[:,:k],k


def data_by_folds(s,splitname):
    keys=(s.genus if splitname=="heldout_genera" else
          s.spatial_holdout_group).astype(str)
    base=np.asarray(s.photographed_white_rate,dtype=float)
    splits=list(GroupKFold(n_splits=FOLDS).split(s,base,keys))
    for train,test in splits:
        Xtrain={}; Xtest={}; axes={}
        a,b,_=fit_block_transform(s.iloc[train][list(GEO)],
                                  s.iloc[test][list(GEO)],None)
        Xtrain["geography"]=a;Xtest["geography"]=b
        for name,features in BLOCKS.items():
            a,b,k=fit_block_transform(s.iloc[train][list(features)],
                                      s.iloc[test][list(features)],features)
            Xtrain[name]=a;Xtest[name]=b;axes[name]=k
        yield train,test,Xtrain,Xtest,axes


def ridge_score(s,train,test,trainX,testX):
    ytrain=s.photographed_white_rate.to_numpy(float)[train]
    ytest=s.photographed_white_rate.to_numpy(float)[test]
    weights=weights_for_genus(s.genus.astype(str).to_numpy()[train])
    wt_test=weights_for_genus(s.genus.astype(str).to_numpy()[test])
    estimator=Ridge(alpha=len(train)*ALPHA_PER_TAXON,fit_intercept=True)
    estimator.fit(trainX,ytrain,sample_weight=weights)
    predicted=estimator.predict(testX)
    reference=float(np.average(ytrain,weights=weights))
    return weighted_scores(ytest,predicted,wt_test,reference)


def cv_abiotic_importance(s):
    rows=[]
    all_blocks=list(BLOCKS)
    for holdout in ("heldout_genera","heldout_geographic_10degree_taxon_centroids"):
        sums={}
        pcs=[]
        for split_id,(train,test,ta,tb,axes) in enumerate(
            data_by_folds(s,holdout),1):
            pcs.append(axes)
            for key,selected in [
                ("geo_only",[]),
                ("ALL_five_environment_groups",all_blocks),
                *[("drop_"+name,[n for n in all_blocks if n!=name])
                  for name in all_blocks],
                *[("only_"+name,[name]) for name in all_blocks],
            ]:
                xa=np.column_stack([ta["geography"],*[
                    ta[name] for name in selected]])
                xb=np.column_stack([tb["geography"],*[
                    tb[name] for name in selected]])
                sse,tss=ridge_score(s,train,test,xa,xb)
                pair=sums.setdefault(key,[])
                pair.append((sse,tss,int(len(test))))
        scores={key:1-sum(q[0] for q in parts)/sum(q[1] for q in parts)
                for key,parts in sums.items()}
        for key,val in scores.items():
            group=key.replace("drop_","").replace("only_","")
            rows.append({
                "heldout_design":holdout,
                "model":key,
                "environment_group":group if group in BLOCKS else "",
                "whole_holdout_group_cv_R2":float(val),
                "delta_R2_full_minus_without_this_group":
                    float(scores["ALL_five_environment_groups"]-val)
                    if key.startswith("drop_") else None,
                "delta_R2_geo_plus_group_minus_geo_only":
                    float(val-scores["geo_only"])
                    if key.startswith("only_") else None,
                "full_model_cv_R2":float(scores["ALL_five_environment_groups"]),
                "geo_only_cv_R2":float(scores["geo_only"]),
                "source_taxa":int(len(s)),
                "source_genera":int(s.genus.nunique()),
                "n_folds":FOLDS,
                "no_training_or_PCA_parameter_leakage_into_test":True,
                "prediction_target":"source taxon photographed white fraction",
                "not_causal_variable_importance":True,
            })
    return pd.DataFrame(rows)


def cv_individual_variable_drop(s):
    # Raw annual/soil/terrain per-feature conditional sensitivity alongside
    # group PCA. Strong correlations mean an isolated variable is redundant,
    # so strongest single raw predictor is especially hard to identify.
    out=[]
    for holdout in ("heldout_genera","heldout_geographic_10degree_taxon_centroids"):
        keys=(s.genus if holdout=="heldout_genera" else
              s.spatial_holdout_group).astype(str)
        score_by_setting={}
        feature_sets={"all_raw_environment":ABIOTIC}
        for variable in ABIOTIC:
            feature_sets["remove_"+variable]=tuple(k for k in ABIOTIC if k!=variable)
        for train,test in GroupKFold(n_splits=FOLDS).split(s,s.photographed_white_rate,keys):
            scaler=StandardScaler()
            base_train=scaler.fit_transform(s.iloc[train][list(GEO)])
            base_test=scaler.transform(s.iloc[test][list(GEO)])
            e_scale=StandardScaler()
            all_train=e_scale.fit_transform(s.iloc[train][list(ABIOTIC)])
            all_test=e_scale.transform(s.iloc[test][list(ABIOTIC)])
            for setting,features in feature_sets.items():
                indices=[ABIOTIC.index(v) for v in features]
                xx=np.column_stack((base_train,all_train[:,indices]))
                tt=np.column_stack((base_test,all_test[:,indices]))
                outcome=ridge_score(s,train,test,xx,tt)
                score_by_setting.setdefault(setting,[]).append(outcome)
        def score(k):
            a=score_by_setting[k]
            return float(1-sum(x[0] for x in a)/sum(x[1] for x in a))
        full=score("all_raw_environment")
        for name in ABIOTIC:
            removed=score("remove_"+name)
            out.append({
                "heldout_design":holdout,"abiotic_variable":name,
                "environment_group":next(g for g,v in BLOCKS.items() if name in v),
                "conditional_removed_variable_cv_R2":removed,
                "all_raw_predictor_cv_R2":full,
                "unique_predictive_delta_R2_full_minus_without":
                    float(full-removed),
                "highly_correlated_environmental_variables_can_mask_unique_effect":True,
                "original_source_predictor_not_intervention":True,
            })
    return pd.DataFrame(out)


def pca_descriptives(s):
    allX=s[list(ABIOTIC)].to_numpy(float)
    z=StandardScaler().fit_transform(allX)
    p=PCA(n_components=min(10,len(ABIOTIC)),svd_solver="full").fit(z)
    axes=[];load=[]
    for i in range(len(p.explained_variance_ratio_)):
        axes.append({
            "pca_block":"all_33_abiotic_full_community_descriptive",
            "axis":f"PC{i+1}",
            "explained_abiotic_variance_fraction":float(p.explained_variance_ratio_[i]),
            "cumulative_variance_fraction":float(np.sum(p.explained_variance_ratio_[:i+1])),
            "full_panel_unsupervised_PCA_not_predictive_cross_validation":True,
        })
        for j,k in enumerate(ABIOTIC):
            load.append({"pca_block":"all_33_abiotic_full_community_descriptive",
                         "axis":f"PC{i+1}","variable":k,
                         "loading":float(p.components_[i,j])})
    for name,block in BLOCKS.items():
        z=StandardScaler().fit_transform(s[list(block)])
        if len(block)==1:
            ev=np.array([1.]);components=np.array([[1.]])
        else:
            obj=PCA(n_components=min(5,len(block)),svd_solver="full").fit(z)
            ev=obj.explained_variance_ratio_;components=obj.components_
        for k in range(len(ev)):
            axes.append({
                "pca_block":name,"axis":f"PC{k+1}",
                "explained_abiotic_variance_fraction":float(ev[k]),
                "cumulative_variance_fraction":float(ev[:k+1].sum()),
                "full_panel_unsupervised_PCA_not_predictive_cross_validation":True,
            })
            for j,feature in enumerate(block):
                load.append({"pca_block":name,"axis":f"PC{k+1}",
                             "variable":feature,"loading":float(components[k,j])})
    corrs=s[list(ABIOTIC)].corr().to_numpy(float)
    pairs=[]
    for i,k in enumerate(ABIOTIC):
        for j in range(i+1,len(ABIOTIC)):
            val=float(corrs[i,j])
            if abs(val)>=.75:
                pairs.append({"variable_A":k,"variable_B":ABIOTIC[j],
                              "pearson_on_matched_taxon_means":val,
                              "strong_collinearity_limits_unique_causal_ranking":True})
    return pd.DataFrame(axes),pd.DataFrame(load),pd.DataFrame(pairs)


def run(release,extended,output):
    require(not output.exists() or not any(output.iterdir()),
            "Do not overwrite a source-signed ecological model")
    photos,p_source,p_extended=read_source(release,extended)
    support=source_support(photos)
    taxon,panel=species_panel(photos)
    scores=cv_abiotic_importance(taxon)
    variables=cv_individual_variable_drop(taxon)
    axis,load,corr=pca_descriptives(taxon)
    output.mkdir(parents=True,exist_ok=True)
    outputs={
        OUT_MODELS:scores,OUT_VARS:variables,OUT_PCA:load,
        OUT_AXES:axis,OUT_CORR:corr,OUT_SUPPORT:support,
    }
    for name,data in outputs.items():
        data.to_csv(output/name,index=False,lineterminator="\n")
    importance=scores.loc[scores.model.str.startswith("drop_")].copy()
    table={}
    for design,rows in importance.groupby("heldout_design"):
        ordered=rows.sort_values(
            "delta_R2_full_minus_without_this_group",ascending=False)
        table[design]=[
            {"group":str(r.environment_group),
             "full_minus_drop_group_cv_R2":
                float(r.delta_R2_full_minus_without_this_group)}
            for r in ordered.itertuples()]
    top={k:v[0]["group"] for k,v in table.items()}
    toppositive=all(table[k][0]["full_minus_drop_group_cv_R2"]>0 for k in table)
    consensus=(len(set(top.values()))==1 and toppositive)
    pca1=axis[(axis.pca_block=="all_33_abiotic_full_community_descriptive") &
              (axis.axis=="PC1")]
    report={
        "schema":"mois_313535_full33_abiotic_PCA_blocked_predictive_ablation_v1",
        "real_source_photo_attempts":prev.N_PHOTOS,
        "measured_hue_labels":prev.N_HUES,
        "white_original_photo_ROIs":prev.WHITE_N,
        "all_19_WorldClim_BIO_and_vapr_solar_and_static_controls_used":True,
        "n_abiotic_variables_considered":len(ABIOTIC),
        "nuisance_geography_and_source_features":list(GEO),
        "source_abiotic_environment_groups":
            {k:list(v) for k,v in BLOCKS.items()},
        "source_panel":panel,
        "eligible_photo_support_by_environment_block":support.to_dict("records"),
        "photo_source_SHA":p_source["full_16_workers_original_pixel_source_SHA256"],
        "extended_abiotic_source_SHA":p_extended["new_environment_parquet_sha256"],
        "rankings_by_heldout_group":table,
        "same_top_positive_group_for_both_heldout_designs":bool(consensus),
        "tentative_strongest_group_only_if_both_holdouts_agree":(
            next(iter(top.values())) if consensus else None),
        "global_PCA_PC1_abiotic_variance_fraction":
            float(pca1.explained_abiotic_variance_fraction.iloc[0]),
        "PCA_does_not_establish_strongest_colour_predictor":True,
        "holdout_group_PCAs_standardized_and_fitted_on_training_only":True,
        "heldout_tests_blocked_by_genus_and_approx_10degree_source_taxon_centroid":True,
        "all_environment_importance_compared_on_exact_same_complete_case_photos":True,
        "soil_missingness_not_imputed":True,
        "full_313535_population_not_equal_to_soil_complete_case_population":True,
        "FCP_intraspecific_colour_polymorphism_not_imported":True,
        "no_verified_open_petals_pigment_fitness_cost_or_evolutionary_50_50":True,
        "annual_1970_2000_climate_normals_only_season_month_SI":True,
        "output_sha256":{name:sha(output/name) for name in outputs},
    }
    (output/OUT_DECISION).write_text(json.dumps(
        report,indent=2,ensure_ascii=False,sort_keys=True,
        allow_nan=False)+"\n")
    print("MOIS_ALL_ENV_SOURCE_CONFUNDING_PCA_AND_BLOCKED_IMPORTANCE",
          json.dumps({
            "source_n":panel,"PC1_abiotic_explained_variance":
                report["global_PCA_PC1_abiotic_variance_fraction"],
            "heldout_block_ranking":table,
            "consensus_strongest":report[
                "tentative_strongest_group_only_if_both_holdouts_agree"],
          },sort_keys=True),flush=True)
    return report


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--existing-photo-feature-release",type=Path,required=True)
    parser.add_argument("--full19-extended-abiotic-source",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    a=parser.parse_args()
    run(a.existing_photo_feature_release,
        a.full19_extended_abiotic_source,a.output)
