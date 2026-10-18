#!/usr/bin/env python3
"""MOIS sensitivity: which climate group best predicts photographed white hue?

Full 19 BIO/VAPR/srad data cover far more original photos than the partially
known source soil. Compare train-only-PCA, held-out-genus and 10-degree geographic
block CV on (1) all climate-supported originals, (2) terrain-complete originals,
(3) all-33-complete originals but using only climate predictors, then compare
the third subset to the previously archived all-33 five-block experiment.
No photo hue is inferred from missing source, no original pixel reprocessing,
no fitness or exact 50:50 adaptive-equilibrium interpretation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold

from scripts.analysis import mois_global313535_all_environment_pca_ablation as base

N=313535
CLIMATE_BLOCKS=("temperature","moisture","solar")
TERRAIN_BLOCKS=(*CLIMATE_BLOCKS,"terrain")
ALL_COVARIATES=base.ABIOTIC
OUT_RANK="mois_313535_climate_vs_soil_sampling_holdout_group_rank.csv"
OUT_COVER="mois_313535_climate_soil_sampling_coverage_and_white_fraction.csv"
OUT_CMP="mois_313535_soil_complete_full33_vs_climate_only_comparison.csv"
OUT_R="mois_313535_climate_source_sample_selection_importance.json"
RUN_REFERENCE="mois-313535-all19-BIO-vapr-PCA-ablation-v1"
NAME_OLD="mois_313535_all_environment_group_ablation_blocked_cv.csv"
RECEIPT_OLD="mois_313535_all_environment_pca_confounded_ranking_receipt.json"


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for part in iter(lambda:f.read(1<<20), b""):
            h.update(part)
    return h.hexdigest()


def build_panel(d,blocks,condition=()):
    fields=tuple(dict.fromkeys((
        *base.GEO,
        *(k for name in blocks for k in base.BLOCKS[name]),
        *condition,
    )))
    source=d.loc[d.classified & d.genus.notna()].dropna(subset=fields).copy()
    source=source.loc[source.inat_taxon_id.astype(str).str.len().gt(0)]
    useful=tuple(dict.fromkeys((
        *base.GEO,
        *(k for name in blocks for k in base.BLOCKS[name]),
    )))
    agg={k:"mean" for k in useful}
    agg.update({"white":"mean","genus":"first","photo_id":"size"})
    s=source.groupby("inat_taxon_id",sort=True).agg(agg).rename(
        columns={"white":"photographed_white_rate",
                 "photo_id":"original_classifiable_photos"})
    counts=s.groupby("genus").size()
    s=s.loc[s.genus.isin(counts[counts>=2].index)].copy()
    require(len(s)>=1000 and s.genus.nunique()>=100,
            "Climate complete photo cohort too small")
    s["latblock"]=np.floor((s.latitude+90)/10).astype(int)
    lon=np.rad2deg(np.arctan2(s.longitude_sin.to_numpy(float),
                             s.longitude_cos.to_numpy(float)))
    s["lonblock"]=np.floor((lon+180)/10).astype(int)
    s["spatial_holdout_group"]=s.latblock.astype(str)+"_"+s.lonblock.astype(str)
    require(s.spatial_holdout_group.nunique()>=base.FOLDS,
            "Not enough distinct geographic photo blocks for crossvalidation")
    return s,{
        "classifiable_original_photo_rows_before_genus_gate":int(len(source)),
        "eligible_source_taxa":int(len(s)),
        "qualifying_genera":int(s.genus.nunique()),
        "geographic_taxon_centroid_blocks":int(s.spatial_holdout_group.nunique()),
        "source_photos_are_NOT_independent_flowering_individuals":True,
    }


def evaluate(s,blocks,design):
    groups=(s.genus if design=="heldout_genera" else
            s.spatial_holdout_group).astype(str)
    y=s.photographed_white_rate.to_numpy(dtype=float)
    cum={}
    keys=["geo_only","all_selected_blocks",*("drop_"+g for g in blocks)]
    component_counts=[]
    for train,test in GroupKFold(n_splits=base.FOLDS).split(s,y,groups):
        xa,xb,_=base.fit_block_transform(s.iloc[train][list(base.GEO)],
                                         s.iloc[test][list(base.GEO)],None)
        tx={"geography":xa};vx={"geography":xb}
        axes={}
        for name in blocks:
            u,v,k=base.fit_block_transform(
                s.iloc[train][list(base.BLOCKS[name])],
                s.iloc[test][list(base.BLOCKS[name])],base.BLOCKS[name])
            tx[name]=u;vx[name]=v;axes[name]=k
        component_counts.append(axes)
        for key in keys:
            selected=([] if key=="geo_only" else
                      list(blocks) if key=="all_selected_blocks" else
                      [b for b in blocks if b!=key[5:]])
            xtrain=np.column_stack([tx["geography"],*(tx[k] for k in selected)])
            xtest=np.column_stack([vx["geography"],*(vx[k] for k in selected)])
            sse,tss=base.ridge_score(s,train,test,xtrain,xtest)
            cum.setdefault(key,[]).append((sse,tss))
    scores={key:float(1-sum(a for a,b in v)/sum(b for a,b in v))
            for key,v in cum.items()}
    rows=[]
    for key in keys:
        feature=key.removeprefix("drop_") if key.startswith("drop_") else ""
        rows.append({
            "heldout_design":design,
            "model":key,
            "excluded_environment_group":feature,
            "cv_weighted_R2":scores[key],
            "full_minus_drop_group_cv_R2":(
                scores["all_selected_blocks"]-scores[key]
                if key.startswith("drop_") else None),
            "full_minus_geography_baseline_R2":(
                scores["all_selected_blocks"]-scores["geo_only"]
                if key=="all_selected_blocks" else None),
            "train_only_standardized_PCA_90pct_up_to_6_per_block":True,
            "individual_predictor_causal_importance_not_identified":True,
        })
    return rows,{"holdout_design":design,
                 "train_fold_environment_axis_counts":component_counts,
                 "full_cv_R2":scores["all_selected_blocks"],
                 "geo_only_R2":scores["geo_only"]}


def source_stats(d,extra_fields):
    source=d.dropna(subset=list(extra_fields))
    photo=source.loc[source.classified]
    species=photo.groupby("inat_taxon_id").white.mean()
    genera=photo.assign(genus=photo.species.map(base.prev.genus))
    gen_means=(genera.groupby(["genus","inat_taxon_id"]).white.mean()
                      .groupby(level=0).mean())
    return {
        "processed_source_photo_attempts":len(source),
        "classifiable_photo_hues":int(len(photo)),
        "photo_white_fraction":float(photo.white.mean()) if len(photo) else None,
        "photo_colour_missingness_rate":float(1-len(photo)/len(source))
            if len(source) else None,
        "source_taxon_equal_white_fraction":float(species.mean())
            if len(species) else None,
        "source_genus_equal_white_fraction":float(gen_means.mean())
            if len(gen_means) else None,
        "percent_original_photo_attempts_retained":100*len(source)/len(d),
    }


def run(release,abiotic,old_result,out):
    require(not out.exists() or not any(out.iterdir()),
            "Never overwrite a previous source-locked environmental analysis")
    d,photorec,abiorec=base.read_source(release,abiotic)
    old_doc=json.loads((old_result/RECEIPT_OLD).read_text())
    require(old_doc["schema"]=="mois_313535_full33_abiotic_PCA_blocked_predictive_ablation_v1"
            and old_doc["extended_abiotic_source_SHA"]==
                abiorec["new_environment_parquet_sha256"]
            and old_doc["photo_source_SHA"]==
                photorec["full_16_workers_original_pixel_source_SHA256"]
            and old_doc["n_abiotic_variables_considered"]==33,
            "Original all-33 actual-photo analysis not source-identical")
    original_abl=old_result/NAME_OLD
    require(sha(original_abl)==old_doc["output_sha256"][NAME_OLD],
            "Original published all-33 predictive ranking CSV SHA changed")
    ref=pd.read_csv(original_abl)
    climates=tuple(dict.fromkeys(
        v for name in CLIMATE_BLOCKS for v in base.BLOCKS[name]))
    panel_defs=[
        ("ALL_climate_measured_source",CLIMATE_BLOCKS,climates),
        ("TERRAIN_covered_climate_and_terrain",TERRAIN_BLOCKS,
         (*climates,*base.BLOCKS["terrain"])),
        ("SOIL_AND_TERRAIN_complete_but_CLIMATE_ONLY",CLIMATE_BLOCKS,
         (*climates,*base.BLOCKS["terrain"],*base.BLOCKS["soil"])),
    ]
    models=[];census=[];panels={};meta={}
    for label,blocks,fields in panel_defs:
        cohort=source_stats(d,fields)
        panel,pop=build_panel(d,blocks,condition=fields)
        meta[label]=pop
        panels[label]=panel
        census.append({"data_population":label,**cohort,**pop})
        for heldout in ("heldout_genera",
                        "heldout_geographic_10degree_taxon_centroids"):
            rows,result=evaluate(panel,blocks,heldout)
            models.extend({"data_population":label,**row,
                           "n_taxa":len(panel),
                           "n_genera":panel.genus.nunique()} for row in rows)
    models=pd.DataFrame(models)
    census=pd.DataFrame(census)
    third="SOIL_AND_TERRAIN_complete_but_CLIMATE_ONLY"
    reference=ref[ref.model.str.startswith("drop_")]
    third_panel=panels[third]
    base_panel=old_doc["source_panel"]
    require(len(third_panel)==base_panel[
        "source_taxa_in_regressions_after_congeneric_gate"]
        and third_panel.genus.nunique()==base_panel[
            "source_nominal_genera_after_gate"]
        and int(meta[third][
            "classifiable_original_photo_rows_before_genus_gate"])==
                base_panel["source_classified_original_photo_rows_exact_complete_case"],
        "The soil/terrain-complete climate subset differs from all-33 original study")
    contrast=[]
    for design in ("heldout_genera",
                   "heldout_geographic_10degree_taxon_centroids"):
        b=models[(models.data_population==third)&
                 (models.heldout_design==design)]
        for group in CLIMATE_BLOCKS:
            h=float(b.loc[b.model=="drop_"+group,
                "full_minus_drop_group_cv_R2"].iloc[0])
            f=float(reference.loc[(reference.heldout_design==design)&
                   (reference.environment_group==group),
                   "delta_R2_full_minus_without_this_group"].iloc[0])
            contrast.append({
                "heldout_design":design,
                "environment_group":group,
                "data_population":third,
                "climate_only_comparator_drop_cv_delta_R2":h,
                "full33_conditional_drop_cv_delta_R2":f,
                "same_taxa_and_photo_support_but_DIFFERENT_predictor_sets":True,
                "not_a_causal_experiment_or_algorithm_independent_confirmation":True,
            })
    ranks={}
    for cohort in models.data_population.unique():
        ranks[cohort]={}
        for dsgn in models.heldout_design.unique():
            subset=models[(models.data_population==cohort)&
                          (models.heldout_design==dsgn)&
                          models.model.str.startswith("drop_")].copy()
            subset=subset.sort_values("full_minus_drop_group_cv_R2",ascending=False)
            ranks[cohort][dsgn]=[
                {"group":row.excluded_environment_group,
                 "delta_R2":float(row.full_minus_drop_group_cv_R2)}
                for row in subset.itertuples()]
    out.mkdir(parents=True,exist_ok=True)
    files={OUT_RANK:models,OUT_COVER:census,OUT_CMP:pd.DataFrame(contrast)}
    for path,data in files.items():data.to_csv(out/path,index=False,lineterminator="\n")
    report={
        "schema":"mois_313535_climate_main_dataset_soil_selection_robustness_v1",
        "full_photo_attempts":base.prev.N_PHOTOS,
        "all_environment_original_BIO19_vapr_soil_terrain_release":RUN_REFERENCE,
        "source_original_photo_sha":old_doc["photo_source_SHA"],
        "source_complete_abiotic_parquet_sha":old_doc["extended_abiotic_source_SHA"],
        "prior_33_abiotic_cv_csv_verified_sha":old_doc["output_sha256"][NAME_OLD],
        "climate_only_predictors":list(climates),
        "same_environment_support_comparison_definitions":
            {x:[*y] for x,y in
             [("climate",climates),
              ("terrain_plus_climate",(*climates,*base.BLOCKS["terrain"])),
              ("full33",(*climates,*base.BLOCKS["terrain"],
                         *base.BLOCKS["soil"]))]},
        "cohort_and_taxon_support":census.to_dict("records"),
        "group_ablation_rankings":ranks,
        "original_all33_on_soil_complete_rankings":
            old_doc["rankings_by_heldout_group"],
        "temperature_is_robust_top_across_three_photographic_sampling_scopes":
            bool(all(
                all(v[0]["group"]=="temperature" and v[0]["delta_R2"]>0
                    for v in item.values()) for item in ranks.values())),
        "soil_missingness_not_imputed_and_only_real_photo_pixels_used":True,
        "blockwise_PCA_and_scaler_fit_only_on_training_folds":True,
        "no_open_petals_UV_pigment_cost_reproductive_fitness_50_50_equilibrium_claim":True,
        "no_fcp_intraspecific_white_darkness_estimator":True,
        "annual_1970_2000_WorldClim_not_weather_in_observation_month":True,
        "output_sha256":{key:sha(out/key) for key in files},
    }
    (out/OUT_R).write_text(json.dumps(report,indent=2,sort_keys=True,
                                       ensure_ascii=False,allow_nan=False)+"\n")
    print("MOIS_CLIMATE_ALL313535_VS_SOIL_COMPLETE_BLOCK_CV",json.dumps({
        "sampling_cohorts":report["cohort_and_taxon_support"],
        "rankings":ranks,
        "temperature_top_all_scope":report[
            "temperature_is_robust_top_across_three_photographic_sampling_scopes"]
    },sort_keys=True),flush=True)
    return report

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--photo-release",type=Path,required=True)
    p.add_argument("--all19-abiotic-release",type=Path,required=True)
    p.add_argument("--previous33-pca-release",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    run(a.photo_release,a.all19_abiotic_release,a.previous33_pca_release,a.output)
