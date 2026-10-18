#!/usr/bin/env python3
"""MOIS: which VISIBLE coloured hue shifts when white is common/rare?

Source-original 313535 RGB attempts; full-19 BIO VAPR solar photo IDs.
Compare 4 automatic hues, not just white vs everything labelled coloured.
All climate terciles are defined on ALL photo attempts before hue filtering.
Predict taxon-level 4-class *photographic* composition on genus- and
geography-heldout source taxa with climate-block PCA trained within folds.
This is not a chemical anthocyanin assay or pigment cost/fitness inference.
"""
from __future__ import annotations

import argparse, hashlib, json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold

from scripts.analysis import mois_global313535_all_environment_pca_ablation as parent
from scripts.analysis import mois_global313535_climate_subset_sampling_robustness as climate

HUES=parent.prev.HUES
CLIMATE=("temperature","moisture","solar")
RAIN=parent.prev.RAIN
SOLAR=parent.prev.SOLAR
FOLDS=5
OUT_COMP="mois_313535_four_visible_hues_water_solar_tertile_composition.csv"
OUT_CV="mois_313535_four_hues_climate_group_predictive_brier_ablation.csv"
OUT_DELTA="mois_313535_coloured_hue_changes_wet_vs_dry_at_solar_tertile.csv"
OUT_RECEIPT="mois_313535_four_hue_climate_heterogeneity_provenance.json"


def require(ok,msg):
    if not ok:raise ValueError(msg)


def sha(path:Path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):h.update(chunk)
    return h.hexdigest()


def classify_climate(d):
    # Build bins without restricting to classifiable photographic hues.
    out=d.copy()
    limits={}
    for field,title in ((RAIN,"rain"),(SOLAR,"solar")):
        value=pd.to_numeric(out[field],errors="coerce")
        valid=value.dropna()
        require(len(valid)>100 and valid.nunique()>=10,
                "Source real photo climate predictor not available")
        lo,hi=np.quantile(valid,[1/3,2/3])
        require(lo<hi,"Source climate tertiles not distinguishable")
        out[title+"_tertile"]=pd.cut(value,
            [-np.inf,lo,hi,np.inf],labels=["low","middle","high"],
            include_lowest=True)
        limits[title]={"low_middle":float(lo),"middle_high":float(hi)}
    return out,limits


def four_hue_source_cells(d):
    rows=[]
    for (w,l),sub in d.dropna(subset=["rain_tertile","solar_tertile"]).groupby(
            ["rain_tertile","solar_tertile"],observed=True):
        known=sub[sub.classified]
        photo=len(sub);n=len(known)
        counts={h:int(known.morph.eq(h).sum()) for h in HUES}
        require(sum(counts.values())==n and n>30,
                "Original photographic four-hue source must partition classified ROI")
        row={
            "rainfall_tertile":str(w),"solar_tertile":str(l),
            "original_processed_photo_attempts":photo,
            "real_original_classifiable_ROI_photos":n,
            "technically_unknown_hue_photo_attempts":photo-n,
            "photo_hue_classification_success":float(n/photo),
            "named_taxa":int(known.inat_taxon_id.nunique()),
            "source_genera":int(known.genus.nunique()),
        }
        for h in HUES:
            row[f"{h}_photo_ROIs"]=counts[h]
            row[f"{h}_among_4class_photos"]=float(counts[h]/n)
        rows.append(row)
    table=pd.DataFrame(rows)
    require(len(table)==9
            and table.original_processed_photo_attempts.sum()<=len(d),
            "Original global water/light 3x3 hue cells are incomplete")
    return table


def wet_dry_changes(cells):
    results=[]
    for solar in ("low","middle","high"):
        a=cells[cells.solar_tertile.eq(solar)].set_index("rainfall_tertile")
        dry,wet=a.loc["low"],a.loc["high"]
        for h in HUES:
            results.append({
                "solar_tertile":solar,"visible_hue":h,
                "dry_photo_fraction":float(dry[f"{h}_among_4class_photos"]),
                "wet_photo_fraction":float(wet[f"{h}_among_4class_photos"]),
                "wet_minus_dry_photo_percentage_points":float(100*(
                    wet[f"{h}_among_4class_photos"]-
                    dry[f"{h}_among_4class_photos"])),
                "not_botanically_confirmed_flower_or_individual":True,
                "different_geographic_taxonomic_solar_cell_compositions":True,
            })
    return pd.DataFrame(results)


def taxon_table(source,fieldnames):
    # A taxon is aggregated only from the exact original photographed rows
    # with ALL predictors present, rather than independently imputing class
    # and climate values from different source photographs.
    x=source[source.classified&source.genus.notna()].dropna(
        subset=list(dict.fromkeys((*parent.GEO,*fieldnames)))).copy()
    x=x[x.inat_taxon_id.astype(str).str.len().gt(0)]
    require(len(x)>12000,"Original photographic class labels missing for all19 source")
    one_hot=[]
    for h in HUES:
        name="photo_"+h
        x[name]=x.morph.eq(h).astype(float)
        one_hot.append(name)
    agg={k:"mean" for k in dict.fromkeys((*parent.GEO,*fieldnames,*one_hot))}
    agg.update({"genus":"first","photo_id":"size"})
    t=x.groupby("inat_taxon_id",sort=True).agg(agg).rename(
        columns={"photo_id":"n_classifiable_source_photos"})
    count=t.groupby("genus").size()
    t=t[t.genus.isin(count[count>=2].index)].copy()
    lon=np.rad2deg(np.arctan2(t.longitude_sin,t.longitude_cos))
    t["spatial_holdout"]=(np.floor((t.latitude+90)/10).astype(int).astype(str)
                          +"_"+np.floor((lon+180)/10).astype(int).astype(str))
    require(t.genus.nunique()>100
            and t.spatial_holdout.nunique()>=FOLDS
            and len(t)>1000
            and np.isclose(t[one_hot].sum(axis=1).to_numpy(float),1).all(),
            "Photo colour fractions or heldout taxon panel invalid")
    return t,one_hot,{"original_classified_photo_rows":len(x),
                       "retained_nominal_taxa":len(t),
                       "eligible_genera":t.genus.nunique(),
                       "geographic_centroid_blocks":t.spatial_holdout.nunique()}


def simplex(y):
    p=np.clip(np.asarray(y,dtype=float),1e-6,1)
    return p/p.sum(axis=1,keepdims=True)


def brier(y,p,weights):
    return float(np.average(((y-p)**2).sum(axis=1),weights=weights))


def predict_fold(s,train,test,ys,blocks):
    yy=ys[train]
    group=s.genus.astype(str).to_numpy()
    train_w=parent.weights_for_genus(group[train])
    test_w=parent.weights_for_genus(group[test])
    train_blocks={};test_blocks={}
    a,b,_=parent.fit_block_transform(s.iloc[train][list(parent.GEO)],
                                      s.iloc[test][list(parent.GEO)],None)
    train_blocks["geography"]=a;test_blocks["geography"]=b
    for name in blocks:
        a,b,_=parent.fit_block_transform(
            s.iloc[train][list(parent.BLOCKS[name])],
            s.iloc[test][list(parent.BLOCKS[name])],parent.BLOCKS[name])
        train_blocks[name]=a;test_blocks[name]=b
    result={}
    for exclude in (None,*blocks,"__geography_only__"):
        keep=([] if exclude=="__geography_only__" else
              [name for name in blocks if name!=exclude])
        tx=np.column_stack([train_blocks["geography"],*[
            train_blocks[name] for name in keep]])
        vx=np.column_stack([test_blocks["geography"],*[
            test_blocks[name] for name in keep]])
        ridge=Ridge(alpha=len(train)*parent.ALPHA_PER_TAXON)
        ridge.fit(tx,yy,sample_weight=train_w)
        prob=simplex(ridge.predict(vx))
        name=("full" if exclude is None else
              "geo_only" if exclude=="__geography_only__" else
              "drop_"+exclude)
        per_hue=(ys[test]-prob)**2
        result[name]={
            "weighted_brier":brier(ys[test],prob,test_w),
            "hue_weighted_brier":[float(np.average(per_hue[:,i],weights=test_w))
                                 for i in range(4)],
            "n_test":int(len(test)),
        }
    return result


def heldout_models(t,one_hot):
    y=t[one_hot].to_numpy(float)
    outputs=[]
    for design in ("heldout_genera","heldout_geographic_10degree_taxon_centroids"):
        groups=t.genus if design=="heldout_genera" else t.spatial_holdout
        records={k:[] for k in ("full","geo_only",
                                *["drop_"+b for b in CLIMATE])}
        for train,test in GroupKFold(n_splits=FOLDS).split(t,y,groups):
            scores=predict_fold(t,train,test,y,CLIMATE)
            for key,val in scores.items():
                records[key].append(val)
        # Test folds hold out disjoint taxa; pool TEST predictions by
        # the true number of held-out taxa, rather than mean five equal folds.
        def pooled(k):
            z=records[k]; weights=np.array([r["n_test"] for r in z],float)
            val=float(np.average([r["weighted_brier"] for r in z],weights=weights))
            hues=np.average(np.array([r["hue_weighted_brier"] for r in z]),
                            weights=weights,axis=0)
            return val,hues
        overall,overall_hues=pooled("full")
        geo,geo_hues=pooled("geo_only")
        for block in ("geo_only",*CLIMATE):
            key=("geo_only" if block=="geo_only" else "drop_"+block)
            minus,minus_hues=pooled(key)
            outputs.append({
                "holdout_design":design,
                "environment_group_removed_or_baseline":block,
                "full_model_brier":overall,
                "compared_model_brier":minus,
                "group_conditional_brier_degradation_when_removed":minus-overall,
                "full_model_improvement_over_geography_baseline":geo-overall,
                "test_taxa":len(t),"test_genera":t.genus.nunique(),
                **{f"{h}_removed_group_brier_increase":
                       float(minus_hues[i]-overall_hues[i])
                   for i,h in enumerate(HUES)},
                "train_only_PCA_standardisation":True,
                "photographic_label_not_biochemical_pigment_measurement":True,
                "no_causal_group_variable_selection_claim":True,
            })
    return pd.DataFrame(outputs)


def run(release,extended,out):
    require(not out.exists() or not any(out.iterdir()),
            "Original four-hue environmental study output already exists")
    p,photo,abiotic=parent.read_source(release,extended)
    cells,cut=classify_climate(p)
    world=four_hue_source_cells(cells)
    delta=wet_dry_changes(world)
    t,names,support=taxon_table(p,tuple(
        feature for block in CLIMATE for feature in parent.BLOCKS[block]))
    cv=heldout_models(t,names)
    require(len(cv)==8 and set(world.rainfall_tertile)==
            {"low","middle","high"},
            "Genuine original class composition and climate heldout failed")
    out.mkdir(parents=True,exist_ok=True)
    frames={OUT_COMP:world,OUT_DELTA:delta,OUT_CV:cv}
    for key,frame in frames.items():
        frame.to_csv(out/key,index=False,lineterminator="\n")
    receipt={
        "schema":"mois_original313535_four_visible_hues_water_solar_heldout_composition_v1",
        "real_original_photo_attempts":parent.prev.N_PHOTOS,
        "original_classified_hues":parent.prev.N_HUES,
        "original_white_hue_ROIs":parent.prev.WHITE_N,
        "four_original_auto_ROI_visible_hue_classes":list(HUES),
        "source_photo_RGB_sha256":photo["full_16_workers_original_pixel_source_SHA256"],
        "original_abiotic_parquet_sha256":abiotic[
            "new_environment_parquet_sha256"],
        "covariate_bins_formed_before_any_photo_hue_filter":True,
        "original_photo_rain_solar_tertile_boundaries":cut,
        "heldout_taxon_photo_source":support,
        "four_hue_photo_labels_NOT_anthocyanin_or_flowering_stage_verified":True,
        "pixel_missingness_explicit_not_nonwhite":True,
        "all_climate_normals_original_1970_2000_not_observation_year_month":True,
        "future_SI_month_by_month_separate_from_main_global_analysis":True,
        "no_adaptive_global_50_50_equilibrium_or_pigment_cost_fitness":True,
        "CSV_source_SHA256":{key:sha(out/key) for key in frames},
    }
    (out/OUT_RECEIPT).write_text(json.dumps(receipt,sort_keys=True,
        indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("MOIS_FOUR_COLOURS_RAIN_SOLAR_CLIMATE_HELDOUT",json.dumps({
        "support":support,
        "wet_dry_hue_changes_pp":delta.to_dict("records"),
        "climate_group_brier":cv.to_dict("records")
    },sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--original-photo-release",type=Path,required=True)
    ap.add_argument("--full19-abiotic-release",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args()
    run(a.original_photo_release,a.full19_abiotic_release,a.output)
