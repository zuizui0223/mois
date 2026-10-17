#!/usr/bin/env python3
"""MOIS expanded-photo global flower-hue ecology, not FCP intraspecific D/W.

Estimate genus-centred, genus-equal photographed four-hue fractions with
geography, BIO climate, landscape topography, and original soil as distinct
covariate blocks on matched source-photo supports. Exploratory observational
cluster CIs / genus-heldout diagnostics; NOT a causal result or a frozen maxT
confirmation. Monthly seasonal source photos are described, not interpreted
as within-individual colour switches or observed-year weather.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

HUES=("white","yellow_orange","red_pink","blue_purple")
GEOG=("latitude","abs_lat","longitude_sin","longitude_cos",
      "wc2_elev_2p5m","new_source_fraction")
CLIM=("wc2_bio01_10m","wc2_bio05_10m","wc2_bio12_10m",
      "wc2_bio15_10m","wc2_annual_srad_normal_10m_Wm2")
TERRAIN=("terrain_slope_2p5m_degree","terrain_TRI_2p5m_m")
SOIL=("historical_original_soil_pH","historical_original_soil_SOC",
      "historical_original_soil_N","historical_original_soil_clay",
      "historical_original_soil_available_water_proxy")
SOURCE_EXPECTED_PHOTOS=248312
SOURCE_EXPECTED_CLASSIFIED=120807
MIN_GENUS_SPECIES=2
FOLDS=5


def require(v:bool,msg:str):
    if not v: raise ValueError(msg)


def digest(path:Path):
    sha=hashlib.sha256()
    with path.open("rb") as f:
        for buf in iter(lambda:f.read(1<<20),b""):
            sha.update(buf)
    return sha.hexdigest()


def genus(name:str)->str|None:
    first=str(name).split()[0] if isinstance(name,str) and name.strip() else ""
    return first if len(first)>=2 and first[0].isupper() and first.isalpha() else None


def make_photo_features(data:pd.DataFrame)->pd.DataFrame:
    f=data.copy()
    f=f.loc[f["morph"].isin(HUES)].copy()
    f["genus"]=f["species"].map(genus)
    f=f.loc[f.genus.notna()].copy()
    f["abs_lat"]=pd.to_numeric(f.latitude,errors="coerce").abs()
    f["longitude_sin"]=np.sin(np.deg2rad(pd.to_numeric(f.longitude,errors="coerce")))
    f["longitude_cos"]=np.cos(np.deg2rad(pd.to_numeric(f.longitude,errors="coerce")))
    f["new_source_fraction"]=(f.original_source!="historical").astype(float)
    for colour in HUES:
        f["y_"+colour]=(f.morph==colour).astype(float)
    for x in (*GEOG,*CLIM,*TERRAIN,*SOIL):
        if x in f:
            f[x]=pd.to_numeric(f[x],errors="coerce")
    return f


def aggregate_species(frame:pd.DataFrame,required:tuple[str,...],
                      *,min_photos:int=1)->tuple[pd.DataFrame,dict]:
    # Crucial: PHOTO rows are selected on the entire block first, so the
    # additional soil/terrain layer and baseline have exactly the same
    # original source photos and within-species spatial supports.
    sub=frame.dropna(subset=list(required)).copy()
    measurements={k:"mean" for k in set((*GEOG,*CLIM,*TERRAIN,*SOIL)) if k in sub}
    measurements.update({"y_"+k:"mean" for k in HUES})
    measurements.update({"genus":"first","morph":"size"})
    g=sub.groupby("inat_taxon_id",sort=True).agg(measurements)
    g=g.rename(columns={"morph":"original_classified_photo_count"})
    g=g.loc[g.original_classified_photo_count>=min_photos]
    counts=g.groupby("genus").size()
    g=g.loc[g.genus.isin(counts[counts>=MIN_GENUS_SPECIES].index)].copy()
    require(len(g)>200 and g.genus.nunique()>50,
            "Too few independently qualifying congeneric species")
    return g,{
        "original_photo_rows_with_all_block_covariates":int(len(sub)),
        "taxa_with_block_photo_support":int(sub.inat_taxon_id.nunique()),
        "taxa_after_min_photo_and_congeneric_gates":int(len(g)),
        "qualified_nominal_genera":int(g.genus.nunique()),
        "min_classified_photos_per_species":min_photos,
        "median_source_photo_count_per_eligible_taxon":
            float(g.original_classified_photo_count.median()),
    }


def design(s:pd.DataFrame,features:tuple[str,...],hue:str)->tuple[np.ndarray,np.ndarray,np.ndarray,list[str]]:
    X=s[list(features)].to_numpy(float)
    mu=X.mean(axis=0)
    sd=X.std(axis=0,ddof=0)
    good=np.isfinite(sd)&(sd>1e-10)
    X=(X[:,good]-mu[good])/sd[good]
    names=[x for x,ok in zip(features,good) if ok]
    require(bool(names),"All covariates are constant")
    for i in range(X.shape[1]):
        X[:,i]-=pd.Series(X[:,i],index=s.index).groupby(s.genus).transform("mean").to_numpy()
    usable=np.isfinite(X).all(axis=0)&(np.std(X,axis=0)>1e-10)
    X=X[:,usable]
    names=[k for k,ok in zip(names,usable) if ok]
    require(len(names)>0, "No among-species contrasts within the qualifying genera")
    y=s["y_"+hue].to_numpy(float)
    y-=pd.Series(y,index=s.index).groupby(s.genus).transform("mean").to_numpy()
    counts=s.groupby("genus").genus.transform("size").to_numpy(float)
    weights=1/counts
    require(np.all(np.isfinite(X)) and np.isfinite(y).all() and
            np.all(np.isfinite(weights)),"Invalid genus-centred outcome/covariate")
    return X,y,weights,names


def genus_cv_r2(s:pd.DataFrame,features:tuple[str,...],hue:str)->float:
    X,y,w,names=design(s,features,hue)
    unique=s.genus.astype(str)
    folds=unique.map(lambda z:int(hashlib.sha256(z.encode()).hexdigest()[:8],16)%FOLDS).to_numpy()
    total_sse=total_sst=0.
    for k in range(FOLDS):
        train=folds!=k
        test=~train
        if train.sum()<50 or test.sum()<20:continue
        coeff=np.linalg.lstsq(X[train]*np.sqrt(w[train])[:,None],
                              y[train]*np.sqrt(w[train]),rcond=1e-9)[0]
        error=y[test]-X[test]@coeff
        total_sse+=float(np.sum(w[test]*error**2))
        total_sst+=float(np.sum(w[test]*y[test]**2))
    return float(1-total_sse/total_sst) if total_sst>0 else float("nan")


def fitted(s:pd.DataFrame,features:tuple[str,...],hue:str):
    X,y,w,names=design(s,features,hue)
    model=sm.WLS(y,X,weights=w).fit(
        cov_type="cluster",cov_kwds={"groups":s.genus.to_numpy()},use_t=False)
    require(np.isfinite(model.params).all() and np.isfinite(model.pvalues).all(),
            "Cluster-robust genus model ill conditioned")
    r2=1-float(np.sum(w*(y-model.fittedvalues)**2)/np.sum(w*y*y))
    z=model.conf_int(alpha=.05)
    effects=[]
    for i,name in enumerate(names):
        effects.append({
          "response_four_hue":hue,"variable":name,
          "coefficient_white_or_other_hue_percentage_points_per_species_source_SD":
              round(float(100*model.params[i]),7),
          "cluster_genus_95pct_CI_low_pp":round(float(100*z[i,0]),7),
          "cluster_genus_95pct_CI_high_pp":round(float(100*z[i,1]),7),
          "nominal_two_sided_cluster_p":float(model.pvalues[i]),
          "BH_q_in_block_hue_family":None,
          "not_a_joint_maxT_confirmatory_interval":True
        })
    return effects,{
        "in_sample_genus_centered_R2":float(r2),
        "heldout_whole_genera_fivefold_R2":genus_cv_r2(s,features,hue),
        "actual_predictor_terms":names,
    }


def block(frame:pd.DataFrame,blockname:str,extra:tuple[str,...],
          *,minimum_photos:int=1)->tuple[list[dict],list[dict],dict]:
    required=(*GEOG,*CLIM,*extra)
    # The same original PHOTO rows qualify for all submodels in this block.
    s,support=aggregate_species(frame,required,min_photos=minimum_photos)
    baseline=GEOG
    final=tuple(dict.fromkeys((*GEOG,*CLIM,*extra)))
    before=[]
    effects=[]
    for h in HUES:
        _,geo=fitted(s,baseline,h)
        _,clim=fitted(s,(*GEOG,*CLIM),h)
        entries,full=fitted(s,final,h)
        for e in entries:
            e["analysis_block"]=blockname
            e["minimum_photos_per_species"]=minimum_photos
            effects.append(e)
        before.append({
            "block":blockname,"hue":h,"min_photos_per_species":minimum_photos,
            "original_photo_support_same_for_each_nested_comparison":True,
            "n_species":int(len(s)),"n_genera":int(s.genus.nunique()),
            "geo_only_genus_centered_R2":geo["in_sample_genus_centered_R2"],
            "geo_plus_climate_genus_centered_R2":clim["in_sample_genus_centered_R2"],
            "geo_climate_extra_genus_centered_R2":full["in_sample_genus_centered_R2"],
            "delta_R2_climate_on_matched_photo_support":clim["in_sample_genus_centered_R2"]-geo["in_sample_genus_centered_R2"],
            "delta_R2_extra_on_matched_photo_support":full["in_sample_genus_centered_R2"]-clim["in_sample_genus_centered_R2"],
            "geo_heldout_genera_fivefold_R2":geo["heldout_whole_genera_fivefold_R2"],
            "geo_climate_heldout_genera_fivefold_R2":clim["heldout_whole_genera_fivefold_R2"],
            "full_heldout_genera_fivefold_R2":full["heldout_whole_genera_fivefold_R2"],
            "delta_heldout_genus_R2_extra":full["heldout_whole_genera_fivefold_R2"]-clim["heldout_whole_genera_fivefold_R2"]
        })
    # Exploratory BH false-discovery adjustment only WITHIN this block's
    # 4 hues × (climate + new soil/terrain variables) family.
    test_names=set((*CLIM,*extra))
    indices=[i for i,e in enumerate(effects) if e["variable"] in test_names]
    if indices:
        p=np.array([effects[i]["nominal_two_sided_cluster_p"] for i in indices])
        q=multipletests(p,alpha=.05,method="fdr_bh")[1]
        for idx,val in zip(indices,q):
            effects[idx]["BH_q_in_block_hue_family"]=float(val)
    support["model_block"]=blockname
    support["response_hue_levels"]=list(HUES)
    support["genus_weight"]="one species equal; each qualifying genus total 1"
    support["source_photo_filter_same_for_nested_geography_climate_extra_models"]=True
    return effects,before,support


def calendar_descriptives(frame:pd.DataFrame)->list[dict]:
    o=frame.loc[frame.observed_on.notna()].copy()
    o["yearmonth"]=pd.to_datetime(o.observed_on,errors="coerce").dt.to_period("M").astype(str)
    o=o.loc[~o.yearmonth.isin(("NaT",""))]
    o["hemisphere"]=np.select([
        o.latitude>=23.5,o.latitude<=-23.5],
        ["northern_temperate","southern_temperate"],default="tropical")
    out=[]
    for (source,month,hemi),p in o.groupby(
        ["original_source","yearmonth","hemisphere"],observed=True):
        if len(p)<3:continue
        source_frac={h:float((p.morph==h).mean()) for h in HUES}
        bytax=p.groupby("inat_taxon_id")
        equal={h:float(bytax["y_"+h].mean().mean()) for h in HUES}
        out.append({
            "original_source":source,"observed_year_month":month,
            "hemisphere":hemi,"classified_source_photos":int(len(p)),
            "taxa_classifiable":int(p.inat_taxon_id.nunique()),
            "genera_classifiable":int(p.genus.nunique()),
            **{"photo_fraction_"+h:round(v,6) for h,v in source_frac.items()},
            **{"species_equal_fraction_"+h:round(v,6) for h,v in equal.items()},
            "weather_is_monthly_normal_not_actual_year":True,
            "no_within_plant_or_FCP_polymorphism_inference":True,
        })
    return out


def main(source:Path,receipt_path:Path,outdir:Path):
    require(not outdir.exists() or not any(outdir.iterdir()),
            "Refuse overwriting a frozen genuine photo/environment analysis")
    rec=json.loads(receipt_path.read_text())
    require(digest(source)==rec["full_photo_static_monthly_normal_csv_sha256"],
            "Source 248k photo-worldclim-environment rows do not match actual GIS receipt")
    require(rec["n_original_photo_ids"]==SOURCE_EXPECTED_PHOTOS
            and rec["n_four_hue_classified"]==SOURCE_EXPECTED_CLASSIFIED,
            "Source actual 248k, 120807 measured photographic hue counts changed")
    cols=["photo_id","inat_taxon_id","species","morph","original_source",
          "observed_on","latitude","longitude","wc2_elev_2p5m",*CLIM,*TERRAIN,*SOIL]
    cols=list(dict.fromkeys(cols))
    raw=pd.read_csv(source,compression="gzip",usecols=cols,
                    low_memory=False,dtype={"photo_id":str,"inat_taxon_id":str})
    require(len(raw)==SOURCE_EXPECTED_PHOTOS and raw.photo_id.nunique()==len(raw),
            "Original photo overlap/dedup violated")
    frame=make_photo_features(raw)
    original_classified=int(raw.morph.isin(HUES).sum())
    require(original_classified==SOURCE_EXPECTED_CLASSIFIED,
            "Original four-state photographed hue source labels changed")
    require(0<len(frame)<=original_classified and
            set(frame.morph.unique()).issubset(HUES),
            "Classifier failures or invalid nominal genus were wrongly admitted")
    tables=[]
    diagnostics=[]
    supports=[]
    for name,extra,threshold in (
        ("climate_over_geography",(),1),
        ("terrain_addition_matching_photo_support",TERRAIN,1),
        ("soil_addition_matching_photo_support",SOIL,1),
        ("climate_two_photos_per_species_sensitivity",(),2),
    ):
        effects,comparisons,support=block(frame,name,extra,minimum_photos=threshold)
        tables.extend(effects); diagnostics.extend(comparisons); supports.append(support)
    # Original historical source-only sensitivity; same photo sample and genus
    # cohort differs from combined sample and is NEVER called replication.
    hist=frame.loc[frame.original_source.eq("historical")].copy()
    eff,comp,sup=block(hist,"historical_only_climate_source_sensitivity",(),minimum_photos=1)
    tables.extend(eff);diagnostics.extend(comp);supports.append(sup)
    periods=calendar_descriptives(frame)
    outdir.mkdir(parents=True,exist_ok=True)
    exports={}
    for name,records in (
        ("mois_global_248k_four_hue_genus_cluster_coefficients.csv",tables),
        ("mois_global_248k_matched_photo_support_block_comparisons.csv",diagnostics),
        ("mois_global_248k_month_hemisphere_source_photo_hues.csv",periods)
    ):
        path=outdir/name
        pd.DataFrame(records).to_csv(path,index=False,lineterminator="\n")
        exports[name]=digest(path)
    result={
        "schema":"mois_frozen_248312_original_photo_colour_global_species_genus_static_abiotic_blocks_v1",
        "status":"EXECUTED_OBSERVATIONAL_GENUS_EQUAL_ABIOTIC_BLOCKS_NOT_FCP_POLYMORPHISM",
        "original_photo_source_rows":len(raw),
        "original_photo_source_four_hue_classified":original_classified,
        "classified_rows_with_valid_nominal_genus_in_models":len(frame),
        "source_four_hue_classified_excluded_nonstandard_genus":original_classified-len(frame),
        "unique_photo_identifiers_verified":True,
        "source_WorldClim_SHAs_verbatim_from_previous_site_assembly":rec["worldclim_source_zip_sha256"],
        "original_site_bio_colours_from_photo_source_different_from_genetic_morphs":True,
        "original_source_environment_join_sha256":digest(source),
        "named_species_aggregation_rule":"species equals one taxon ID, photo hue means among classified supported originals",
        "nominal_genus":"first word from source taxon name, genus≥2 named taxa",
        "nuisance_location_features":list(GEOG),
        "long_term_climate_features":list(CLIM),
        "terrain_2p5_arcminute_grain_only":list(TERRAIN),
        "soil_features_original_historical_source_only":list(SOIL),
        "support_by_block":supports,
        "multiple_comparisons":"Exploratory Benjamini-Hochberg q across four hue categories × climate+additional block variables. Nominal cluster CIs are not simultaneous maxT confirmation.",
        "cluster_inference":"Nominal genus-robust sandwich 95% CI; observational linear-probability species photo fractions",
        "heldout_validation":"Deterministic five-fold whole-genus-heldout deviations relative to observed within-genus mean reference; no causal inference",
        "source_sampled_calendar_months_desc_only":len(periods),
        "photo_time_of_day_verified":False,
        "actual_observed_year_month_weather_fitted":False,
        "static_longterm_and_seasonal_predictors_separated":True,
        "FCP_within_species_D_W_or_local_50km_statistic_calculated":False,
        "files_sha256":exports,
        "analysis_not_formal_original_MOIS_maxT_protocol_confirmation":True,
    }
    (outdir/"mois_248k_global_abiotic_block_model_receipt.json").write_text(
        json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("MOIS_248K_GLOBAL_GENUS_EQUAL_ABIOTIC_MODELS",json.dumps({
        "supports":supports,"monthly_source_groups":len(periods)},sort_keys=True))
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--original-photo-site-worldclim",type=Path,required=True)
    parser.add_argument("--original-site-source-receipt",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    a=parser.parse_args()
    main(a.original_photo_site_worldclim,a.original_site_source_receipt,a.output)
