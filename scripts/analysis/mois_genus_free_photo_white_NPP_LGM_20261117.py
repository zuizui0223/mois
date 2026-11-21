#!/usr/bin/env python3
"""MOIS primary floral-colour analysis WITHOUT ANY GENUS adjustment.

Unit is one accepted species (37,104 source accepted species) with ONE
species-modal photographic 9-colour label. Source independent photo metadata
are used as quality/sampling covariates, not to infer true species flower petals.
GBIF supplies current range, MOD17 23-yr NPP supplies annual carbon flux,
WorldClim v1.4 supplies current and LGM paleoclimate, ESA woody AGB/NDVI
supplies current habitat metrics.

NO genus variable, genus FE, genus clusters, genus weights or genus prior in
estimation. A separate genuine GBIF 100-km regional/spatial-block model
without genus is already available in the signed primary source artifact.
The species HC3 CI DOES NOT substitute for phylogenetic or spatial
independence, and positive slopes are not climate adaptation or refugia.
"""
from __future__ import annotations

import argparse,hashlib,io,json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.special import expit
from statsmodels.stats.multitest import multipletests

from scripts.analysis.mois_lgm_climate_history_species_20261116 import HUES,require,sha
from scripts.analysis.mois_flower_colour_LGM_NPP_GBIF_20261117 import (
    NPP_NAME,NPP_COVER,LGM_TEMP,LGM_RAIN,CTRL
)

SOURCE_MODEL_SHA="1da364ff63c506f9e2c06d1d71c9609a2c48fec7a428a21ba2c3c4c6ff8f4800"
SOURCE_MODEL_MEMBER="GBIF37104_photo9_species_NPP23years_LGM_AGB_NDVI.parquet"
SOURCE_TRAITS_SHA="abcd3aee766c85a92b86a8f93f8ed696bbe49d0f6cf2f64e5f75b457abf78e6a"
SOURCE_TRAITS_MEMBER="current_exact_accepted_COLXR_numeric9_traits.csv.gz"
SOURCE_TRAITS_INNER_SHA="b06d8980860ab11d1e9cff1af7dd1d92b3bfa9de4f2d227bc7e58809f931663a"
CONFIDENCE="sampling_consensus_single"

def read_source(main_archive, traits_archive):
    import pyarrow.parquet as pq
    require(sha(main_archive)==SOURCE_MODEL_SHA,"source species photo NPP LGM ZIP SHA differs")
    with ZipFile(main_archive) as z:
        source=z.read(SOURCE_MODEL_MEMBER)
        receipt=json.loads(z.read("receipt.json"))
    require(sha(source)==receipt["output_SHA256"][SOURCE_MODEL_MEMBER]
            and receipt["source_37k_species"]==37104
            and receipt["source_5493931_GBIF_species_cells"]==5493931,
            "signed accepted species traits are not intact")
    species=pq.read_table(io.BytesIO(source)).to_pandas()
    species["key"]=species.key.astype(str)
    require(len(species)==37104 and species.key.is_unique
            and species.source_hue9.isin(HUES).all(),
            "bad source 37104 species labels")

    require(sha(traits_archive)==SOURCE_TRAITS_SHA,
            "signed original nine-colour trait ZIP changed")
    with ZipFile(traits_archive) as z:
        b=z.read(SOURCE_TRAITS_MEMBER)
    require(sha(b)==SOURCE_TRAITS_INNER_SHA,
            "photo mean9 original numeric data SHA differs")
    traits=pd.read_csv(io.BytesIO(b),compression="gzip",
                       dtype={"gbif_colxr_accepted_species_key":str},
                       low_memory=False)
    names=["observed_mean_"+x for x in HUES]
    require(len(traits)==37104
            and traits.gbif_colxr_accepted_species_key.is_unique
            and set(names).issubset(traits.columns),"photo mean9 source drift")
    m=traits[names].to_numpy(dtype=float)
    require(np.isfinite(m).all() and (m>=0).all()
            and np.allclose(m.sum(axis=1),1,atol=5e-5),
            "source mean nine-class vectors are invalid")
    traits["mean9_hue"]=np.asarray(HUES)[m.argmax(axis=1)]
    traits["mean9_white_frac"]=traits["observed_mean_white"].astype(float)
    traits=traits[["gbif_colxr_accepted_species_key","mean9_hue",
                   "mean9_white_frac"]]
    species=species.merge(traits,left_on="key",
                          right_on="gbif_colxr_accepted_species_key",
                          how="left",validate="one_to_one")
    require(len(species)==37104 and species.mean9_hue.notna().all(),
            "cannot map mean photo 9 hues to GBIF accepted species")
    species["photo_modal_vs_mean9_agree"]=(species.source_hue9
                                           ==species.mean9_hue)
    require((~species.photo_modal_vs_mean9_agree).sum()==2466,
            "original 2466 mean versus modal nine hue mismatches changed")
    return species

def sample_frame(s,min_npp=.75,min_lgm=.75):
    numerical=[*CTRL,NPP_NAME,NPP_COVER,LGM_TEMP,LGM_RAIN,
               "lgm_climate_coverage","occupied_grid_cells"]
    good=(s[numerical].notna().all(axis=1)
          &(s[NPP_COVER]>=min_npp)
          &(s.lgm_climate_coverage>=min_lgm)
          &(s.BIO12>=0)
          &(s.ESA_AGB_woody_Mg_ha_2020>=0)
          &(s.records_passed_spatial_QA>=0)
          &s.hue_status.notna())
    d=s.loc[good].copy()
    require(d.key.is_unique,"duplicate species, violates one-species unit")
    return d

CONT=("NPP","LGM_warming","LGM_rain","AGB","NDVI",
      "current_temp","current_log_precip","latitude_abs",
      "log_GBIF_effort","log_GBIF_occupied_cells")
FOCUS=("NPP","LGM_warming","LGM_rain","AGB","NDVI")

def design(d,status_controls=True):
    """Standardized full species-level design, no taxon rank used."""
    v={
      "NPP":d[NPP_NAME].astype(float),
      "LGM_warming":d[LGM_TEMP].astype(float),
      "LGM_rain":d[LGM_RAIN].astype(float),
      "AGB":np.log1p(d.ESA_AGB_woody_Mg_ha_2020.astype(float)),
      "NDVI":d.MODIS_NDVI_2001_2021_mean.astype(float),
      "current_temp":d.BIO01.astype(float),
      "current_log_precip":np.log1p(d.BIO12.astype(float)),
      "latitude_abs":d.mean_abs_latitude_over_occupied_cells.astype(float),
      "log_GBIF_effort":np.log1p(d.records_passed_spatial_QA.astype(float)),
      "log_GBIF_occupied_cells":np.log1p(d.occupied_grid_cells.astype(float))
    }
    x=pd.DataFrame(v,index=d.index)
    scale=x.std()
    require((scale>1e-9).all() and np.isfinite(scale).all()
            and np.isfinite(x.to_numpy()).all(),
            "no variance or invalid environmental exposure")
    x=(x-x.mean())/scale
    if status_controls:
        dummy=pd.get_dummies(d.hue_status.astype(str),prefix="photo_status",
                             drop_first=True,dtype=float)
        x=pd.concat([x,dummy],axis=1)
    x=sm.add_constant(x,has_constant="add")
    require(x.columns.is_unique and not any("genus" in c.lower() for c in x),
            "no genus information should enter genus-free model")
    return x,{k:float(v) for k,v in scale.items()}

def estimate(d,color="white",label="modal",adjustment="environment",
             logistic=False):
    label_col="source_hue9" if label=="modal" else "mean9_hue"
    require(label_col in d and color in HUES,
            "requested colour label not in source traits")
    y=d[label_col].eq(color).astype(float)
    if adjustment=="unadjusted":
        z=d[NPP_NAME].astype(float)
        sd=float(z.std())
        require(sd>0,"NPP ecological exposure does not vary")
        x=sm.add_constant(pd.DataFrame({"NPP":(z-z.mean())/sd},index=d.index),
                         has_constant="add")
    elif adjustment=="environment":
        x,all_sd=design(d,status_controls=True)
    elif adjustment=="no_photo_status":
        x,all_sd=design(d,status_controls=False)
    else:
        raise ValueError("unknown genus-free model")
    if len(d)<120 or y.nunique()<2:
        return {"status":"insufficient_white_nonwhite_species",
                "n_species":len(d)}
    if logistic:
        model=sm.GLM(y,x,family=sm.families.Binomial()).fit(
            cov_type="HC3",maxiter=110)
        effect={}
        for f in FOCUS:
            if f not in x:continue
            lo,hi=model.conf_int().loc[f].tolist()
            effect[f]={"odds_ratio_per_species_SD":float(np.exp(model.params[f])),
                       "odds_ratio_95CI":[float(np.exp(lo)),float(np.exp(hi))],
                       "p":float(model.pvalues[f])}
        return {"status":"fitted","n_species":len(d),
                "model":"genus_free_logistic_binomial_HC3",
                "n_white":int(y.sum()),"white_fraction":float(y.mean()),
                "effects":effect}
    model=sm.OLS(y,x).fit(cov_type="HC3",use_t=False)
    effect={}
    for f in FOCUS:
        if f not in x:continue
        lo,hi=model.conf_int().loc[f].tolist()
        effect[f]={"beta_white_probability_per_species_SD":float(model.params[f]),
                   "95CI":[float(lo),float(hi)],
                   "p":float(model.pvalues[f])}
    return {"status":"fitted","n_species":len(d),
            "n_white":int(y.sum()),"n_nonwhite":int(len(y)-y.sum()),
            "white_fraction":float(y.mean()),
            "model":"genus_free_one_accepted_species_one_vote_OLS_HC3",
            "n_model_columns":len(x.columns),
            "effects":effect,
            "r_squared_in_sample":float(model.rsquared),
            "spatial_or_phylogenetic_independence_established":False}

def run(main,traits,output):
    out=Path(output)
    require(not out.exists(),"do not overwrite genus-free original run")
    species=read_source(main,traits)
    d=sample_frame(species)
    require(len(d)>30000,"photo species coverage changed")
    matches=d[d.photo_modal_vs_mean9_agree].copy()
    high=d[d.hue_status.eq(CONFIDENCE)].copy()
    low=d[d.mean_abs_latitude_over_occupied_cells<30].copy()
    upper=d[d.mean_abs_latitude_over_occupied_cells>=30].copy()
    cohorts={
      "all_source_modal":d,
      "all_mean9_top":d,
      "modal_mean9_agree_only":matches,
      "consensus_3plus_independent_observations":high,
      "absolute_latitude_under30":low,
      "absolute_latitude_30plus":upper,
      "NPP_100percent_recorded_cells":sample_frame(species,min_npp=1.0),
    }
    report={"source_GBIF_accepted_species":37104,
        "source_GBIF_species_x_cells":5493931,
        "n_complete_main":len(d),
        "n_source_modal_mean9_disagreements":2466,
        "photo_hue_is_true_independent_petal_annotation":False,
        "all_models_use_no_genus_FEs_clusters_weights":True,
        "species_coverage_nonrandom_and_geospatial_correlated":True,
        "NPP_2001_2023_is_not_floral_biomass":True,
        "LGM_climate_at_present_range_is_not_refugia":True,
        "base":{},"nine_hue_exploratory":{}}
    for name,sub in cohorts.items():
        which="mean9" if name=="all_mean9_top" else "modal"
        report["base"][name]={
          "n_species":len(sub),
          "unadjusted_NPP":estimate(sub,label=which,adjustment="unadjusted"),
          "full_10ecology_photo_status":estimate(sub,label=which),
          "full_10ecology_no_photo_status":estimate(
                sub,label=which,adjustment="no_photo_status")
        }
        if name in ("all_source_modal","all_mean9_top"):
            report["base"][name]["binomial_logistic_HC3"]=estimate(
                sub,label=which,logistic=True)
        fit=report["base"][name]["full_10ecology_photo_status"]
        print("MOIS_GENUS_FREE_WHITE_FLOWER",json.dumps({
            "cohort":name,"n_species":len(sub),
            "white":fit.get("n_white"),
            "NPP_beta":fit.get("effects",{}).get("NPP",{}).get(
                "beta_white_probability_per_species_SD"),
            "NPP_p":fit.get("effects",{}).get("NPP",{}).get("p"),
            "LGM_beta":fit.get("effects",{}).get("LGM_warming",{}).get(
                "beta_white_probability_per_species_SD"),
            "LGM_p":fit.get("effects",{}).get("LGM_warming",{}).get("p")
        },sort_keys=True),flush=True)
    # One-vs-rest exploratory nine colour model, same species, same X.
    for col in HUES:
        report["nine_hue_exploratory"][col]=estimate(d,color=col,label="modal")
    for feature in ("NPP","LGM_warming","AGB"):
        keys=[h for h,v in report["nine_hue_exploratory"].items()
              if v.get("status")=="fitted" and feature in v["effects"]]
        q=multipletests([
            report["nine_hue_exploratory"][h]["effects"][feature]["p"]
            for h in keys],method="fdr_bh")[1]
        for h,qq in zip(keys,q):
            report["nine_hue_exploratory"][h].setdefault("q_BH_nine",{})[
                feature]=float(qq)
    out.mkdir(parents=True)
    brief=out/"GBIF_photo_white_genus_free_environment_models.json"
    brief.write_text(json.dumps(report,indent=2,sort_keys=True,default=str)+"\n")
    tax=out/"GBIF37104_photo9_mean9_modal_NPP_LGM_species_quality.parquet"
    keep=["key","species","source_hue9","mean9_hue","mean9_white_frac",
          "photo_modal_vs_mean9_agree","hue_status","BIO01","BIO12",
          "MODIS_NDVI_2001_2021_mean","ESA_AGB_woody_Mg_ha_2020",
          "mean_abs_latitude_over_occupied_cells",NPP_NAME,NPP_COVER,
          LGM_TEMP,LGM_RAIN,"records_passed_spatial_QA",
          "occupied_grid_cells","lgm_climate_coverage"]
    species[keep].to_parquet(tax,index=False)
    receipt={"source_main_sha256":SOURCE_MODEL_SHA,
       "source_mean9_trait_sha256":SOURCE_TRAITS_SHA,
       "source_37104_accepted_species":37104,
       "source_GBIF_species_cells":5493931,
       "response_is_accepted_species_photo9_hue":True,
       "original_photo_modal_vs_mean9_discrepancies":2466,
       "no_genus_in_models":True,
       "HC3_does_not_correct_shared_phylogeny_or_GBIF_geography":True,
       "results":report,
       "output_sha256":{f.name:sha(f) for f in (brief,tax)}}
    (out/"receipt.json").write_text(json.dumps(receipt,indent=2,sort_keys=True,
                                              default=str)+"\n")
    return receipt

if __name__=="__main__":
    p=argparse.ArgumentParser()
    for k in ("main","traits","output"):
        p.add_argument("--"+k,required=True)
    a=p.parse_args()
    run(a.main,a.traits,a.output)
