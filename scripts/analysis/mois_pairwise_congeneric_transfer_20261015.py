#!/usr/bin/env python3
"""Out-of-genus transfer of BETWEEN-SPECIES PHOTO-COLOUR DIFFERENCES.

A pair of previously unseen congeners has environmental predictor values;
predict their *difference* in classified-photo white fractions. Prediction
needs neither the test-genus mean white phenotype nor any test labels.
Evaluating pairwise differences naturally uses labels AFTER prediction.
This is an original-source photo trait, NOT biological fitness or genetics.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np

PHOTO_SHA = "2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0"
ENV_SHA = "6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca"
BLOCKS = {
    "geography":["abs_lat","lon_sin","lon_cos"],
    "moisture":["wc_bio12","wc_vapr_annual_mean"],
    "temp_radiation":["wc_bio5","wc_srad_annual_mean"],
    "elevation":["wc_elevation_m"],
    "soil5":["soil_pH","soil_SOC","soil_N","soil_clay","soil_available_water_proxy"],
}
MODELS = {
    "null":[],
    "geography":["geography"],
    "geography_moisture":["geography","moisture"],
    "geography_fullclimate":["geography","moisture","temp_radiation","elevation"],
    "geography_fullclimate_soil5":["geography","moisture","temp_radiation","elevation","soil5"],
    "geography_soil5":["geography","soil5"],
}
SALTS = (None,202610151,202610152,202610153)
BOOT = 1999
BOOT_SEED = 202610152041
ALPHA = .01

def _loader():
    path = Path(__file__).with_name("mois_full_biogeography_climate_soil_audit_20261015.py")
    spec = importlib.util.spec_from_file_location("mois_original_environment_owner",path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m

def folds(genera,salt,original_fold):
    if salt is None:
        return np.array([original_fold(g) for g in genera],np.int8)
    return np.array([
        int.from_bytes(hashlib.sha256((str(salt)+"|"+str(g)).encode()).digest()[:8],"big")%5
        for g in genera],np.int8)

def direct_pairwise_squared_error(truth,pred):
    """Actual unordered pairs; test fixture only, O(n^2)."""
    n=len(truth)
    if n<2:raise ValueError("Need two sampled congeners")
    return float(np.mean([
        ((truth[i]-truth[j])-(pred[i]-pred[j]))**2
        for i in range(n) for j in range(i+1,n)
    ]))

def pairwise_group_squared_error(truth,pred):
    """Exact O(n) identity; requires NO test-genus target mean to predict.

    All unordered pairs of species within a test genus have equal weight.
    This agrees with a prospective pair-of-taxa environmental contrast task.
    """
    n=len(truth)
    if n<2:raise ValueError("Need two sampled congeners")
    delta=np.asarray(truth,float)-np.asarray(pred,float)
    return float(2*n/(n-1)*np.mean((delta-delta.mean())**2))

def source_species(photo_zip,env_zip):
    loader=_loader()
    if (loader.PHOTO_SHA,loader.ENV_SHA)!=(PHOTO_SHA,ENV_SHA):
        raise ValueError("Original FCP source SHA contract changed")
    _,soil,_,species,_=loader.original(photo_zip,env_zip)
    if (len(soil),len(species),species.genus.nunique(),int(species.n_photos.sum()))!=(31925,16801,1665,26163):
        raise ValueError("Original photo/source genus complete-case denominator changed")
    return loader,species

def evaluate(species,fold_of_genera,nboot=BOOT):
    genus,gi=np.unique(species.genus.to_numpy(str),return_inverse=True)
    size=np.bincount(gi).astype(float)
    if size.min()<3:raise ValueError("Source genus must have at least 3 observed taxa")
    weights=1/size[gi]
    feature_names=[x for v in BLOCKS.values() for x in v]
    x=species[feature_names].to_numpy(float)
    y=species.white.to_numpy(float)
    # These centring operations use ONLY predictor/site values and training
    # labels. The test-genus label mean is used exclusively to EVALUATE the
    # contrast, never to form a pairwise prediction.
    X=x-np.stack([np.bincount(gi,weights=x[:,j])/size for j in range(x.shape[1])],axis=1)[gi]
    Y=y-(np.bincount(gi,weights=y)/size)[gi]
    idx={name:[feature_names.index(x) for block in blocks for x in BLOCKS[block]]
         for name,blocks in MODELS.items()}
    fold=fold_of_genera[gi]
    if not np.array_equal(np.unique(fold),np.arange(5)):
        raise ValueError("Every source genus fold must appear")
    residual_errors=np.empty((len(y),len(MODELS)),float)
    # Full test-species predicted photocolour differences can be formed
    # from beta · (predictors_species_i - predictors_species_j), independent
    # of the unknown white average of the newly encountered genus.
    for k in range(5):
        tr=fold!=k
        te=~tr
        w=weights[tr]/weights[tr].sum()
        for col,indices in enumerate(idx.values()):
            if not indices:
                pred=np.zeros(te.sum(),float)
            else:
                sub=X[:,indices]
                means=(sub[tr]*w[:,None]).sum(axis=0)
                sd=np.maximum(np.sqrt((w[:,None]*np.square(sub[tr]-means)).sum(axis=0)),1e-9)
                z=(sub[tr]-means)/sd
                z_test=(sub[te]-means)/sd
                beta=np.linalg.solve((z.T*w)@z+ALPHA*np.eye(len(indices)),
                                     z.T@(w*Y[tr]))
                pred=z_test@beta
            residual_errors[te,col]=Y[te]-pred
    per_genus=np.stack([
        np.bincount(gi,weights=np.square(residual_errors[:,j]))/size
        for j in range(len(MODELS))
    ],axis=1)
    # Center only for numerical stability in case fold-normalization produces
    # a tiny nonzero genus prediction mean; true genus label means cancel in
    # actual difference (y_i-y_j) and ARE NOT required as predictor inputs.
    for j in range(len(MODELS)):
        means=np.bincount(gi,weights=residual_errors[:,j])/size
        centered=residual_errors[:,j]-means[gi]
        per_genus[:,j]=np.bincount(gi,weights=centered**2)/size
    pairwise=per_genus*(2*size/(size-1))[:,None]
    indices={name:j for j,name in enumerate(MODELS)}
    points=pairwise.mean(axis=0)
    return genus,fold,points,pairwise,indices

def summarize(genus,fold,points,pairwise,idx,seed):
    draw=np.random.default_rng(seed).integers(0,len(genus),size=(BOOT,len(genus)),dtype=np.int32)
    boot=pairwise[draw].mean(axis=1)
    def contrast(a,b):
        i,j=idx[a],idx[b]
        diff=boot[:,i]-boot[:,j]
        return {"point":round(float(points[i]-points[j]),12),
                "whole_genus_fixed_prediction_bootstrap95":np.round(np.quantile(diff,[.025,.975]),12).tolist(),
                "fraction_bootstrap_improved":round(float(np.mean(diff>0)),9)}
    return {
        "heldout_fold_species_counts":np.bincount(fold,minlength=5).tolist(),
        "pairwise_white_difference_squared_error":{
            k:round(float(points[j]),12) for k,j in idx.items()
        },
        "moisture_increment_over_geography":contrast("geography","geography_moisture"),
        "soil5_increment_over_fullclimate":contrast("geography_fullclimate","geography_fullclimate_soil5"),
        "soil5_vs_geography":contrast("geography","geography_soil5"),
    }

def build(photo_zip,env_zip):
    loader,species=source_species(photo_zip,env_zip)
    genus=np.unique(species.genus.to_numpy(str))
    cases={}
    for index,salt in enumerate(SALTS):
        split=folds(genus,salt,loader.fold)
        gen,fold,point,pairwise,ix=evaluate(species,split)
        if not np.array_equal(gen,genus):
            raise AssertionError("Genus order drift")
        label="original" if salt is None else "salt"+str(salt)
        cases[label]=summarize(gen,fold,point,pairwise,ix,BOOT_SEED+index)
    if abs(cases["original"]["pairwise_white_difference_squared_error"]["geography"]-.3162849480328433)>1e-9:
        raise AssertionError("Independent original-source pairwise task did not reproduce support")
    return {
        "schema":"mois_original_photo_white_congeneric_pairwise_transfer_v1",
        "source_sha256":{"photo":PHOTO_SHA,"environment":ENV_SHA},
        "population":{"qualified_photo_white_taxa":len(species),
                      "qualifying_genera":len(genus),"classified_source_photo_rows":int(species.n_photos.sum()),
                      "soil36_classifiable_photo_opportunities":31925},
        "folds":5,"n_source_genus_assignments":len(SALTS),
        "ridge_penalty":ALPHA,"bootstrap_draws":BOOT,"bootstrap_seed_base":BOOT_SEED,
        "target":"For an unseen nominal genus: photographically white species fraction of taxon i minus taxon j; requires predictor-site covariates for i and j but NO phenotypic mean of the test genus at prediction time.",
        "loss":"Mean of squared difference errors among ALL unordered original eligible species pairs WITHIN each test genus; each eligible genus contributes one equal weight.",
        "pairwise_identity":"mean_{i<j}((y_i-y_j)-(yhat_i-yhat_j))² = [2*n_g/(n_g-1)]*mean_i[(e_i-mean_g(e))²]",
        "training":"Within-genus-centred species predictor and photo-white fractions in training taxa; genus-equal ridge0.01 no intercept; five SHA256 genus-disjoint folds; predictors scaled on training-fold genera only.",
        "uncertainty":"1999 resamples of fixed per-genus heldout pairwise losses, NOT model refitting, full research selection correction or external botanical replication.",
        "comparisons":cases,
        "finding":"Across all four genus assignments, genus-mean-centred environmental information has small positive POINT contribution to pairwise congeneric difference prediction from adding rainfall+vapour to broad geography, but EVERY original-source fold-specific 95% interval crosses zero. Increment of Soil5 on full climate changes sign across genus assignments and is not clearly different from zero.",
        "limits":[
            "Pairwise comparisons are among observed classifiable-photo taxon fractions; a source photograph is not a genetic morph, flower-pigment assay or fitness outcome.",
            "Predictor values of two test-congeneric species must be known; no test-genus white phenotype labels or mean white phenotype are needed to FORM their predicted difference.",
            "The test pair's real photo-colour labels ARE used to EVALUATE prediction error; this is a retrospective source-only ecological test.",
            "Genus source labels are not a phylogenetic covariance model; broad site averaging and photo classifiability can bias all outcomes.",
            "Invariance of a small positive point sign across four nonindependent fold assignments is NOT independent statistical confirmation.",
            "No real DEM terrain slope/aspect/TPI/TRI, reproductive success, metabolic pigment costs or natural selection were measured.",
            "No source species names, individual genus names, photo IDs, photographs, site coordinates or attributed media are exported."
        ]
    }

def compare(a,b,key="root"):
    if isinstance(a,dict):
        if not isinstance(b,dict) or set(a)!=set(b):raise AssertionError("Changed keys: "+key)
        for k,val in a.items():compare(val,b[k],key+"."+k)
    elif isinstance(a,list):
        if not isinstance(b,list) or len(a)!=len(b):raise AssertionError("Length drift "+key)
        for i,z in enumerate(a):compare(z,b[i],key+f"[{i}]")
    elif isinstance(a,(float,int)) and not isinstance(a,bool):
        if abs(a-b)>1e-8:raise AssertionError((key,a,b))
    elif a!=b:raise AssertionError((key,a,b))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--photo-zip",type=Path,required=True)
    ap.add_argument("--env-zip",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    mode=ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write",action="store_true")
    mode.add_argument("--check",action="store_true")
    args=ap.parse_args()
    results=build(args.photo_zip,args.env_zip)
    if args.write:
        args.output.write_text(json.dumps(results,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
        print("WROTE source-safe original-photo pairwise congeneric transfer")
    else:
        compare(results,json.loads(args.output.read_text(encoding="utf8")))
        print("PASS 16801 original taxa, 1665 genera, four genus fold assignments and pairwise white-difference comparisons")

if __name__=="__main__":
    main()
