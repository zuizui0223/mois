#!/usr/bin/env python3
"""Source-verified soil-versus-climate block comparison on identical taxa.

Requires owner-held SHA256 original FCP photo/environment ZIPs. Uses the
source cohort loader in mois_full_biogeography_climate_soil_audit_20261015.py.
Only aggregate losses/intervals are exported. No geolocations, taxon names,
flower photographs or source photographer rights are redistributed.
"""
from __future__ import annotations
import argparse,importlib.util,json
from pathlib import Path
import numpy as np

ALPHAS=(.001,.01,.1)
DRAW_COUNT=1999
SEED=202610151515

def source_module():
    path=Path(__file__).with_name("mois_full_biogeography_climate_soil_audit_20261015.py")
    spec=importlib.util.spec_from_file_location("mois_soil_original_loader",path)
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m

def calculate(photo_zip:Path,environment_zip:Path):
    source=source_module()
    _,soil,_,taxa,_=source.original(photo_zip,environment_zip)
    if (len(soil),len(taxa),taxa.genus.nunique(),int(taxa.n_photos.sum()))!=(31925,16801,1665,26163):
        raise ValueError("Original soil36 full-source cohort drift")
    G,C,S,W=source.G,source.C,source.S,source.W
    other24=[x for x in C if x not in W]
    models={
        "geography":G,"geography_soil5":G+S[:5],"geography_soil10":G+S,
        "geography_water":G+W,"geography_other24climate":G+other24,
        "geography_full26climate":G+C,
        "geography_full26climate_soil5":G+C+S[:5],
        "geography_full26climate_soil10":G+C+S,
    }
    genus,gi=np.unique(taxa.genus.to_numpy(str),return_inverse=True)
    n=np.bincount(gi); weights=1/n[gi]
    fold=np.array([source.fold(x) for x in genus])[gi]
    X={name:taxa[columns].to_numpy(float) for name,columns in models.items()}
    y=taxa.white.to_numpy(float)
    columns=list(models)
    result={}
    for alpha in ALPHAS:
        err=np.empty((len(y),len(columns)),float)
        for k in range(5):
            train=fold!=k;test=~train
            w=weights[train]/weights[train].sum()
            intercept=float(np.sum(w*y[train]))
            for j,name in enumerate(columns):
                array=X[name]
                mu=np.sum(array[train]*w[:,None],axis=0)
                sd=np.maximum(np.sqrt(np.sum(w[:,None]*(array[train]-mu)**2,axis=0)),1e-9)
                centered=(array[train]-mu)/sd
                beta=np.linalg.solve((centered.T*w)@centered+alpha*np.eye(centered.shape[1]),
                                    (centered.T*w)@(y[train]-intercept))
                pred=np.clip(intercept+((array[test]-mu)/sd)@beta,0,1)
                err[test,j]=(y[test]-pred)**2
        per_genus=np.stack([np.bincount(gi,weights=err[:,j])/n
                            for j in range(len(columns))],axis=1)
        observed=per_genus.mean(axis=0)
        if alpha==.01:
            expected={"geography":.22981006389209194,
                      "geography_water":.22745535180555312,
                      "geography_full26climate":.22813288536129273,
                      "geography_full26climate_soil5":.22853734681634522,
                      "geography_full26climate_soil10":.2286610993702406}
            for key,value in expected.items():
                if not np.isclose(observed[columns.index(key)],value,rtol=0,atol=1e-10):
                    raise AssertionError("Independent replayer contradicts frozen original-source CV "+key)
        rng=np.random.default_rng(SEED+int(alpha*1e6))
        sample=rng.integers(0,len(genus),size=(DRAW_COUNT,len(genus)),dtype=np.int32)
        boot=per_genus[sample].mean(axis=1)
        geography=columns.index("geography")
        rows={}
        for j,name in enumerate(columns):
            rows[name]={
                "loss":round(float(observed[j]),9),
                "gain_vs_geography":round(float(observed[geography]-observed[j]),9),
                "gain_ci95":np.round(np.quantile(boot[:,geography]-boot[:,j],[.025,.975]),9).tolist()}
        contrasts={}
        comparisons={
            "water_vs_soil5":("geography_water","geography_soil5"),
            "water_vs_soil10":("geography_water","geography_soil10"),
            "increment_water_given_other24climate":("geography_full26climate","geography_other24climate"),
            "increment_soil5_given_full26":("geography_full26climate_soil5","geography_full26climate"),
            "increment_soil10_given_full26":("geography_full26climate_soil10","geography_full26climate")}
        for key,(left,right) in comparisons.items():
            i,j=columns.index(left),columns.index(right)
            d=boot[:,j]-boot[:,i]
            contrasts[key]={
                "gain":round(float(observed[j]-observed[i]),9),
                "genus_bootstrap95":np.round(np.quantile(d,[.025,.975]),9).tolist()}
        result[str(alpha)]={"models":rows,"contrasts":contrasts}
    return {
        "schema":"mois_soil_crossdomain_block_audit_v1",
        "date_jst":"2026-10-15",
        "source_sha256":{"photo":source.PHOTO_SHA,"environment":source.ENV_SHA},
        "population":{"n_species":len(taxa),"n_genera":len(genus),
                      "n_original_classifiable_photos":int(taxa.n_photos.sum())},
        "protocol":{"folds":5,"fold_species_counts":np.bincount(fold,minlength=5).tolist(),
                    "same_taxa":True,"genus_holdout":True,
                    "ridge_alpha_grid":list(ALPHAS),
                    "training_only_feature_scaling":True,
                    "target":"taxon mean of originally classified photo-white fractions",
                    "loss":"equal-genus mean squared prediction error",
                    "bootstrap_draws":DRAW_COUNT,"bootstrap_seed_base":SEED,
                    "bootstrap_is_fixed_prediction_resampling_not_nested_refitting":True},
        "by_alpha":result,
        "scope":"retrospective source association and prediction, not independent ecological replication or fitness selection",
        "terrain":"No original DEM neighbourhoods; local slope/aspect/TPI/TRI still not tested.",
        "no_raw_rows_exported":True}

def check_numeric(a,b,where="root"):
    if isinstance(a,dict):
        assert isinstance(b,dict) and set(a)==set(b),where
        for key in a:check_numeric(a[key],b[key],where+"."+key)
    elif isinstance(a,list):
        assert isinstance(b,list) and len(a)==len(b),where
        for j,(x,y) in enumerate(zip(a,b)):check_numeric(x,y,where+"["+str(j)+"]")
    elif isinstance(a,(int,float)) and not isinstance(a,bool):
        assert isinstance(b,(float,int)) and abs(a-b)<1e-9,(where,a,b)
    else:assert a==b,(where,a,b)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--photo-zip",required=True,type=Path)
    p.add_argument("--env-zip",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    sub=p.add_mutually_exclusive_group(required=True)
    sub.add_argument("--write",action="store_true")
    sub.add_argument("--check",action="store_true")
    a=p.parse_args()
    receipt=calculate(a.photo_zip,a.env_zip)
    if a.write:a.output.write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    else:check_numeric(receipt,json.loads(a.output.read_text(encoding="utf8")))
    print("PASS matched original-soil36 water-vs-soil and unique-climate block analysis")

if __name__=="__main__":main()
