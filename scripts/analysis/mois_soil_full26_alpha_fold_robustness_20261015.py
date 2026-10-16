#!/usr/bin/env python3
"""Source-only soil-block prediction robustness across ridge strengths and genus folds.

Requires owner-held SHA256-verified original FCP photo+environment ZIPs;
imports their exact source cohort loader from the committed macro ecology audit.
Emits ONLY global aggregate loss differences and grouped uncertainty, not
species, genera, georeferences, photo IDs or photographer attribution.
"""
from __future__ import annotations
import argparse, hashlib, importlib.util, json
from pathlib import Path
import numpy as np

ALPHAS=(.0003,.001,.003,.01,.03,.1,.3,1.)
SALTS=(None,202610151,202610152,202610153)
SEED=202610151888
NBOOT=1999

def source():
    p=Path(__file__).with_name("mois_full_biogeography_climate_soil_audit_20261015.py")
    spec=importlib.util.spec_from_file_location("mois_source_full_climate_soil",p)
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m

def split(genus,fold_fn,salt):
    if salt is None:return np.array([fold_fn(g) for g in genus],dtype=np.int8)
    return np.array([
        int.from_bytes(hashlib.sha256((str(salt)+"|"+str(g)).encode()).digest()[:8],"big")%5
        for g in genus],dtype=np.int8)

def model_losses(x,y,gi,ct,weight,fold,alpha):
    losses=np.empty((len(y),3),float)
    for k in range(5):
        train=fold!=k;test=~train
        ww=weight[train]/weight[train].sum()
        ybar=float(np.sum(y[train]*ww))
        for j,Z in enumerate(x):
            mu=np.sum(Z[train]*ww[:,None],axis=0)
            sd=np.maximum(np.sqrt(np.sum(ww[:,None]*(Z[train]-mu)**2,axis=0)),1e-9)
            xx=(Z[train]-mu)/sd
            coef=np.linalg.solve((xx.T*ww)@xx+alpha*np.eye(Z.shape[1]),(xx.T*ww)@(y[train]-ybar))
            pred=np.clip(ybar+(Z[test]-mu)/sd@coef,0,1)
            losses[test,j]=(y[test]-pred)**2
    genus_loss=np.stack([np.bincount(gi,weights=losses[:,j])/ct for j in range(3)],axis=1)
    return genus_loss.mean(axis=0),genus_loss

def calculate(photo_zip,env_zip):
    m=source()
    _,soil_photos,_,taxa,_=m.original(photo_zip,env_zip)
    if (len(soil_photos),len(taxa),taxa.genus.nunique(),int(taxa.n_photos.sum()))!=(31925,16801,1665,26163):
        raise ValueError("Exact original source soil36 support changed")
    fixed=json.loads((Path(__file__).resolve().parents[2]/
       "results/mois_full26_climate_soil10_transfer_20261015.json").read_text(encoding="utf8"))
    genus,gi=np.unique(taxa.genus.to_numpy(str),return_inverse=True)
    ct=np.bincount(gi);weights=1/ct[gi]
    x=[taxa[cols].to_numpy(float) for cols in (m.G+m.C,m.G+m.C+m.SOIL5,m.G+m.C+m.S)]
    y=taxa.white.to_numpy(float)
    group={}
    negatives=[0,0]
    per_alpha={str(a):{"soil5_negative":0,"soil10_negative":0} for a in ALPHAS}
    for rank,salt in enumerate(SALTS):
        key="original" if salt is None else "salt"+str(salt)
        fold=split(genus,m.fold,salt)[gi]
        if not np.array_equal(np.unique(fold),np.arange(5)):
            raise ValueError("Unoccupied original-genus heldout fold")
        table={"soil5":[],"soil10":[]}
        original_alpha=None
        for alpha in ALPHAS:
            avg,g_loss=model_losses(x,y,gi,ct,weights,fold,alpha)
            gains=avg[0]-avg[1:]
            for j,nm in enumerate(("soil5","soil10")):
                table[nm].append(round(float(gains[j]),9))
                if gains[j]<0:
                    negatives[j]+=1
                    per_alpha[str(alpha)][nm+"_negative"]+=1
            if alpha==.01:
                if salt is None:
                    for j,label in enumerate((
                         "geo_plus_all26_climate","geo_plus_all26_climate_plus_soil5",
                         "geo_plus_all26_climate_plus_soil10")):
                        if abs(float(avg[j])-fixed["model_results"][label]["heldout_species_photo_white_squared_error"])>1e-11:
                            raise AssertionError("Source original full26/soil10 genus-heldout model mismatch")
                rng=np.random.default_rng(SEED+rank)
                draw=rng.integers(0,len(genus),size=(NBOOT,len(genus)),dtype=np.int32)
                boot=g_loss[draw].mean(axis=1)
                inc=boot[:,0,None]-boot[:,1:]
                original_alpha={
                    "climate26_loss":round(float(avg[0]),9),
                    "soil5_increment":round(float(gains[0]),9),
                    "soil10_increment":round(float(gains[1]),9),
                    "soil5_bootstrap95":np.round(np.quantile(inc[:,0],[.025,.975]),9).tolist(),
                    "soil10_bootstrap95":np.round(np.quantile(inc[:,1],[.025,.975]),9).tolist(),
                    "bootstrap_seed":SEED+rank}
        group[key]={"source_species_per_heldout_fold":np.bincount(fold,minlength=5).tolist(),
                    "soil_block_increment_by_alpha":table,"at_original_alpha_001":original_alpha}
    return {
      "schema":"mois_soil_full26_alpha_fold_predictive_robustness_v1",
      "source_photo_sha256":m.PHOTO_SHA,"source_environment_sha256":m.ENV_SHA,
      "population":{"original_classifiable_soil36_photos":31925,
                    "qualified_source_photos":26163,"qualified_species":16801,"qualified_genera":1665},
      "same_species_across_all_models":True,
      "ridge_alpha_grid":list(ALPHAS),
      "genus_split_salts":["original",*[("salt"+str(s)) for s in SALTS[1:]]],
      "n_settings":len(ALPHAS)*len(SALTS),
      "summary":{"soil5_negative_gain_settings":negatives[0],
                 "soil10_negative_gain_settings":negatives[1],
                 "soil5_positive_gain_settings":len(ALPHAS)*len(SALTS)-negatives[0],
                 "soil10_positive_gain_settings":len(ALPHAS)*len(SALTS)-negatives[1],
                 "by_alpha_signs":per_alpha},
      "results":group,
      "prediction_protocol":{
        "target":"original classifiable photographed-white fraction per named source taxon",
        "training_test_split":"5 heldout SHA256 genus folds for every scheme, no shared genus",
        "training_feature_standardization":"equal-genus training means and SD only",
        "ridge_penalty":"fixed for every model within one setting, intercept unpenalized",
        "prediction_range":"clipped [0,1]",
        "loss":"heldout squared error averaged by named genus, not fitness",
        "bootstrap":"1999 resamples of fixed whole-genus heldout losses at alpha=.01, no model refitting"},
      "limitations":[
        "Four heldout assignments and eight penalties are 32 nonindependent checks on the identical original photo-derived flora.",
        "Changing penalty and fold changes prediction tasks; tiny positive single-settings can occur. No one setting demonstrates a soil adaptation mechanism.",
        "Even persistently negative soil predictive increments are not proof that soil is biologically irrelevant or cannot matter at local scale.",
        "Soil36 photo/site eligibility is nonrandom by elevation and region; photo hue labels are incomplete and classifier-dependent.",
        "Actual slope, aspect, TPI, TRI and terrain moisture accumulation are unmeasured without an independent original site-linked DEM.",
        "The source contains no relative fitness, flower pigment chemistry or genetically validated white/coloured morph frequencies.",
        "Original attributed photo/site, original individual named taxon/genus and positions are not exported."
      ]}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--photo-zip",type=Path,required=True)
    p.add_argument("--env-zip",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    flags=p.add_mutually_exclusive_group(required=True)
    flags.add_argument("--write",action="store_true")
    flags.add_argument("--check",action="store_true")
    opt=p.parse_args()
    result=calculate(opt.photo_zip,opt.env_zip)
    if opt.write:opt.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    else:assert result==json.loads(opt.output.read_text(encoding="utf8")),"Original-source alpha/fold check differs"
    print("PASS original-source full26 climate + soil model robustness, 4 folds × 8 ridge values")

if __name__=="__main__":main()
