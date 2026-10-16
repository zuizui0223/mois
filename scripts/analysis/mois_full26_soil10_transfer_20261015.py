#!/usr/bin/env python3
"""Full climate26 vs climate26+soil5/soil10 predictive block test.

Exact original FCP photo/environment source SHA checks come from
mois_full_biogeography_climate_soil_audit_20261015.py. All three predictors
share IDENTICAL 16,801 taxa and 1,665 original nominal genera. The training
folds hide entire genera, scaling is fitted only on training source species,
the ridge alpha is fixed and numeric source results contain no location data.
"""
from __future__ import annotations
import argparse, importlib.util, json
from pathlib import Path
import hashlib
import numpy as np

SEED=202610151710
N_BOOT=1999

def audit_module():
 path=Path(__file__).with_name("mois_full_biogeography_climate_soil_audit_20261015.py")
 spec=importlib.util.spec_from_file_location("mois_original_full_environment",path)
 m=importlib.util.module_from_spec(spec)
 spec.loader.exec_module(m)
 return m

def calculate(photo_zip:Path,env_zip:Path):
 m=audit_module()
 _,_,_,species,_=m.original(photo_zip,env_zip)
 if (len(species),int(species.genus.nunique()))!=(16801,1665):
  raise ValueError("Original soil36 species/genus full predictor support drift")
 cols={
  "geo_plus_all26_climate":m.G+m.C,
  "geo_plus_all26_climate_plus_soil5":m.G+m.C+m.SOIL5,
  "geo_plus_all26_climate_plus_soil10":m.G+m.C+m.S,
 }
 fold=species.genus.map(m.fold).to_numpy(int)
 sizes=species.genus.value_counts()
 weight=species.genus.map(1/sizes).to_numpy(float)
 y=species.white.to_numpy(float)
 losses=np.empty((len(species),len(cols)))
 for k in range(5):
  train=fold!=k;test=~train
  w=weight[train]/weight[train].sum()
  mu_y=float(np.sum(y[train]*w))
  for j,(label,features) in enumerate(cols.items()):
   X=species[features].to_numpy(float)
   avg=np.sum(X[train]*w[:,None],axis=0)
   sd=np.maximum(np.sqrt(np.sum((X[train]-avg)**2*w[:,None],axis=0)),1e-9)
   Z=(X[train]-avg)/sd;T=(X[test]-avg)/sd
   coef=np.linalg.solve((Z.T*w)@Z+.01*np.eye(len(features)),(Z.T*w)@(y[train]-mu_y))
   pred=np.clip(mu_y+T@coef,0,1)
   losses[test,j]=(y[test]-pred)**2
 g,gi=np.unique(species.genus.to_numpy(),return_inverse=True)
 n=np.bincount(gi)
 loss_g=np.stack([np.bincount(gi,weights=losses[:,j])/n for j in range(len(cols))],axis=1)
 fitted=loss_g.mean(axis=0)
 rng=np.random.default_rng(SEED)
 draws=rng.integers(0,len(g),size=(N_BOOT,len(g)))
 boot=loss_g[draws].mean(axis=1)
 return {
  "schema":"mois_original_full26_climate_soil10_same_cohort_unseen_genus_v1",
  "photo_sha256":m.PHOTO_SHA,"environment_sha256":m.ENV_SHA,
  "eligible_species":16801,"eligible_genera":1665,
  "cross_validation":"five disjoint SHA256 hashed genus folds; no source genus is shared between train/test within a fold",
  "n_folds":5,"ridge_alpha":.01,
  "scaling":"training-fold equal-genus mean/SD only; prediction clipped to [0,1]",
  "target":"classified photo white fraction per source named taxon; heldout equal genus squared error, not genotype or lifetime fitness",
  "bootstrap":{"unit":"heldout observed whole genera, no model refitting","draws":N_BOOT,"seed":SEED},
  "model_results":{key:{
    "heldout_species_photo_white_squared_error":float(fitted[j]),
    "incremental_soil_prediction_gain_over_geo_climate26":float(fitted[0]-fitted[j]),
    "incremental_gain_genus_bootstrap95":np.quantile(boot[:,0]-boot[:,j],[.025,.975]).tolist()
   } for j,key in enumerate(cols)},
  "interpretation":"On the same original soil36-complete taxon support, all26 climate plus geography did not gain from additional 5 or 10 SoilGrids predictor fields under this fixed ridge/unseen-genus loss setup. Negative heldout gain does NOT show soil has no ecological or causal effect.",
  "limits":"Climate26 includes site elevation; the added soil5/soil10 fields are often correlated, sand-silt-clay texture is compositional; unlicensed photo/site data and real DEM terrain remain absent from public GitHub."
 }

def main():
 p=argparse.ArgumentParser()
 p.add_argument("--photo-zip",required=True,type=Path)
 p.add_argument("--env-zip",required=True,type=Path)
 p.add_argument("--output",required=True,type=Path)
 mode=p.add_mutually_exclusive_group(required=True)
 mode.add_argument("--write",action="store_true")
 mode.add_argument("--check",action="store_true")
 args=p.parse_args()
 out=calculate(args.photo_zip,args.env_zip)
 if args.write:
  args.output.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
  print("WROTE exact-source full26 versus plus-soil36 genus-heldout model audit")
 else:
  old=json.loads(args.output.read_text(encoding="utf8"))
  def compare(a,b,path="root"):
   if isinstance(a,dict):
    assert isinstance(b,dict) and set(a)==set(b),path
    for k,v in a.items():compare(v,b[k],path+"."+k)
   elif isinstance(a,list):
    assert isinstance(b,list) and len(a)==len(b),path
    for i,v in enumerate(a):compare(v,b[i],path+"["+str(i)+"]")
   elif isinstance(a,(float,int)) and not isinstance(a,bool):
    assert isinstance(b,(float,int)) and abs(a-b)<1e-9,(path,a,b)
   else:assert a==b,(path,a,b)
  compare(out,old)
  print("PASS original species-match full26 climate and full soil10 heldout tests")

if __name__=="__main__":main()
