#!/usr/bin/env python3
"""Source-locked within-genus photo-white environment block audit.

Two deliberately different targets:
(1) in-sample weighted WLS variance partition of genus-centred white-photo
    fraction, with five blocks of original soil36-complete species predictors;
(2) heldout-GENUS prediction of taxon deviations from their own sampled
    genus mean, not ordinary prospective absolute species-white forecasting.

The raw FCP photograph/environment ZIPs stay private. Uses full original
SHA256/identity/colour/soil36 eligibility checks in companion source loader.
Only original aggregate numeric receipts are publicly written.
"""
from __future__ import annotations
import argparse,hashlib,importlib.util,json
from pathlib import Path
from math import factorial
import numpy as np

SHAPLEY_DRAW=999
SHAPLEY_SEED=202610151610
CV_DRAW=1999
CV_SEED=202610151710
BLOCKS={
 "Geography":["abs_lat","lon_sin","lon_cos"],
 "Moisture":["wc_bio12","wc_vapr_annual_mean"],
 "Temperature and radiation":["wc_bio5","wc_srad_annual_mean"],
 "Elevation":["wc_elevation_m"],
 "Soil5":["soil_pH","soil_SOC","soil_N","soil_clay","soil_available_water_proxy"]
}
CV_MODELS={
 "null_mean_genus_fraction":[],
 "geography":["Geography"],
 "geography_elevation":["Geography","Elevation"],
 "geography_moisture":["Geography","Moisture"],
 "geography_moisture_temperature_elevation":["Geography","Moisture","Temperature and radiation","Elevation"],
 "geography_moisture_temperature_elevation_soil5":["Geography","Moisture","Temperature and radiation","Elevation","Soil5"],
 "geography_soil5":["Geography","Soil5"],
 "soil5_only":["Soil5"],
 "moisture_only":["Moisture"]
}

def source_module():
    path=Path(__file__).with_name("mois_full_biogeography_climate_soil_audit_20261015.py")
    spec=importlib.util.spec_from_file_location("fcp_original_environment_audit",path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def load(photo_zip:Path,env_zip:Path):
    source=source_module()
    _,soil,_,species,norm=source.original(photo_zip,env_zip)
    if (len(species),int(species.genus.nunique()),int(species.n_photos.sum()))!=(16801,1665,26163):
        raise ValueError("Soil36 source taxon/genus/photo cohort drift")
    names=[x for v in BLOCKS.values() for x in v]
    if len(set(names))!=13:raise ValueError("Duplicated source environmental block variables")
    genus,gi=np.unique(species.genus.to_numpy(),return_inverse=True)
    genus_n=np.bincount(gi)
    w=1/genus_n[gi]
    # Source photo-site standardisation is algebraically irrelevant to WLS R²,
    # but train-fold standardisation below uses observed genus-weighted SD.
    X=species[names].to_numpy(float)/np.array([norm[k] for k in names])
    yy=species.white.to_numpy(float)
    ybar=np.bincount(gi,weights=yy)/genus_n
    Y=yy-ybar[gi]
    means=np.stack([np.bincount(gi,weights=X[:,j])/genus_n for j in range(len(names))],axis=1)
    X=X-means[gi]
    if not np.all(np.isfinite(X)) or not np.all(np.isfinite(Y)):
        raise ValueError("Original source environmental or photo-white centering invalid")
    return source,soil,species,names,genus,gi,genus_n,w,X,Y

def moments(gi,w,X,Y,n_genera:int):
    n,d=X.shape
    aa=np.zeros((n_genera,d,d),float)
    bb=np.zeros((n_genera,d),float)
    sst=np.zeros(n_genera,float)
    np.add.at(aa,gi,X[:,:,None]*X[:,None,:]*w[:,None,None])
    np.add.at(bb,gi,X*Y[:,None]*w[:,None])
    np.add.at(sst,gi,w*Y**2)
    return aa,bb,sst

def fitted_r2_by_subset(A,b,sst,indices):
    scores=np.zeros(32,float)
    if sst<=0:raise ValueError("No genus-centred photographic-white variance")
    for mask in range(1,32):
        cols=indices[mask]
        mat=A[np.ix_(cols,cols)]
        rhs=b[cols]
        scores[mask]=float(rhs@np.linalg.solve(mat,rhs)/sst)
    return scores

def shapley(scores):
    weight=[factorial(i)*factorial(4-i)/factorial(5) for i in range(5)]
    contrib=np.zeros(5,float)
    for k in range(5):
        for mask in range(32):
            if mask & (1<<k):continue
            contrib[k]+=weight[mask.bit_count()]*(scores[mask|(1<<k)]-scores[mask])
    if not np.isclose(contrib.sum(),scores[31],rtol=0,atol=1e-12):
        raise AssertionError("Shapley allocations do not sum to in-sample R²")
    return contrib

def in_sample(a,b,sst,n_genera:int):
    keys=list(BLOCKS)
    col_indices={key:np.array([list(sum((v for v in BLOCKS.values()),[])).index(x) for x in val])
                 for key,val in BLOCKS.items()}
    masks={mask:np.concatenate([col_indices[keys[i]] for i in range(5) if mask&(1<<i)]).astype(int)
           for mask in range(1,32)}
    score=fitted_r2_by_subset(a.sum(axis=0),b.sum(axis=0),sst.sum(),masks)
    values=shapley(score)
    partial=np.array([score[31]-score[31^(1<<i)] for i in range(5)])
    rng=np.random.default_rng(SHAPLEY_SEED)
    boot=np.empty((SHAPLEY_DRAW,5))
    boot_partial=np.empty((SHAPLEY_DRAW,5))
    for j in range(SHAPLEY_DRAW):
        draws=rng.multinomial(n_genera,np.full(n_genera,1/n_genera))
        m=np.einsum("g,gij->ij",draws,a,optimize=True)
        rhs=draws@b
        denom=float(draws@sst)
        ss=fitted_r2_by_subset(m,rhs,denom,masks)
        boot[j]=shapley(ss)
        boot_partial[j]=[ss[31]-ss[31^(1<<i)] for i in range(5)]
    return {
      "genus_centered_photo_white_SST":float(sst.sum()),
      "in_sample_total_R2":float(score[31]),
      "predictor_blocks":BLOCKS,
      "shapley_bootstrap":{"unit":"original whole genus","seed":SHAPLEY_SEED,"draws":SHAPLEY_DRAW,
         "method":"genus-resampled WLS moment reweighting; positive Shapley values are nonnegative by construction, not P-values"},
      "block_attribution":{name:{
         "shapley_fraction_of_total_within_genus_white_variance":float(values[k]),
         "shapley_genus_bootstrap95":np.quantile(boot[:,k],[.025,.975]).tolist(),
         "partial_R2_given_all_other_four_blocks":float(partial[k]),
         "partial_R2_genus_bootstrap95":np.quantile(boot_partial[:,k],[.025,.975]).tolist()}
          for k,name in enumerate(keys)}
    }

def residual_heldout_cv(source,species,genus,gi,n,w,X,Y,features):
    # Uses sampled genotype-free source genus means in each TEST genus to
    # define the heldout target Y. It is only a retrospective diagnostic of
    # transfer of within-genus slope, NOT prospective absolute predictions.
    feature_idx={name:i for i,name in enumerate(features)}
    defs={name:[feature_idx[x] for b in blocks for x in BLOCKS[b]]
          for name,blocks in CV_MODELS.items()}
    folds=species.genus.map(source.fold).to_numpy(int)
    if np.bincount(folds,minlength=5).tolist()!=[3212,3505,3803,3306,2975]:
        raise AssertionError("Original source SHA256 genus-fold allocation drift")
    model_names=list(CV_MODELS)
    err=np.empty((len(species),len(model_names)))
    for k in range(5):
        train=folds!=k
        test=~train
        ww=w[train]/w[train].sum()
        for j,name in enumerate(model_names):
            ix=defs[name]
            if not ix:
                pred=np.zeros(test.sum())
            else:
                features=X[:,ix]
                centre=np.sum(features[train]*ww[:,None],axis=0)
                sd=np.maximum(np.sqrt(np.sum(ww[:,None]*(features[train]-centre)**2,axis=0)),1e-9)
                Z=(features[train]-centre)/sd
                V=(features[test]-centre)/sd
                beta=np.linalg.solve((Z.T*ww)@Z+.01*np.eye(len(ix)),(Z.T*ww)@Y[train])
                pred=V@beta
            err[test,j]=(Y[test]-pred)**2
    per_genus=np.stack([np.bincount(gi,weights=err[:,j])/n for j in range(len(model_names))],axis=1)
    point=per_genus.mean(axis=0)
    rng=np.random.default_rng(CV_SEED)
    boot=per_genus[rng.integers(0,len(genus),size=(CV_DRAW,len(genus)))].mean(axis=1)
    ix_g=model_names.index("geography")
    ix_climate=model_names.index("geography_moisture_temperature_elevation")
    return {
      "target":"Sample-mean-centred photo-white fraction among species within a heldout genus. The test-genus mean is retrospectively used to define label/predictors: this is NOT prospective unknown-genus absolute-white prediction.",
      "loss":"squared error of genus-centred species photograph-white fractions, weighted equally by test genus",
      "source_folds":"five SHA256 original genus-disjoint folds",
      "fold_species_counts":np.bincount(folds,minlength=5).tolist(),
      "model_setup":"training-only equal-genus standardization, zero intercept on genus-centred target, fixed ridge alpha=0.01, no residual prediction clipping",
      "bootstrap":{"seed":CV_SEED,"draws":CV_DRAW,"unit":"resample already-heldout genus mean losses; no model refitting"},
      "model_results":{name:{
        "heldout_residual_squared_error":float(point[j]),
        "gain_over_geography":float(point[ix_g]-point[j]),
        "gain_over_geography_genus_bootstrap95":np.quantile(boot[:,ix_g]-boot[:,j],[.025,.975]).tolist(),
        "increment_vs_geo_climate_elevation":float(point[ix_climate]-point[j]),
        "increment_vs_geo_climate_elevation_genus_bootstrap95":
            np.quantile(boot[:,ix_climate]-boot[:,j],[.025,.975]).tolist()
      } for j,name in enumerate(model_names)}
    }

def build(photo_zip,env_zip):
    source,soil,species,features,genus,gi,n,w,X,Y=load(photo_zip,env_zip)
    A,b,sst=moments(gi,w,X,Y,len(genus))
    ins=in_sample(A,b,sst,len(genus))
    if abs(ins["in_sample_total_R2"]-.0062023660640665)>1e-9:
        raise ValueError("Soil36 source thirteen-predictor in-sample genus WLS R² did not match original")
    cv=residual_heldout_cv(source,species,genus,gi,n,w,X,Y,features)
    if abs(cv["model_results"]["geography_moisture"]["gain_over_geography"]-.00023194)>1e-7:
        raise ValueError("Exact source within-genus residual CV moisture gain drift")
    return {
      "schema":"mois_original_16801_environmental_block_partition_plus_residual_transfer_v1",
      "date":"2026-10-15",
      "original_photo_sha256":source.PHOTO_SHA,
      "original_environment_sha256":source.ENV_SHA,
      "eligible_species":len(species),"eligible_genera":len(genus),
      "source_photo_hue_observations_of_eligible_species":int(species.n_photos.sum()),
      "same_support_soil36_photo_opportunities":len(soil),
      "shapley_explanatory_in_sample":ins,
      "within_genus_residual_out_of_genus":cv,
      "interpretation":"Soil5 can have measurable descriptive within-genus explained-photo-white variance in the source while failing to incrementally improve absolute white fraction prediction for new genera. A second genus-residual prediction diagnostic is not a forward prediction and also does not establish a reproducible incremental soil signal.",
      "limitations":[
        "Source colour outcomes are digitally classified photographs, not independently validated plant pigments, genetically segregating colours or reproductive fitness.",
        "Source soil36 is a highly selected subset with geographic and altitude-dependent predictor availability.",
        "All13-model total within-genus in-sample R2 is tiny. Shapley shares are noncausal and nonnegative by mathematical construction, not discovery P-values.",
        "Genus-heldout residual diagnostic uses the TEST genus's own observed species mean to define the target; comparing it with prospective absolute unseen-genus predictions requires keeping this target distinction explicit.",
        "Resampling already-generated heldout predictions does not resample training or fix nonrandom photographic source selection.",
        "This block set uses five selected soil variables, not all ten redundant soil texture terms.",
        "Elevation is point altitude; no actual independently sampled DEM slope/aspect/TPI/TRI terrain effect was estimated."]
    }

def compare(a,b,path="root"):
    if isinstance(a,dict):
        if not isinstance(b,dict) or set(a)!=set(b):raise AssertionError("JSON keys changed: "+path)
        for k,x in a.items():compare(x,b[k],path+"."+k)
    elif isinstance(a,list):
        if not isinstance(b,list) or len(a)!=len(b):raise AssertionError("JSON length changed: "+path)
        for i,x in enumerate(a):compare(x,b[i],path+f"[{i}]")
    elif isinstance(a,(float,int)) and not isinstance(a,bool):
        if not isinstance(b,(float,int)) or abs(a-b)>1e-7:
            raise AssertionError((path,a,b))
    elif a!=b:
        raise AssertionError((path,a,b))

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--photo-zip",required=True,type=Path)
    p.add_argument("--env-zip",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    act=p.add_mutually_exclusive_group(required=True)
    act.add_argument("--write",action="store_true")
    act.add_argument("--check",action="store_true")
    a=p.parse_args()
    r=build(a.photo_zip,a.env_zip)
    if a.write:
        a.output.write_text(json.dumps(r,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
        print("WROTE private-source replayed environmental block decomposition and residual transfer")
    else:
        compare(r,json.loads(a.output.read_text(encoding="utf8")))
        print("PASS original soil36 photo-source Shapley and within-genus residual heldout CV")

if __name__=="__main__":main()
