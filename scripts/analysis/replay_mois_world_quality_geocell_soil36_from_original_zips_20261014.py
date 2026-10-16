#!/usr/bin/env python3
"""Replay frozen original global photo-colour maxT under three sensitivity scopes.

Requires the *two exact original*, SHA256-frozen tabular-photo and all36 site
ZIPs plus the original 1280-row coefficient CSV.  Never writes photo IDs,
genus strings, taxon labels, original geographical coordinates or raw tables.

The source-study units and 999 genus bootstrap are preserved. This audit is
retrospective source-numerical re-estimation, NOT independent biology or a
programme-wide correction for hypothesis/region/model selection.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import zipfile
from pathlib import Path
import numpy as np
import pandas as pd
from replay_mois_original_matched_spatial_576_20261014 import (
    original_source, species_frame, read_member, KEYS, CLIMATE, GEO, HUES,
)

SOIL=["soil_pH","soil_SOC","soil_N","soil_clay","soil_available_water_proxy",
      "soil_cec_0_30cm_source_raw","soil_sand_0_30cm_source_raw",
      "soil_silt_0_30cm_source_raw","soil_bdod_0_30cm_source_raw",
      "soil_cfvo_0_30cm_source_raw"]
AXES=["wc_bio5","wc_bio12","wc_bio14","wc_bio15","wc_srad_annual_mean",
      "wc_vapr_annual_mean"]
SCOPE_PARAMS={
    "strict_photo_climate26":{
        "source_mode":"climate","quality":"strict","adjust":"genus_geo",
        "n_species":5848,"n_genera":727,"n_photos":8538,
        "n_subgroups":727,"seed":2026101455},
    "within_genus_geocell_climate26":{
        "source_mode":"climate","quality":"all","adjust":"genus_geocell",
        "n_species":9370,"n_genera":845,"n_photos":14268,
        "n_subgroups":1847,"seed":2026102253},
    "global_complete_soil36":{
        "source_mode":"soil","quality":"all","adjust":"genus_geo",
        "n_species":16801,"n_genera":1665,"n_photos":31925,
        "n_subgroups":1665,"seed":2026101397},
}


def species_for_scope(photo_zip:Path, env_zip:Path, scope:str):
    # Independent raw input source reader is from the mois companion source
    # verifier, not the original FCP regression/CI function.
    images, clmean, clsd=original_source(photo_zip,env_zip)
    pars=SCOPE_PARAMS[scope]
    if scope=="strict_photo_climate26":
        cols=["photo_id","flower_effective_pixels","nuisance_pixel_fraction"]
        with zipfile.ZipFile(photo_zip) as z:
            breadth=read_member(z,"fcp42111_real_site_abiotic/breadth_flower_colour_climate_soil.csv.gz",cols)
            regional=read_member(z,"fcp42111_real_site_abiotic/taxon_cell_flower_colour_climate_soil.csv.gz",cols)
        quality=pd.concat([breadth,regional],ignore_index=True).drop_duplicates("photo_id")
        assert len(quality)==100543 and quality.photo_id.is_unique
        quality["is_strict"]=quality.flower_effective_pixels.ge(23826)&quality.nuisance_pixel_fraction.le(.1375717742195135)
        tagged=images.merge(quality[["photo_id","is_strict"]],on="photo_id",validate="one_to_one")
        images=tagged.loc[tagged.is_strict].copy()
        return species_frame(images),clmean,clsd
    if scope=="within_genus_geocell_climate26":
        return species_frame(images),clmean,clsd
    assert scope=="global_complete_soil36"
    with zipfile.ZipFile(env_zip) as z:
        env=read_member(z,"fcp42111_original_photo_all_climate_soil/all_original_100543_photo_sites_expanded_environment.csv.gz",KEYS+SOIL)
    assert env.photo_id.nunique()==100543
    expanded=images.merge(env,on=KEYS,how="left",validate="one_to_one")
    complete=expanded.loc[expanded[SOIL].notna().all(axis=1)].copy()
    assert (len(complete),complete.inat_taxon_id.nunique())==(31925,20786)
    columns=list(GEO)+list(CLIMATE)+SOIL
    means=complete[columns].mean()
    scales=complete[columns].std(ddof=0).where(lambda v:v>=1e-12,1.0)
    grouped=complete.groupby("inat_taxon_id",sort=True)
    species=grouped[columns+["hue_"+h for h in HUES]].mean()
    species["genus"]=grouped.genus.first()
    species["n_photos"]=grouped.photo_id.size()
    return species,means,scales


def normal_equations(spp:pd.DataFrame,means,scales,scope:str):
    p=SCOPE_PARAMS[scope]
    bycell=scope=="within_genus_geocell_climate26"
    axes=AXES+([] if scope!="global_complete_soil36" else ["soil_available_water_proxy"])
    predictors=axes if bycell else list(GEO)+axes
    x=(spp[predictors]-means[predictors])/scales[predictors]
    y=spp[["hue_"+h for h in HUES]]
    groups=spp.genus.astype(str) + ("__"+spp.cell.astype(str) if bycell else "")
    size=groups.value_counts()
    selected=groups.isin(size[size>=3].index)
    spp=spp.loc[selected].copy()
    groups=groups.loc[selected]
    assert (len(spp),spp.genus.nunique(),int(spp.n_photos.sum()))==(
        p["n_species"],p["n_genera"],p["n_photos"])
    assert groups.nunique()==p["n_subgroups"]
    x=x.loc[selected].to_numpy(float)
    y=y.loc[selected].to_numpy(float)
    namelist,ix=np.unique(groups.to_numpy(str),return_inverse=True)
    ngrp=np.bincount(ix)
    xmeans=np.stack([np.bincount(ix,weights=x[:,j])/ngrp for j in range(x.shape[1])],axis=1)
    ymeans=np.stack([np.bincount(ix,weights=y[:,j])/ngrp for j in range(y.shape[1])],axis=1)
    x-=xmeans[ix]
    y-=ymeans[ix]
    genera=sorted(spp.genus.unique())
    gmap={v:i for i,v in enumerate(genera)}
    ngen=spp.genus.value_counts()
    gm=np.zeros((len(genera),x.shape[1],x.shape[1]))
    gy=np.zeros((len(genera),x.shape[1],4))
    for j,gr in enumerate(namelist):
        active=ix==j
        genus=str(spp.loc[groups.eq(gr),"genus"].iloc[0])
        w=1.0/int(ngen[genus])
        gm[gmap[genus]]+=(x[active].T@x[active])*w
        gy[gmap[genus]]+=(x[active].T@y[active])*w
    return gm,gy,axes


def solve_bootstrap(gm,gy,axes,scope):
    seed=SCOPE_PARAMS[scope]["seed"]
    n=len(gm)
    draws=np.random.default_rng(seed).multinomial(
        n,np.repeat(1/n,n),size=999)
    mat=gm.sum(axis=0)
    rhs=gy.sum(axis=0)
    bm=np.einsum("bg,gij->bij",draws,gm,optimize=True)
    br=np.einsum("bg,gik->bik",draws,gy,optimize=True)
    nenv=len(axes)
    p=np.empty((nenv,4))
    bs=np.empty((999,nenv,4))
    for j in range(nenv):
        if scope=="within_genus_geocell_climate26":
            p[j]=rhs[j]/mat[j,j]
            bs[:,j]=br[:,j]/bm[:,j,j,None]
        else:
            cols=np.array([0,1,2,j+3])
            p[j]=np.linalg.solve(mat[np.ix_(cols,cols)],rhs[cols])[-1]
            bs[:,j]=np.linalg.solve(bm[:,cols,:][:,:,cols],br[:,cols])[:,-1]
    deviation=bs-p[None,:,:]
    sd=deviation.std(axis=0,ddof=1)
    assert (sd>1e-12).all()
    maxT=float(np.quantile(np.max(np.abs(deviation/sd),axis=(1,2)),.95))
    lo=p-maxT*sd
    hi=p+maxT*sd
    qlo,qhi=np.quantile(bs,[.025,.975],axis=0)
    assert np.allclose(p.sum(axis=1),0,atol=1e-10)
    return p,qlo,qhi,lo,hi,maxT


def replay(photo:Path,environment:Path,source_csv:Path,scope:str):
    if scope not in SCOPE_PARAMS:raise ValueError(scope)
    spp,means,sd=species_for_scope(photo,environment,scope)
    gm,gy,axes=normal_equations(spp,means,sd,scope)
    points,ql,qh,lo,hi,crit=solve_bootstrap(gm,gy,axes,scope)
    ledger=pd.read_csv(source_csv)
    assert len(ledger)==1280
    pars=SCOPE_PARAMS[scope]
    error=0.0
    verified=0
    supported=[]
    for j,axis in enumerate(axes):
        for k,hue in enumerate(HUES):
            rec=ledger.loc[(ledger["mode"]==pars["source_mode"]) &
              (ledger["quality"]==pars["quality"]) & (ledger["region"]=="world") &
              (ledger["adjustment"]==pars["adjust"]) &
              (ledger["weight"]=="genus_equal") & (ledger["anchor"]==axis) &
              (ledger["hue"]==hue)]
            assert len(rec)==1,(scope,axis,hue)
            row=rec.iloc[0]
            for val,col in (
              (points[j,k]*100,"direction_pp_per_env_global_sd"),
              (ql[j,k]*100,"cluster_ci_lo_pp"),(qh[j,k]*100,"cluster_ci_hi_pp"),
              (lo[j,k]*100,"maxT_ci_lo_pp"),(hi[j,k]*100,"maxT_ci_hi_pp")):
                diff=abs(float(row[col])-float(val))
                error=max(error,diff)
                if diff>1e-9:raise AssertionError((scope,axis,hue,col,diff))
                verified+=1
            sig=bool(lo[j,k]>0 or hi[j,k]<0)
            assert sig==bool(row["maxT_excludes_zero"])
            if sig:supported.append(axis+":"+hue)
    assert verified==(24 if scope!="global_complete_soil36" else 28)*5
    return {"schema":"mois_original_world_quality_spatial_soil36_independent_replay_v1",
      "status":"PASS_ALL_ORIGINAL_COEFFICIENTS_PERCENTILE_AND_MAXT_INTERVALS",
      "source_scope":scope,"photo_archive_sha256":"2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0",
      "environment_archive_sha256":"6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca",
      "eligible_named_species":pars["n_species"],"eligible_genera":pars["n_genera"],
      "source_photos":pars["n_photos"],"source_groups":pars["n_subgroups"],
      "bootstrap_draws":999,"original_bootstrap_seed":pars["seed"],
      "within_scope_maxT_critical":crit,"checked_numeric_point_and_CI_values":verified,
      "max_abs_error_pp":error,"significant_after_within_scope_maxT":supported,
      "not_independent_biology":True,"no_original_source_individual_ID_or_coordinates_exported":True,
      "source_rights_and_stable_reviewer_data_deposit_not_completed":True}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--photo",type=Path,required=True)
    p.add_argument("--environment",type=Path,required=True)
    p.add_argument("--ledger",type=Path,required=True)
    p.add_argument("--scope",choices=list(SCOPE_PARAMS),required=True)
    p.add_argument("--out",type=Path)
    args=p.parse_args()
    result=replay(args.photo,args.environment,args.ledger,args.scope)
    if args.out:args.out.write_text(json.dumps(result,indent=2)+"\n",encoding="utf8")
    print(json.dumps(result,indent=2))


if __name__=="__main__":
    main()
