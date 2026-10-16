#!/usr/bin/env python3
"""Compare between-genus versus within-genus flower-white/moisture slopes.

Owner-only original measured-photo/environment ZIP input; exports ONLY
aggregate results, no taxon names, photo IDs, coordinates or attribution.
The result is retrospective descriptive biogeography, not fitness selection.
"""
from __future__ import annotations
import argparse, hashlib, io, json, zipfile
from pathlib import Path
import numpy as np
import pandas as pd

PHOTO_SHA = '2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0'
ENV_SHA = '6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca'
HUES = ('white', 'yellow_orange', 'red_pink', 'blue_purple')
AXES = ('wc_bio12', 'wc_vapr_annual_mean')
GEO = ('abs_lat', 'lon_sin', 'lon_cos')
FEATURES = [*(f'wc_bio{i}' for i in range(1,20)), 'wc_srad_annual_mean','wc_srad_monthly_cv',
            'wc_wind_annual_mean','wc_wind_monthly_cv','wc_vapr_annual_mean','wc_vapr_monthly_cv',
            'wc_elevation_m']
SEED, BOOT = 2026101515, 1999


def sha(path: Path) -> str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1 << 20), b''):h.update(block)
    return h.hexdigest()


def read(zip_path: Path, member: str, cols: list[str]) -> pd.DataFrame:
    with zipfile.ZipFile(zip_path) as z:
        with z.open(member) as f:
            return pd.read_csv(f,compression='gzip',usecols=cols,low_memory=False)


def load_species(photo_zip: Path, env_zip: Path):
    if (sha(photo_zip),sha(env_zip)) != (PHOTO_SHA,ENV_SHA):raise ValueError('Original input SHA mismatch')
    for f in (photo_zip,env_zip):
        with zipfile.ZipFile(f) as z:
            if z.testzip() is not None:raise ValueError('Original source ZIP corrupted')
    cols=['inat_taxon_id','observation_id','photo_id','species','morph','measurement_status']
    br=read(photo_zip,'fcp42111_real_site_abiotic/breadth_flower_colour_climate_soil.csv.gz',cols)
    ce=read(photo_zip,'fcp42111_real_site_abiotic/taxon_cell_flower_colour_climate_soil.csv.gz',cols)
    rows=pd.concat([br,ce],ignore_index=True)
    if (len(br),len(ce),int(rows.photo_id.duplicated().sum()))!=(42111,85337,26905):raise ValueError('Photo source denominator drift')
    if rows.groupby('photo_id').inat_taxon_id.nunique().max()!=1:raise ValueError('Different taxa assigned same photo')
    rows=rows.drop_duplicates('photo_id')
    if len(rows)!=100543:raise ValueError('Photo source unique count drift')
    eco=read(env_zip,'fcp42111_original_photo_all_climate_soil/all_original_100543_photo_sites_expanded_environment.csv.gz',
             ['inat_taxon_id','observation_id','photo_id','latitude','longitude',*FEATURES])
    if len(eco)!=100543:raise ValueError('Climate source denominator drift')
    m=rows.merge(eco,on=['inat_taxon_id','observation_id','photo_id'],validate='one_to_one')
    eligible=(m.measurement_status.eq('classified_four_state_morph')&m.morph.isin(HUES)&
              m.latitude.between(-90,90)&m.longitude.between(-180,180)&m[FEATURES].notna().all(axis=1))
    m=m.loc[eligible].copy()
    if len(m)!=45836:raise ValueError('Climate26 source cohort drift')
    m['genus']=m.species.str.split().str[0]
    m['abs_lat']=m.latitude.abs()
    m['lon_sin']=np.sin(np.deg2rad(m.longitude))
    m['lon_cos']=np.cos(np.deg2rad(m.longitude))
    for hue in HUES:m['hue_'+hue]=m.morph.eq(hue).astype(float)
    norm={a:{'mean':float(m[a].mean()),'sd':float(m[a].std(ddof=0))} for a in [*GEO,*AXES]}
    fields=[*GEO,*AXES,*('hue_'+h for h in HUES)]
    species=m.groupby('inat_taxon_id',sort=True).agg({**{x:'mean' for x in fields},'genus':'first','photo_id':'size'}).rename(columns={'photo_id':'n_photo'}).reset_index()
    n=species.genus.value_counts()
    species=species.loc[species.genus.isin(n.index[n>=3])].copy().reset_index(drop=True)
    if (len(species),int(species.genus.nunique()))!=(20694,1934):raise ValueError('Qualified species/genus cohort drift')
    return species,norm,{'n_unique_original_photos':100543,'n_photo_climate26':len(m),'n_climate26_species':int(m.inat_taxon_id.nunique()),
                         'n_primary_taxa':len(species),'n_primary_genera':int(species.genus.nunique()),'n_primary_photos':int(species.n_photo.sum())}


def genus_moments(s: pd.DataFrame, norm: dict, axis: str):
    names=[*GEO,axis]
    gids,gi=np.unique(s.genus.to_numpy(),return_inverse=True)
    sizes=np.bincount(gi).astype(float)
    xx=s[names].to_numpy(float)
    means=np.array([norm[a]['mean'] for a in names]);scales=np.array([norm[a]['sd'] for a in names])
    X=(xx-means)/scales
    y=s['hue_white'].to_numpy(float)
    xm=np.stack([np.bincount(gi,weights=X[:,j])/sizes for j in range(len(names))],axis=1)
    ym=np.bincount(gi,weights=y)/sizes
    Xw=X-xm[gi];yw=y-ym[gi];w=1/sizes[gi]
    # Each genus contributes one total weight to within-genus contrasts.
    cross=(Xw[:,:,None]*Xw[:,None,:])*w[:,None,None]
    xy=(Xw*yw[:,None])*w[:,None]
    AW=np.zeros((len(gids),4,4));bW=np.zeros((len(gids),4));np.add.at(AW,gi,cross);np.add.at(bW,gi,xy)
    # Genus-equal between-genus slope with genus-mean geographic covariates.
    B=np.column_stack((np.ones(len(gids)),xm))
    AB=B[:,:,None]*B[:,None,:];bB=B*ym[:,None]
    # Pooled, genus-equal species slope with geography; point comparator only.
    P=np.column_stack((np.ones(len(s)),X))
    pooled=np.linalg.solve((P.T*w)@P,(P.T*w)@y)[-1]
    return AW,bW,AB,bB,pooled


def regress(aw,bw,ab,bb):
    betaw=np.linalg.solve(aw,bw)[-1]
    betab=np.linalg.solve(ab,bb)[-1]
    return np.array([betaw,betab,betab-betaw])


def analyze(s:pd.DataFrame,norm:dict,axes=AXES,boot=BOOT,seed=SEED):
    moments={};points={};total_g=int(s.genus.nunique())
    for axis in axes:
        AW,bW,AB,bB,pooled=genus_moments(s,norm,axis)
        point=regress(AW.sum(axis=0),bW.sum(axis=0),AB.sum(axis=0),bB.sum(axis=0))
        points[axis]=point*100
        moments[axis]=(AW,bW,AB,bB,float(pooled*100))
    # Same whole-genus multinomial cluster draws for both axes and scales.
    rng=np.random.default_rng(seed)
    draws=rng.multinomial(total_g,[1/total_g]*total_g,size=boot)
    estimates={}
    for axis in axes:
        AW,bW,AB,bB,pool=moments[axis]
        boot_aw=np.einsum('bg,gij->bij',draws,AW,optimize=True)
        boot_bw=np.einsum('bg,gi->bi',draws,bW,optimize=True)
        boot_ab=np.einsum('bg,gij->bij',draws,AB,optimize=True)
        boot_bb=np.einsum('bg,gi->bi',draws,bB,optimize=True)
        wb=np.linalg.solve(boot_aw,boot_bw[...,None])[:,-1,0]
        bb=np.linalg.solve(boot_ab,boot_bb[...,None])[:,-1,0]
        sample=np.column_stack([wb,bb,bb-wb])*100
        estimate={}
        for j,term in enumerate(('within_nominal_genus','between_genus_means','between_minus_within')):
            ci=np.quantile(sample[:,j],[.025,.975])
            estimate[term]={'point_pp_per_original_photo_sd':float(points[axis][j]),'whole_genus_percentile_CI95_pp':[float(ci[0]),float(ci[1])],
                            'bootstrap_positive_fraction':float(np.mean(sample[:,j]>0))}
        estimate['pooled_genus_equal_geographic_adjusted_pp_per_original_photo_sd']=pool
        estimate['photo_environment_source_sd']=norm[axis]['sd']
        estimates[axis]=estimate
    return {'n_bootstrap':boot,'bootstrap_seed':seed,'bootstrap_genus_unit':total_g,
            'axes':estimates,
            'note':'Within and between are different observational covariate contrasts on identical species. Between-minus-within has paired same-genus bootstrap but is not a causal mixture decomposition; pooled slope is not algebraic sum.'}


def build(photo_zip:Path,env_zip:Path) -> dict:
    s,norm,support=load_species(photo_zip,env_zip)
    result=analyze(s,norm)
    support_sensitivity=[]
    for limit in (3,5,10,20):
        n=s.genus.value_counts()
        subset=s.loc[s.genus.isin(n.index[n>=limit])]
        row={'min_classified_climate_species_per_genus':limit,'n_species':int(len(subset)),'n_genera':int(subset.genus.nunique()),'point_coefficients_pp':{}}
        for axis in AXES:
            aw,bw,ab,bb,_=genus_moments(subset,norm,axis)
            values=100*regress(aw.sum(axis=0),bw.sum(axis=0),ab.sum(axis=0),bb.sum(axis=0))
            row['point_coefficients_pp'][axis]={'within':float(values[0]),'between':float(values[1]),'between_minus_within':float(values[2])}
        support_sensitivity.append(row)
    if not (abs(result['axes']['wc_bio12']['within_nominal_genus']['point_pp_per_original_photo_sd']-2.254545235123587)<1e-9 and
            abs(result['axes']['wc_vapr_annual_mean']['within_nominal_genus']['point_pp_per_original_photo_sd']-3.951850325152891)<1e-9):
        raise AssertionError('Failed to reproduce exact historical 26-climate world genus-centred white slopes')
    return {'schema':'mois_same_source_genus_between_within_moisture_v1','date':'2026-10-15',
            'source_sha256':{'photo':PHOTO_SHA,'environment':ENV_SHA},'source_cohort':support,
            'design':'Exact original climate26 photo-eligible species with >=3 source species/nominal genus. Four-hue source white photo fraction per species. For each moisture axis, regress genus-centred white versus genus-centred axis controlling genus-centred broad photo-site abs latitude and longitude sine/cosine; a genus-equal model. Separately regress genus-mean white against genus-mean axis and the same three genus-mean geographic nuisance covariates, with each genus weight one. One source photo-site environmental SD anchors both slopes. Pooled genus-equal ordinary slope recorded for comparison only. Genus-level paired bootstrap 1999 samples same genus draws for within, between, their difference and both axes.',
            'estimates':result,'descriptive_genus_support_thresholds_without_new_bootstrap_intervals':support_sensitivity,'limitations':['Source-exposed post hoc distinction; within and between observational contrasts are not identical biological estimands or a causal decomposition.',
            'Genus means reflect named and classified photograph opportunities, unequal taxon sampling and broad photo-site geography; genus labels are not a phylogeny.',
            'Confidence intervals are ordinary within-design whole-genus bootstrap percentiles; not whole-programme model-selection-corrected.',
            'Neither flowering-plant relative abundance nor within-population morph genetic frequencies, pollen/seed fitness, or pigment costs were measured.',
            'Raw photo labels/environment measured in archived original data; no new photo pixels remeasured. Source archive is not publicly licensed for sharing and stays outside GitHub.']}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--photo-zip',type=Path,required=True)
    ap.add_argument('--env-zip',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    action=ap.add_mutually_exclusive_group(required=True)
    action.add_argument('--write',action='store_true');action.add_argument('--check',action='store_true')
    args=ap.parse_args()
    result=build(args.photo_zip,args.env_zip)
    if args.write:
        args.output.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
        print('WROTE aggregate-only source genus between/within humidity result',args.output)
    else:
        frozen=json.loads(args.output.read_text(encoding='utf8'))
        assert frozen==result,'Original source numeric replay differs from committed aggregate'
        print('PASS source-exact paired genus/moisture decomposition and 1999 genus bootstrap')

if __name__=='__main__':main()
