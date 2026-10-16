#!/usr/bin/env python3
"""Original FCP photographed-white geography/vapour omitted-variable identity.

Uses exact SHA256-verified original measured-photo and site-environment ZIPs.
Outputs only aggregate covariances, coefficients and paired genus-bootstrap
intervals, never photo IDs, species/genus names or source coordinates.
Statistical coefficient changes do NOT establish causal ecological mediation.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import zipfile
import numpy as np
import pandas as pd

PHOTO_SHA='2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0'
ENV_SHA='6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca'
PHOTO_DIR='fcp42111_real_site_abiotic/'
ENV_MEMBER='fcp42111_original_photo_all_climate_soil/all_original_100543_photo_sites_expanded_environment.csv.gz'
PHOTO_FIELDS=['inat_taxon_id','observation_id','photo_id','species','morph','measurement_status']
CLIM=[*(f'wc_bio{i}' for i in range(1,20)),'wc_srad_annual_mean','wc_srad_monthly_cv','wc_wind_annual_mean','wc_wind_monthly_cv','wc_vapr_annual_mean','wc_vapr_monthly_cv','wc_elevation_m']
VARS=('abs_lat','lon_sin','lon_cos','wc_elevation_m','wc_bio12','wc_vapr_annual_mean')
MODELS={'geography_only':VARS[:4],'geography_plus_vapour':(*VARS[:4],VARS[5]),'geography_plus_rain':(*VARS[:4],VARS[4]),'geography_both_moisture':VARS}
N_BOOT=1999
SEED=202610151216

def sha256(f):
    h=hashlib.sha256()
    with open(f,'rb') as stream:
        for chunk in iter(lambda:stream.read(1<<20),b''):h.update(chunk)
    return h.hexdigest()

def read(zfile,member,columns):
    with zipfile.ZipFile(zfile) as archive:
        with archive.open(member) as stream:
            return pd.read_csv(stream,compression='gzip',usecols=columns,low_memory=False)

def source_species(photo_zip,env_zip):
    for path,expected in ((photo_zip,PHOTO_SHA),(env_zip,ENV_SHA)):
        if sha256(path)!=expected:raise ValueError('Source archive SHA256 mismatch')
        with zipfile.ZipFile(path) as archive:
            if archive.testzip() is not None:raise ValueError('Original ZIP integrity failure')
    breadth=read(photo_zip,PHOTO_DIR+'breadth_flower_colour_climate_soil.csv.gz',PHOTO_FIELDS)
    cells=read(photo_zip,PHOTO_DIR+'taxon_cell_flower_colour_climate_soil.csv.gz',PHOTO_FIELDS)
    union=pd.concat([breadth,cells],ignore_index=True)
    if (len(breadth),len(cells),int(union.photo_id.duplicated().sum()))!=(42111,85337,26905):
        raise ValueError('Historical source-photo denominators changed')
    if union.groupby('photo_id').inat_taxon_id.nunique().max()!=1:
        raise ValueError('Original photo ID has multiple taxa')
    union=union.drop_duplicates('photo_id')
    if len(union)!=100543:raise ValueError('Original unique source-photo count changed')
    env=read(env_zip,ENV_MEMBER,['inat_taxon_id','observation_id','photo_id','latitude','longitude',*CLIM])
    data=union.merge(env,on=['inat_taxon_id','observation_id','photo_id'],validate='one_to_one')
    if len(data)!=100543:raise ValueError('Original site-photo identity join changed')
    valid=(data.measurement_status.eq('classified_four_state_morph') &
           data.morph.isin(('white','yellow_orange','red_pink','blue_purple')) &
           data.latitude.between(-90,90) & data.longitude.between(-180,180) &
           data[CLIM].notna().all(axis=1))
    data=data.loc[valid].copy()
    if len(data)!=45836:raise ValueError('Historical complete climate26 photo panel changed')
    data['abs_lat']=data.latitude.abs()
    data['lon_sin']=np.sin(np.deg2rad(data.longitude))
    data['lon_cos']=np.cos(np.deg2rad(data.longitude))
    data['white']=data.morph.eq('white').astype(float)
    data['genus']=data.species.str.split().str[0]
    norms={v:{'mean':float(data[v].mean()),'sd':float(data[v].std(ddof=0))} for v in VARS}
    aggreg={v:'mean' for v in VARS}
    aggreg.update({'genus':'first','white':'mean','photo_id':'size'})
    species=data.groupby('inat_taxon_id',sort=True).agg(aggreg).rename(columns={'photo_id':'n_photo'}).reset_index()
    genus_size=species.genus.value_counts()
    species=species.loc[species.genus.isin(genus_size.index[genus_size>=3])].reset_index(drop=True)
    if (len(species),species.genus.nunique(),int(species.n_photo.sum()))!=(20694,1934,38498):
        raise ValueError('Original congeneric species panel changed')
    return species,norms

def genus_moments(species,norms):
    genus,gi=np.unique(species.genus.to_numpy(),return_inverse=True)
    counts=np.bincount(gi).astype(float)
    w=1/counts[gi]
    x=np.stack([(species[key].to_numpy(float)-norms[key]['mean'])/norms[key]['sd'] for key in VARS],axis=1)
    y=species.white.to_numpy(float)
    xm=np.stack([np.bincount(gi,weights=x[:,j])/counts for j in range(len(VARS))],axis=1)
    ym=np.bincount(gi,weights=y)/counts
    z=x-xm[gi];yc=y-ym[gi]
    A=np.zeros((len(genus),len(VARS),len(VARS)))
    b=np.zeros((len(genus),len(VARS)))
    np.add.at(A,gi,z[:,:,None]*z[:,None,:]*w[:,None,None])
    np.add.at(b,gi,z*yc[:,None]*w[:,None])
    return A,b

def calculate(photo_zip,env_zip,nboot=N_BOOT,seed=SEED):
    species,norms=source_species(photo_zip,env_zip)
    ga,gb=genus_moments(species,norms)
    A=ga.sum(axis=0);b=gb.sum(axis=0)
    indices={name:np.array([VARS.index(x) for x in cols]) for name,cols in MODELS.items()}
    def solve(matrix,outcome,index):
        return 100*np.linalg.solve(matrix[np.ix_(index,index)],outcome[index])
    points={name:solve(A,b,index) for name,index in indices.items()}
    original=float(solve(A,b,np.array([0,1,2,4]))[-1])
    if abs(original-2.254545235123587)>1e-9:
        raise ValueError('Main original white/BIO12 coefficient no longer replays')
    rng=np.random.default_rng(seed)
    draws=rng.multinomial(len(ga),np.ones(len(ga))/len(ga),size=nboot)
    bsA=np.einsum('bg,gij->bij',draws,ga,optimize=True)
    bsb=np.einsum('bg,gi->bi',draws,gb,optimize=True)
    sampled={}
    for name,idx in indices.items():
        sub=bsA[:,idx,:][:,:,idx]
        sampled[name]=100*np.linalg.solve(sub,bsb[:,idx,None])[:,:,0]
    comparisons={}
    pairs=(('add_vapour_to_geography','geography_only','geography_plus_vapour'),
           ('add_vapour_to_geography_and_rain','geography_plus_rain','geography_both_moisture'))
    for label,reduced,full in pairs:
        idx=indices[reduced];vap=VARS.index('wc_vapr_annual_mean')
        gamma=np.linalg.solve(A[np.ix_(idx,idx)],A[idx,vap])
        vap_beta=float(points[full][-1])
        omitted=gamma*vap_beta
        p0=points[reduced];p1=points[full][:-1]
        if not np.allclose(p0-p1,omitted,atol=1e-9):
            raise ValueError('Weighted omitted-variable identity violated')
        vap_r2=float((A[idx,vap]@np.linalg.solve(A[np.ix_(idx,idx)],A[idx,vap]))/A[vap,vap])
        changes={}
        for key,units,phys in (('abs_lat','10_deg_absolute_latitude',10),
                               ('wc_elevation_m','1000_m_elevation',1000)):
            col=list(MODELS[reduced]).index(key)
            scale=phys/norms[key]['sd']
            shift=(sampled[full][:,col]-sampled[reduced][:,col])*scale
            changes[key]={
                'unit':units,
                'reduced_white_coef_pp':round(float(p0[col]*scale),9),
                'with_vapour_white_coef_pp':round(float(p1[col]*scale),9),
                'increase_after_adding_vapour_pp':round(float((p1[col]-p0[col])*scale),9),
                'exact_omitted_variable_identity_contribution_pp':round(float(-omitted[col]*scale),9),
                'source_conditional_paired_genus_bootstrap95_for_coef_increase_pp':
                    np.round(np.quantile(shift,[.025,.975]),6).tolist(),
                'reduced_vs_full_same_taxa':True,
            }
        comparisons[label]={
            'reduced':reduced,'full':full,
            'vapour_model_coefficient_pp_per_source_sd':round(vap_beta,9),
            'vapour_variation_explained_by_reduced_geo_and_climate_predictors_R2':round(vap_r2,9),
            'vapour_collinearity_vif':round(float(1/(1-vap_r2)),9),
            'regression_of_vapour_on_reduced_predictors_standardized_gamma':
                {name:round(float(x),9) for name,x in zip(MODELS[reduced],gamma)},
            'geographic_coef_changes':changes,
        }
    def correlation(i,j):return round(float(A[i,j]/np.sqrt(A[i,i]*A[j,j])),9)
    return {
        'schema':'mois_within_genus_geography_vapour_omitted_variable_identity_v1',
        'source_photo_sha256':PHOTO_SHA,'source_environment_sha256':ENV_SHA,
        'source_support':{'photo_opportunities':100543,'photo_climate26_classified':45836,
                          'taxa':20694,'genera':1934,'original_classified_photos_in_eligible_taxa':38498},
        'original_genus_equal_BIO12_pp_per_sd_reproduced':round(original,9),
        'normalization':{k:round(norms[k]['sd'],9) for k in VARS},
        'within_genus_genus_equal_correlations':{
            'absolute_latitude_vs_vapour':correlation(0,5),
            'elevation_vs_vapour':correlation(3,5),
            'precipitation_vs_vapour':correlation(4,5)},
        'paired_genus_bootstrap':{
            'repetitions':nboot,'seed':seed,
            'uncertainty':'ordinary paired percentile 95% within-design, NOT selection-, multiple-model- or phylogeny-corrected'},
        'identity_comparisons':comparisons,
        'interpretation':'On identical source species/genus weights, geography-coefficient changes after adding vapour are exactly explained algebraically by the fitted vapour-white association times the projection of vapour on included geography/rain variables. This is omitted-variable algebra, not ecological causal mediation, independent altitude selection or floral-pigment economics.',
        'limits':[
            'Genus-centred photographic species associations: not a species-tree or genetically identified floral morph outcome.',
            'Omitted-variable identity holds in observational weighted least squares and does NOT estimate causal mediation.',
            'VIF is source-conditional shared support, not a standalone biological trait effect.',
            'Bootstrap sensitivity intervals are post hoc and lack study-wide hypothesis/model selection correction.',
            'Source photo-classification missingness and geography-specific opportunistic sampling remain unaddressed.',
            'No original photo IDs, taxon names, genotype frequencies, photo media or geographic source rows are exported.']
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--photo-zip',required=True,type=Path)
    p.add_argument('--env-zip',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path)
    mode=p.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write',action='store_true')
    mode.add_argument('--check',action='store_true')
    args=p.parse_args()
    result=calculate(args.photo_zip,args.env_zip)
    if args.write:
        args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
        print('WROTE aggregate-only original-source geographic coefficient shifts')
    else:
        assert json.loads(args.output.read_text(encoding='utf8'))==result
        print('PASS SHA256-original geographic coefficient identity and paired genus bootstrap')

if __name__=='__main__':main()
