#!/usr/bin/env python3
"""Pooled vs within-genus source-white variation with 13 environment covariates.

Reconstructs the *exact same* 16,801 species / 1,665 nominal genera from
SHA256-verified original FCP photo/soil ZIPs as the earlier 0.6202% within-genus
Shapley receipt. Compares species-equal and genus-equal pooled regressions.
Original source photo IDs, species names and coordinates NEVER enter output.
Nothing here is a causal flower pigment or pollinator selection estimate.
"""
from __future__ import annotations
import argparse
import importlib.util
import json
from pathlib import Path
from math import factorial
import numpy as np

BLOCKS={
 "Geography":["abs_lat","lon_sin","lon_cos"],
 "Moisture":["wc_bio12","wc_vapr_annual_mean"],
 "Temperature and radiation":["wc_bio5","wc_srad_annual_mean"],
 "Elevation":["wc_elevation_m"],
 "Soil5":["soil_pH","soil_SOC","soil_N","soil_clay","soil_available_water_proxy"],
}
BOOT_DRAWS=999
BOOT_SEED=202610161210
NULL_DRAWS=1999
NULL_SEED=202610161214
FEATURES=[f for features in BLOCKS.values() for f in features]
if len(FEATURES)!=13 or len(set(FEATURES))!=13: raise ValueError("Expected 13 distinct source covariates")
INDICES={mask:np.array([j for k,b in enumerate(BLOCKS) if mask&(1<<k)
                        for j,f in enumerate(FEATURES) if f in BLOCKS[b]], dtype=int)
         for mask in range(1,32)}


def source_module():
    path=Path(__file__).with_name('mois_full_biogeography_climate_soil_audit_20261015.py')
    spec=importlib.util.spec_from_file_location('mois_source_soil36_rebuild_for_global13',path)
    module=importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def load_source(photo_zip:Path,environment_zip:Path):
    source=source_module()
    _,soil,_,species,norm=source.original(photo_zip,environment_zip)
    if (len(species),species.genus.nunique(),int(species.n_photos.sum()),len(soil))!=(16801,1665,26163,31925):
        raise ValueError('Original soil36 cohort mismatch')
    labels,gi=np.unique(species.genus.to_numpy(),return_inverse=True)
    n=np.bincount(gi).astype(int)
    x=species[FEATURES].to_numpy(float)/np.array([norm[key] for key in FEATURES])
    y=species.white.to_numpy(float)
    if not np.isfinite(x).all() or not np.isfinite(y).all() or not np.all((0<=y)&(y<=1)):
        raise ValueError('Invalid original photographed species colours or environments')
    return source,labels,gi,n,x,y,norm,soil


def block_shapley(scores):
    out=np.zeros(5,float)
    for k in range(5):
        for mask in range(32):
            if mask & (1<<k): continue
            c=mask.bit_count()
            out[k]+=factorial(c)*factorial(4-c)/factorial(5)*(scores[mask|(1<<k)]-scores[mask])
    if not np.isclose(out.sum(),scores[31],atol=1e-11,rtol=0):
        raise AssertionError('Shapley shares do not sum to pooled total R2')
    return out


def pooled_fit(gi,n,x,y,weights):
    weights=np.asarray(weights,dtype=float)
    ybar=float(np.average(y,weights=weights))
    xbar=np.average(x,weights=weights,axis=0)
    xc=x-xbar;yc=y-ybar
    a=(xc.T*weights)@xc; b=xc.T@(weights*yc)
    sst=float(np.sum(weights*yc**2))
    scores=np.zeros(32,float)
    for mask,idx in INDICES.items():
        mat=a[np.ix_(idx,idx)];rhs=b[idx]
        scores[mask]=float(rhs@np.linalg.solve(mat,rhs)/sst)
    betas=100*np.linalg.solve(a,b)
    gmean=np.bincount(gi,weights=weights*y)/np.bincount(gi,weights=weights)
    genus_bss=float(np.sum(weights*(gmean[gi]-ybar)**2))
    within_y=y-gmean[gi]
    groupweights=np.bincount(gi,weights=weights)
    within_x=x-np.stack([np.bincount(gi,weights=weights*x[:,j])/groupweights
                         for j in range(x.shape[1])],axis=1)[gi]
    wa=(within_x.T*weights)@within_x; wb=within_x.T@(weights*within_y)
    wsst=float(np.sum(weights*within_y**2))
    w_r2=float(wb@np.linalg.solve(wa,wb)/wsst)
    gr2=genus_bss/sst
    share=block_shapley(scores)
    return {
      'weighting': 'one vote per source genus' if np.allclose(groupweights,groupweights[0]) else 'one vote per sampled source species',
      'photo_white_mean_of_species_proportions':ybar,
      'pooled_all13_R2':float(scores[31]),
      'pooled_all13_coefficient_pp_per_original_soil36_photo_site_SD':{key:float(betas[i]) for i,key in enumerate(FEATURES)},
      'pooled_block_only_R2':{k:float(scores[1<<i]) for i,k in enumerate(BLOCKS)},
      'pooled_block_unique_partial_R2':{k:float(scores[31]-scores[31^(1<<i)]) for i,k in enumerate(BLOCKS)},
      'pooled_block_shapley_fraction_total_white_variance':{k:float(share[i]) for i,k in enumerate(BLOCKS)},
      'same_cohort_genus_identity_only_R2':gr2,
      'same_cohort_within_genus_all13_R2_of_remaining_within_white_variance':w_r2,
      'same_cohort_genus_identity_plus_all13_R2':float(gr2+(1-gr2)*w_r2),
      'added_env_after_genus_R2_of_total_white_variance':float((1-gr2)*w_r2),
      'genus_identity_R2_is_in_sample_many_group_means_not_species_tree':True,
    }


def pooled_genus_moments(gi,n,x,y):
    m=len(n);d=x.shape[1];w=1/n[gi]
    ax=np.zeros((m,d,d)); sx=np.zeros((m,d));xy=np.zeros((m,d))
    sy=np.zeros(m);y2=np.zeros(m);within=np.zeros(m)
    np.add.at(ax,gi,w[:,None,None]*x[:,:,None]*x[:,None,:])
    np.add.at(sx,gi,w[:,None]*x)
    np.add.at(xy,gi,w[:,None]*x*y[:,None])
    np.add.at(sy,gi,w*y)
    np.add.at(y2,gi,w*y**2)
    np.add.at(within,gi,w*(y-sy[gi])**2)
    return ax,sx,xy,sy,y2,within


def bootstrap_genus_pooled(gi,n,x,y,point):
    ng=len(n);k=x.shape[1]
    xx,sx,xy,sy,yy2,within=pooled_genus_moments(gi,n,x,y)
    rng=np.random.default_rng(BOOT_SEED)
    draws=rng.multinomial(ng,np.repeat(1/ng,ng),size=BOOT_DRAWS)
    xx_t=(draws@xx.reshape(ng,-1)).reshape(BOOT_DRAWS,k,k)
    sx_t=draws@sx;xy_t=draws@xy;sy_t=draws@sy;y2_t=draws@yy2
    a=xx_t - sx_t[:,:,None]*sx_t[:,None,:]/ng
    b=xy_t-sx_t*sy_t[:,None]/ng
    sst=y2_t-sy_t**2/ng
    betas=np.linalg.solve(a,b[...,None])[...,0]
    scores=np.zeros((BOOT_DRAWS,32))
    for mask,idx in INDICES.items():
        mat=a[:,idx[:,None],idx[None,:]];rhs=b[:,idx]
        co=np.linalg.solve(mat,rhs[...,None])[...,0]
        scores[:,mask]=np.sum(rhs*co,axis=1)/sst
    R2=scores[:,31]
    shares=np.zeros((BOOT_DRAWS,5))
    for j in range(BOOT_DRAWS):shares[j]=block_shapley(scores[j])
    between=1-(draws@within)/sst
    pointbeta=np.array([point['pooled_all13_coefficient_pp_per_original_soil36_photo_site_SD'][name]
                        for name in FEATURES])
    bpp=100*betas;se=bpp.std(axis=0,ddof=1)
    maxT=np.max(abs((bpp-pointbeta)/se),axis=1)
    crit=float(np.quantile(maxT,.95))
    co={name:{'point_pp_per_original_soil36_photo_site_SD':float(pointbeta[i]),
              'genus_bootstrap_percentile95':np.quantile(bpp[:,i],[.025,.975]).tolist(),
              'family_13_variable_maxT95':[float(pointbeta[i]-crit*se[i]),float(pointbeta[i]+crit*se[i])]}
        for i,name in enumerate(FEATURES)}
    return {
      'draws':BOOT_DRAWS,'seed':BOOT_SEED,
      'unit':'resample full nominal genera; refit weighted pooled WLS including intercept per draw',
      'pooled_all13_R2_bootstrap_percentile95':np.quantile(R2,[.025,.975]).tolist(),
      'genus_identity_only_R2_bootstrap_percentile95':np.quantile(between,[.025,.975]).tolist(),
      'block_shapley_R2_bootstrap_percentile95':{name:np.quantile(shares[:,i],[.025,.975]).tolist() for i,name in enumerate(BLOCKS)},
      'simultaneous_13_conditional_coefficients_critical_maxT95':crit,
      'all13_coefficients':co,
      'scope':'retrospective, conditional 13-variable model, ordinary genus-bootstrap and within-13-family simultaneous limits; not study-wide model-selection correction',
    }


def genus_shuffle(gi,n,y,observed):
    rng=np.random.default_rng(NULL_SEED)
    null_g=[];null_s=[]
    wg=1/n[gi]
    for _ in range(NULL_DRAWS):
        yy=rng.permutation(y);mu_g=np.bincount(gi,weights=yy)/n
        m_g=float(np.average(yy,weights=wg))
        null_g.append(float(np.sum(wg*(mu_g[gi]-m_g)**2)/np.sum(wg*(yy-m_g)**2)))
        m_s=float(np.mean(yy))
        null_s.append(float(np.sum((mu_g[gi]-m_s)**2)/np.sum((yy-m_s)**2)))
    def record(arr,pt):
        values=np.array(arr)
        return {'observed':float(pt),'mean':float(values.mean()),
                'central_95_envelope':np.quantile(values,[.025,.975]).tolist(),
                'one_sided_p_plus_one':float((1+np.count_nonzero(values>=pt))/(NULL_DRAWS+1))}
    return {'draws':NULL_DRAWS,'seed':NULL_SEED,
            'shuffle':'permute observed species photographed-white fractions across fixed genus membership / genus sizes',
            'genus_equal':record(null_g,observed['genus_equal']),
            'species_equal':record(null_s,observed['species_equal']),
            'interpretation':'many fitted genus group means mechanically explain positive variance even under shuffle; conditional null is not phylogenetic signal or fitness'}


def make(photo_zip,env_zip):
    source,genera,gi,n,x,y,site_norm,soil=load_source(photo_zip,env_zip)
    gw=1/n[gi]
    equal_genus=pooled_fit(gi,n,x,y,gw)
    equal_species=pooled_fit(gi,n,x,y,np.ones(len(y)))
    # Additional *global* species-equal sensitivity: no >=3 species/genus gate.
    # Retains 20,786 unique original soil36-complete species and 31,925 photos.
    col={name:'mean' for name in FEATURES}
    col.update({'genus':'first','white':'mean','photo_id':'size'})
    all_species=soil.groupby('inat_taxon_id',sort=True).agg(col).reset_index()
    if (len(all_species),all_species.genus.nunique(),int(all_species.photo_id.sum()))!=(20786,4853,31925):
        raise ValueError('All soil36 original global taxon/genus/photo denominator drift')
    all_genera,all_id=np.unique(all_species.genus.to_numpy(),return_inverse=True)
    all_n=np.bincount(all_id)
    all_X=all_species[FEATURES].to_numpy(float)/np.array([site_norm[k] for k in FEATURES])
    all_y=all_species.white.to_numpy(float)
    all_global=pooled_fit(all_id,all_n,all_X,all_y,np.ones(len(all_y)))
    old_within=.006202366064066511
    if abs(equal_genus['same_cohort_within_genus_all13_R2_of_remaining_within_white_variance']-old_within)>1e-10:
        raise AssertionError('Original source validated within-genus all13 R2 does not replay')
    if abs(all_global['pooled_all13_R2']-.012806646762199388)>1e-9:
        raise AssertionError('Original unqualified 20,786-species pooled environmental R2 drift')
    if abs(equal_genus['pooled_all13_R2']-.01999939742187049)>1e-9:
        raise AssertionError('Original pooled whole-colour WLS R2 drift')
    report={
      'schema':'mois_source_identical_soil36_overall_white_13predictor_taxonomic_grain_v1',
      'date':'2026-10-16',
      'photo_source_archive_SHA256':source.PHOTO_SHA,'environment_source_archive_SHA256':source.ENV_SHA,
      'original_source_cohort':{'eligible_classifiable_photo_rows':int(26163),
                                'eligible_species':len(y),'eligible_genera':len(n),
                                'soil36_complete_photo_rows_before_genus_gate':31925},
      'white_target':'Per named species, observed proportion of original four-colour-classifiable flower photos labelled white; unclassifiable photos excluded; not a global population true-white census',
      'predictor_blocks_13_distinct':BLOCKS,
      'environment_scale':'Original soil36-complete classified photo-site SD of each source predictor; irrelevant to R2 but applies to pp coefficients',
      'source_photo_site_sd_for_all13':{key:site_norm[key] for key in FEATURES},
      'pooled_overall_white_genus_equal':equal_genus,
      'pooled_overall_white_species_equal':equal_species,
      'all_soil36_classifiable_species_without_congeneric_gate_species_equal':{
        'eligible_species':len(all_y),'observed_named_genera':len(all_genera),
        'source_classifiable_photo_rows':int(all_species.photo_id.sum()),
        'species_equal_observed_photo_white_mean':all_global['photo_white_mean_of_species_proportions'],
        'all13_in_sample_pooled_white_R2':all_global['pooled_all13_R2'],
        'all13_five_block_shapley':all_global['pooled_block_shapley_fraction_total_white_variance'],
        'BIO12_full13_slope_pp_per_original_source_soil_site_SD':all_global['pooled_all13_coefficient_pp_per_original_soil36_photo_site_SD']['wc_bio12'],
        'warning':'Only species-equal descriptive overall environmental fit; no >=3 congeneric species condition; no bootstrap family or independent occurrence-site replication for this sensitivity.'
      },
      'genus_equal_999_cluster_bootstrap':bootstrap_genus_pooled(gi,n,x,y,equal_genus),
      'fixed_genus_size_species_white_permutation_null':genus_shuffle(gi,n,y,{
        'genus_equal':equal_genus['same_cohort_genus_identity_only_R2'],
        'species_equal':equal_species['same_cohort_genus_identity_only_R2']}),
      'source_replay_consistency':{'original_13_variable_genus_centered_R2':old_within,
        'same_pop_genus_equal_genus_only_R2':equal_genus['same_cohort_genus_identity_only_R2'],
        'different_original_20725_species_genus_clustering_R2':0.343762761,
        'cannot_subtract_two_different_population_R2s':True},
      'limits':['Pooled 13-variable explanation is in-sample observational association and not prospective predictive validation.',
                'Between-genus dummy R2 uses 1,665 fitted group means and is not molecular phylogenetic heritability or causal species-pool effect.',
                'Genus-dummy R2, genus-centred 13-variable R2 and environment-only pooled R2 answer different but here same-source-target questions; genus+13 total R2 is a valid same-source weighted decomposition.',
                'Five-block Shapley increments are nonnegative allocations of overlapping predictive information, not causal effects or significance tests.',
                'Source photo labels and complete soil36 data are selected nonrandomly; many source photographs have unclassifiable floral hues.',
                'No source linked actual terrain slope/aspect/TPI/TRI, pigment production cost, pollinator behaviour or relative reproductive fitness.',
                'Some original photo/site species names, exact photo IDs and private locations are held only in the rights-controlled source ZIPs; never stored in public JSON.']
    }
    return report


def compare(a,b,path='root'):
    if isinstance(a,dict):
        if not isinstance(b,dict) or set(a)!=set(b):raise AssertionError('different JSON keys at '+path)
        for k,v in a.items():compare(v,b[k],path+'.'+k)
    elif isinstance(a,list):
        if not isinstance(b,list) or len(a)!=len(b):raise AssertionError('different array lengths '+path)
        for j,(v,z) in enumerate(zip(a,b)):compare(v,z,path+f'[{j}]')
    elif isinstance(a,(int,float)) and not isinstance(a,bool):
        if not isinstance(b,(int,float)) or abs(a-b)>1e-8:raise AssertionError((path,a,b))
    elif a!=b:raise AssertionError((path,a,b))


def self_test():
    rng=np.random.default_rng(53013)
    g=np.repeat(np.arange(8),12)
    n=np.bincount(g)
    X=rng.normal(size=(len(g),13))
    y=np.array([(j%4==0) for j in range(len(g))],float)
    x=pooled_fit(g,n,X,y,1/n[g])
    xc=np.column_stack([np.ones(len(g)),X]);ww=np.sqrt(1/n[g]);
    beta=np.linalg.lstsq(xc*ww[:,None],y*ww,rcond=None)[0]
    pred=xc@beta
    r2=1-np.sum((pred-y)**2/n[g])/np.sum((y-x['photo_white_mean_of_species_proportions'])**2/n[g])
    assert abs(r2-x['pooled_all13_R2'])<1e-10
    assert abs(sum(x['pooled_block_shapley_fraction_total_white_variance'].values())-r2)<1e-10
    assert abs(x['same_cohort_genus_identity_only_R2']+(1-x['same_cohort_genus_identity_only_R2'])*x['same_cohort_within_genus_all13_R2_of_remaining_within_white_variance']-x['same_cohort_genus_identity_plus_all13_R2'])<1e-10
    print('PASS synthetic source-free pooled 13-variable WLS/Shapley and between/within-genus identity; NOT source ecology')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--photo-zip',type=Path)
    p.add_argument('--env-zip',type=Path)
    p.add_argument('--output',type=Path)
    p.add_argument('--self-test',action='store_true')
    mode=p.add_mutually_exclusive_group()
    mode.add_argument('--write',action='store_true');mode.add_argument('--check',action='store_true')
    args=p.parse_args()
    if args.self_test:self_test();return
    if not args.photo_zip or not args.env_zip or not args.output or not (args.write or args.check):
        p.error('Original photo and environment ZIP paths, --output, --write/--check required')
    report=make(args.photo_zip,args.env_zip)
    if args.write:
        args.output.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
        print('WROTE source-verified pooled floral-white 13-variable global/genus grain result')
    else:
        compare(report,json.loads(args.output.read_text(encoding='utf8')))
        print('PASS original-source pooled and genus-controlled white 13-variable analysis')

if __name__=='__main__':main()