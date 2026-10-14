#!/usr/bin/env python3
"""FCP M1 ecological photo-hue direction: hold species constant before changing spatial control.

Three stages for each scope on a common named-genus bootstrap:
 A = genera which have sufficient within-cell species, but ALL their sampled species;
 B = only genus x coarse-cell groups with >=3 species, global-coordinate nuisance;
 C = EXACTLY THE SAME species as B, replacing broad coordinate nuisance with
     within-(genus x coarse-cell) conditioning.

The resulting A -> B change describes within-genus opportunity selection, and
B -> C contrasts model conditioning ON IDENTICAL species/hue/environment rows.
Neither difference identifies climate causality, microhabitat adaptation or selection.
Parent PR193 original data and results are not changed. Post-outcome diagnostic.
"""
from __future__ import annotations
import argparse
import itertools
import json
from pathlib import Path
import numpy as np
import pandas as pd
import fcp_abiotic_hue_direction_macro_20261013 as src

SCHEMA='fcp_macro_hue_same_species_spatial_conditioning_20261014_v1'
ANCHORS=src.ANCHOR_CLIMATE
N_HUE=len(src.CLASSES)
N_ENV=len(ANCHORS)
SEED=2026101422
REGIONS=('world','E90_E180')
STAGES=('full_photo_in_matched_genera','same_species_with_broad_geography','same_species_with_genus_cell')


def prepare_parent(photo_zip:Path,env_zip:Path):
    d,pools,_=src.original_source(photo_zip,env_zip)
    feats=tuple(v for k,b in src.BLOCKS.items() if k!='soil' for v in b)
    full=d.loc[pools['climate']].copy()
    normcol=(*src.GEO,*feats)
    mean=full[list(normcol)].mean()
    scale=full[list(normcol)].std(ddof=0).where(lambda z:z>=1e-12,1.0)
    return full,mean,scale,feats


def three_stages(species:pd.DataFrame,means:pd.Series,scales:pd.Series):
    """Source phenotype/environment rows cannot change from stage B to C."""
    x=species.copy()
    key=x.genus.astype(str)+'|'+x.cell.astype(str)
    nc=key.value_counts()
    chosen=x.loc[key.map(nc).ge(3)].copy()
    genuses=sorted(chosen.genus.unique())
    if len(genuses)<30:
        return None,{'status':'HOLD_FEWER_THAN_30_MATCHED_GENERA','n':len(genuses)}
    broader=x.loc[x.genus.isin(genuses)].copy()
    if not chosen.inat_taxon_id.isin(broader.inat_taxon_id).all():
        raise AssertionError('matched species absent from parent genus population')
    if chosen.inat_taxon_id.duplicated().any():
        raise AssertionError('species must be unique in each sampled sector')
    ca,ma=src.cluster_suff(broader,ANCHORS,means,scales,geos=True)
    cb,mb=src.cluster_suff(chosen,ANCHORS,means,scales,geos=True)
    cc,mc=src.cluster_suff(chosen,ANCHORS,means,scales,geos=False,cell_conditioned=True)
    if any(q is None for q in (ca,cb,cc)):
        return None,{'status':'HOLD_SOURCE_SUPPORT_DISAPPEARED','details':[ma,mb,mc]}
    gensets=[set(t[3] for t in q.values()) for q in (ca,cb,cc)]
    if gensets[0]!=gensets[1] or gensets[1]!=gensets[2]:
        raise AssertionError('stage models do not contain identical named genera')
    if ma['n_genera']!=mb['n_genera'] or mb['n_genera']!=mc['n_genera']:
        raise AssertionError('stage genus count drift')
    if mb['n_species']!=mc['n_species'] or mb['n_photo_rows']!=mc['n_photo_rows']:
        raise AssertionError('matched species or photo opportunity changed with conditioning')
    return (ca,cb,cc),{
        'status':'OK','n_matched_named_genera':len(genuses),
        'n_original_species_all_genera':int(ma['n_species']),
        'n_matched_species':int(mb['n_species']),
        'n_original_photos_all_genera':int(ma['n_photo_rows']),
        'n_matched_original_photos':int(mb['n_photo_rows']),
        'n_geographic_source_cells':int(mc['n_geo_cells']),
        'n_nested_cell_groups':int(mc['n_groups']),
        'only_cell_adjustment_changes_between_last_two_stages':True,
        'original_colour_hues_unchanged':True}


def stats_for_genus(stage:dict,genus_names:list[str],weight:str,cell:bool)->tuple[np.ndarray,np.ndarray]:
    """Stage sufficient stats per genus/anchor, preserving exactly the parent PR193 estimator.
    Matrix shape=(genera, axes, 4,4), response=(genera,axes,4,4) for broad geo;
    for cell-conditioned stage, packed in the [3,3] coefficient and last row of RHS.
    """
    ng=len(genus_names)
    gm=np.zeros((ng,N_ENV,4,4))
    gr=np.zeros((ng,N_ENV,4,N_HUE))
    ix={genus:i for i,genus in enumerate(genus_names)}
    genus_species={g:0 for g in genus_names}
    for _,(m,r,n,g) in stage.items():
        if g not in ix:raise ValueError('unexpected genus')
        genus_species[g]+=n
    for _,(m,r,n,g) in stage.items():
        k=ix[g]
        ww=1.0 if weight=='species_equal' else 1.0/genus_species[g]
        for j in range(N_ENV):
            if cell:
                gm[k,j,3,3]+=ww*m[j,j]
                gr[k,j,3,:]+=ww*r[j,:]
            else:
                cols=[0,1,2,3+j]
                gm[k,j]+=ww*m[np.ix_(cols,cols)]
                gr[k,j]+=ww*r[cols,:]
    if any(n<3 for n in genus_species.values()):raise ValueError('genus with <3 plants')
    return gm,gr


def fit(m:np.ndarray,r:np.ndarray,draws:np.ndarray,cell:bool):
    M=m.sum(axis=0);R=r.sum(axis=0)
    mm=np.einsum('bg,gjkl->bjkl',draws,m,optimize=True)
    rr=np.einsum('bg,gjkl->bjkl',draws,r,optimize=True)
    if cell:
        if np.any(M[:,3,3]<=1e-10) or np.any(mm[:,: ,3,3]<=1e-10):
            raise ValueError('no within-cell environmental variation')
        p=R[:,3,:]/M[:,3,3,None]
        bs=rr[:,:,3,:]/mm[:,:,3,3,None]
    else:
        if np.max([np.linalg.cond(z) for z in M])>1e10:raise ValueError('ill-conditioned broad geo')
        p=np.linalg.solve(M,R)[:,3,:]
        bs=np.linalg.solve(mm,rr)[:,:,3,:]
    assert np.isfinite(p).all() and np.isfinite(bs).all()
    if not np.allclose(p.sum(axis=-1),0,atol=1e-9):raise ValueError('four categorical slopes not sum 0')
    return p,bs


def max_t_envelopes(p:np.ndarray,boot:np.ndarray)->tuple[np.ndarray,float]:
    dev=boot-p[None,:,:]
    sd=dev.std(axis=0,ddof=1)
    valid=sd>1e-12
    if not valid.any():
        return np.stack((p,p),axis=-1),0.0  # an exactly zero contrast is not significant
    t=np.max(np.abs(dev[:,valid]/sd[valid]),axis=1)
    crit=float(np.quantile(t,.95))
    return np.stack((p-crit*sd,p+crit*sd),axis=-1),crit


def paired_compare(stages:tuple[dict,dict,dict],nboot:int,seed:int,weight:str)->dict:
    names=sorted({v[3] for v in stages[1].values()})
    if len(names)<30:return {'status':'HOLD_INSUFFICIENT_GENERA'}
    draws=np.random.default_rng(seed).multinomial(len(names),np.repeat(1/len(names),len(names)),size=nboot)
    vals=[]
    for i,stage in enumerate(stages):
        m,r=stats_for_genus(stage,names,weight,cell=i==2)
        try:p,boot=fit(m,r,draws,cell=i==2)
        except (ValueError,np.linalg.LinAlgError) as exc:return {'status':'HOLD_ILL_CONDITIONED','reason':str(exc)[:150]}
        vals.append((p,boot))
    estimates={}
    for i,label in enumerate(STAGES):
        point,boot=vals[i];lo,hi=np.quantile(boot,[.025,.975],axis=0)
        estimates[label]={'point':point.tolist(),'CI95_percentile':np.stack((lo,hi),axis=-1).tolist()}
    contrasts={}
    for label,i,j in [('opportunity_within_matched_genera',0,1),('spatial_adjustment_on_same_species',1,2),('combined_coverage_and_spatial_adjustment',0,2)]:
        x=vals[i][0]-vals[j][0]
        xx=vals[i][1]-vals[j][1]
        ci=np.quantile(xx,[.025,.975],axis=0)
        maxci,max_t=max_t_envelopes(x,xx)
        contrasts[label]={
            'definition':STAGES[i]+' MINUS '+STAGES[j],
            'point':x.tolist(),
            'paired_genus_bootstrap_CI95':np.stack(ci,axis=-1).tolist(),
            'family_maxT_95_CI':maxci.tolist(),
            'family_maxT_95_critical':max_t,
            'n_pairs_in_family':N_ENV*N_HUE,
            'n_maxT_significant':int(((maxci[:,:,0]>0)|(maxci[:,:,1]<0)).sum())}
    return {'status':'OK','n_bootstrap':nboot,'n_named_genera':len(names),
            'same_genus_draws_all_three_stages':True,'species_and_outcome_identical_in_last_two_stages':True,
            'estimate':estimates,'stage_differences':contrasts}


def as_rows(results):
    rows=[]
    for scope,rec in results['regions'].items():
        for weight,fit in rec.get('fit',{}).items():
            if fit['status']!='OK':continue
            for label,st in fit['estimate'].items():
                for j,anchor in enumerate(ANCHORS):
                    for k,hue in enumerate(src.CLASSES):
                        rows.append({'region':scope,'weight':weight,'type':'stage','comparison':label,
                            'anchor':anchor,'hue':hue,'point_pp':100*st['point'][j][k],
                            'CI95_lo_pp':100*st['CI95_percentile'][j][k][0],
                            'CI95_hi_pp':100*st['CI95_percentile'][j][k][1],
                            'maxT95_lo_pp':np.nan,'maxT95_hi_pp':np.nan})
            for label,st in fit['stage_differences'].items():
                for j,anchor in enumerate(ANCHORS):
                    for k,hue in enumerate(src.CLASSES):
                        rows.append({'region':scope,'weight':weight,'type':'paired_difference','comparison':label,
                            'anchor':anchor,'hue':hue,'point_pp':100*st['point'][j][k],
                            'CI95_lo_pp':100*st['paired_genus_bootstrap_CI95'][j][k][0],
                            'CI95_hi_pp':100*st['paired_genus_bootstrap_CI95'][j][k][1],
                            'maxT95_lo_pp':100*st['family_maxT_95_CI'][j][k][0],
                            'maxT95_hi_pp':100*st['family_maxT_95_CI'][j][k][1]})
    return pd.DataFrame(rows)


def run(photo_archive:Path,environment_archive:Path,nboot=999):
    data,mu,sd,features=prepare_parent(photo_archive,environment_archive)
    out={'schema':SCHEMA,'date':'2026-10-14','role':'retrospective existing-photo sampling and spatial support decomposition; NOT new confirmatory biological evidence',
         'photo_source_sha256':src.SOURCE_PHOTO_SHA,'abiotic_source_sha256':src.SOURCE_ENV_SHA,
         'full_photo_opportunities':100543,'full_photo_colour_classifiable':45953,
         'climate26_eligible_source_photos':int(len(data)),
         'climate26_eligible_species':int(data.inat_taxon_id.nunique()),
         'frozen_anchors':list(ANCHORS),'frozen_four_hues':list(src.CLASSES),
         'study_scope':list(REGIONS),'study_cohorts_have_distinct_taxa':True,
         'regions':{},'hard_nonclaims':[
            'The matched B/C species, photo hues and environment measurements are identical; its contrast isolates MODEL AND SPATIAL CONDITIONING differences, not causal ecology.',
            'Matching discards species lacking congeneric same-coarse-cell replication, and reduces environmental range; bootstrap intervals cannot repair MNAR.',
            'Coarse equal-area cell controls are not 50km geodesic local controls nor a Gaussian process or dated phylogeny.',
            'Photo classifier categories are not confirmed petal anthocyanins, within-population inherited colour morphs or experimental pollinator fitness.',
            'Scopes, axes, weights and posthoc interpretation were chosen after original PR193 outcome; no independent prospective biology.'
          ]}
    for region in REGIONS:
        sub=data if region=='world' else data.loc[data.sector_index.eq(3)]
        sp=src.species_frame(sub,features)
        tr,meta=three_stages(sp,mu,sd)
        rec={'status':meta['status'],'support':meta,'fit':{}}
        if tr is not None:
            for weight in ('species_equal','genus_equal'):
                rec['fit'][weight]=paired_compare(tr,nboot,SEED+len(region)*11,weight)
        out['regions'][region]=rec
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--photo-archive',type=Path,required=True)
    p.add_argument('--environment-archive',type=Path,required=True)
    p.add_argument('--bootstrap',type=int,default=999)
    p.add_argument('--outdir',type=Path,required=True)
    a=p.parse_args()
    r=run(a.photo_archive,a.environment_archive,a.bootstrap)
    target=a.outdir;target.mkdir(exist_ok=True,parents=True)
    (target/'result.json').write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    rows=as_rows(r);rows.to_csv(target/'matched_spatial_support_differences.csv',index=False)
    print('RESULT',json.dumps({'support':{k:v['support'] for k,v in r['regions'].items()},
      'n_full_table_rows':len(rows)},ensure_ascii=False))
if __name__=='__main__':main()
