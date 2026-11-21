#!/usr/bin/env python3
"""Test whether nonwhite 8-photo-hue clustering survives fixing white IDs
and six broad present-day geographic species-pool affinity zones.

Reproduces published 199-permutation white-fixed null (p=.53) identically.
Isolates exchangeable-species support from added region conditioning.
Source hues are NOT independently verified petal colours.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from scripts.analysis.mois_global_hue9_community_assembly_genus_null_20261118 import (
    signed_original,baseline,empirical_entropy,SEED,GBIF_SHA,ENV_SHA,check,sha
)
from scripts.analysis.mois_global_hue9_regional_pool_conditional_null_20261121 import (
    geographic_regions,genus_breadth_codes,species_affinity,varying_groups,
    ZONE_NAMES,THRESH
)
from scripts.analysis.mois_nonwhite8_assembly_white_fixed_null_20261118 import count_chromatic

PERMS=199
SCHEMA='mois_global_photo8_white_fixed_georegion_conditioned_20261121_v1'

def calculate_groups(x,sp,zone):
    hue=x['hue']
    primary,occupied,spg,cov=genus_breadth_codes(x,sp)
    affinity,share,cross=species_affinity(x,zone)
    # Original study's exact identical within-genus×occupied-grid breadth hue
    # interchange among only nonwhite species: the reference null.
    old=[]
    for group in x['variable_groups']:
        chrom=group[hue[group]>0]
        if len(chrom)>1 and len(set(hue[chrom]))>1:old.append(chrom)
    chrom_spp=np.flatnonzero(occupied & (hue>0))
    reg=varying_groups(primary*7+affinity,chrom_spp,hue)
    mobile=np.unique(np.concatenate(reg))
    matched=varying_groups(primary,mobile,hue)
    check(set(np.concatenate(matched))==set(mobile),
          'Matched chromatic8 group mobility does not cover same species IDs')
    check(len(old)>300 and sum(map(len,old))>1000 and
          len(reg)>100 and sum(map(len,reg))>500,
          'No independent nonwhite 8-hue interchanges in source')
    return {
       'original_genus_breadth_whitefixed8':old,
       'matched_exchangeable_whitefixed8':matched,
       'sixzone_genus_breadth_whitefixed8':reg
    },spg,affinity


def per_model(x,zone,model_groups,spg):
    hue=x['hue']
    nonwhite=hue[x['species']]>0
    cell=x['cell'][nonwhite]
    species=x['species'][nonwhite]
    n=np.bincount(cell,minlength=x['ngrids'])
    check(int(n.sum())>1500000,'Original nonwhite source richness unexpectedly low')
    original=count_chromatic(x,hue,cell,species)
    obs=empirical_entropy(original)
    masks={'nonwhite_richness_ge30':n>=30,
           'nonwhite_richness_ge100':n>=100}
    check(int(masks['nonwhite_richness_ge30'].sum())==10388 and
          int(masks['nonwhite_richness_ge100'].sum())==6187,
          'Earlier white fixed 8-hue denominator not reproduced')
    for k in range(6):masks[ZONE_NAMES[k]+'_nonwhite_ge30']=(n>=30)&(zone==k)
    idx={k:np.flatnonzero(v) for k,v in masks.items()}
    for k in ZONE_NAMES:
        check(len(idx[k+'_nonwhite_ge30'])>=100,
              'One source geography zone has insufficient nonwhite support')
    results={};spatial={}
    for key,groups in model_groups.items():
        rng=np.random.default_rng(SEED+803)
        nums={name:np.empty(PERMS,float) for name in masks}
        exp=np.zeros(x['ngrids'],float)
        for t in range(PERMS):
            perm=hue.copy()
            for g in groups:
                perm[g]=hue[g][rng.permutation(len(g))]
            check(np.array_equal(perm==0,hue==0),'White species IDs moved in a white-fixed null')
            hb=empirical_entropy(count_chromatic(x,perm,cell,species))
            exp+=hb
            for name,index in idx.items():nums[name][t]=float(hb[index].mean())
            if (t+1)%50==0:print('MOIS_WHITE_FIXED_EIGHT_REGIONAL_PERM',key,t+1,flush=True)
        exp/=PERMS
        stat={}
        for name,ii in idx.items():
            v=nums[name]
            real=float(obs[ii].mean());mu=float(v.mean());delta=real-mu
            stat[name]={
                'n_real_cells':len(ii),
                'observed_chromatic_H':real,
                'null_expected_chromatic_H':mu,
                'observed_minus_expected_H':delta,
                'global_whole_species_permutation_p_two_sided':
                    float((1+int(np.count_nonzero(np.abs(v-mu)>=abs(delta)-1e-12)))/(PERMS+1)),
                'global_whole_species_permutation_p_lower_tail':
                    float((1+int(np.count_nonzero(v<=real+1e-12)))/(PERMS+1)),
                'null_95pct_interval':np.percentile(v,[2.5,97.5]).tolist()
            }
        results[key]={
            'n_exchangeable_nonwhite_species':int(sum(map(len,groups))),
            'n_exchangeable_strata':len(groups),
            'n_genera_with_variable_nonwhite_strata':len(set(int(spg[g[0]]) for g in groups)),
            'stats':stat}
        spatial[key]=exp
        print('MOIS_WHITEFIXED_REGION_HUE8_SUMMARY',key,
              json.dumps({'support':{k:v for k,v in results[key].items() if k!='stats'},
                     'global':stat['nonwhite_richness_ge30']}),flush=True)
    # Both source and seed reproduce older full photo8 status result.
    old=results['original_genus_breadth_whitefixed8']['stats']['nonwhite_richness_ge30']
    check(abs(old['observed_minus_expected_H']-(-.005570))<.000003
          and abs(old['global_whole_species_permutation_p_two_sided']-.53)<1e-12,
          'Previously signed white fixed genus-breadth 8-hue source null differs')
    return results,obs,n,spatial


def run(args):
    dest=Path(args.out)
    check(not dest.exists(),'Refuse overwriting source backed floral QC output')
    pair,sp=signed_original(args.gbif,args.env)
    x=baseline(pair,sp,null_model='genus_breadth')
    zone,lon,lat=geographic_regions(x)
    models,spg,affinity=calculate_groups(x,sp,zone)
    res,obs,counts,exp=per_model(x,zone,models,spg)
    dest.mkdir(parents=True)
    coords=np.asarray(x['cell_key'].tolist(),dtype=int)
    c=pd.DataFrame({
       'grid_x':coords[:,0],'grid_y':coords[:,1],
       'macrozone':np.asarray(ZONE_NAMES)[zone],
       'n_nonwhite_source_photo9_species':counts,
       'observed_nonwhite8_Shannon':obs
    })
    for key in models:
        c['null_expected_nonwhite8_H_'+key]=exp[key]
        c['observed_minus_null_nonwhite8_H_'+key]=obs-exp[key]
    file=dest/'global_100km_nonwhite8_whitefixed_geography_conditioned_cells.csv.gz'
    c.to_csv(file,index=False,compression='gzip')
    receipt={
      'schema':SCHEMA,'source_GBIF_zip_SHA256':GBIF_SHA,
      'source_ENV_zip_SHA256':ENV_SHA,
      'species_source':len(sp), 'real_speciesx100km_source':len(pair),
      'n_global_cells':x['ngrids'],
      'null_permutations':PERMS,'seed':SEED+803,
      'white_species_photo9_identity_and_occupied_range_fixed':True,
      'nonwhite_species_geographic_range_fixed':True,
      'six_macrozone_strata_are_not_biogeographical_realms':True,
      'conditional_region_affinity_min_fraction':THRESH,
      'null_model':res,
      'no_original_true_petal_colour_validation':True,
      'no_species_wide_colour_adaptation_inference':True,
      'no_global_hue_registry_change':True,
      'source_photo9_file_sha256':sha(file)
    }
    summary=dest/'global_nonwhite8_species_region_whitefixed_null_summary.json'
    summary.write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
    (dest/'receipt.json').write_text(json.dumps({
        'schema':SCHEMA,'analysis_summary_SHA256':sha(summary),
        'source_global_8chromatic_cells_SHA256':sha(file),
        'source_GBIF_SHA256':GBIF_SHA,'source_ENV_SHA256':ENV_SHA
    },indent=2)+'\n')
    print('MOIS_GLOBAL_WHITEFIXED_8_REGIONAL_CONDITIONAL_DONE',json.dumps({
       k:{'mobilized_species':v['n_exchangeable_nonwhite_species'],
          'global_30':v['stats']['nonwhite_richness_ge30']}
       for k,v in res.items()},sort_keys=True),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser()
    for k in ('gbif','env','out'):p.add_argument('--'+k,required=True)
    run(p.parse_args())