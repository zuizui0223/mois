#!/usr/bin/env python3
"""Falsification of global photo9 aggregation by coarse *species affinity* geography.

Same source 37,104 species, 5,493,931 real species×100km GBIF grid pairs,
same seed, same global species-hue shuffle conditional on genus×range-breadth.
A geography-stratified null holds broad centroid-like range affinity fixed;
matched-mobility null separates loss of exchangeable species from geography.

NEITHER a phylogenetic null NOR evidence of natural selection/petal truth.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from pyproj import Transformer
from scripts.analysis.mois_global_hue9_community_assembly_genus_null_20261118 import (
    signed_original,baseline,counts_by_cell,empirical_entropy,HUES,SEED,
    GBIF_SHA,ENV_SHA,check,sha
)

SCHEMA='mois_photo9_world_species_affinity_region_conditioned_null_20261121_v1'
# Published pre-analysis operations: six coarse zones, not biological biomes.
ZONE_NAMES=('Americas_N','Americas_S','Africa_Eurasia_W_N',
            'Africa_Eurasia_W_S','Asia_Oceania_E_N','Asia_Oceania_E_S')
REGION_STRATA_NAMES=ZONE_NAMES+('multi_zone_35pct_plus',)
THRESH=.65
PERMS=199


def geographic_regions(x):
    z=np.asarray(x['cell_key'].tolist(),dtype=float)
    check(z.shape==(18577,2),'Unexpected actual GBIF grid centers')
    # Source grid_x/grid_y are floor((EPSG:6933 metres)/100000).
    xx=(z[:,0]+.5)*100000
    yy=(z[:,1]+.5)*100000
    lon,lat=Transformer.from_crs(6933,4326,always_xy=True).transform(xx,yy)
    # The signed original 18,577 grid includes exactly 27 extreme polar
    # EPSG:6933 Y edge cells (-74 or +73); longitude is still defined but
    # inverse latitude is undefined outside the equal-area projection domain.
    # Keep their geographic latitude genuinely missing; the sign of the
    # original equal-area northing is sufficient only to classify hemisphere.
    invalid_lat=~np.isfinite(lat)
    check(int(invalid_lat.sum())==27 and np.isfinite(lon).all()
          and (np.abs(lon)<=181).all()
          and (np.abs(lat[~invalid_lat])<=90).all()
          and np.isin(z[invalid_lat,1],[-74,73]).all(),
          'Unexpected projection invalid geography beyond 27 polar edge cells')
    # 3 longitudinal partitions crossed with hemisphere by signed northing.
    lonpart=np.where(lon < -30,0,np.where(lon<60,1,2)).astype(np.int8)
    hemi=(yy<0).astype(np.int8)
    check(np.array_equal(hemi[~invalid_lat],(lat[~invalid_lat]<0).astype(np.int8)),
          'Signed northing hemisphere disagrees with known geographic latitude')
    zone=(lonpart*2+hemi).astype(np.int8)
    check(set(np.unique(zone))==set(range(6)),
          'Six coarse source zones lack any actual grid cells')
    return zone,lon,lat


def genus_breadth_codes(x,species):
    spg,_=pd.factorize(species.genus.astype(str),sort=True)
    cov=np.bincount(x['species'],minlength=len(species))
    occupied=cov>0
    broad=np.digitize(cov,[2,6,21,101]).astype(np.int16)
    primary=spg.astype(np.int64)*5+broad
    check(int(occupied.sum())==36758 and
          int((~occupied).sum())==346,
          'Original species without GBIF footprint changed')
    return primary,occupied,spg,cov


def species_affinity(x,zone):
    n=len(x['hue'])
    # Every species contributes one vote per occupied grid (NOT GBIF photos).
    ix=x['species'].astype(np.int64)
    cross=np.bincount(ix*6+zone[x['cell']].astype(np.int64),
                      minlength=n*6).reshape(n,6)
    ngrid=cross.sum(axis=1)
    top=cross.argmax(axis=1)
    share=cross.max(axis=1)/np.maximum(ngrid,1)
    region=np.where(share>=THRESH,top,6).astype(np.int8)
    region=np.where(ngrid>0,region,7)
    check(np.array_equal(ngrid,np.bincount(ix,minlength=n)),'Footprint counts altered')
    return region,share,cross


def varying_groups(codes,available,hue):
    ix=np.argsort(codes[available],kind='stable')
    ordered=available[ix]
    c=codes[ordered]
    offsets=np.r_[0,1+np.flatnonzero(np.diff(c))]
    ends=np.r_[offsets[1:],len(ordered)]
    variable=[ordered[a:b] for a,b in zip(offsets,ends)
              if b-a>=2 and len(set(hue[ordered[a:b]]))>=2]
    return variable


def stats_for_groups(groups,spg):
    return {
       'n_variable_strata':len(groups),
       'n_exchangeable_species':int(sum(len(a) for a in groups)),
       'n_genera_in_variable_strata':int(len(set(int(spg[a[0]]) for a in groups))),
       'n_min_group_size':int(min(len(a) for a in groups)),
       'n_max_group_size':int(max(len(a) for a in groups))
    }


def nulls(x,species,model,zone,region,base,occupied,spg,affinity):
    available=np.flatnonzero(occupied)
    hue=x['hue']
    if model=='genus_breadth':
        groups=x['variable_groups']
        check(sum(map(len,groups))==25325 and len(groups)==3430,
              'Original strict genus×occupancy null exchangeable IDs changed')
    elif model=='genus_breadth_hemisphere':
        # A separate less restrictive geographic influence test.
        hemi=np.where(np.isin(region,[0,2,4]),0,
                      np.where(np.isin(region,[1,3,5]),1,2))
        groups=varying_groups(base*3+hemi,available,hue)
    elif model in ('genus_breadth_sixzones','matched_exchangeable_genus_breadth'):
        geo_groups=varying_groups(base*7+region,available,hue)
        if model=='genus_breadth_sixzones':
            groups=geo_groups
        else:
            # For a fair baseline with exactly same movable species as geo model,
            # group them by genus×breadth but DO NOT condition on macrozone.
            mobility=np.unique(np.concatenate(geo_groups))
            groups=varying_groups(base,mobility,hue)
            check(sum(map(len,groups))==len(mobility),
                  'Matched geographic mobility is not identical')
            check(set(np.concatenate(groups))==set(mobility),
                  'Matched geographic null changes movable species')
    else: raise ValueError('Unregistered model')
    check(len(groups)>100 and sum(map(len,groups))>500,
          'Coarse geographic pool null has inadequate shuffling support')
    return groups


def test_models(x,spg,model_groups,zone,permutations):
    original=counts_by_cell(x,x['hue'])
    h_obs=empirical_entropy(original)
    ncell=x['cell_n']
    masks={'all_richness_ge30':ncell>=30,'all_richness_ge100':ncell>=100}
    for k in range(6):masks[ZONE_NAMES[k]+'_richness_ge30']=(ncell>=30)&(zone==k)
    mids={k:np.flatnonzero(v) for k,v in masks.items()}
    check(int(masks['all_richness_ge30'].sum())==12728
          and int(masks['all_richness_ge100'].sum())==9827,
          'Previously signed global inclusion masks changed')
    # Equal macrozone weight is a second pre-registered summary; does not
    # let a high-cell-count continent dominate the null evidence.
    check(min(len(mids[k+'_richness_ge30']) for k in ZONE_NAMES)>100,
          'One macrozone lacks enough real source cells')
    observations={key:float(np.mean(h_obs[ii])) for key,ii in mids.items()}
    observations['six_macrozone_equal_weight_mean']=float(np.mean([
        observations[k+'_richness_ge30'] for k in ZONE_NAMES]))
    cellstats={}
    summaries={}
    for label,groups in model_groups.items():
        rng=np.random.default_rng(SEED)
        sums=np.zeros(x['ngrids'],dtype=np.float64)
        vals={name:np.empty(permutations,dtype=float) for name in observations}
        for i in range(permutations):
            hh=x['hue'].copy()
            # One global shuffled species label shared by every occupied grid.
            for g in groups:hh[g]=x['hue'][g][rng.permutation(len(g))]
            hp=empirical_entropy(counts_by_cell(x,hh))
            sums+=hp
            for key,ii in mids.items():vals[key][i]=float(np.mean(hp[ii]))
            vals['six_macrozone_equal_weight_mean'][i]=np.mean([
                vals[k+'_richness_ge30'][i] for k in ZONE_NAMES])
            if (i+1)%50==0:
                print('MOIS_REGIONAL_NULL_PROGRESS',label,i+1,flush=True)
        avg=sums/permutations
        check(np.isfinite(avg).all(),'Nonfinite expected Shannon')
        records={}
        for key,obs in observations.items():
            null=vals[key]
            mu=float(np.mean(null))
            diff=obs-mu
            p=(1+np.count_nonzero(np.abs(null-mu)>=abs(diff)-1e-13))/(permutations+1)
            p_left=(1+np.count_nonzero(null<=obs+1e-12))/(permutations+1)
            records[key]={
                'n_real_source_cells':int(len(mids[key])) if key in mids else int(masks['all_richness_ge30'].sum()),
                'observed_mean_H':obs,'null_mean_H':mu,
                'observed_minus_null_H':diff,
                'null_empirical_95pct_interval':np.percentile(null,[2.5,97.5]).tolist(),
                'global_species_label_two_sided_permutation_p':float(p),
                'global_species_label_lower_tail_p':float(p_left)
            }
        summaries[label]={
            'group_support':stats_for_groups(groups,spg),
            'statistics':records
        }
        cellstats[label]=avg
        print('MOIS_REGIONAL_HUE9_REAL_NULL_SUMMARY',label,
              json.dumps({'support':summaries[label]['group_support'],
                  'ge30':records['all_richness_ge30'],
                  'equal_macrozone':records['six_macrozone_equal_weight_mean']},
                  sort_keys=True),flush=True)
    # Use exact published 199 permutation design and reference as a source gate.
    base=summaries['genus_breadth']['statistics']['all_richness_ge30']
    check(abs(base['observed_mean_H']-1.199973)<2e-6
          and abs(base['null_mean_H']-1.230773)<2e-6
          and abs(base['observed_minus_null_H']-(-.030799))<2e-6
          and abs(base['global_species_label_two_sided_permutation_p']-.010)<1e-12,
          'Previously signed 199-permutation source result was not reproduced')
    return summaries,original,cellstats


def main(args):
    root=Path(args.out)
    check(not root.exists(),'No rewriting frozen signed results allowed')
    pairs,species=signed_original(args.gbif,args.env)
    x=baseline(pairs,species,null_model='genus_breadth')
    zone,lon,lat=geographic_regions(x)
    base,occupied,spg,cov=genus_breadth_codes(x,species)
    region,share,cross=species_affinity(x,zone)
    models=('genus_breadth','genus_breadth_hemisphere',
            'matched_exchangeable_genus_breadth','genus_breadth_sixzones')
    groups={m:nulls(x,species,m,zone,region,base,occupied,spg,share)
            for m in models}
    # Strict mobility identity for geo vs matched baseline.
    check(set(np.concatenate(groups['genus_breadth_sixzones']))==
          set(np.concatenate(groups['matched_exchangeable_genus_breadth'])),
          'Mobility matched groups not sharing identical species support')
    summary,observed,mu=test_models(x,spg,groups,zone,args.permutations)
    root.mkdir(parents=True)
    coords=np.asarray(x['cell_key'].tolist(),dtype=int)
    out=pd.DataFrame({
         'grid_x':coords[:,0], 'grid_y':coords[:,1],
         'center_lat':lat, 'center_lon':lon,
         'macrozone':np.asarray(ZONE_NAMES)[zone],
         'n_species':x['cell_n'],
         'observed_hue9_Shannon':empirical_entropy(observed)
    })
    for name in models:
        out['expected_Shannon_'+name]=mu[name]
        out['observed_minus_null_'+name]=out['observed_hue9_Shannon']-mu[name]
    output=root/'MOIS_global_photo9_100km_conditional_geography_18577_cells.csv.gz'
    out.to_csv(output,index=False,compression='gzip')
    affinity_out=root/'MOIS_region_affinity_accepted_species_summary.csv.gz'
    pd.DataFrame({
        'species_key':species.key.astype(str),
        'genus':species.genus.astype(str),
        'source_photo9_hue':species.source_hue9.astype(str),
        'GBIF_100km_occupied_cells':cov,
        'species_macroregion_affinity':np.array(list(REGION_STRATA_NAMES)+['absent_no_GBIF'])[region],
        'fraction_occurrences_in_dominant_macrozone':share
    }).to_csv(affinity_out,index=False,compression='gzip')
    receipt={
       'schema':SCHEMA,'source_GBIF_original_ZIP_SHA256':GBIF_SHA,
       'source_species_ENV_original_ZIP_SHA256':ENV_SHA,
       'source_real_accepted_species':len(species),
       'source_speciesx100km_cells':len(pairs),
       'source_distinct_100km_cells':x['ngrids'],
       'null_permutations':args.permutations,
       'random_seed':SEED,
       'zone_definition':'EPSG6933 100km center longitude; bands <-30, -30..60, >=60; N/S by source projected northing sign. 27 polar edge-cell center latitudes are undefined and stay missing. Not ecoregions',
       'n_source_EPSG6933_polar_edge_cells_with_undefined_latitude':int(np.count_nonzero(~np.isfinite(lat))),
       'species_affinity_stratum_definition':'>=65% of occupied present GBIF grid cells in one of six zones, else multizone; group by species genus, five occupancy breadth bins, and macrozone affinity',
       'region_counts_species':{
           str(k):int((region==k).sum()) for k in range(8)},
       'region_counts_cells':{
           ZONE_NAMES[k]:int((zone==k).sum()) for k in range(6)},
       'null_models':summary,
       'source_photo9_is_not_validated_petal_colour':True,
       'GBIF_occurrence_is_not_native_wild_presence_or_abundance':True,
       'regional_affinity_from_present_GBIF_not_historical_phylogeography':True,
       'species_full_GBIF_occupancy_footprint_fixed_in_every_null':True,
       'matched_exchangeable_species_set_identical_between_geo_and_unrestricted':True,
       'weaker_null_can_lose_significance_due_to_fixed_group_diversity':True,
       'biological_adaptation_or_mechanism_not_identified':True,
       'old_MOIS_global_hue_label_registry_not_changed':True,
       'outputs_SHA256':{
           output.name:sha(output), affinity_out.name:sha(affinity_out)
       }
    }
    report=root/'ecological_species_affinity_conditional_photo9_permutation_summary.json'
    report.write_text(json.dumps(receipt,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    (root/'receipt.json').write_text(json.dumps({
        'schema':SCHEMA,'summary_SHA256':sha(report),
        'outputs_SHA256':receipt['outputs_SHA256'],
        'source_GBIF_SHA256':GBIF_SHA,'source_ENV_SHA256':ENV_SHA
    },indent=2,allow_nan=False)+'\n')
    print('MOIS_REGIONAL_NULL_ALL_MODELS_DONE',json.dumps({
       k:{'n_movable':v['group_support']['n_exchangeable_species'],
          'global':v['statistics']['all_richness_ge30'],
          'equal_regions':v['statistics']['six_macrozone_equal_weight_mean']}
       for k,v in summary.items()},sort_keys=True),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--gbif',required=True)
    parser.add_argument('--env',required=True)
    parser.add_argument('--out',required=True)
    parser.add_argument('--permutations',type=int,default=PERMS)
    args=parser.parse_args()
    check(args.permutations==PERMS,'Must use 199 frozen whole-catalog permutations')
    main(args)