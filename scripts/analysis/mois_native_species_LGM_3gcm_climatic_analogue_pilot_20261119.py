#!/usr/bin/env python3
"""Audited species-level LGM climatic analogues in modern WCVP-native regions.

This is *potential climatic matching*, NOT historical occupancy, refugia, a
species-distribution model with validated absences, or evolutionary inference.
"""
from __future__ import annotations
import argparse,hashlib,io,json
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold
from scripts.analysis.mois_global_hue9_community_assembly_genus_null_20261118 import (
    signed_original,baseline,HUES,sha,check,GBIF_SHA,ENV_SHA
)
from scripts.analysis.mois_wcvp_native_hue9_geographical_null_20261118 import (
    load_palaeoclimate_grid,assign_cells_broad_region,get_wcvp_source,
    match_gbif_to_wcvp,regional_species_cell_rows,LGM_SHA,TDWG_GEO_SHA
)

ECO_SHA='4f5d8fc3a308482fc3645d796a779d1a4725c51ed2afbf6a06e61fd6213f3ee6'
ECO_CSV='all_current_GBIF_100km_cells_genus_breadth_entropy_ecology.csv.gz'
GCM=('CCSM4','MIROC_ESM','MPI_ESM_P')
MIN_NATIVE_GRIDS=50
MAX_PER_HUE=100


def source_environment(path):
    check(sha(path)==ECO_SHA,'Original modern environment ZIP SHA drift')
    with ZipFile(path) as z:
        rec=json.loads(z.read('receipt.json'))
        for p,digest in rec['numeric_outputs_sha256'].items():
            check(sha(z.read(p))==digest,'Modern environment member digest drift '+p)
        d=pd.read_csv(io.BytesIO(z.read(ECO_CSV)),compression='gzip')
    check(len(d)==12295 and not d[['grid_x','grid_y']].duplicated().any(),
          'Original modern climate sample deviated from signed 12295 cells')
    check(rec['source_cohort']['strict_complete_30plus_cells']==12295,
          'Signed original environmental gate changed')
    return d


def climate_sources(econ, lgm,regions):
    fields=['grid_x','grid_y','threeGCM_valid_all_models']
    for g in GCM:
        fields.extend([f'{g}_delta_temperature_C',f'{g}_delta_log_precip'])
    check(set(fields).issubset(lgm),'3GCM source delta column missing')
    d=econ[['grid_x','grid_y','present_temperature_C','present_precip_mm']].merge(
        lgm[fields],on=['grid_x','grid_y'],validate='one_to_one')
    d=d.merge(regions[['grid_x','grid_y','code','region_strict']],
        on=['grid_x','grid_y'],validate='one_to_one')
    check(len(d)==12295,'Modern climate/LGM/TDWG source grids differ')
    d=d.loc[d.region_strict & d.threeGCM_valid_all_models & d.code.notna()].copy()
    values=['present_temperature_C','present_precip_mm'] + [f'{g}_delta_{v}' for g in GCM for v in ('temperature_C','log_precip')]
    d=d.loc[d[values].notna().all(axis=1) & d.present_precip_mm.ge(0)].copy()
    check(len(d)>6500 and d.code.nunique()>100,'Too few actual LGM and native-region grid cells')
    # WorldClim1.4 current and LGM use 1960-90 baseline (source receipt).
    # Delta definitions are present-minus-LGM on temperature and log1p precip.
    d['modern_log1p_rain']=np.log1p(d.present_precip_mm)
    for g in GCM:
        d[f'{g}_LGM_temperature_C']=d.present_temperature_C-d[f'{g}_delta_temperature_C']
        d[f'{g}_LGM_log1p_rain']=d.modern_log1p_rain-d[f'{g}_delta_log_precip']
    past=[f'{g}_LGM_{v}' for g in GCM for v in ('temperature_C','log1p_rain')]
    check(np.isfinite(d[past].to_numpy(float)).all(),'Nonfinite source-verified palaeoclimate')
    # Some interpolated coastal cells can yield a negative back-transformed
    # log1p precipitation. Exclude rather than silently set precipitation zero.
    good=d[[f'{g}_LGM_log1p_rain' for g in GCM]].ge(-1e-8).all(axis=1)
    exclusion=int((~good).sum());d=d.loc[good].copy()
    check(len(d)>6400,'Native palaeoclimate test cells exhausted by physical QA')
    d['cell_index']=np.arange(len(d),dtype=np.int32)
    d['block500']=(d.grid_x//5).astype(str)+'|'+(d.grid_y//5).astype(str)
    return d,exclusion


def envelope_stats(y,grid):
    q=np.quantile(y,[.05,.95],axis=0)
    check(np.isfinite(q).all() and (q[1]>=q[0]).all(),'Invalid empirical native climate envelope')
    return q[0],q[1]


def within(matrix,lo,hi):
    return (matrix[:,0]>=lo[0])&(matrix[:,0]<=hi[0])&(
            matrix[:,1]>=lo[1])&(matrix[:,1]<=hi[1])


def study_species(k,subset,grid,modern,past_by_gcm,idx_by_code):
    n=len(subset)
    occupied=subset.cell_index.to_numpy(np.int32)
    codes=subset.code.unique().tolist()
    possible=np.unique(np.concatenate([idx_by_code[v] for v in codes]))
    check(np.isin(occupied,possible).all(), 'Native occupancy outside allowed candidate TDWG3 regions')
    block=subset.block500.to_numpy(str)
    nb=len(np.unique(block))
    lo,hi=envelope_stats(modern[occupied],grid)
    region_modern=modern[possible]
    now=within(region_modern,lo,hi)
    nnow=int(now.sum())
    counts={}
    stacked=[]
    for g in GCM:
        a=within(past_by_gcm[g][possible],lo,hi)
        stacked.append(a)
        counts[g]=int(a.sum())
    same=np.logical_and.reduce(stacked)
    any_gcm=np.logical_or.reduce(stacked)
    occupied_three=0
    for g in GCM:
        occupied_three += within(past_by_gcm[g][occupied],lo,hi).astype(np.int8)
    cv_recall=[];cv_lift=[];folds=0
    if nb>=5:
        splitter=GroupKFold(n_splits=5)
        for train,test in splitter.split(occupied,groups=block):
            if len(train)<25 or len(test)<3:continue
            ll,hh=envelope_stats(modern[occupied[train]],grid)
            rr=within(modern[occupied[test]],ll,hh).mean()
            bg=within(region_modern,ll,hh).mean()
            cv_recall.append(float(rr));cv_lift.append(float(rr/max(bg,0.01)))
            folds+=1
    rec={
        'source_specieskey':k,
        'n_modern_native_occupied_cells':n,
        'n_modern_native_500km_blocks':nb,
        'n_WCVP_TDWG3_native_regions':len(codes),
        'n_candidate_cells_in_current_native_regions':len(possible),
        'climate_envelope_T5_C':float(lo[0]),
        'climate_envelope_T95_C':float(hi[0]),
        'climate_envelope_logP5':float(lo[1]),
        'climate_envelope_logP95':float(hi[1]),
        'modern_climate_analog_cell_count':nnow,
        'modern_climate_analog_fraction':float(nnow/len(possible)),
        'LGM_threeGCM_consensus_climate_analog_cells':int(same.sum()),
        'LGM_anyGCM_climate_analog_cells':int(any_gcm.sum()),
        'LGM_threeGCM_consensus_fraction_region':float(same.mean()),
        'LGM_anyGCM_fraction_region':float(any_gcm.mean()),
        'LGM_consensus_vs_modern_analog_cell_ratio':float((int(same.sum())+.5)/(nnow+.5)),
        'LGM_occupied_sites_threeGCM_consensus_fraction':float(np.mean(occupied_three==3)),
        'n_spatial_500km_cv_folds':folds,
        'modern_heldout_native_recall':float(np.mean(cv_recall)) if cv_recall else np.nan,
        'modern_heldout_native_to_region_background_lift':float(np.mean(cv_lift)) if cv_lift else np.nan,
    }
    for g in GCM:
        rec[f'{g}_LGM_analog_cells']=counts[g]
        rec[f'{g}_LGM_vs_modern_analog_ratio']=float((counts[g]+.5)/(nnow+.5))
    return rec


def run(args):
    root=Path(args.out)
    check(not root.exists(),'Refuse overwriting source-verified output')
    pairs,species=signed_original(args.GBIF,args.ENV)
    source=baseline(pairs,species,null_model='genus_breadth')
    check(source['ngrids']==18577 and source['n_accepted_species_with_at_least_one_GBIF_cell']==36758,
          'GBIF accepted source coverage changed')
    lgm=load_palaeoclimate_grid(args.LGM)
    regions=assign_cells_broad_region(lgm,args.TDWG)
    check(regions.region_strict.sum()==9288 and regions.region_moderate.sum()==10391,
          'Original native classification regional sample has drifted')
    names,dist,manifest=get_wcvp_source(args.WCVP)
    mapped=match_gbif_to_wcvp(species,names)
    check(mapped.wcvp_match_exact_binomial.sum()==36762,'WCVP exact accepted matching changed')
    ids=regions.loc[regions.region_moderate,['grid_x','grid_y']].drop_duplicates()
    src=regional_species_cell_rows(pairs,species,names,dist,regions,ids)
    check(int((src.region_strict&src.status.eq('native_at_TDWG')).sum())==2271852,
          'Original strict region-native 2.27M species×grid source did not reproduce')
    check(int((src.region_strict&src.status.eq('introduced_at_TDWG')).sum())==436422,
          'Original strict introduced 436K species×grid source did not reproduce')
    climate,invalid_past_rain=climate_sources(source_environment(args.ECO),lgm,regions)
    src=src.loc[src.region_strict & src.status.eq('native_at_TDWG')].copy()
    src=src.merge(climate[['grid_x','grid_y','code','cell_index','block500']],
                  on=['grid_x','grid_y'],how='inner',validate='many_to_one')
    check(not src[['specieskey','cell_index']].duplicated().any(),
          'One WCVP-native source species replicated in grid')
    source_sp=species[['key','species','genus','source_hue9','hue_status']].copy()
    agg=src.groupby('specieskey').agg(
        n_grid=('cell_index','nunique'),n_spatial_blocks=('block500','nunique'))
    eligible=source_sp.merge(agg,left_on='key',right_index=True,how='inner')
    eligible=eligible.loc[eligible.n_grid.ge(MIN_NATIVE_GRIDS)&eligible.n_spatial_blocks.ge(5)].copy()
    check(len(eligible)>1000 and eligible.source_hue9.nunique()==9,
          'Insufficient nine-hue source-native species with wide modern range')
    eligible['deterministic_hash']=eligible.key.map(lambda x:hashlib.sha256(
        ('MOIS-20261119-native-paleoclimate:'+str(x)).encode()).hexdigest())
    eligible=eligible.sort_values('deterministic_hash',kind='stable')
    selected=eligible.groupby('source_hue9',sort=False).head(MAX_PER_HUE).copy()
    check(selected.source_hue9.nunique()==9 and len(selected)>250,
          'All original nine hue classes must be represented in pilot')
    native=src.loc[src.specieskey.isin(selected.key)].copy()
    modern=climate[['present_temperature_C','modern_log1p_rain']].to_numpy(float)
    past={g:climate[[f'{g}_LGM_temperature_C',f'{g}_LGM_log1p_rain']].to_numpy(float) for g in GCM}
    idx_by_code={k:g.cell_index.to_numpy(np.int32) for k,g in climate.groupby('code')}
    records=[]
    grp=native.groupby('specieskey',sort=False)
    for i,k in enumerate(selected.key):
        info=study_species(str(k),grp.get_group(str(k)),climate,modern,past,idx_by_code)
        records.append(info)
        if (i+1)%50==0:print('MOIS_LGM_NATIVE_CLIMATE_ANALOG_SPECIES_COMPLETED',i+1,flush=True)
    results=pd.DataFrame(records).merge(selected[['key','species','genus','source_hue9','hue_status','n_grid']],
                                        left_on='source_specieskey',right_on='key',validate='one_to_one')
    check(len(results)==len(selected) and results.source_specieskey.is_unique,
          'Pilot species climate results incomplete')
    group=results.groupby('source_hue9').agg(
        n_species=('source_specieskey','count'),
        median_current_native_sites=('n_modern_native_occupied_cells','median'),
        median_current_region_analogue_cells=('modern_climate_analog_cell_count','median'),
        median_LGM_threeGCM_consensus_cells=('LGM_threeGCM_consensus_climate_analog_cells','median'),
        median_LGM_ensemble_consensus_to_modern_ratio=('LGM_consensus_vs_modern_analog_cell_ratio','median'),
        median_500km_heldout_native_recall=('modern_heldout_native_recall','median'),
        median_500km_heldout_background_lift=('modern_heldout_native_to_region_background_lift','median'),
    ).reset_index()
    root.mkdir(parents=True)
    t=root/'WCVP_strict_native_species_LGM_3GCM_climatic_analogue_pilot.csv.gz'
    results.to_csv(t,index=False,compression='gzip')
    s=root/'nine_hue_pilot_climate_analogues_summary.csv'
    group.to_csv(s,index=False)
    source_receipt={
      'source_accepted_photo_hue_taxa':37104,
      'full_original_GBIF_species_grid_pairs':5493931,
      'original_GBIF_occurrence_records':396142825,
      'strict_WCVP_native_source_species_grid_pairs':2271852,
      'WCVP_TDWG3_regional_status_not_individual_wild_provenance':True,
      'source_WCVP_unique_exact_matched_species':36762,
      'all_original_global_GBIF_grid_cells':18577,
      'shared_1960_1990_current_3GCM_LGM_climate_not_1970_2000_WC21':True,
      'original_complete_environment_cells':12295,
      'WCVP_native_and_3GCM_climate_complete_candidate_grid_cells':int(len(climate)),
      'excluded_negative_backcalculated_log1p_LGM_rain_cells':int(invalid_past_rain),
      'n_species_eligible_native50_5spatialblocks':int(len(eligible)),
      'n_species_in_stratified_deterministic_pilot':int(len(selected)),
      'n_pilot_by_hue':{k:int(v) for k,v in selected.source_hue9.value_counts().items()},
      'source_zip_SHA256':{'GBIF':GBIF_SHA,'ENVoriginal':ENV_SHA,
            'LGM':LGM_SHA,'ENVcurrent':ECO_SHA,'TDWG3_geojson':TDWG_GEO_SHA},
      'source_WCVP_manifest_SHA256':sha(Path(args.WCVP)/'mois350399_wcvp_TDWG3_source_only_taxon_native_receipt.json'),
      'current_native_climate_90pct_rectangle_not_mechanistic_SDM':True,
      'candidate_domain_current_WCVP_TDWG3_regions_with_3GCM_land_climate_only':True,
      'no_glacial_fossil_or_pollen_or_phylogeography_evidence':True,
      'LGM_analogue_is_not_LGM_species_occupancy_or_refugium':True,
      'photo9_modal_roi_not_independently_validated_petal_colour':True,
      'stratified_hue_pilot_does_not_estimate_global_colour_prevalence':True,
      'selection_deterministic_hash_ignores_LGM_outcome':True,
      'source_data_had_no_new_image_labeling':True,
    }
    full=root/'source_screening_and_3GCM_PILOT.json'
    full.write_text(json.dumps(source_receipt,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    source_receipt['outputs_sha256']={p.name:sha(p) for p in (t,s,full)}
    (root/'receipt.json').write_text(json.dumps(source_receipt,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    print('MOIS_REAL_WCVP_NATIVE_SPECIES_LGM_3GCM_PILOT_SUCCESS',json.dumps({
        'eligible':len(eligible),'pilot':len(selected),
        'climate_cells':len(climate),
        'by_hue':group.to_dict('records'),
        'native_3GCM_consensus_nonzero':int(results.LGM_threeGCM_consensus_climate_analog_cells.gt(0).sum()),
        'n_with_holdout_lift_gt1':int(results.modern_heldout_native_to_region_background_lift.gt(1).sum()),
        'scientific_limit':'Current climatic niche analogues, not actual LGM populations or refugia'
    },ensure_ascii=False,allow_nan=False),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser()
    for v in ('GBIF','ENV','LGM','WCVP','TDWG','ECO','out'):
        p.add_argument('--'+v,required=True)
    run(p.parse_args())