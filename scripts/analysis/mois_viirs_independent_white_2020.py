#!/usr/bin/env python3
"""2020 MOIS original VIIRS/WorldPop bright-dark pairs vs independent Island direct white records.

Read third-party species traits only in a private runner; publish AGGREGATE JSON only.
No raw per-species colour, individual source coordinates or native/extinction claims.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd

PINS={
 'canonical':'f033661d76e07898a5c13199a55c5d0ce6e55353daa48b0058c456979e90644c',
 'wcvp_lookup':'dc5bb0e3686fde0cfd327e4336245d0ceae3c7657b3bd12841bb814cc94d1eb2',
 'wcvp_names':'9cd6aa7ca29f00c64283bcc08dcb778e212b644b837e64070051339b11d8b7b3',
 'island_direct':'7e86170218814abe2845da9a0282aecfdc930819b9d0c439bae8df4a0679c577',
 'candidate_cells':'d8d0b7d5327a9add2107b2c6a1fa0e5aa6358d0813ef90da6f70e478e7a4f336',
 'matched_pairs':'32b53d23d4beac39846e3c618f7a4f5d8b24164a126cf63c4cc530b2e720aae6',
}
SINGLE={'white':1,'yellow_orange':0,'red_pink':0,'blue_purple':0,'green_brown_inconspicuous':0}

def require(ok,msg):
    if not bool(ok):raise ValueError('MOIS_ALAN_DIRECT_COLOUR: '+msg)

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda:f.read(4<<20),b''):h.update(block)
    return h.hexdigest()

def exact_direct_species_labels(frame):
    required={'accepted_species','trait_name','normalized_value','quality','resolution_status','classification'}
    require(required<=set(frame),'Direct trait source schema changed')
    x=frame.loc[frame.trait_name.eq('flower_primary_color') &
        frame.quality.isin(['high','medium']) & frame.resolution_status.eq('resolved') &
        ~frame.classification.isin(['cultivar_contamination','partial_self_compatibility'])].copy()
    x['white_direct']=x.normalized_value.map(SINGLE)
    x=x.dropna(subset=['white_direct'])[['accepted_species','white_direct']]
    require(x.accepted_species.is_unique,'Duplicate direct species colour')
    require(x.white_direct.isin([0,1]).all(),'Mixed or unknown colour leaked into known denominator')
    return x

def source_taxonomy(names,accepted,direct):
    require({'inat_taxon_id','wcvp_match_state','wcvp_accepted_plant_name_id'}<=set(names),'WCVP lookup schema changed')
    require({'wcvp_accepted_plant_name_id','wcvp_accepted_name'}<=set(accepted),'WCVP accepted species changed')
    names=names.loc[names.wcvp_match_state.eq('unique_exact_name_to_accepted_id')].copy()
    names.inat_taxon_id=names.inat_taxon_id.astype(str)
    accepted=accepted[['wcvp_accepted_plant_name_id','wcvp_accepted_name']].drop_duplicates()
    require(accepted.groupby('wcvp_accepted_plant_name_id').wcvp_accepted_name.nunique().le(1).all(),
            'WCVP accepted species ID has contradictory accepted names')
    j=names[['inat_taxon_id','wcvp_accepted_plant_name_id']].merge(
        accepted,on='wcvp_accepted_plant_name_id',how='left',validate='many_to_one')
    j=j.merge(direct,left_on='wcvp_accepted_name',right_on='accepted_species',how='left',validate='many_to_one')
    require(j.inat_taxon_id.is_unique,'Original source taxon maps more than once')
    return j[['inat_taxon_id','wcvp_accepted_plant_name_id','white_direct']]

def aggregate_known_photo_taxa(canonical,mapping,selected_cells):
    needed={'source_subarm','observed_on','latitude','longitude','climate_complete',
            'photo_id','observation_id','inat_taxon_id'}
    require(needed<=set(canonical),'Canonical source schema changed')
    require(canonical.photo_id.is_unique and canonical.observation_id.is_unique,
            'Photo / observation source identity drift')
    x=canonical.loc[canonical.source_subarm.eq('historical_main')&
        pd.to_numeric(canonical.observed_on.astype(str).str[:4],errors='coerce').eq(2020)&
        canonical.climate_complete&canonical.latitude.between(-60,60)&
        canonical.longitude.between(-180,180)].copy()
    x['cell_id']=(np.floor((x.longitude+180)/.05).astype(int).astype(str)+'_'+
                  np.floor((x.latitude+90)/.05).astype(int).astype(str))
    x=x.loc[x.cell_id.isin(selected_cells.cell_id)].copy()
    require(len(x)>5000,'Frozen selected 2020 source cells unexpectedly sparse')
    x.inat_taxon_id=x.inat_taxon_id.astype(str)
    x=x.merge(mapping,on='inat_taxon_id',how='left',validate='many_to_one')
    good=x.loc[x.white_direct.notna()&x.wcvp_accepted_plant_name_id.notna(),
               ['cell_id','wcvp_accepted_plant_name_id','white_direct']].drop_duplicates()
    require(not good.duplicated(['cell_id','wcvp_accepted_plant_name_id']).any(),
            'Accepted species duplicated within a location with different labels')
    cell=good.groupby('cell_id',observed=True).agg(
        direct_white_share=('white_direct','mean'),
        n_documented_species=('white_direct','size'),
        n_white=('white_direct','sum')).reset_index()
    return cell,{'source_attempts_with_selected_cells':len(x),
        'source_attempts_with_direct_colour':int(x.white_direct.notna().sum()),
        'source_distinct_direct_taxa':int(x.loc[x.white_direct.notna(),'inat_taxon_id'].nunique()),
        'cells_with_direct_3plus':int((cell.n_documented_species>=3).sum()),
        'cells_with_direct_5plus':int((cell.n_documented_species>=5).sum())}

def paired_direct(pairs,cell,*,minimum=3):
    require(pairs.bright_cell_id.is_unique and pairs.dark_cell_id.is_unique,
            'Original source cell reused across independent matched pairs')
    a=cell[['cell_id','n_documented_species','direct_white_share']].rename(
        columns={'cell_id':'bright_cell_id','n_documented_species':'bright_direct_n',
                 'direct_white_share':'bright_white'})
    b=cell[['cell_id','n_documented_species','direct_white_share']].rename(
        columns={'cell_id':'dark_cell_id','n_documented_species':'dark_direct_n',
                 'direct_white_share':'dark_white'})
    d=pairs.merge(a,on='bright_cell_id',how='left',validate='many_to_one').merge(
        b,on='dark_cell_id',how='left',validate='many_to_one')
    d=d.loc[d.bright_direct_n.ge(minimum)&d.dark_direct_n.ge(minimum)].copy()
    d['direct_white_difference_pp']=100*(d.bright_white-d.dark_white)
    require(d.direct_white_difference_pp.between(-100,100).all(),'Invalid matched direct colour percentage')
    return d

def block_bootstrap(pairs,col,seed=20261105,iters=7000):
    if len(pairs)<20:return None
    groups=[f[col].to_numpy(float) for _,f in pairs.groupby('match_block',sort=True)]
    rng=np.random.default_rng(seed)
    vals=[]
    for _ in range(iters):
        ids=rng.integers(0,len(groups),len(groups))
        vals.append(float(np.concatenate([groups[j] for j in ids]).mean()))
    return list(map(float,np.percentile(vals,[2.5,97.5])))

def run(canonical,wcvp_lookup,wcvp_names,island_direct,candidate_cells,matched_pairs,out):
    files={'canonical':Path(canonical),'wcvp_lookup':Path(wcvp_lookup),
           'wcvp_names':Path(wcvp_names),'island_direct':Path(island_direct),
           'candidate_cells':Path(candidate_cells),'matched_pairs':Path(matched_pairs)}
    checks={k:sha(p) for k,p in files.items()}
    require(checks==PINS,'Signed scientific original source SHA changed')
    require(not Path(out).exists(),'Refuse to overwrite an existing scientific result')
    orig=pd.read_csv(files['canonical'],compression='gzip',low_memory=False,dtype={'inat_taxon_id':str,'photo_id':str})
    n=pd.read_csv(files['wcvp_lookup'],compression='gzip',dtype=str)
    a=pd.read_csv(files['wcvp_names'],compression='gzip',
                  usecols=['wcvp_accepted_plant_name_id','wcvp_accepted_name'],dtype=str)
    direct=pd.read_csv(files['island_direct'],compression='gzip',dtype=str)
    c=pd.read_csv(files['candidate_cells'])
    p=pd.read_csv(files['matched_pairs'],dtype={'bright_cell_id':str,'dark_cell_id':str})
    require(len(c)==1059 and len(p)==128,'External MOIS/VIIRS source match baseline changed')
    m=source_taxonomy(n,a,exact_direct_species_labels(direct))
    cell,stats=aggregate_known_photo_taxa(orig,m,c)
    f=paired_direct(p,cell)
    report={'schema':'mois_2020_real_worldpop_viirs_independent_documented_white_taxa_aggregate_v1',
        'signed_inputs_SHA256':checks,'candidate_cells':len(c),
        'population_matched_pairs':len(p),'direct_colour_supported_pairs':len(f),
        'nominal_taxon_to_WCVP_accepted_unique_source':True,
        'independent_flower_colour_from_Island_high_medium_direct_only':True,
        'multistate_or_unresolved_never_nonwhite':True,
        'source':stats,
        'direct_white_bright_minus_dark_mean_pp':float(f.direct_white_difference_pp.mean()) if len(f)>=20 else None,
        'direct_white_median_pp':float(f.direct_white_difference_pp.median()) if len(f)>=20 else None,
        'direct_white_95pct_spatial_block_bootstrap_pp':block_bootstrap(f,'direct_white_difference_pp'),
        'direct_white_bright_higher_pairs':int((f.direct_white_difference_pp>0).sum()),
        'direct_white_bright_lower_pairs':int((f.direct_white_difference_pp<0).sum()),
        'direct_white_equal_pairs':int((f.direct_white_difference_pp==0).sum()),
        'same_supported_pair_original_photo_ROI_white_mean_pp':float(f.white_diff_bright_minus_dark_pp.mean()) if len(f)>=20 else None,
        'source_white_category_means_literature_single_state_white_NOT_moth_guild':True,
        'cross_sectional_recorded_species_NOT_native_flora_extirpation':True,
        'rights_restricted_individual_species_labels_published':False,
        'source_photo_images_downloaded':0,'new_roi_inference':0,
        'MNAR_sensitivity_scenarios_run':False}
    out=Path(out);out.mkdir(parents=True)
    (out/'independent_white_population_matched_aggregate_only.json').write_text(
        json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n')
    print('MOIS_ALAN_DIRECT_WHITE_2020_SOURCE_VERIFIED',
          json.dumps({k:report[k] for k in ('population_matched_pairs',
          'direct_colour_supported_pairs','direct_white_bright_minus_dark_mean_pp',
          'direct_white_95pct_spatial_block_bootstrap_pp')},sort_keys=True),flush=True)
    return report

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    for name in PINS:ap.add_argument('--'+name.replace('_','-'),dest=name,required=True)
    ap.add_argument('--out',required=True)
    ar=ap.parse_args()
    run(**vars(ar))
