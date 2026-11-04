#!/usr/bin/env python3
"""Original-image genus identity and regional heldout strict9 palette (numeric only).

WCVP accepted genus does not establish native-island status, pollination guild,
flower size, floral dullness, calibrated chroma, or plant fitness.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import pandas as pd
import numpy as np

METRICS=('strict9_white','strict9_blue','strict9_purple','strict9_bronze','palette_entropy')
PAIRS=(('arctic','north_subpolar'),('south_temperate','north_temperate'),('island_cached','mainland_cached'))
EXPECTED_REGION_ROWS=30634
PINS={
 'regional':'5f9dee3b0ab41e42d0b1b1dd06650eda2b954b89cb118ddc8bcdd1a813da1479',
 'taxa':'dc5bb0e3686fde0cfd327e4336245d0ceae3c7657b3bd12841bb814cc94d1eb2',
 'wcvp':'9cd6aa7ca29f00c64283bcc08dcb778e212b644b837e64070051339b11d8b7b3',
}

def require(condition,msg):
    if not bool(condition):raise ValueError('MOIS_GENUS_HOLDOUT: '+msg)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as fh:
        for part in iter(lambda:fh.read(4*1024*1024),b''):h.update(part)
    return h.hexdigest()

def accepted_taxonomy(names,distribution):
    require({'inat_taxon_id','wcvp_match_state','wcvp_accepted_plant_name_id'}<=set(names),'Invalid source taxon mapping')
    require({'wcvp_accepted_plant_name_id','wcvp_accepted_name'}<=set(distribution),'Invalid accepted-name source')
    names=names.copy();distribution=distribution.copy()
    for c in ('inat_taxon_id','wcvp_accepted_plant_name_id'):
        names[c]=names[c].fillna('').astype(str).str.strip()
    distribution['wcvp_accepted_plant_name_id']=distribution.wcvp_accepted_plant_name_id.fillna('').astype(str).str.strip()
    distribution['wcvp_accepted_name']=distribution.wcvp_accepted_name.fillna('').astype(str).str.strip()
    require(names.inat_taxon_id.is_unique and names.inat_taxon_id.ne('').all(),'Duplicate original source taxon ID')
    require(distribution.wcvp_accepted_plant_name_id.ne('').all(),'Missing accepted species ID')
    need=distribution.groupby('wcvp_accepted_plant_name_id').wcvp_accepted_name.nunique()
    require(need.le(1).all(),'Conflicting accepted names for same species ID')
    valid=names.wcvp_match_state.eq('unique_exact_name_to_accepted_id')
    key=distribution[['wcvp_accepted_plant_name_id','wcvp_accepted_name']].drop_duplicates()
    selected=names.loc[valid,['inat_taxon_id','wcvp_accepted_plant_name_id']].merge(
        key,on='wcvp_accepted_plant_name_id',how='left',validate='many_to_one')
    acceptable=selected.wcvp_accepted_name.fillna('').str.fullmatch(r'[A-Z][A-Za-z-]+ [a-z][A-Za-z-]+')
    selected=selected.loc[acceptable].copy()
    selected['genus']=selected.wcvp_accepted_name.str.split().str[0]
    selected=selected.rename(columns={'wcvp_accepted_plant_name_id':'accepted_id'})
    require(selected.inat_taxon_id.is_unique and selected.accepted_id.ne('').all(),'Invalid accepted botanical identity')
    diagnostics={'source_taxon_ids':len(names),'unique_accepted_id_status':int(valid.sum()),
                 'resolved_accepted_binomial_taxa':len(selected),
                 'rejected_ambiguous_names':int(names.wcvp_match_state.eq('multiple_accepted_ids_ambiguous').sum()),
                 'rejected_unresolved_or_unverified':int((~valid).sum()),
                 'accepted_id_without_reliable_binomial':int(valid.sum()-len(selected)),
                 'accepted_genera':int(selected.genus.nunique()),
                 'accepted_species_ids':int(selected.accepted_id.nunique()),
                 'WCVP_region_distribution_not_used_for_photo_native_status':True}
    return selected,diagnostics

def collapse_to_accepted_species(table,mapping,metrics):
    require({'sector','inat_taxon_id','n',*metrics}<=set(table),'Regional taxon palette columns missing')
    d=table.copy();d.inat_taxon_id=d.inat_taxon_id.fillna('').astype(str)
    require(not d.duplicated(['sector','inat_taxon_id']).any(),'Duplicated regional source taxon')
    require(d.n.gt(0).all(),'Non-positive original strict9 photo support')
    require(np.isfinite(d[list(metrics)].to_numpy(float)).all(),'Nonfinite original measured ROI traits')
    require(d.inat_taxon_id.ne('').all(),'Missing original taxon IDs')
    j=d.merge(mapping[['inat_taxon_id','accepted_id','genus']],how='left',on='inat_taxon_id',validate='many_to_one')
    supported=j.loc[j.accepted_id.notna()].copy()
    for c in metrics:supported[c]=supported[c].astype(float)*supported.n.to_numpy()
    sums=supported.groupby(['sector','accepted_id','genus'],as_index=False,observed=True).agg(
        {'n':'sum',**{m:'sum' for m in metrics}})
    for c in metrics:sums[c]=sums[c]/sums.n
    require(not sums.duplicated(['sector','accepted_id']).any(),'Duplicate accepted species in regional source')
    stats={'input_nominal_taxon_sector_rows':len(j),'resolved_taxon_sector_rows':len(supported),
           'ambiguous_or_unmapped_input_rows':len(j)-len(supported),
           'accepted_species_sector_rows':len(sums),'source_synonym_row_collapses':len(supported)-len(sums),
           'resolved_photo_support':int(supported.n.sum()),'total_photo_support':int(j.n.sum()),
           'unmapped_original_photo_support':int(j.loc[j.accepted_id.isna(),'n'].sum())}
    return sums,stats

def genus_composition(species,focal,reference,metric):
    """Exact four-term equal-genus assemblage identity; not evolutionary causation."""
    a=species.loc[species.sector.eq(focal)].groupby('genus')[metric].mean()
    b=species.loc[species.sector.eq(reference)].groupby('genus')[metric].mean()
    require(len(a)>0 and len(b)>0,'No accepted genera in requested region')
    common=a.index.intersection(b.index);one=a.index.difference(b.index);two=b.index.difference(a.index)
    fa=1/len(a);fb=1/len(b)
    terms={
        'common_genus_trait_pp':100*float(((a.loc[common]-b.loc[common])*(fa+fb)/2).sum()),
        'common_genus_weight_pp':100*float(((a.loc[common]+b.loc[common])*(fa-fb)/2).sum()),
        'focal_only_genera_pp':100*float(a.loc[one].sum()*fa),
        'reference_only_genera_pp':-100*float(b.loc[two].sum()*fb),
    }
    delta=100*float(a.mean()-b.mean())
    require(abs(sum(terms.values())-delta)<1e-9,'Genus composition decomposition does not close')
    return {'focal':focal,'reference':reference,'metric':metric,
            'focal_accepted_genera':len(a),'reference_accepted_genera':len(b),
            'shared_genera':len(common),'shared_genus_fraction_focal':len(common)/len(a),
            'shared_genus_fraction_reference':len(common)/len(b),
            'genus_equal_difference_pp':delta,
            'common_genus_only_difference_pp':100*float((a.loc[common]-b.loc[common]).mean()) if len(common) else None,**terms}

def score_holdout(species,focal,reference,metrics=METRICS,min_species=3):
    """Only novel focal accepted species predicted from reference-only genus means."""
    require(min_species>=2,'No singleton genus extrapolation')
    a=species.loc[species.sector.eq(focal)].copy()
    b=species.loc[species.sector.eq(reference)].copy()
    require(not a.accepted_id.duplicated().any() and not b.accepted_id.duplicated().any(),'Duplicated accepted species')
    train_species=set(b.accepted_id)
    novel=a.loc[~a.accepted_id.isin(train_species)].copy()
    means=b.groupby('genus',as_index=False,observed=True).agg(
        {'accepted_id':'nunique',**{metric:'mean' for metric in metrics}})
    means=means.rename(columns={'accepted_id':'train_species_in_genus',
                              **{metric:'pred_'+metric for metric in metrics}})
    means=means.loc[means.train_species_in_genus.ge(min_species)]
    scored=novel.merge(means,on='genus',how='inner',validate='many_to_one').sort_values('accepted_id').reset_index(drop=True)
    require(scored.accepted_id.is_unique and not scored.accepted_id.isin(train_species).any(),
            'Target-to-train accepted species leakage')
    result={'focal':focal,'reference':reference,
            'min_distinct_training_species_per_genus':min_species,
            'focal_accepted_species':len(a),'reference_accepted_species':len(b),
            'shared_accepted_species':len(a)-len(novel),'novel_target_species':len(novel),
            'supported_novel_species':len(scored),
            'novel_support_fraction':len(scored)/len(novel) if len(novel) else None,
            'all_focal_species_support_fraction':len(scored)/len(a) if len(a) else None,
            'supported_target_genera':int(scored.genus.nunique()),'metrics':{}}
    for metric in metrics:
        if not len(scored):
            result['metrics'][metric]={'mae_genus':None,'mae_reference_global_mean':None,
                'mse_genus':None,'mse_reference_global_mean':None,
                'mse_skill_vs_reference_global_mean':None,'observed_minus_genus_prediction_pp':None}
            continue
        y=scored[metric].to_numpy(float);pred=scored['pred_'+metric].to_numpy(float)
        baseline=float(b[metric].mean())
        denom=float(np.mean((y-baseline)**2));numer=float(np.mean((y-pred)**2))
        result['metrics'][metric]={
            'mae_genus':float(np.mean(np.abs(y-pred))),
            'mae_reference_global_mean':float(np.mean(np.abs(y-baseline))),
            'mse_genus':numer,'mse_reference_global_mean':denom,
            'mse_skill_vs_reference_global_mean':float(1-numer/denom) if denom>1e-15 else None,
            'observed_minus_genus_prediction_pp':100*float(np.mean(y-pred))}
    scored.insert(0,'focal',focal);scored.insert(1,'reference',reference)
    return result,scored

def run(regional,names,distribution,prior_receipt,output):
    paths={'regional':Path(regional),'taxa':Path(names),'wcvp':Path(distribution)}
    dest=Path(output);require(not dest.exists(),'Refuse to overwrite a previous genus result')
    checks={k:sha(p) for k,p in paths.items()}
    require(checks==PINS,'Input source SHA mismatch')
    prior=json.loads(Path(prior_receipt).read_text())
    require(prior['schema']=='mois_original_canonical_taxon_palette_environment_20261030_v1'
            and prior['output_sha256'].get('taxon_region_palette.csv.gz')==checks['regional']
            and prior['canonical_observations']==1032882 and prior['strict9_success']==401703
            and prior['MNAR_scenarios_run'] is False,'Not the previously audited numeric-only region result')
    source=pd.read_csv(paths['regional'],compression='gzip',dtype={'inat_taxon_id':str})
    require(len(source)==EXPECTED_REGION_ROWS,'Unexpected original sector taxon counts')
    names_df=pd.read_csv(paths['taxa'],compression='gzip',dtype={'inat_taxon_id':str,
        'wcvp_accepted_plant_name_id':str},low_memory=False)
    dist=pd.read_csv(paths['wcvp'],compression='gzip',
        dtype={'wcvp_accepted_plant_name_id':str},
        usecols=['wcvp_accepted_plant_name_id','wcvp_accepted_name'],
        low_memory=False).drop_duplicates()
    mapping,taxdiag=accepted_taxonomy(names_df,dist)
    species,coverage=collapse_to_accepted_species(source,mapping,METRICS)
    require(len(species)>0 and source.sector.nunique()==6 and
            coverage['resolved_taxon_sector_rows']>28000,'Insufficient verified taxon mapping')
    genus=[];scores=[];predictions=[]
    for a,b in PAIRS:
        for metric in METRICS:genus.append(genus_composition(species,a,b,metric))
        for focal,ref in ((a,b),(b,a)):
            record,scored=score_holdout(species,focal,ref,METRICS,3)
            scores.append({'focal':focal,'reference':ref,
                **{k:v for k,v in record.items() if k not in ('focal','reference','metrics')},
                **{('white_'+key):value for key,value in record['metrics']['strict9_white'].items()},
                **{('bronze_'+key):value for key,value in record['metrics']['strict9_bronze'].items()}})
            predictions.append(scored)
    dest.mkdir(parents=True)
    outputs={
        'accepted_species_region.csv.gz':species,
        'accepted_genus_composition.csv':pd.DataFrame(genus),
        'region_heldout_genus_performance.csv':pd.DataFrame(scores),
        'region_heldout_predictions.csv.gz':pd.concat(predictions,ignore_index=True),
    }
    for filename,df in outputs.items():
        if filename.endswith('.gz'):
            df.to_csv(dest/filename,index=False,compression={'method':'gzip','mtime':0})
        else:df.to_csv(dest/filename,index=False)
    regional_support=species.groupby('sector',as_index=False).agg(
        accepted_species=('accepted_id','size'),genera=('genus','nunique'),
        strict9_photo_support=('n','sum'))
    regional_support.to_csv(dest/'accepted_region_coverage.csv',index=False)
    report={
      'schema':'mois_existing_9palette_wcvp_genus_composition_crossregion_holdout_20261030_v1',
      'prior_source_receipt_SHA256':sha(prior_receipt),'input_sha256':checks,
      'source_10_month_2026_WCVP_mapping':taxdiag,'matched_source_coverage':coverage,
      'region_coverage':regional_support.to_dict(orient='records'),
      'heldout_summary':scores,'out_of_focal_sector_training':True,
      'heldout_novel_accepted_species_only':True,'photo_pixels_downloaded':0,
      'original_image_roi_inference':0,'MNAR_scenarios_run':False,
      'technical_unknown_imputed':False,'family_taxonomy_in_WCVP_derived_source':False,
      'background_contrast_or_dull_flower_measured':False,
      'regional_assignments_are_not_a_species_level_native_island_validation':True,
      'predictions_are_descriptive_not_independent_botanical_validation':True,
      'output_sha256':{f.name:sha(f) for f in sorted(dest.iterdir())}}
    (dest/'genus_holdout_receipt.json').write_text(
        json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n')
    print('MOIS_WCVP_GENUS_HOLDOUT_COMPLETE',json.dumps({
        'sector_accepted_species_rows':len(species),
        'mapped_legacy_taxon_sector_rows':coverage['resolved_taxon_sector_rows'],
        'heldout_directional_comparisons':len(scores),
        'genera_metrics':len(genus),'all_image_downloads':0,
        'all_MNAR_scenarios':0},sort_keys=True),flush=True)
    for row in scores:
        print('HELDOUT_REGION',row['focal'],'reference',row['reference'],
              'new_species',row['novel_target_species'],
              'supported',row['supported_novel_species'],
              'genus_vs_global_mse_skill',
              row['white_mse_skill_vs_reference_global_mean'],flush=True)
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('regional','names','distribution','prior_receipt','output'):
        p.add_argument('--'+k.replace('_','-'),dest=k,required=True)
    args=p.parse_args()
    run(args.regional,args.names,args.distribution,args.prior_receipt,args.output)
