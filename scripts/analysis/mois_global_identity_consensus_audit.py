"""Audit existing numerical ROI9 observations; never fetch images or public API records.

Exact-photo monthly joins are conservative. Unmapped numeric IDs remain unresolved.
A single-colour consensus survives cross-namespace aliasing if either namespace
alone contains >=3 clear agreeing observations. This does not prove petal identity.
"""
from __future__ import annotations
import argparse,hashlib,io,json,zipfile
from pathlib import Path
import numpy as np
import pandas as pd
from scripts.analysis.mois_species_first_palette9 import build_catalogue,PALETTE

UUID=r'[0-9a-fA-F]{8}-(?:[0-9a-fA-F]{4}-){3}[0-9a-fA-F]{12}'
ZIP_SHA={
 'photos':'fd7a9462f1832c2d6405611bbc71d45e49fac1f8e098069f37945eb19fee4247',
 'initial':'6cc7689938dcf13baf6498b954286924aa97f2b63b3418c41181b18e346c0d2b',
 'current':'ef6578a8934280de654d9ff3adf9c8458a041feca463dc52b3cd87f501ba527e',
 'snapshot':'63cb48ff89839bb4d031da450c6e23f7260c83f793b52b41add5f954f56dc172'}

def require(ok,message):
    if not bool(ok):raise ValueError('MOIS_GLOBAL_IDENTITY: '+message)

def canonicalize(photos,source):
    result=photos.copy().reset_index(drop=True);src=source.copy()
    cols=['original_photo_id','original_numeric_observation_id','original_taxon_id','snapshot_observation_uuid','snapshot_taxon_id']
    for col in ('photo_id','observation_id','inat_taxon_id'):
        result[col]=result[col].fillna('').astype(str).str.strip()
    for col in cols:src[col]=src[col].fillna('').astype(str).str.strip()
    require(result.photo_id.is_unique and src.original_photo_id.is_unique,'Duplicate photo source')
    result['source_observation_id']=result.observation_id
    result['identity_status']='unresolved_other'
    numeric=result.observation_id.str.fullmatch(r'[1-9][0-9]*')
    result.loc[numeric,'identity_status']='numeric_unresolved'
    result.loc[result.observation_id.str.fullmatch(UUID),'identity_status']='source_uuid'
    result.loc[result.observation_id.eq(''),'identity_status']='missing'
    valid=src.snapshot_observation_uuid.str.fullmatch(UUID)
    associations=src.loc[valid].groupby('original_numeric_observation_id').snapshot_observation_uuid.nunique()
    conflicts=set(associations.index[associations.gt(1)])
    src=src.set_index('original_photo_id').reindex(result.photo_id).reset_index(drop=True).fillna('')
    match=(numeric & src.original_numeric_observation_id.eq(result.observation_id)
        & src.original_taxon_id.eq(result.inat_taxon_id) & src.snapshot_taxon_id.eq(result.inat_taxon_id)
        & src.snapshot_observation_uuid.str.fullmatch(UUID) & ~result.observation_id.isin(conflicts))
    result.loc[match,'observation_id']=src.loc[match,'snapshot_observation_uuid'].str.lower().to_numpy()
    result.loc[match,'identity_status']='historical_exact_photo_uuid'
    result.loc[numeric & result.source_observation_id.isin(conflicts),'identity_status']='numeric_conflicting_snapshot_held'
    return result

def namespace_bounds(observations):
    d=observations[['inat_taxon_id','observation_id']].drop_duplicates().copy()
    numeric=d.observation_id.str.fullmatch(r'[1-9][0-9]*')
    require((numeric|d.observation_id.str.fullmatch(UUID)).all(),'Unrecognized successful observation ID')
    d['namespace']=np.where(numeric,'numeric','uuid')
    n=d.groupby(['inat_taxon_id','namespace']).size().unstack(fill_value=0).reindex(columns=['numeric','uuid'],fill_value=0)
    n['independent_lower_bound']=n[['numeric','uuid']].max(axis=1)
    n['independent_upper_bound']=n[['numeric','uuid']].sum(axis=1)
    return n

def read_member(folder,kind,name):
    raw=(Path(folder)/(kind+'.zip')).read_bytes()
    require(hashlib.sha256(raw).hexdigest()==ZIP_SHA[kind],kind+' ZIP SHA changed')
    with zipfile.ZipFile(io.BytesIO(raw)) as z:return z.read(name)

def csv(raw):
    return pd.read_csv(io.BytesIO(raw),compression='gzip',keep_default_na=False,dtype={c:str for c in ['inat_taxon_id','original_photo_id','original_numeric_observation_id','original_taxon_id','snapshot_taxon_id']})

def run(folder,output):
    out=Path(output);require(not out.exists(),'Refuse overwrite')
    initial=csv(read_member(folder,'initial','mois_species_first_42111_nine_palette_catalogue.csv.gz'))
    current=csv(read_member(folder,'current','mois_42111_taxa_corrected_strict9_counts.csv.gz'))
    accepted=csv(read_member(folder,'current','mois_37103_COLXR_corrected_strict9_counts.csv.gz'))
    snapshot=csv(read_member(folder,'snapshot','mois_historical_100543_original_photo_id_s3_uuid_date_real_snapshot.csv.gz'))
    raw=read_member(folder,'photos','mois_1181146_unique_original_photo_continuous9_by_arm.parquet')
    require(hashlib.sha256(raw).hexdigest()=='3fbd0f0ec7eb347b3bbc7ab509ec9312ec6d797ed932dcf7bb425959898c3835','Parquet SHA changed')
    photos=pd.read_parquet(io.BytesIO(raw))
    require(len(photos)==1181146 and photos.photo_id.is_unique and int(photos.strict9.sum())==451643,'Original counts changed')
    require(len(initial)==len(current)==42111 and len(accepted)==37103,'Full taxonomy not retained')
    identities=canonicalize(photos[['photo_id','observation_id','inat_taxon_id','strict9']],snapshot)
    vec=photos.loc[photos.strict9].copy()
    before=build_catalogue(vec,initial,initial).set_index('inat_taxon_id')
    initial=initial.set_index('inat_taxon_id').loc[before.index]
    current=current.set_index('inat_taxon_id').loc[before.index]
    fields=['hue_status','dominant_hue9','observed_hue9_set','independent_roi9_observations','distinct_palette_votes','ambiguous_observation_rois']
    means=['observed_mean_'+x for x in PALETTE]
    for col in fields:require(before[col].equals(initial[col]),'Original classifier replay differs: '+col)
    require(np.allclose(before[means].to_numpy(float),initial[means].replace('',np.nan).to_numpy(float),atol=5.01e-7,rtol=0,equal_nan=True),'Original means differ beyond saved rounding')
    newer=vec.copy();lookup=identities.set_index('photo_id').observation_id
    newer['observation_id']=lookup.loc[newer.photo_id].to_numpy()
    after=build_catalogue(newer,initial.reset_index(),initial.reset_index()).set_index('inat_taxon_id').loc[before.index]
    require(before.hue_status.equals(after.hue_status) and before.dominant_hue9.equals(after.dominant_hue9),'Status changes require separate reviewed application')
    changed=before.independent_roi9_observations.ne(after.independent_roi9_observations)
    changed_ids=before.index[changed]
    require((current.loc[changed_ids,'after_strict9']==before.loc[changed_ids,'after_strict9']).all(),'Cannot transfer baseline corrections over later successful pixels')
    fixed=current.copy()
    for col in fields+means:fixed.loc[changed_ids,col]=after.loc[changed_ids,col]
    fixed['identity_snapshot_patch_applied']=fixed.index.isin(changed_ids)
    require(fixed.after_strict9.equals(current.after_strict9),'Numerical photo success count changed')
    consensus=current.hue_status.eq('sampling_consensus_single')
    require((current.loc[consensus,'after_strict9']==before.loc[consensus,'after_strict9']).all(),'Current consensus has additional photos not replayed')
    bounds=namespace_bounds(newer.loc[newer.inat_taxon_id.isin(current.index[consensus])]).reindex(current.index[consensus])
    bounds=bounds.join(current.loc[consensus,['species','dominant_hue9']])
    bounds['status_after_snapshot_join']=after.loc[bounds.index,'hue_status']
    bounds['robust_to_cross_namespace_aliases']=bounds.independent_lower_bound.ge(3)&bounds.status_after_snapshot_join.eq('sampling_consensus_single')
    fixed['consensus_alias_bound_verified']=fixed.index.isin(bounds.index[bounds.robust_to_cross_namespace_aliases])
    a=accepted.set_index('inat_taxon_id');hits=a.index.intersection(changed_ids)
    require((a.loc[hits,'independent_roi9_observations']==current.loc[hits,'independent_roi9_observations']).all(),'COL XR bridge source mismatch')
    for col in means+['independent_roi9_observations','distinct_palette_votes','ambiguous_observation_rois']:a.loc[hits,col]=fixed.loc[hits,col]
    a['identity_snapshot_patch_applied']=a.index.isin(changed_ids)
    diff=pd.DataFrame(index=changed_ids);diff['species']=current.loc[changed_ids,'species']
    for col in ['independent_roi9_observations','hue_status','dominant_hue9']:
        diff['before_'+col]=current.loc[changed_ids,col];diff['after_'+col]=fixed.loc[changed_ids,col]
    group_alias=identities.groupby(['inat_taxon_id','observation_id']).source_observation_id.nunique()
    out.mkdir(parents=True)
    files={'current_42111_identity_snapshot_patch.csv.gz':fixed.reset_index(),
           'accepted_37103_identity_snapshot_patch.csv.gz':a.reset_index(),
           'consensus_2876_alias_lower_bounds.csv':bounds.reset_index(),
           'six_taxa_observation_corrections.csv':diff.reset_index()}
    hashes={}
    for name,frame in files.items():
        p=out/name
        frame.to_csv(p,index=False,compression={'method':'gzip','mtime':0} if name.endswith('.gz') else None)
        hashes[name]=hashlib.sha256(p.read_bytes()).hexdigest()
    profile=identities.groupby('identity_status').agg(photo_rows=('photo_id','size'),strict9=('strict9','sum')).to_dict(orient='index')
    report={'schema':'mois_global_identity_consensus_audit_v1','input_zip_sha256':ZIP_SHA,
        'original_fraction_dtypes':{p:str(photos['strict9_'+p].dtype) for p in PALETTE},
        'baseline_photo_attempts':len(photos),'baseline_ROI9':len(vec),'current_total_ROI9_unchanged':467565,
        'current_taxa_retained':len(fixed),'current_accepted_species_retained':len(a),
        'identity_profile':profile,'unique_observation_ids_before':int(photos.observation_id.nunique()),
        'unique_observation_ids_after':int(identities.observation_id.nunique()),
        'within_taxon_alias_groups':int(group_alias.gt(1).sum()),
        'baseline_success_observations_before':int(before.independent_roi9_observations.sum()),
        'baseline_success_observations_after':int(after.independent_roi9_observations.sum()),
        'current_success_observations_before':int(current.independent_roi9_observations.sum()),
        'current_success_observations_after_partial_patch':int(fixed.independent_roi9_observations.sum()),
        'taxa_with_observation_count_corrected':len(changed_ids),'accepted_species_corrected':len(hits),
        'baseline_status_changes':int(before.hue_status.ne(after.hue_status).sum()),
        'baseline_dominant_hue_changes':int(before.dominant_hue9.ne(after.dominant_hue9).sum()),
        'current_dominant_hue_taxa':int(fixed.dominant_hue9.ne('').sum()),
        'current_consensus_taxa':int(consensus.sum()),'consensus_alias_robust':int(bounds.robust_to_cross_namespace_aliases.sum()),
        'consensus_uuid_only':int(bounds.numeric.eq(0).sum()),'consensus_numeric_only':int(bounds.uuid.eq(0).sum()),
        'consensus_mixed_namespace':int((bounds.numeric.gt(0)&bounds.uuid.gt(0)).sum()),
        'unresolved_numeric_success_photos':int(newer.observation_id.str.fullmatch(r'[1-9][0-9]*').sum()),
        'new_images':0,'new_ROI_inferences':0,'new_iNaturalist_API_requests':0,
        'all_catalogue_observations_fully_canonicalized':False,'species_wide_fixed_petal_colour_proven':False,
        'output_sha256':hashes}
    (out/'receipt.json').write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
    print(json.dumps(report,sort_keys=True,indent=2))
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--inputs',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();run(a.inputs,a.output)
