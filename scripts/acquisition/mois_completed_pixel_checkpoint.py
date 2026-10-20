#!/usr/bin/env python3
"""Recover only sealed image outcomes from active waves; never acquire images.

A checkpoint is a partial compute-completion subset, not a completed probability
cohort. Preserve original unknowns, source identities and baseline hue labels.
The original full-cohort collectors and their releases remain untouched.
"""
from __future__ import annotations
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import re

HUES=('white','yellow_orange','red_pink','blue_purple')
BASE_FILE='mois_source_verified_FCP447500_and_93m_p0002_combined_unique_photo_ROIs.parquet'
BASE_RECEIPT='mois_original_photo_RGB_FCP447500_plus_first93m_id_union_receipt.json'
COUNTERS=('attempted','downloaded','classifiable','acquisition_failed','unclassifiable_acquired')
WAVES={
    'p0004':(39353343403,18119,'c9e2cb8b80e696bf34c0e5f1cdba2bd3579e64c5126781c9aa75989222a437eb'),
    'p0008':(39371199127,36375,'0038024ee40adfbc71aa427a112948ec5ce53cb53a7aae9879c1b32bc06eabde'),
    'shadow':(39369517433,569,'285e9f86e13f4d9330c71ee3719898268338069e3db59f639c99bd88b7e84fd1'),
    'p0016':(39374811084,73872,'16baf617396dd7b388debc65fd25eb8f9041f49ca5e6f78d0b3e2898f6dfeed5'),
    'p0032':(39401484641,148071,'3628794fbb70f9e1515dabea6db731862265c0ad6c251732cfd9ab4df4badce8'),
    'geo05':(39417990675,22278,'3b80efe07287c026a34eec7636038e61fa5aa4e16d49aeb139a533bb46c8feda'),
}

def require(ok,message):
    if not ok:raise ValueError(message)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for part in iter(lambda:f.read(1<<20),b''):h.update(part)
    return h.hexdigest()

def read_rows(path):
    with Path(path).open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f))

def identity(rows,key,label):
    result={r[key]:r for r in rows}
    require(len(result)==len(rows) and all(result),label+' identity missing or duplicated')
    return result

def collect_wave(source,workers,method):
    source,workers=Path(source),Path(workers)
    info=json.loads((source/'manifest.json').read_text())
    require(info['schema']=='mois_full_source_observation_sampling_cohort_v1'
            and info['colour_measured'] is False and info['full_frame_photo_ids']==93209351,
            'Unknown original source plan')
    units=info['worker_ids']
    require(units==list(range(info['work_units'])), 'Invalid source unit sequence')
    targets={};all_source=[]
    for u in units:
        f=source/'shards'/f'shard-{u:03d}.csv'
        require(sha(f)==info['photo_shard_sha256'][str(u)],'Source SHA mismatch')
        rows=read_rows(f)
        require(0<len(rows)<=256 and all(r['shard']==str(u) for r in rows),'Invalid source count/unit')
        identity(rows,'measurement_id','Source');targets[u]=rows;all_source.extend(rows)
    require(len(all_source)==info['new_photo_attempts_planned'],'Original source count mismatch')
    identity(all_source,'photo_id','Source photo');identity(all_source,'measurement_id','Source measurement')
    found={}
    for p in workers.rglob('receipt.json'):
        if not re.fullmatch(r'shard_[0-9]+',p.parent.name):continue
        u=int(p.parent.name.split('_')[1])
        require(u in targets,'Unrecognized completed source unit')
        require(u not in found,'Duplicate completed unit')
        found[u]=p
    complete=[];counters=Counter({k:0 for k in COUNTERS});audits=[]
    for u in sorted(found):
        rp=found[u];r=json.loads(rp.read_text());folder=rp.parent
        require(r['schema']=='mois_full_original93m_real_RGB_photo_unit_v1'
                and r['shard']==u and r['photo_colours_measured'] is True
                and r['raw_licensed_image_bytes_persisted'] is False,'Not a sealed image outcome')
        require(r['method']==method,'Frozen original method changed')
        require(r['source_sha256']==info['photo_shard_sha256'][str(u)]
                and sha(folder/'photos.csv')==r['photos_sha256'],'Completed original photo SHA mismatch')
        terminal=folder/'sealed'/f'terminal_p{u:03d}_b0000000.csv'
        require(terminal.is_file(),'Missing completed image terminal')
        require(sha(terminal)==r['terminal_sha256'],'Original terminal SHA mismatch')
        photos=read_rows(folder/'photos.csv');ends=identity(read_rows(terminal),'measurement_id','Terminal')
        by_mid=identity(photos,'measurement_id','Completed photo')
        expected=identity(targets[u],'measurement_id','Original source')
        require(set(by_mid)==set(expected)==set(ends),'Completed source identity coverage changed')
        c=Counter({k:0 for k in COUNTERS})
        for mid,original in expected.items():
            row=by_mid[mid];end=ends[mid]
            require(all(row.get(k)==v for k,v in original.items()),'Original photo/source identity changed')
            cols=['morph','measurement_status','image_sha256']+[k for k in end if k.startswith('flower_fraction_')]
            require(all(row.get(k)==end.get(k) for k in cols),'Original terminal hue or image changed')
            known=row['morph'] in HUES
            failed=row['measurement_status']=='image_acquisition_failed'
            require(row['morph'] in (*HUES,'mixed_uncertain')
                    and known==(row['measurement_status']=='classified_four_state_morph')
                    and row['classified_four_hue'].lower()==str(known).lower()
                    and not (known and failed),'Unknown image cannot become a coloured flower')
            require(failed or bool(re.fullmatch(r'[0-9a-f]{64}',row['image_sha256'])),
                    'Acquired original image SHA missing')
            c['attempted']+=1;c['downloaded']+=not failed;c['classifiable']+=known
            c['acquisition_failed']+=failed;c['unclassifiable_acquired']+=not known and not failed
            complete.append(row)
        require(all(c[k]==r[k] for k in COUNTERS),'Original measured outcome count mismatch')
        counters.update(c)
        audits.append({'unit':u,'receipt_sha256':sha(rp),'photos_sha256':r['photos_sha256'],
                       'terminal_sha256':r['terminal_sha256']})
    identity(complete,'photo_id','Checkpoint photo') if complete else None
    missing=[u for u in units if u not in found]
    return complete,{'planned_photo_attempts':len(all_source),**dict(counters),
        'pending_photo_attempts':len(all_source)-len(complete),
        'completed_units':sorted(found),'missing_units':missing,'wave_complete':not missing,
        'source_manifest_sha256':sha(source/'manifest.json'),'completed_unit_audit':audits}

PROJECTION=('photo_id','measurement_id','observation_id','inat_taxon_id','species',
    'observed_on','latitude','longitude','positional_accuracy','morph','measurement_status','image_sha256')

def union_projection(old,new):
    import pandas as pd
    require({'photo_id','morph','classified','white','arm'}<=set(old),'Baseline columns missing')
    require({'photo_id','morph','measurement_status','checkpoint_wave'}<=set(new),'New measured columns missing')
    require(old.photo_id.astype(str).is_unique and new.photo_id.astype(str).is_unique,'Duplicate measured photo ID')
    require(not set(old.photo_id.astype(str)) & set(new.photo_id.astype(str)), 'Measured photo overlap')
    require(old.classified.eq(old.morph.isin(HUES)).all() and old.white.eq(old.morph.eq('white')).all(),
            'Baseline original hue flags changed')
    frames=[]
    for frame,is_old in ((old,True),(new,False)):
        d=pd.DataFrame(index=frame.index)
        for col in PROJECTION:
            d[col]=frame[col].fillna('').astype(str) if col in frame else ''
        require(d.photo_id.str.fullmatch(r'\d+').all(),'Invalid actual original photo ID')
        # SHA-verified legacy sources have other original unclassified labels.
        # Preserve them verbatim; only NEW ROI-v4 rows use the five-state schema.
        if not is_old:
            require(d.morph.isin((*HUES,'mixed_uncertain')).all(),'Unknown source hue state')
        d['classified']=d.morph.isin(HUES);d['white']=d.morph.eq('white')
        d['technical_failure']=frame['technical_failure'].astype(bool) if is_old and 'technical_failure' in frame else d.measurement_status.eq('image_acquisition_failed')
        d['source_arm']=frame['arm'].fillna('').astype(str) if is_old else frame.get('source_arm',pd.Series('',index=frame.index)).fillna('').astype(str)
        d['checkpoint_wave']='baseline' if is_old else frame.checkpoint_wave.astype(str)
        frames.append(d.reset_index(drop=True))
    union=pd.concat(frames,ignore_index=True)
    require(union.iloc[:len(old)].reset_index(drop=True).equals(frames[0]),'Baseline projected values changed')
    require(union.photo_id.is_unique,'Duplicate checkpoint photo identity')
    return union,{'baseline_photo_attempts':len(old),'new_completed_photo_attempts':len(new),
        'combined_unique_photo_attempts':len(union),'baseline_four_hue_classified':int(old.classified.sum()),
        'new_actual_four_hue_classified':int(frames[1].classified.sum()),
        'combined_four_hue_classified':int(union.classified.sum()),
        'new_acquisition_failed':int(frames[1].technical_failure.sum()),
        'baseline_projected_rows_preserved':True}

def run(inputs,output):
    import pandas as pd
    from scripts.analysis.mois_reassemble_original_new4096_measurements import METHOD
    inputs,output=Path(inputs),Path(output)
    require(not output.exists() or not any(output.iterdir()),'Refuse to overwrite a checkpoint')
    base=inputs/'baseline'/BASE_FILE;doc=json.loads((inputs/'baseline'/BASE_RECEIPT).read_text())
    require(doc['schema']=='mois_original_photo_ROI_FCP447500_plus_first_full93m_p0002_ID_union_v1'
            and doc['all_genuinely_attempted_unique_original_photo_IDs']==612483
            and doc['all_original_photo_four_hue_classified']==295801
            and doc['duplicate_source_photo_IDs_in_final_union']==0
            and sha(base)==doc['completed_output_sha256'][base.name], 'Baseline source receipt/SHA mismatch')
    old=pd.read_parquet(base)
    require(len(old)==612483 and int(old.classified.sum())==295801,'Baseline actual row count changed')
    all_new=[];waves={};snapshots={}
    for wave,(run_id,planned,source_sha) in WAVES.items():
        root=inputs/'waves'/wave
        require(sha(root/'source/manifest.json')==source_sha,'Pinned wave source manifest SHA mismatch')
        rows,receipt=collect_wave(root/'source',root/'workers',METHOD)
        require(receipt['planned_photo_attempts']==planned,'Frozen wave planned count changed')
        for row in rows:row['checkpoint_wave']=wave
        all_new.extend(rows);receipt['source_run_id']=run_id;waves[wave]=receipt
        snapshots[wave]=json.loads((root/'artifact_snapshot.json').read_text())
    require(bool(all_new),'No completed genuine image units available')
    new=pd.DataFrame(all_new);combined,counts=union_projection(old,new)
    output.mkdir(parents=True,exist_ok=True)
    new.to_parquet(output/'new_completed_original_photo_ROIs.parquet',index=False,compression='zstd')
    combined.to_parquet(output/'combined_actual_photo_checkpoint.parquet',index=False,compression='zstd')
    summary=combined.groupby(['checkpoint_wave','morph'],dropna=False).size().rename('photo_attempts').reset_index()
    summary.to_csv(output/'completed_photo_hue_counts.csv',index=False)
    report={'schema':'mois_sealed_actual_pixel_checkpoint_v1',**counts,'waves':waves,
        'status':'COMPLETE_CAPTURE_OF_AVAILABLE_SEALED_UNITS_NOT_COMPLETE_ACQUISITION_COHORT',
        'source_artifact_snapshots':snapshots,'baseline_parquet_sha256':sha(base),
        'baseline_receipt_sha256':sha(inputs/'baseline'/BASE_RECEIPT),
        'image_downloads_or_colour_reinference_performed':False,
        'partial_completion_is_not_representative_probability_sample':True,
        'ecological_analysis_run':False,'original_full_wave_releases_untouched':True,
        'combined_table_is_metadata_projection_full_originals_retained_in_source_releases':True,
        'output_sha256':{p.name:sha(p) for p in sorted(output.iterdir()) if p.is_file()}}
    (output/'checkpoint_receipt.json').write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n')
    print('MOIS_ACTUAL_SEALED_PIXELS_RECOVERED_WITHOUT_IMAGE_REDOWNLOAD',json.dumps(counts,sort_keys=True),flush=True)
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inputs',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.inputs,a.output)
