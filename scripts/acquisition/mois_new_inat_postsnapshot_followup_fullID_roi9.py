#!/usr/bin/env python3
"""Signed hash-disjoint second 64-taxon new iNaturalist photo wave, no old pixels.

Reuses original full 1,231,244-ID registry from the first independently
signed post-snapshot source, and excludes the first new wave's two IDs and
their measured raw image SHA. Prohibits auto-appending old atlas.
"""
from __future__ import annotations
import argparse, hashlib, json, shutil
from pathlib import Path
from datetime import datetime, timezone
from urllib.parse import urlsplit
import pandas as pd
from scripts.acquisition.mois_new_inat_postsnapshot_verified_original_roi9 import (
    check, write_json, UUID, ARM, EXPECTED_ATTEMPTED, WORKER_SCHEMA,
    audit as base_audit,
)
from scripts.acquisition.mois_zero9_recent_inat_new_created_photo_source import (
    ATLAS_SHA, source_sha, SCHEMA as API_SCHEMA, CC,
)
from scripts.acquisition.mois_full_index_real_original_photo_measure import sha as sha256

SECOND_OFFSET=64
FIRST_NEW_COUNT=2
def prep(meta, receipt, atlas, earlier, earlier_audit, output):
    meta=Path(meta); receipt=Path(receipt); atlas=Path(atlas)
    earlier=Path(earlier); earlier_audit=Path(earlier_audit); output=Path(output)
    check(not output.exists(), 'Refuse to reuse or overwrite an already frozen source SHA')
    now=json.loads(receipt.read_text())
    check(now.get('schema')==API_SCHEMA
          and now.get('sampled_zero_ROI9_source_taxa')==64
          and now.get('source_deterministic_hash_offset')==SECOND_OFFSET
          and now.get('original_source_zero_ROI9_taxa')==2229
          and now.get('original_photo_pixels_downloaded')==0
          and now.get('genuine_nine_fraction_ROI_measured')==0
          and now.get('old_photo_ID_1_231_244_ledger_disjointness_verified') is False
          and now.get('source_taxa_API_errors')==0,
          'Second hash-disjoint taxon source not signed or has API errors')
    check(source_sha(meta)==now['output_sha256'][meta.name], 'API metadata changed')
    check(source_sha(atlas)==ATLAS_SHA, 'Old 42,111 atlas changed')
    census=pd.read_csv(atlas,dtype={'inat_taxon_id':str})
    check(len(census)==42111 and census.inat_taxon_id.is_unique
          and int(census.after_strict9.eq(0).sum())==2229,
          'Old atlas census/identity changed')
    prevm=json.loads((earlier/'manifest.json').read_text())
    prevg=json.loads((earlier/'gate_receipt.json').read_text())
    reg=earlier/'historical_1231244_original_photo_observation_taxon_ID_registry.csv.gz'
    check(reg.is_file()
          and prevm.get('schema')==WORKER_SCHEMA
          and prevm.get('source_arm')==ARM
          and prevm.get('source_previous_attempts')==EXPECTED_ATTEMPTED
          and prevm.get('source_prior_no_reused_original_photo_or_observation') is True
          and prevm.get('new_photo_attempts_planned')==FIRST_NEW_COUNT
          and prevg.get('historical_verified_photo_ID_count')==EXPECTED_ATTEMPTED
          and prevg.get('new_photo_IDs_historical_collision_count')==0
          and prevg.get('new_observation_UUID_historical_collision_count')==0
          and prevg.get('historical_full_registry_sha256')==sha256(reg)
          and prevg.get('source_manifest_sha256')==sha256(earlier/'manifest.json')
          and prevm.get('source_previous_ID_registry_SHA256')==sha256(reg)
          and sha256(earlier/'shards'/'shard-000.csv')==prevm['photo_shard_sha256']['0'],
          'Cannot validate prior signed 1,231,244 + first 2 source gate')
    first=pd.read_csv(earlier/'shards'/'shard-000.csv',dtype=str,keep_default_na=False)
    check(len(first)==FIRST_NEW_COUNT and first.photo_id.is_unique
          and first.observation_id.is_unique and first.inat_taxon_id.is_unique,
          'Earlier two photo/observation/taxon IDs incomplete')
    pa=json.loads((earlier_audit/'receipt.json').read_text())
    audited=earlier_audit/'mois_new_inat_post20261028_actual_original_nine_roi.csv.gz'
    check(pa.get('original_licensed_photo_attempts')==FIRST_NEW_COUNT
          and pa.get('actual_image_downloaded')==FIRST_NEW_COUNT
          and pa.get('actual_genuine_strict9_numeric_success')==1
          and pa.get('actual_technical_unknowns')==1
          and pa.get('photo_observation_taxon_collision_count')==0
          and pa.get('source_photo_manifest_sha256')==sha256(earlier/'manifest.json')
          and pa.get('historical_full_ID_registry_sha256')==sha256(reg)
          and pa.get('audited_derivative_sha256')==sha256(audited)
          and pa.get('official_iNat_42111_original_frame_auto_appended') is False,
          'Earlier real pixel/audit receipt not independently complete')
    firstbytes=pa.get('new_photo_image_SHA256',{})
    check(set(firstbytes)==set(first.photo_id)
          and len(set(firstbytes.values()))==FIRST_NEW_COUNT
          and all(len(x)==64 and int(x,16)>=0 for x in firstbytes.values()),
          'Original image SHA absent from prior two photo audit')
    history=pd.read_csv(reg,dtype=str,keep_default_na=False,
                        usecols=['photo_id','observation_id','inat_taxon_id'])
    check(len(history)==EXPECTED_ATTEMPTED and history.photo_id.is_unique
          and history.photo_id.ne('').all(),
          'Historic 1,231,244 photo IDs missing or colliding')
    current=pd.read_csv(meta,dtype=str,keep_default_na=False)
    check(len(current)==now['source_taxa_with_new_post_snapshot_API_photo_opportunities']
          and len(current)==now['new_independent_original_iNaturalist_photo_IDs_seen']
          and 0<=len(current)<=64,'Second source counts mismatch')
    check(current.photo_id.is_unique and current.observation_id.is_unique
          and current.inat_taxon_id.is_unique, 'Second source self-overlap')
    check(current.inat_taxon_id.isin(census.loc[census.after_strict9.eq(0),'inat_taxon_id']).all(),
          'Second photo not a still-unmeasured source taxon')
    check(not set(current.inat_taxon_id)&set(first.inat_taxon_id),
          'Hash-disjoint 64-taxon followup somehow selected prior taxon again')
    check(not set(current.photo_id)&(set(history.photo_id)|set(first.photo_id))
          and not set(current.observation_id)&(set(history.observation_id)|set(first.observation_id)),
          'Old 1,231,244 or first new 2 photo/observation IDs reused')
    for t in current.itertuples(index=False):
        u=urlsplit(t.photo_url_large)
        check(t.photo_id.isdigit() and UUID.fullmatch(t.observation_id) is not None
              and t.photo_license in CC
              and u.scheme=='https' and u.hostname=='inaturalist-open-data.s3.amazonaws.com'
              and u.path.startswith('/photos/'+t.photo_id+'/')
              and u.path.lower().endswith(('.jpg','.jpeg','.png')),
              'Bad CC license, raw photo URL, photo ID or observation UUID')
        created=datetime.fromisoformat(t.source_observation_created_at.replace('Z','+00:00'))
        check(created.tzinfo is not None and created.astimezone(timezone.utc)>=datetime(2026,10,28,tzinfo=timezone.utc),
              'Not a genuinely post-snapshot new source observation')
    output.mkdir(parents=True)
    filename='historical_1231244_original_photo_observation_taxon_ID_registry.csv.gz'
    shutil.copyfile(reg,output/filename)
    (output/'shards').mkdir()
    if len(current):
        q=current.copy()
        q['measurement_id']='MOIS-POST20261028-OFFSET64-'+q.photo_id
        q['image_filename']=q.measurement_id+'.jpg'
        q['shard']='0'
        shard=output/'shards'/'shard-000.csv'
        q.to_csv(shard,index=False,lineterminator='\n')
        hashes={'0':sha256(shard)}
    else: hashes={}
    m={
        'schema':WORKER_SCHEMA,'source_arm':ARM,'colour_measured':False,
        'source_previous_attempts':EXPECTED_ATTEMPTED,
        'source_previous_ID_registry_SHA256':sha256(output/filename),
        'source_prior_no_reused_original_photo_or_observation':True,
        'source_first_new_photo_IDs':first.photo_id.tolist(),
        'source_first_new_observation_UUIDs':first.observation_id.tolist(),
        'source_first_new_image_SHA256':firstbytes,
        'source_first_new_real_audit_SHA256':sha256(earlier_audit/'receipt.json'),
        'source_post_snapshot_API_receipt_SHA256':source_sha(receipt),
        'source_post_snapshot_API_metadata_SHA256':source_sha(meta),
        'source_fixed_hash_offset':SECOND_OFFSET,
        'new_photo_attempts_planned':len(current),
        'work_units':int(len(current)>0),
        'worker_ids':[0] if len(current) else [],
        'photo_shard_sha256':hashes,
        'raw_image_pixels_downloaded':False,'new_nineROI_measured':False,
        'official_iNat_42111_original_frame_auto_appended':False,
        'species_wide_fixed_colour_proven':False
    }
    write_json(output/'manifest.json',m)
    g={
        'schema':'mois_post_snapshot_second_64taxa_firstplusallhistoricID_v1',
        'historical_verified_photo_ID_count':len(history),
        'new_photo_IDs_historical_collision_count':0,
        'new_observation_UUID_historical_collision_count':0,
        'previous_new_photo_source_count':FIRST_NEW_COUNT,
        'prior_new_original_raw_image_SHA_count':len(firstbytes),
        'new_CC_photo_opportunities':len(current),
        'source_manifest_sha256':sha256(output/'manifest.json'),
        'historical_full_registry_sha256':sha256(output/filename),
        'source_only':True,'strict9_measured':0,'new_photo_pixels_downloaded':0
    }
    write_json(output/'gate_receipt.json',g)
    print('MOIS_HASH_DISJOINT_SECOND64_REAL_1231244_PLUS2_ID_GATE',
          json.dumps({k:g[k] for k in ('historical_verified_photo_ID_count',
                 'previous_new_photo_source_count','new_CC_photo_opportunities',
                 'new_photo_IDs_historical_collision_count',
                 'new_observation_UUID_historical_collision_count')},sort_keys=True),flush=True)

def audit_second(source, worker, output):
    source=Path(source);worker=Path(worker);output=Path(output)
    m=json.loads((source/'manifest.json').read_text())
    g=json.loads((source/'gate_receipt.json').read_text())
    check(g.get('schema')=='mois_post_snapshot_second_64taxa_firstplusallhistoricID_v1'
          and m.get('source_fixed_hash_offset')==SECOND_OFFSET
          and m.get('source_first_new_photo_IDs') is not None
          and g.get('source_manifest_sha256')==sha256(source/'manifest.json')
          and len(m['source_first_new_image_SHA256'])==FIRST_NEW_COUNT,
          'Missing signed second-wave first+history ID evidence')
    src=pd.read_csv(source/'shards'/'shard-000.csv',dtype=str,keep_default_na=False)
    oldids=set(m['source_first_new_photo_IDs'])
    oldobs=set(m['source_first_new_observation_UUIDs'])
    check(not (set(src.photo_id)&oldids) and not(set(src.observation_id)&oldobs),
          'Second source duplicates first independent photo or observation')
    raw=worker/'shard_000'/'photos.csv'
    receipt=worker/'shard_000'/'receipt.json'
    check(raw.is_file() and receipt.is_file(),'Second wave original pixels missing')
    photos=pd.read_csv(raw,dtype=str,keep_default_na=False)
    prior_raw_sha=set(m['source_first_new_image_SHA256'].values())
    check(not set(photos.image_sha256[photos.image_sha256.ne('')]) & prior_raw_sha,
          'Raw pixel SHA collision with first two new actual originals')
    base_audit(source,worker,output)
    ar=json.loads((output/'receipt.json').read_text())
    ar.update({
        'additional_first_new_photo_IDs_excluded':FIRST_NEW_COUNT,
        'additional_first_new_raw_image_SHA_collisions':0,
        'additional_first_new_observation_UUID_collisions':0,
        'full_historical_raw_image_SHA_collision_check_possible':False,
        'official_iNat_42111_original_frame_auto_appended':False
    })
    write_json(output/'receipt.json',ar)
    print('MOIS_POSTSNAPSHOT_SECOND64_AUTHENTIC_ROI9_AUDIT',
          json.dumps({k:ar[k] for k in ('original_licensed_photo_attempts',
                                      'actual_image_downloaded','actual_genuine_strict9_numeric_success',
                                      'actual_technical_unknowns')},sort_keys=True),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser()
    sub=p.add_subparsers(dest='command',required=True)
    a=sub.add_parser('gate')
    for k in ('meta','receipt','atlas','earlier','earlier-audit','output'):a.add_argument('--'+k,required=True)
    b=sub.add_parser('audit')
    for k in ('source','worker','output'):b.add_argument('--'+k,required=True)
    v=p.parse_args()
    if v.command=='gate':prep(v.meta,v.receipt,v.atlas,v.earlier,v.earlier_audit,v.output)
    else:audit_second(v.source,v.worker,v.output)
