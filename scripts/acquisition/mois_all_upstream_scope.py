"""MOIS full upstream scope; counts are never substituted for measured hues.

Use the complete audited 2026-09-27 source, not the 32768/photo-species
subsample. The existing FCP bulk importer is reused only for source IDs.
No FCP within-species statistical analyses are run.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import shutil
import sqlite3

COUNT_SHA='ef26c16252dd785f2d08bde3b18938c7624d73a54386fa3d3974dceb09ed97a9'
ARCHIVE_SHA='bfc1a97c436bf5989fccd43f12392953128489b76b92cea24da335ad73da4d8e'
FRAME_SHA='5d871bd1f1190c68fd513bb527e6c93de1606dc35ff8cfff5f6187850dc382cc'
TOTAL=93209351
TAXA=42111
COUNTS='all_42111_species_photo_opportunity_counts.csv'

def require(ok, message):
    if not ok: raise ValueError(message)

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8<<20),b''):h.update(b)
    return h.hexdigest()

def read_counts(path):
    rows=[]; seen=set(); total=0
    with Path(path).open(newline='',encoding='utf-8') as f:
        for row in csv.DictReader(f):
            tid=int(row['inat_taxon_id']); n=int(row['n_unique_source_photo_ids'])
            require(tid not in seen,'Duplicate source taxon')
            require(tid>0 and n>=0 and row['species'] and row['snapshot_date']=='2026-09-27','Invalid upstream count')
            seen.add(tid); total+=n
            rows.append({**row,'inat_taxon_id':tid,'n_unique_source_photo_ids':n})
    return rows,total

def only(root,name):
    found=list(Path(root).rglob(name))
    require(len(found)==1,'Expected one source file: '+name)
    return found[0]

def database_counts(path):
    with sqlite3.connect(Path(path).resolve().as_uri()+'?mode=ro',uri=True) as db:
        counts={int(t):int(n) for t,n in db.execute('SELECT taxon,COUNT(*) FROM photos WHERE ambiguous=0 GROUP BY taxon')}
        quarantined=db.execute('SELECT COUNT(*) FROM photos WHERE ambiguous!=0').fetchone()[0]
    return counts,int(quarantined)

def verify_counts(root,out):
    path=only(root,COUNTS)
    r=json.loads(only(root,'photo_opportunity_counts_only_receipt.json').read_text())
    require(digest(path)==COUNT_SHA==r['counts_file_sha256'],'Upstream count SHA differs')
    rows,n=read_counts(path)
    require(n==TOTAL==r['n_target_taxon_licensed_photos'] and len(rows)==TAXA==r['source_species'],'Not the full upstream population')
    require(r['source_species_frame_sha256']==FRAME_SHA and r['source_archive_attestation']['source_archive_sha256']==ARCHIVE_SHA and r['photo_cap'] is None,'Source archive or scope differs')
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    record={
        'schema':'mois_all_upstream_scope_v1',
        'source_snapshot_date':'2026-09-27',
        'source_photo_candidates':n,'source_frame_taxa':len(rows),
        'taxa_with_photos':sum(x['n_unique_source_photo_ids']>0 for x in rows),
        'photo_cap':None,'per_species_photo_cap':None,'per_observation_photo_cap':None,
        'new_genus_only':False,'new_species_only':False,
        'retain_existing_other_source_measurements':True,
        'same_taxon_and_observation_multiple_photos_allowed':True,
        'source_counts_sha256':COUNT_SHA,'source_archive_sha256':ARCHIVE_SHA,
        'source_frame_sha256':FRAME_SHA,
        'source_filters':'Use the audited upstream licence, research-grade and coordinate eligibility unchanged; conflicting photo IDs quarantined',
        'reported_conflicting_photo_ids':r['n_conflicting_photo_ids_quarantined'],
        'all_photo_ids_indexed':False,'all_image_pixels_measured':False,
        'ecological_analysis_run':False,
    }
    shutil.copyfile(path,out/COUNTS)
    (out/'scope.json').write_text(json.dumps(record,indent=2,sort_keys=True)+'\n')
    print('MOIS_ALL_UPSTREAM_SCOPE',json.dumps(record),flush=True)
    return record

def admit_stage(stage_root,counts_root,out):
    scope=verify_counts(counts_root,out)
    rp=only(stage_root,'fcp_staged_sqlite_source_receipt.json')
    r=json.loads(rp.read_text());db=rp.parent/r['sqlite_basename']
    require(r['schema']=='fcp_staged_sqlite_photo_universe_v1' and r['status']=='SOURCE_SQLITE_STAGED_SHARDS_NOT_YET_EMITTED','Source is not a complete staged photo-ID table')
    require(r['source_species']==TAXA and r['n_target_taxon_licensed_photos']==TOTAL and r['source_species_frame_sha256']==FRAME_SHA,'Wrong upstream subset')
    require(r['verified_original_dated_archive']['source_archive_sha256']==ARCHIVE_SHA and r['original_photo_ids_deduplicated_globally'],'Different source or duplicate photo IDs')
    require(digest(db)==r['sqlite_sha256'],'Staged database checksum mismatch')
    got,quarantine=database_counts(db)
    rows,_=read_counts(only(counts_root,COUNTS))
    expected={x['inat_taxon_id']:x['n_unique_source_photo_ids'] for x in rows if x['n_unique_source_photo_ids']}
    require(got==expected and quarantine==scope['reported_conflicting_photo_ids'],'Full photo IDs disagree with original per-taxon counts')
    scope.update(all_photo_ids_indexed=True,source_sqlite_sha256=r['sqlite_sha256'],source_sqlite_bytes=db.stat().st_size,
                 classification_state='PENDING_OR_REUSE_PRIOR_MATCHING_PHOTO_RESULT; not newly classified',
                 frozen_model_results_must_be_reused_before_new_inference=True)
    out=Path(out);shutil.copyfile(rp,out/rp.name)
    (out/'scope.json').write_text(json.dumps(scope,indent=2,sort_keys=True)+'\n')
    print('MOIS_ALL_UPSTREAM_IDS_ADMITTED',json.dumps(scope),flush=True)
    return scope

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--counts',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--stage',type=Path)
    a=p.parse_args()
    if a.stage:admit_stage(a.stage,a.counts,a.output)
    else:verify_counts(a.counts,a.output)
