#!/usr/bin/env python3
"""Export every eligible original MOIS photo ID and observation from verified SQLite.

Read each table only once. Two taxon-modulo-partitioned tables are exported
separately rather than doing a source-wide 93m-row join. There is NO per-taxon,
per-observation, per-genus or total photo cap. No original photo pixels or hues
are inferred here. Join photos and observations on observation_uuid inside
the SAME partition during subsequent source/photo acquisition.
"""
from __future__ import annotations

import argparse
from collections import Counter
from contextlib import ExitStack
import csv
import gzip
import hashlib
import json
from pathlib import Path
import sqlite3
from datetime import datetime, timezone

PHOTO_COLS = ('photo_id','observation_uuid','inat_taxon_id',
              'photo_license','image_extension','photographer_id')
OBS_COLS = ('observation_uuid','inat_taxon_id','observer_id',
            'observed_on','latitude','longitude','positional_accuracy')
EXPECTED_TOTAL = 93209351
EXPECTED_TAXA = 42111
SNAPSHOT_DATE = '2026-09-27'
COUNT_SHA = 'ef26c16252dd785f2d08bde3b18938c7624d73a54386fa3d3974dceb09ed97a9'

def require(ok,msg):
    if not ok: raise ValueError(msg)

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(8<<20),b''):
            h.update(chunk)
    return h.hexdigest()

def expected_counts(path):
    counts={}
    with Path(path).open(newline='',encoding='utf-8') as f:
        reader=csv.DictReader(f)
        require({'inat_taxon_id','species','n_unique_source_photo_ids','snapshot_date'}.issubset(reader.fieldnames or []),
                'Count table schema does not match original 42111 source')
        for row in reader:
            tax=int(row['inat_taxon_id'])
            n=int(row['n_unique_source_photo_ids'])
            require(tax>0 and tax not in counts and n>=0 and bool(row['species'])
                    and row['snapshot_date']==SNAPSHOT_DATE,
                    'Original species-level photo count row invalid or duplicated')
            counts[tax]=n
    return counts

def _writer_files(root,kind,parts,columns,stack):
    base=Path(root)/kind
    base.mkdir(parents=True,exist_ok=True)
    sinks=[]
    writers=[]
    for i in range(parts):
        target=base/f'part-{i:03d}.csv.gz'
        handle=stack.enter_context(gzip.open(target,'wt',encoding='utf-8',newline='',compresslevel=1))
        writer=csv.writer(handle,lineterminator='\n')
        writer.writerow(columns)
        sinks.append(target)
        writers.append(writer)
    return sinks,writers

def export_index(db_path,count_csv,out_dir,*,parts=128,
                 expected_photos=EXPECTED_TOTAL,expected_taxa=EXPECTED_TAXA):
    db_path=Path(db_path)
    out_dir=Path(out_dir)
    count_csv=Path(count_csv)
    require(parts>=1 and parts<=256,'Invalid number of source partitions')
    require(db_path.is_file() and count_csv.is_file(),'Missing full SQLite photo-ID frame/counts')
    require(not out_dir.exists() or not any(out_dir.iterdir()),
            'Refusing to overwrite a previous source-index result')
    reference=expected_counts(count_csv)
    require(len(reference)==expected_taxa and sum(reference.values())==expected_photos,
            'Full original photo source count mismatch')
    out_dir.mkdir(parents=True,exist_ok=True)
    reference_digest=digest(count_csv)
    photo_counts=Counter()
    photo_totals=[0]*parts
    observation_totals=[0]*parts
    with sqlite3.connect(db_path.resolve().as_uri()+'?mode=ro',uri=True) as db:
        quarantined=int(db.execute('SELECT COUNT(*) FROM photos WHERE ambiguous!=0').fetchone()[0])
        count=int(db.execute('SELECT COUNT(*) FROM photos WHERE ambiguous=0').fetchone()[0])
        require(count==expected_photos,'Full photo table count differs from upstream original counts')
        with ExitStack() as stack:
            paths,w=_writer_files(out_dir,'photos',parts,PHOTO_COLS,stack)
            it=db.execute('SELECT photo_id,uuid,taxon,license,extension,photographer '
                          'FROM photos WHERE ambiguous=0')
            while True:
                batch=it.fetchmany(10000)
                if not batch:break
                for pid,uuid,taxon,license,ext,photographer in batch:
                    tax=int(taxon)
                    require(tax in reference and str(pid).isdecimal() and bool(uuid),
                            'Photo ID/taxon/observation foreign key missing')
                    shard=tax%parts
                    w[shard].writerow((pid,uuid,tax,license,ext,photographer))
                    photo_counts[tax]+=1
                    photo_totals[shard]+=1
        require(dict(photo_counts)=={k:v for k,v in reference.items() if v},
                'Photo-ID partitioning lost per-taxon original source counts')
        with ExitStack() as stack:
            obs_paths,w=_writer_files(out_dir,'observations',parts,OBS_COLS,stack)
            it=db.execute('SELECT uuid,taxon,observer,observed_on,lat,lon,accuracy FROM observations')
            while True:
                batch=it.fetchmany(10000)
                if not batch:break
                for uuid,taxon,observer,day,lat,lon,accuracy in batch:
                    tax=int(taxon)
                    require(tax in reference and bool(uuid),
                            'Observation taxon not present in original frame')
                    shard=tax%parts
                    w[shard].writerow((uuid,tax,observer,day,lat,lon,accuracy))
                    observation_totals[shard]+=1
    require(sum(photo_totals)==expected_photos and
            len(photo_counts)==sum(1 for n in reference.values() if n>0),
            'Original source photo total or taxon coverage differs')
    shards=[]
    for i in range(parts):
        pp=out_dir/'photos'/f'part-{i:03d}.csv.gz'
        po=out_dir/'observations'/f'part-{i:03d}.csv.gz'
        shards.append({'part':i,'photo_count':photo_totals[i],
                       'observation_count':observation_totals[i],
                       'photos_sha256':digest(pp),'observations_sha256':digest(po),
                       'photos_bytes':pp.stat().st_size,
                       'observations_bytes':po.stat().st_size})
    receipt={
        'schema':'mois_full_42111_93209351_uncapped_original_photo_observation_index_v1',
        'source_snapshot_date':SNAPSHOT_DATE,
        'photo_caps':None,
        'per_species_photo_cap':None,
        'per_observation_photo_cap':None,
        'same_observation_multiple_photos_allowed':True,
        'photos_emitted':sum(photo_totals),
        'taxa_in_frame':len(reference),
        'taxa_with_original_photos':len(photo_counts),
        'observations_emitted':sum(observation_totals),
        'conflicting_photo_ids_excluded':quarantined,
        'original_per_taxon_counts_sha256':reference_digest,
        'original_photos_are_globally_unique_by_source_SQLite_primary_key':True,
        'photo_colours_measured':False,
        'raw_licensed_photo_pixels_downloaded_or_redistributed':False,
        'source_photo_ids_full_without_hue_or_botanical_flowering_assertions':True,
        'shards':shards,
        'finished_utc':datetime.now(timezone.utc).isoformat(),
    }
    (out_dir/'manifest.json').write_text(json.dumps(receipt,ensure_ascii=False,sort_keys=True,indent=2)+'\n')
    print('MOIS_FULL_SOURCE_ALL_PHOTO_ID_AND_OBSERVATION_INDEX_COMPLETE',json.dumps({
        k:receipt[k] for k in ['photos_emitted','taxa_in_frame','taxa_with_original_photos',
                             'observations_emitted','conflicting_photo_ids_excluded']
    },sort_keys=True),flush=True)
    return receipt

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stage-db',type=Path,required=True)
    p.add_argument('--counts',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--partitions',type=int,default=128)
    a=p.parse_args()
    require(digest(a.counts)==COUNT_SHA,'Pinned original source count SHA changed')
    export_index(a.stage_db,a.counts,a.output,parts=a.partitions)

if __name__=='__main__':main()
