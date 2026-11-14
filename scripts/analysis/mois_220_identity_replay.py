"""Canonicalize the existing 220-photo QC queue using cached official identities only."""
import argparse,csv,hashlib,io,json,zipfile
from pathlib import Path
from collections import Counter,defaultdict
from scripts.analysis.mois_220_identity_annotations import kind,link,SOURCE_SHA,CSV_SHA
CACHE_SHA='34e940fd2c08b465048532f84168b5c1600b2bdad431213106db944a766ad770'

def canonicalize(rows):
    result=[dict(r) for r in rows];groups=defaultdict(list)
    for r in result:
        if r['record_status']=='resolved':r['canonical_identity_evidence']='official_api_identity'
        elif kind(r['observation_id'])=='uuid':
            r['canonical_observation_uuid']=r['observation_id'].lower()
            r['canonical_identity_evidence']='source_uuid_only'
            r['record_status']='not_fetched_uuid_batch_rejected'
        else:raise ValueError('Unresolved numeric observation ID cannot be canonicalized')
        groups[r['canonical_observation_uuid']].append(r)
    for uid,group in groups.items():
        if len({r['inat_taxon_id'] for r in group})!=1:raise ValueError('Canonical observation crosses source taxa')
        first=min(group,key=lambda r:(int(r.get('review_order',1)),int(r['photo_id'])))
        for r in group:
            r['same_observation_photo_count']=len(group)
            r['represented_by_photo_id']=first['photo_id']
            r['independent_representative']=r is first
    return result

def replay(source,metadata,output):
    out=Path(output)
    if out.exists():raise ValueError('Refuse overwrite')
    assert hashlib.sha256(Path(source).read_bytes()).hexdigest()==SOURCE_SHA
    assert hashlib.sha256(Path(metadata).read_bytes()).hexdigest()==CACHE_SHA
    with zipfile.ZipFile(source) as z:raw=z.read('mois_96_review_original_photo_IDs.csv')
    assert hashlib.sha256(raw).hexdigest()==CSV_SHA
    rows=list(csv.DictReader(io.StringIO(raw.decode())))
    assert len(rows)==220 and len({r['photo_id'] for r in rows})==220
    terms={};observations=[];cache_receipts={}
    with zipfile.ZipFile(metadata) as z:
        assert len(z.namelist())==6
        for name in sorted(z.namelist()):
            data=z.read(name);cache_receipts[name]=hashlib.sha256(data).hexdigest();data=json.loads(data)
            if name.endswith('controlled_terms.json'):
                for a in data['results']:
                    for v in a['values']:terms[(a['id'],v['id'])]=(a.get('label',str(a['id'])),v.get('label',str(v['id'])))
            else:
                assert '/numeric_id_' in name
                assert data['total_results']==len(data['results'])
                observations.extend(data['results'])
    assert len(observations)==len({o['uuid'] for o in observations})==110
    result=canonicalize(link(rows,observations,terms))
    reps=[r for r in result if r['independent_representative']]
    aliases=[r for r in result if r['same_observation_photo_count']>1]
    resolved=[r for r in result if r['record_status']=='resolved']
    annotated=[r for r in resolved if any(a['attribute']=='Flowers and Fruits' and a['value']=='Flowers' and (a['vote_score'] is None or a['vote_score']>=0) for a in json.loads(r['observation_annotations']))]
    assert len(reps)==211 and len(aliases)==18 and len(resolved)==119
    assert all(r['photo_membership']=='matched' and r['taxon_agreement']=='exact' for r in resolved)
    assert len({r['inat_taxon_id'] for r in reps})==96
    assert all(r['photo_flower_stage']=='unreviewed' and r['strict9']=='False' for r in result)
    out.mkdir(parents=True)
    files={'all_220_photos_canonical_ids.csv':result,'independent_211_observation_review.csv':reps,
           'nine_alias_pairs_18_photos.csv':aliases,'119_photo_links_to_110_flower_annotated_observations.csv':annotated}
    hashes={}
    for name,records in files.items():
        p=out/name
        with p.open('w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=list(result[0]));w.writeheader();w.writerows(records)
        hashes[name]=hashlib.sha256(p.read_bytes()).hexdigest()
    receipt={'schema':'mois_220_mixed_id_cached_canonicalization_v1',
        'source_zip_sha256':SOURCE_SHA,'source_csv_sha256':CSV_SHA,'numeric_api_cache_zip_sha256':CACHE_SHA,
        'api_capture_run':43593568216,'api_capture_result':'partial_success_numeric110_UUID_batch_HTTP422',
        'numeric_api_records':110,'source_id_kind_counts':dict(Counter(kind(r['observation_id']) for r in rows)),
        'all_source_photo_rows':220,'source_taxa_retained':96,'independent_observation_units':len(reps),
        'aliased_duplicate_observation_groups':9,'aliased_excess_photo_rows':9,
        'live_photo_membership_and_taxon_matched_rows':len(resolved),'live_verified_unique_observations':110,
        'original_uuid_only_rows':101,'flower_annotated_unique_observations':len({r['canonical_observation_uuid'] for r in annotated}),
        'flower_annotated_photo_links':len(annotated),'taxa_with_flower_annotated_observation':len({r['inat_taxon_id'] for r in annotated}),
        'photo_stage_directly_reviewed':0,'new_ROI9_successes':0,'new_original_images':0,'API_calls_during_replay':0,
        'source_metadata_sha256':cache_receipts,'output_sha256':hashes,
        'inference_limit':'Observation Flowers annotation is not photo-local flower or petal truth; source UUID-only records are not live-verified. Full atlas independence has not been re-audited.'}
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(json.dumps(receipt,indent=2,sort_keys=True))
    return receipt
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--metadata',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();replay(a.source,a.metadata,a.output)
