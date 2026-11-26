"""Read-only identity/annotation audit: no image bytes or inferred flower labels."""
import argparse,csv,datetime,hashlib,io,json,re,time,urllib.request,urllib.error,zipfile
from pathlib import Path
from collections import Counter
SOURCE_SHA='eb2a8fa4dd15655b9f1cd17910b45c00408f871ca5db44cf559cf37e23011a94'
CSV_SHA='b94432e97f97b467dd012d7f31d4a0db250fbae7b518d91395a9bb557c8a8c8c'

def kind(value):
    if re.fullmatch(r'[1-9][0-9]*',str(value)): return 'numeric_id'
    if re.fullmatch(r'[0-9a-f]{8}-(?:[0-9a-f]{4}-){3}[0-9a-f]{12}',str(value),re.I): return 'uuid'
    raise ValueError('Malformed source observation identifier: '+str(value))

def link(rows,observations,terms):
    index={}
    for obs in observations:
        for key in (str(obs['id']),str(obs['uuid'])):
            if key in index and index[key]['uuid']!=obs['uuid']:raise ValueError('API identity collision')
            index[key]=obs
    result=[]
    for row in rows:
        r=dict(row,source_observation_id_kind=kind(row['observation_id']),
            canonical_observation_id='',canonical_observation_uuid='',current_taxon_id='',
            current_taxon_name='',taxon_agreement='unverified',record_status='not_returned',
            photo_membership='unverified',photo_license_code='',photo_url_from_api='',
            observation_annotations='',photo_flower_stage='unreviewed')
        obs=index.get(row['observation_id'])
        if obs:
            tax=obs.get('taxon') or {};photos={str(p['id']):p for p in obs.get('photos',[])}
            for op in obs.get('observation_photos',[]):
                p=op.get('photo') or {}
                if p.get('id') is not None:photos[str(p['id'])]=p
            r.update(canonical_observation_id=str(obs['id']),canonical_observation_uuid=str(obs['uuid']),
                current_taxon_id=str(tax.get('id','')),current_taxon_name=tax.get('name',''),
                taxon_agreement='exact' if str(tax.get('id',''))==row['inat_taxon_id'] else 'changed_or_unresolved',
                record_status='resolved',photo_membership='matched' if row['photo_id'] in photos else 'not_in_current_record')
            photo=photos.get(row['photo_id'],{})
            r['photo_license_code']=photo.get('license_code') or ''
            r['photo_url_from_api']=photo.get('url') or ''
            anns=[]
            for a in obs.get('annotations',[]):
                key=(a.get('controlled_attribute_id'),a.get('controlled_value_id'))
                attr,value=terms.get(key,(str(key[0]),str(key[1])))
                anns.append({'attribute':attr,'value':value,'vote_score':a.get('vote_score')})
            r['observation_annotations']=json.dumps(anns,ensure_ascii=False,sort_keys=True)
        result.append(r)
    return result

def request_json(url,path):
    # Fail closed on throttling or network error; do not convert outage into absence.
    req=urllib.request.Request(url,headers={'User-Agent':'MOIS-read-only-identity-audit/1.0','Accept':'application/json'})
    with urllib.request.urlopen(req,timeout=45) as response:raw=response.read()
    parsed=json.loads(raw);path.write_bytes(raw)
    return parsed,{'url':url,'sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw)}

def run(source,out):
    source=Path(source);out=Path(out)
    if out.exists():raise ValueError('Refuse overwrite')
    assert hashlib.sha256(source.read_bytes()).hexdigest()==SOURCE_SHA
    with zipfile.ZipFile(source) as z:raw=z.read('mois_96_review_original_photo_IDs.csv')
    assert hashlib.sha256(raw).hexdigest()==CSV_SHA
    rows=list(csv.DictReader(io.StringIO(raw.decode())))
    assert len(rows)==220 and len({r['photo_id'] for r in rows})==220
    kinds=Counter(kind(r['observation_id']) for r in rows)
    assert kinds=={'uuid':110,'numeric_id':110}
    out.mkdir(parents=True);cache=out/'metadata_cache';cache.mkdir()
    termdata,receipt=request_json('https://api.inaturalist.org/v1/controlled_terms',cache/'controlled_terms.json')
    requests=[receipt];terms={}
    for a in termdata.get('results',[]):
        for v in a.get('values',[]):terms[(a['id'],v['id'])]=(a.get('label',str(a['id'])),v.get('label',str(v['id'])))
    observations=[]
    # Batch homogeneous namespaces, as requested by iNaturalist API recommended practices.
    for namespace in ('numeric_id','uuid'):
        keys=sorted({r['observation_id'] for r in rows if kind(r['observation_id'])==namespace})
        for start in range(0,len(keys),25):
            time.sleep(1.2)
            url='https://api.inaturalist.org/v1/observations/'+','.join(keys[start:start+25])+'?per_page=200'
            data,receipt=request_json(url,cache/f'{namespace}_{start:03d}.json')
            requests.append(receipt);observations.extend(data.get('results',[]))
    results=link(rows,observations,terms)
    p=out/'photo_observation_identity_annotations.csv'
    with p.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(results[0]));w.writeheader();w.writerows(results)
    canonical=Counter(r['canonical_observation_uuid'] for r in results if r['record_status']=='resolved')
    collisions={k:v for k,v in canonical.items() if v>1}
    summary={'schema':'mois_220_source_identity_annotations_v1','source_artifact_SHA256':SOURCE_SHA,
        'source_csv_SHA256':CSV_SHA,'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'source_rows':220,'source_id_kind_counts':dict(kinds),
        'record_status_counts':dict(Counter(r['record_status'] for r in results)),
        'photo_membership_counts':dict(Counter(r['photo_membership'] for r in results)),
        'taxon_agreement_counts':dict(Counter(r['taxon_agreement'] for r in results)),
        'canonical_unique_observations':len(canonical),'canonical_duplicates':collisions,
        'observation_annotation_patterns':dict(Counter(r['observation_annotations'] for r in results)),
        'new_image_downloads':0,'new_roi9':0,'photo_flower_stage_reviewed':0,
        'requests':requests,'output_csv_SHA256':hashlib.sha256(p.read_bytes()).hexdigest()}
    (out/'receipt.json').write_text(json.dumps(summary,indent=2,sort_keys=True,ensure_ascii=False)+'\n')
    print('MOIS_IDENTITY_ANNOTATION_AUDIT',json.dumps(summary,sort_keys=True,ensure_ascii=False))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();run(a.source,a.output)
