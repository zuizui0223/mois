"""One licensed existing photo per unresolved taxon for bounded visual QA.
No new photo IDs, no original-resolution fetch, no ROI inference or colour labels.
All image bytes remain in short-retention artifact, never committed to repository.
"""
import argparse,csv,hashlib,io,json,re,time,urllib.request,zipfile
from collections import Counter
from pathlib import Path
from urllib.parse import urlsplit,urlunsplit
SOURCE_SHA='48224a2f74a56d3c91f370e1a402e6e730d6294dbab623a5cfd447f7eb37054f'
CACHE_SHA='34e940fd2c08b465048532f84168b5c1600b2bdad431213106db944a766ad770'

def image_url(row):
    if any(row.get(k)!=v for k,v in [('record_status','resolved'),('photo_membership','matched'),('taxon_agreement','exact')]):
        raise ValueError('Unverified photo membership or taxonomy')
    if row.get('photo_license_code') not in {'cc0','cc-by','cc-by-sa','cc-by-nc','cc-by-nc-sa','cc-by-nd','cc-by-nc-nd'}:
        raise ValueError('No supported licence')
    u=urlsplit(row['photo_url_from_api'])
    if u.scheme!='https' or u.netloc!='inaturalist-open-data.s3.amazonaws.com' or u.query or u.fragment:
        raise ValueError('Unapproved photo URL')
    m=re.fullmatch(r'/photos/([1-9][0-9]*)/square\.(jpg|jpeg|png)',u.path)
    if not m or m.group(1)!=row['photo_id']:raise ValueError('Photo identity mismatch')
    return urlunsplit((u.scheme,u.netloc,u.path.replace('/square.','/medium.'),'',''))

def run(source,cache,output):
    out=Path(output)
    if out.exists():raise ValueError('Refuse overwrite')
    assert hashlib.sha256(Path(source).read_bytes()).hexdigest()==SOURCE_SHA
    assert hashlib.sha256(Path(cache).read_bytes()).hexdigest()==CACHE_SHA
    with zipfile.ZipFile(source) as z:
        rows=list(csv.DictReader(io.StringIO(z.read('all_220_photos_canonical_ids.csv').decode())))
    metadata={}
    with zipfile.ZipFile(cache) as z:
        for name in z.namelist():
            if 'numeric_id_' not in name:continue
            for obs in json.loads(z.read(name))['results']:
                for p in obs.get('photos',[]):metadata[str(p['id'])]=p
    eligible=[]
    for r in rows:
        try:u=image_url(r)
        except ValueError:continue
        p=metadata[r['photo_id']]
        assert p['url']==r['photo_url_from_api'] and p['license_code']==r['photo_license_code']
        if p.get('hidden') or p.get('flags') or not p.get('attribution'):continue
        eligible.append(dict(r,review_url=u,attribution=p['attribution']))
    eligible.sort(key=lambda r:(-int(r['after_photo_attempts']),r['inat_taxon_id'],int(r['review_order']),int(r['photo_id'])))
    chosen=[];seen=set()
    for r in eligible:
        if r['inat_taxon_id'] in seen:continue
        seen.add(r['inat_taxon_id']);chosen.append(r)
    assert len(chosen)==90 and len({r['canonical_observation_uuid'] for r in chosen})==90
    out.mkdir();(out/'images').mkdir()
    results=[]
    for n,row in enumerate(chosen,1):
        rec={k:row[k] for k in ['photo_id','inat_taxon_id','species','canonical_observation_id','canonical_observation_uuid','after_photo_attempts','photo_page_url','photo_license_code','attribution','review_url']}
        rec.update(review_sequence=n,retrieval_status='unattempted',image_sha256='',byte_count=0,error='',visual_stage='unreviewed',strict9_promoted=False)
        try:
            req=urllib.request.Request(row['review_url'],headers={'User-Agent':'MOIS-bounded-noncommercial-image-QA/1.0'})
            with urllib.request.urlopen(req,timeout=20) as resp:
                if not resp.headers.get('Content-Type','').startswith('image/'):raise ValueError('Not an image response')
                data=resp.read(8000001)
                if len(data)>8000000:raise ValueError('Image exceeds bounded derivative byte cap')
            if not (data.startswith(b'\xff\xd8\xff') or data.startswith(b'\x89PNG')):raise ValueError('Unknown image signature')
            ext=row['review_url'].rsplit('.',1)[-1]
            name='images/'+row['photo_id']+'.'+ext
            (out/name).write_bytes(data)
            rec.update(retrieval_status='downloaded',image_sha256=hashlib.sha256(data).hexdigest(),byte_count=len(data),image_path=name)
        except Exception as e:
            rec.update(retrieval_status='request_failed',error=str(e)[:400],image_path='')
        results.append(rec)
        print('MOIS_VISUAL_QA_SOURCE',n,rec['photo_id'],rec['retrieval_status'],flush=True)
        time.sleep(1.05)
    p=out/'image_source_manifest.csv'
    with p.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(results[0]));w.writeheader();w.writerows(results)
    (out/'LICENSE_AND_SCOPE.txt').write_text('Existing iNaturalist photographs, native medium derivatives. Individual attribution, CC licence and photo page are in image_source_manifest.csv. No pixels committed to Git; artifact retention 1 day. Research visual QA only. No numerical ROI9, new independent image or species-fixed colour claim. Review results are AI-assisted triage, not an expert botanical ground-truth audit.\n')
    receipt=dict(schema='mois_bounded_existing90_visual_QA_sources_v1',source_zip_sha256=SOURCE_SHA,metadata_zip_sha256=CACHE_SHA,taxa_selected=len(chosen),retrieval_status=dict(Counter(r['retrieval_status'] for r in results)),new_original_photo_IDs=0,original_resolution_images_fetched=0,derivative_requests=len(results),new_ROI9_successes=0,direct_visual_reviews_completed=0,manifest_sha256=hashlib.sha256(p.read_bytes()).hexdigest())
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(json.dumps(receipt,sort_keys=True),flush=True)
    if not any(r['retrieval_status']=='downloaded' for r in results):raise RuntimeError('No image pixels obtained')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--metadata',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();run(a.source,a.metadata,a.output)
