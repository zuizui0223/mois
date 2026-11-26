"""Bounded remeasurement of 38 EXISTING photo IDs; fixed ROI-v4 on provider large.
Never count repeated image IDs as newly acquired photographs. No main-atlas write.
"""
import argparse, csv, hashlib, io, json, math, re, zipfile
from collections import Counter
from pathlib import Path
from urllib.parse import urlsplit
import pandas as pd

PALETTE=('white','yellow','orange','red','pink','magenta','purple','blue','bronze')
HUES={'white','yellow_orange','red_pink','blue_purple'}
ZIP_SHA={'qa':'eb2a8fa4dd15655b9f1cd17910b45c00408f871ca5db44cf559cf37e23011a94',
 'metadata':'34e940fd2c08b465048532f84168b5c1600b2bdad431213106db944a766ad770',
 'atlas':'d26c524ca50e8d4912285980e5d1d958ce64a3cfa7428e6164d3f89ee0d09fcb'}
TRIAGE_BLOB='f44e054bf6b6ebc60ba0ac69b80abd5aa67a9d34'

def require(ok,msg):
    if not bool(ok):raise ValueError('MOIS_EXISTING38: '+msg)

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def large_url(url,photo_id):
    u=urlsplit(url)
    m=re.fullmatch(r'/photos/([0-9]+)/(?:square|small|medium|large|original)\.(jpg|jpeg|png)',u.path)
    require(u.scheme=='https' and u.netloc=='inaturalist-open-data.s3.amazonaws.com' and not u.query and not u.fragment and m is not None,'unapproved image URL')
    require(m.group(1)==str(photo_id),'photo ID differs from provider path')
    return f'https://{u.netloc}/photos/{photo_id}/large.{m.group(2)}'

def classify(morph,fractions):
    if morph not in HUES:return False,''
    v=[float(x) for x in fractions]
    require(len(v)==9 and all(math.isfinite(x) and x>=0 for x in v) and abs(sum(v)-1)<=.001,'classified result has invalid nine fractions')
    order=sorted(range(9),key=lambda i:v[i],reverse=True)
    return True,(PALETTE[order[0]] if v[order[0]]>=.6 and v[order[0]]-v[order[1]]>=.2 else '')

def signed_zip(folder,key):
    p=Path(folder)/(key+'.zip');require(sha(p)==ZIP_SHA[key],key+' ZIP SHA drift')
    return zipfile.ZipFile(p)

def prepare(inputs,triage,out):
    out=Path(out);require(not out.exists(),'refuse overwrite source')
    raw=Path(triage).read_bytes()
    require(hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()==TRIAGE_BLOB,'visual triage changed')
    t=pd.read_csv(io.BytesIO(raw),dtype=str,keep_default_na=False)
    t=t.loc[t.route.eq('petaloid_roi_recheck')]
    require(len(t)==38 and t.photo_id.is_unique and t.inat_taxon_id.is_unique,'38 fixed photos/taxa missing')
    with signed_zip(inputs,'qa') as z:
        q=pd.read_csv(z.open('mois_96_review_original_photo_IDs.csv'),dtype=str,keep_default_na=False).set_index('photo_id')
    with signed_zip(inputs,'atlas') as z:
        a=pd.read_csv(z.open('current_42111_identity_snapshot_patch.csv.gz'),compression='gzip',dtype={'inat_taxon_id':str}).set_index('inat_taxon_id')
    require(len(a)==42111 and a.index.is_unique and a.loc[t.inat_taxon_id,'after_strict9'].eq(0).all(),'taxa already have numeric ROI9')
    meta={}
    with signed_zip(inputs,'metadata') as z:
        for n in z.namelist():
            if '/numeric_id_' not in n:continue
            for obs in json.loads(z.read(n))['results']:
                for photo in obs.get('photos',[]):
                    pid=str(photo['id']); pair=(obs,photo)
                    if pid in meta:require(meta[pid][0]['uuid']==obs['uuid'],'provider photo reassigned')
                    meta[pid]=pair
    rows=[]
    for r in t.itertuples(index=False):
        old=q.loc[r.photo_id];obs,p=meta[r.photo_id]
        require(old.strict9=='False' and old.inat_taxon_id==r.inat_taxon_id,'not a prior failed photo')
        require(str(obs['taxon']['id'])==r.inat_taxon_id and old.observation_id in (str(obs['id']),obs['uuid']),'source identity mismatch')
        license=p.get('license_code') or ''
        require(license in ('cc0','cc-by','cc-by-sa','cc-by-nc','cc-by-nc-sa','cc-by-nd','cc-by-nc-nd'),'no CC licence')
        mid='MOIS-EXISTING38-'+r.photo_id
        rows.append(dict(photo_id=r.photo_id,inat_taxon_id=r.inat_taxon_id,species=a.loc[r.inat_taxon_id,'species'],observation_id=obs['uuid'],measurement_id=mid,image_filename=mid+'.jpg',photo_license=license,photo_url_large=large_url(p['url'],r.photo_id),shard=0,prior_strict9=False,image_version='provider_large_derivative',prior_unique_photo_id=True))
    d=pd.DataFrame(rows);require(d.observation_id.is_unique,'repeated observation')
    (out/'shards').mkdir(parents=True);dest=out/'shards/shard-000.csv';d.to_csv(dest,index=False,lineterminator='\n')
    m={'schema':'mois_full_source_observation_sampling_cohort_v1','source_arm':'MOIS_EXISTING38_PETALOID_REMEASUREMENT_LARGE_V1','colour_measured':False,'worker_ids':[0],'work_units':1,'photo_shard_sha256':{'0':sha(dest)},'measurement_attempts_planned':38,'new_unique_photo_ids':0,'new_unique_observations':0,'repeat_existing_ids':True,'image_version':'provider_large_derivative_NOT_original','method_change':False,'main_atlas_auto_append':False,'input_ZIP_SHA256':ZIP_SHA,'triage_git_blob':TRIAGE_BLOB}
    (out/'manifest.json').write_text(json.dumps(m,indent=2,sort_keys=True)+'\n')
    print('MOIS_EXISTING38_SOURCE_FROZEN',json.dumps(m,sort_keys=True))

def audit(source,workers,out):
    source=Path(source);workers=Path(workers);out=Path(out)
    require(not out.exists(),'refuse overwrite audit')
    m=json.loads((source/'manifest.json').read_text());p=source/'shards/shard-000.csv'
    require(sha(p)==m['photo_shard_sha256']['0'],'source modified')
    s=pd.read_csv(p,dtype=str,keep_default_na=False)
    r=json.loads((workers/'shard_000/receipt.json').read_text())
    f=workers/'shard_000/photos.csv';require(sha(f)==r['photos_sha256'] and r['source_sha256']==sha(p) and r['attempted']==38,'worker receipt mismatch')
    terminal=workers/'shard_000/sealed/terminal_p000_b0000000.csv'
    require(sha(terminal)==r['terminal_sha256'],'raw numerical terminal mismatch')
    d=pd.read_csv(f,dtype=str,keep_default_na=False)
    require(len(d)==38 and d.photo_id.is_unique and set(d.photo_id)==set(s.photo_id),'worker ID census mismatch')
    for col in ('inat_taxon_id','observation_id','photo_url_large'):
        require(d.set_index('photo_id')[col].sort_index().equals(s.set_index('photo_id')[col].sort_index()),'worker source changed: '+col)
    states=[classify(x.morph,[getattr(x,'flower_fraction_'+c) for c in PALETTE]) for x in d.itertuples(index=False)]
    d['strict9_remeasured']=[x[0] for x in states];d['clear_hue9_remeasured']=[x[1] for x in states]
    for c in PALETTE:d.loc[~d.strict9_remeasured,'flower_fraction_'+c]=''
    raw=pd.read_csv(terminal,dtype=str,keep_default_na=False)
    extras=[c for c in raw.columns if c not in d.columns or c=='measurement_id']
    d=d.merge(raw[extras],on='measurement_id',validate='one_to_one')
    require(int(d.strict9_remeasured.sum())==r['classifiable'],'worker success count differs')
    require(d.loc[d.measurement_status.ne('image_acquisition_failed'),'image_sha256'].str.fullmatch(r'[0-9a-f]{64}').all(),'acquired image lacks SHA')
    out.mkdir(parents=True);dest=out/'existing38_large_ROI9_remeasurement.csv';d.to_csv(dest,index=False)
    summary={'schema':'mois_existing38_large_frozen_roi9_remeasurement_audit_v1','repeat_existing_photo_ids':38,'distinct_taxa':38,'new_unique_photos':0,'new_unique_observations':0,'large_images_acquired':r['downloaded'],'image_version':'provider_large_derivative_NOT_original','numerical_ROI9_successes':int(d.strict9_remeasured.sum()),'clear_sample_photo_hues':int(d.clear_hue9_remeasured.ne('').sum()),'measurement_status_counts':dict(Counter(d.measurement_status)),'morph_counts':dict(Counter(d.morph)),'failure_reasons_counts':dict(Counter(d['failure_reasons'])) if 'failure_reasons' in d else {},'main_catalogue_modified':False,'main_unique_photo_attempts_unchanged':1231337,'main_numeric_ROI9_unchanged_pending_review':467565,'true_petal_masks_independently_verified':False,'source_csv_SHA256':sha(p),'worker_receipt_SHA256':sha(workers/'shard_000/receipt.json'),'terminal_SHA256':sha(terminal),'output_SHA256':{dest.name:sha(dest)},'method':r['method']}
    (out/'receipt.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    print('MOIS_EXISTING38_REMEASUREMENT_AUDITED',json.dumps(summary,sort_keys=True))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','audit']);p.add_argument('--inputs');p.add_argument('--triage');p.add_argument('--source');p.add_argument('--workers');p.add_argument('--output',required=True);a=p.parse_args()
    if a.stage=='prepare':prepare(a.inputs,a.triage,a.output)
    else:audit(a.source,a.workers,a.output)
