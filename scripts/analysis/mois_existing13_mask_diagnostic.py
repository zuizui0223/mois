"""One diagnostic replay of the 13 recovered numerical successes, not new data.
Fail if source image bytes change. Preserve actual masks and credited small panels.
"""
import argparse,hashlib,io,json,zipfile
from pathlib import Path
import numpy as np
import pandas as pd
from PIL import Image,ImageOps,ImageDraw
from scripts.acquisition.run_fcp_fresh_image_partition import image_loader,verify_frozen_method
from fcp_pipeline.flower_roi_v4_runtime import FrozenFlowerColourEstimator
from fcp_pipeline.photo_first_measurement import nearest_palette_counts,flower_only_fractions

AUDIT_SHA='7f4ce68441cf13b7c4f6ee0cced47d4dd7da4e7a87ade2909e3897b1e9dca7e7'
META_SHA='34e940fd2c08b465048532f84168b5c1600b2bdad431213106db944a766ad770'

def require(ok,msg):
    if not bool(ok):raise ValueError('MOIS_MASK_DIAG: '+msg)

def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def run(audit,metadata,runtime,sam,output):
    require(digest(audit)==AUDIT_SHA and digest(metadata)==META_SHA,'source ZIP changed')
    with zipfile.ZipFile(audit) as z:
        raw=z.read('existing38_large_ROI9_remeasurement.csv');r=json.loads(z.read('receipt.json'))
    require(hashlib.sha256(raw).hexdigest()==r['output_SHA256']['existing38_large_ROI9_remeasurement.csv'],'numerical source changed')
    d=pd.read_csv(io.BytesIO(raw),dtype={'photo_id':str,'inat_taxon_id':str})
    d=d.loc[d.strict9_remeasured].copy();require(len(d)==13 and d.photo_id.is_unique,'fixed 13 success IDs changed')
    attribution={}
    with zipfile.ZipFile(metadata) as z:
        for name in z.namelist():
            if '/numeric_id_' not in name:continue
            for obs in json.loads(z.read(name))['results']:
                for p in obs.get('photos',[]):attribution[str(p['id'])]=p.get('attribution','')
    root=Path(runtime);out=Path(output);require(not out.exists(),'refuse overwrite')
    detector=root/'data/atlas/qualification/roi_v4_training/jrc_yolo11n_last_v4.pt'
    method=verify_frozen_method(root/'scripts/analysis/measure_random_photo_first_blind_partition.py',detector,Path(sam))
    require(method==r['method'],'frozen method differs')
    contract=json.loads((root/'docs/supporting/jbi_atlas_roi_estimator_contract_v4.json').read_text())
    execution=json.loads((root/'docs/supporting/random_photo_first_measurement_execution_v1.json').read_text())
    estimator=FrozenFlowerColourEstimator(detector,Path(sam),contract,torch_threads=int(execution['runtime']['torch_threads_per_worker']))
    out.mkdir();(out/'masks').mkdir();(out/'review').mkdir();rows=[]
    for row in d.itertuples(index=False):
        payload=image_loader(row.photo_url_large)
        require(hashlib.sha256(payload).hexdigest()==row.image_sha256,'large bytes changed; refuse replay '+row.photo_id)
        with Image.open(io.BytesIO(payload)) as im:
            measured=estimator.measure(im)
            oriented=ImageOps.exif_transpose(im).convert('RGB');rgb=np.asarray(oriented)
        require(measured.get('automated_colour_state_status')=='automated_colour_state_admitted','runtime admission drift '+row.photo_id)
        mask=np.asarray(measured['flower_mask'],dtype=bool)
        require(mask.shape==rgb.shape[:2] and mask.any(),'invalid actual flower mask')
        counts=nearest_palette_counts(rgb[mask]);fractions=flower_only_fractions(counts)
        require(all(abs(fractions[k]-float(getattr(row,'flower_fraction_'+k)))<=1e-10 for k in fractions),'numeric fractions drift '+row.photo_id)
        y,x=np.where(mask);bounds=[int(x.min()),int(y.min()),int(x.max()+1),int(y.max()+1)]
        mp=out/'masks'/f'{row.photo_id}.png';Image.fromarray(mask.astype('uint8')*255).save(mp)
        dark=rgb.copy();dark[~mask]=(dark[~mask].astype(float)*.12).astype('uint8')
        frames=[oriented,Image.fromarray(dark),Image.fromarray(dark).crop(bounds)]
        panel=Image.new('RGB',(1500,565),'white');draw=ImageDraw.Draw(panel)
        for i,img in enumerate(frames):panel.paste(ImageOps.contain(img,(490,490)),(i*500,55))
        draw.text((5,3),f'{row.species} | photo {row.photo_id} | full / retained-mask / mask bounds',fill='black')
        draw.text((5,20),attribution[row.photo_id],fill='black')
        draw.text((5,37),f'https://www.inaturalist.org/photos/{row.photo_id} | {row.photo_license} | large derivative',fill='black')
        pp=out/'review'/f'{row.photo_id}.jpg';panel.save(pp,quality=92)
        rows.append({'photo_id':row.photo_id,'inat_taxon_id':row.inat_taxon_id,'species':row.species,'image_sha256':row.image_sha256,'width':rgb.shape[1],'height':rgb.shape[0],'mask_pixels':int(mask.sum()),'mask_bbox':json.dumps(bounds),'mask_SHA256':digest(mp),'panel_SHA256':digest(pp),'photo_license':row.photo_license,'attribution':attribution[row.photo_id],'prior_numeric9_exactly_reproduced':True,'source_image_bytes_exact_match':True,'botanical_mask_review':'pending'})
        print('MASK_DIAGNOSTIC_COMPLETED',row.photo_id,flush=True)
    df=pd.DataFrame(rows);df.to_csv(out/'mask_manifest.csv',index=False)
    summary={'schema':'mois_existing13_same_bytes_frozen_mask_diagnostic_v1','diagnostic_replayed_photos':13,'source_audit_zip_SHA256':AUDIT_SHA,'image_bytes_match_prior_run':True,'all13_numeric_vectors_reproduced':True,'new_unique_photos':0,'new_numeric9_successes_added':0,'main_catalogue_modified':False,'original_image_bytes_retained':False,'credited_review_derivatives_and_masks_retained':True,'method':method,'output_SHA256':{str(p.relative_to(out)):digest(p) for p in out.rglob('*') if p.is_file()}}
    (out/'receipt.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    print(json.dumps(summary,sort_keys=True))

if __name__=='__main__':
    p=argparse.ArgumentParser()
    for k in ('audit','metadata','runtime','sam','output'):p.add_argument('--'+k,required=True)
    a=p.parse_args();run(a.audit,a.metadata,a.runtime,a.sam,a.output)
