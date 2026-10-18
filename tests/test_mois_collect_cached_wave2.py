import hashlib
import json
from pathlib import Path
import pandas as pd
import pytest
from scripts.acquisition import mois_collect_cached_wave2 as m

PLAN='mois_global_angiosperm_second4096_source_only_original_photoids.csv'
RECEIPT='mois_global_angiosperm_cached_second4096_source_only_receipt.json'

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

@pytest.fixture
def source(tmp_path):
    root=tmp_path/'source';root.mkdir()
    rows=pd.DataFrame([{'photo_id':str(10000+i),'measurement_id':f'MOIS-{i:024d}',
        'image_filename':f'MOIS-{i:024d}.jpg','inat_taxon_id':str(20000+i),
        'species':f'Flora species{i}','observation_id':f'observation-{i}',
        'observer_id':str(i+1),'observed_on':'2025-06-10','source_snapshot_date':'2026-10-18',
        'latitude':'35.0','longitude':'135.0','photo_license':'cc-by',
        'photo_url_large':f'https://inaturalist-open-data.s3.amazonaws.com/photos/{10000+i}/large.jpg',
        'source_arm':'mois_broad_angiosperm_cached_wave2_new_species','shard':str(i%16)} for i in range(4096)])
    rows.to_csv(root/PLAN,index=False)
    hashes={}
    for i in range(16):
        f=root/f'mois_broad_secondwave_source_shard_{i:02d}_256.csv'
        rows[rows.shard.eq(str(i))].to_csv(f,index=False);hashes[str(i)]=digest(f)
    doc={'schema':'mois_angiosperm_cached_second4096_new_species_source_only_v1',
        'status':'SOURCE_ONLY_4096_RESERVED_NO_NEW_PIXEL_HUES','new_photo_pixels_actually_measured':False,
        'second4096_photo_candidates':4096,'parent_original_metadata_snapshot_date':'2026-10-18',
        'first4096_and_new4096_disjoint_photo_IDs':True,
        'original313535_measured_photo_IDs_excluded_via_parent_signed_source':True,
        'other_fifth4096_and_SI4096_reservations_excluded_via_parent':True,
        'source_only_second4096_csv_sha256':digest(root/PLAN),'source_only_shard_SHA256':hashes}
    (root/RECEIPT).write_text(json.dumps(doc));return root

def test_signed_original_4096_and_256_source_rows(source):
    rows,doc=m.verify_source(source,3)
    assert len(rows)==256 and rows.photo_id.is_unique
    assert rows.shard.eq('3').all() and 'morph' not in rows
    assert doc['new_photo_pixels_actually_measured'] is False

def test_source_tamper_rejected(source):
    f=source/'mois_broad_secondwave_source_shard_00_256.csv'
    f.write_text(f.read_text().replace('cc-by','cc0',1))
    with pytest.raises(ValueError):m.verify_source(source,0)

def fake_runner(calls):
    def run(_root,part,batch,size,out,**kwargs):
        calls.append((part,batch,size))
        worker,acq,n=kwargs['indexed_source'].get_batch(part,batch,size)
        assert n==1 and len(acq)==256
        out.mkdir(parents=True,exist_ok=True)
        labels=pd.DataFrame({'measurement_id':worker.measurement_id,
            'morph':['white']*80+['red_pink']*20+['mixed_uncertain']*156,
            'measurement_status':['classified_four_state_morph']*100+['not_evaluable_roi_or_flip_gate']*153+['image_acquisition_failed']*3,
            'image_sha256':['a'*64]*253+['']*3})
        f=out/f'terminal_p{part:03d}_b{batch:07d}.csv';labels.to_csv(f,index=False)
        return {'status':'COMPLETE_WORK_UNIT','n_source_photo_ids':256,'n_terminal_photo_ids':256,
            'image_bytes_persisted':False,'reused_historical_measurement_csv':False,
            'n_classifiable':100,'n_acquisition_failed':3,'number_batches_for_this_partition':1,
            'terminal_csv_sha256':digest(f),'method':{'source_commit':'fixed-fixture-only'}}
    return run

def test_256_photos_use_one_evaluator_and_resume_without_repeating(source,tmp_path):
    calls=[];out=tmp_path/'results'
    r=m.measure(source,0,tmp_path/'fcp',tmp_path/'sam',out,runner=fake_runner(calls))
    assert calls==[(0,0,256)]
    assert (r['attempted'],r['classified'],r['acquisition_failed'],r['downloaded'])==(256,100,3,253)
    assert r['unclassified']==156 and r['ecological_models_run'] is False
    again=m.measure(source,0,tmp_path/'fcp',tmp_path/'sam',out,runner=fake_runner(calls))
    assert len(calls)==1 and r==again
    rows=pd.read_csv(out/'shard_00'/'photos.csv',dtype=str,keep_default_na=False)
    assert len(rows)==256 and rows.photo_id.is_unique
    assert (rows.morph=='mixed_uncertain').sum()==156
    assert 'photo_license' in rows and 'observed_on' in rows

def test_incomplete_collection_never_claims4096(tmp_path):
    with pytest.raises(ValueError):m.collect(tmp_path,tmp_path/'combined')

def test_corrupt_completed_cache_not_silently_reused(source,tmp_path):
    calls=[];out=tmp_path/'results'
    m.measure(source,1,tmp_path/'fcp',tmp_path/'sam',out,runner=fake_runner(calls))
    p=out/'shard_01'/'photos.csv';p.write_text(p.read_text().replace('white','blue_purple',1))
    with pytest.raises(ValueError):m.measure(source,1,tmp_path/'fcp',tmp_path/'sam',out,runner=fake_runner(calls))
    assert len(calls)==1

def test_full16_completed_shards_save_only_data_and_counts(source,tmp_path):
    pytest.importorskip('pyarrow')  # Mandatory in CI, optional in a minimal local runtime.
    out=tmp_path/'results';calls=[]
    for i in range(16):
        m.measure(source,i,tmp_path/'fcp',tmp_path/'sam',out,runner=fake_runner(calls))
    r=m.collect(out,tmp_path/'published')
    assert len(calls)==16 and r['attempted']==4096 and r['classified']==1600
    assert r['acquisition_failed']==48 and r['downloaded']==4048
    assert not r['ecological_models_run']
    rows=pd.read_parquet(tmp_path/'published'/'wave2_measured_photos.parquet')
    assert len(rows)==4096 and rows.photo_id.is_unique
    assert digest(tmp_path/'published'/'wave2_measured_photos.parquet')==r['photo_data_sha256']
