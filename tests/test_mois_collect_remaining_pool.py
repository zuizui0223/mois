import json
import hashlib
from pathlib import Path
import pandas as pd
import pytest
from scripts.acquisition import mois_collect_remaining_pool as m


def rows(n, start=100):
    return pd.DataFrame([{
        'photo_id': str(start+i), 'measurement_id': f'MOIS-{start+i:024d}',
        'image_filename': f'MOIS-{start+i:024d}.jpg', 'inat_taxon_id': str(i%3+1),
        'species': f'Campanula species{i%3}', 'observation_id': f'uuid-{start+i}',
        'observer_id': str(start+i), 'observed_on': '2025-06-12',
        'latitude': '36.1', 'longitude': '139.5', 'photo_license': 'cc-by',
        'photo_url_large': f'https://inaturalist-open-data.s3.amazonaws.com/photos/{start+i}/large.jpg',
        'source_snapshot_date': '2026-10-18', 'source_arm': 'cached_original',
        'shard': '0'} for i in range(n)])


def fixture_source(root):
    root.mkdir()
    photo = rows(256)
    path=root/'candidates.csv'; photo.to_csv(path,index=False,lineterminator='\n')
    digest=hashlib.sha256(path.read_bytes()).hexdigest()
    (root/'source.json').write_text(json.dumps({'schema':'mois_cached_remaining_source_v1',
        'status':'SOURCE_ONLY', 'count':256,'shards':1,'data_sha256':digest,
        'no_ecological_analysis':True, 'parents': {'pool':'a'*64, 'second_measurement':'b'*64}}))
    return photo


def test_selects_unprocessed_photos_not_one_photo_per_species():
    pool=rows(1024)
    selected=m.select_remaining(pool,pool.iloc[:256],pool.iloc[256:512],512)
    assert len(selected)==512 and selected.photo_id.is_unique
    assert selected.photo_id.tolist()==pool.iloc[512:].photo_id.tolist()
    assert selected.inat_taxon_id.nunique()==3
    assert selected.shard.value_counts().to_dict()=={'0':256,'1':256}


def test_insufficient_or_duplicate_pool_rejected():
    pool=rows(1024)
    with pytest.raises(ValueError):m.select_remaining(pool,pool.iloc[:256],pool.iloc[256:512],768)
    with pytest.raises(ValueError):m.select_remaining(pd.concat([pool,pool.iloc[:1]]),pool.iloc[:256],pool.iloc[256:512],512)


def test_source_identity_and_unsafe_photo_url(tmp_path):
    source=tmp_path/'source'; photo=fixture_source(source)
    got, receipt=m.source_paths(source,0)
    assert got.equals(photo)
    path=source/'candidates.csv'
    path.write_text(path.read_text().replace('/photos/100/','/photos/999/'))
    with pytest.raises(ValueError):m.source_paths(source,0)


def test_real_result_checkpoint_reused_without_second_inference(tmp_path):
    root=tmp_path/'source'; source=fixture_source(root); calls=[]
    def fake_runner(_,index,batch,size,sealed,**kwargs):
        calls.append(index)
        worker, urls, units=kwargs['indexed_source'].get_batch(index,batch,size)
        assert size==256 and units==1
        sealed.mkdir(parents=True,exist_ok=True)
        table=pd.DataFrame({'measurement_id':source.measurement_id,
            'morph':['white']*70+['red_pink']*30+['mixed_uncertain']*156,
            'measurement_status':['classified_four_state_morph']*100+
                ['not_evaluable_roi_or_flip_gate']*150+['image_acquisition_failed']*6,
            'image_sha256':['e'*64]*250+['']*6})
        path=sealed/f'terminal_p{index:03d}_b0000000.csv';table.to_csv(path,index=False)
        return {'status':'COMPLETE_WORK_UNIT','n_source_photo_ids':256,'n_terminal_photo_ids':256,
            'n_classifiable':100,'n_acquisition_failed':6,'image_bytes_persisted':False,
            'reused_historical_measurement_csv':False,'number_batches_for_this_partition':1,
            'terminal_csv_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'method':{'test':True}}
    output=tmp_path/'measurements'
    r=m.measure(root,0,tmp_path,tmp_path,output,runner=fake_runner)
    assert (r['attempted'],r['downloaded'],r['classified'],r['acquisition_failed'])==(256,250,100,6)
    assert r['unclassified']==156 and r['acquired_unclassified']==150
    assert r['ecological_models_run'] is False
    assert m.measure(root,0,tmp_path,tmp_path,output,runner=fake_runner)==r
    assert calls==[0]
    p=output/'shard_00'/'photos.csv';p.write_text(p.read_text()+'\n')
    with pytest.raises(ValueError):m.measure(root,0,tmp_path,tmp_path,output,runner=fake_runner)
    assert calls==[0]
