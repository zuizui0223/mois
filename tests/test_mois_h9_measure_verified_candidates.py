import importlib.util
import json
from pathlib import Path
import pandas as pd
import pytest

MODULE = Path(__file__).parents[1] / 'scripts/analysis/mois_h9_measure_verified_candidates.py'

def load_module():
    assert MODULE.exists(), 'The H9 single-process measurement executor is not implemented'
    spec = importlib.util.spec_from_file_location('h9_measure', MODULE)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m

@pytest.fixture
def api():
    return load_module()

def source(tmp, m):
    rows = []
    for i in range(3):
        mid = f'MOIS-{i:024d}'
        rows.append(dict(measurement_id=mid, image_filename=mid+'.jpg',
            photo_id=str(500+i), inat_taxon_id=str(100+i), species=f'Novelus s{i}',
            observation_id=str(900+i), observer_id=str(200+i), observed_on='2025-05-18',
            latitude='36.1', longitude='137.2', photo_license='cc-by',
            photo_url_large=f'https://inaturalist-open-data.s3.amazonaws.com/photos/{500+i}/large.jpg',
            source_sampling_arm='fcp_h9_flowering_annotated_rescue_after_unknown',
            source_annotation='flowers_and_fruits:flowers_(12:13)'))
    frame = pd.DataFrame(rows)
    tmp.mkdir(exist_ok=True)
    frame.to_csv(tmp / m.SOURCE, index=False)
    receipt = dict(schema=m.SOURCE_SCHEMA, new_photo_pixels_actually_measured=False,
        rescued_flowering_annotated_new_photo_SOURCE_ONLY_IDs=len(frame),
        source_new_photo_candidates_SHA256=m.digest(tmp/m.SOURCE))
    (tmp / m.SOURCE_RECEIPT).write_text(json.dumps(receipt))
    return frame


def fake_runner(m, calls, failure=False):
    def runner(root, part, batch, size, out, **kwargs):
        worker, key, n = kwargs['indexed_source'].get_batch(part,batch,size)
        calls.append(size)
        assert list(worker) == ['measurement_id','image_filename','photo_license']
        assert set(key) == {*worker.columns, 'photo_url_large'} and n == 1
        if failure:
            raise TimeoutError('Provider temporarily unavailable')
        out.mkdir(parents=True, exist_ok=True)
        frame = pd.DataFrame({'measurement_id':worker.measurement_id,
            'morph':['white','red_pink','mixed_uncertain'],
            'measurement_status':['classified_four_state_morph']*2+['image_acquisition_failed']})
        for hue in m.PALETTE:
            frame['flower_fraction_'+hue] = [float(hue=='white'), float(hue=='red'), 0.0]
        terminal = out/'terminal_p000_b0000000.csv'
        frame.to_csv(terminal, index=False)
        r = dict(schema='fcp_fresh_blind_terminal_partition_v1', status='COMPLETE_WORK_UNIT',
            n_source_photo_ids=size,n_terminal_photo_ids=size,n_classifiable=2,n_acquisition_failed=1,
            compute_batch_size_not_photo_cap=size, number_batches_for_this_partition=1,
            terminal_csv_sha256=m.digest(terminal), method={'synthetic_test_only':True},
            selected_source_worker_sha256=m.digest_bytes(worker.to_csv(index=False).encode()),
            image_bytes_persisted=False,reused_historical_measurement_csv=False)
        (out/'terminal_p000_b0000000.json').write_text(json.dumps(r))
        return r
    return runner


def test_one_process_blind_source_and_failure_accounting(tmp_path):
    api = load_module()
    frame = source(tmp_path/'source',api); calls=[]
    r=api.run(tmp_path/'source',tmp_path/'runtime',tmp_path/'sam',tmp_path/'result',
        executor=fake_runner(api,calls))
    assert calls == [3]
    assert r['attempted_photos']==3 and r['classified_photos']==2
    assert r['download_failures']==1 and r['hue_unknown_after_download']==0
    saved=pd.read_csv(tmp_path/'result'/api.RESULT,dtype=str)
    assert saved.photo_id.tolist()==frame.photo_id.tolist()
    assert saved.classified.tolist()==['True','True','False']
    assert r['white_count']==1 and r['model_changed'] is False
    assert r['annotation_stratum_not_unbiased_prevalence'] is True


def test_resume_valid_seal_without_redownloading(tmp_path):
    api = load_module()
    source(tmp_path/'source',api); calls=[]
    original=api.run(tmp_path/'source',tmp_path/'runtime',tmp_path/'sam',tmp_path/'result',
        executor=fake_runner(api,calls))
    (tmp_path/'result'/api.RECEIPT).unlink()
    (tmp_path/'result'/api.RESULT).unlink()
    def must_not_run(*a,**k):
        raise AssertionError('A valid sealed pixel unit was needlessly recomputed')
    resumed=api.run(tmp_path/'source',tmp_path/'runtime',tmp_path/'sam',tmp_path/'result',
        executor=must_not_run)
    assert resumed['classified_photos']==original['classified_photos']
    assert resumed['reused_sealed_work_unit'] is True


def test_source_or_pixel_tampering_aborts(tmp_path):
    api = load_module()
    source(tmp_path/'source',api)
    api.run(tmp_path/'source',tmp_path/'runtime',tmp_path/'sam',tmp_path/'result',executor=fake_runner(api,[]))
    terminal=tmp_path/'result/work/terminal_p000_b0000000.csv'
    terminal.write_text(terminal.read_text().replace('red_pink','white'))
    with pytest.raises(ValueError,match='digest'):
        api.run(tmp_path/'source',tmp_path/'runtime',tmp_path/'sam',tmp_path/'result',executor=fake_runner(api,[]))


def test_transient_failure_does_not_create_colour_summary(tmp_path):
    api = load_module()
    source(tmp_path/'source',api)
    with pytest.raises(TimeoutError):
        api.run(tmp_path/'source',tmp_path/'runtime',tmp_path/'sam',tmp_path/'result',executor=fake_runner(api,[],failure=True))
    assert not (tmp_path/'result'/api.RECEIPT).exists()


def test_source_validation_and_safe_file_names(tmp_path):
    api = load_module()
    frame=source(tmp_path/'source',api)
    frame.loc[0,'image_filename']='../../bad.jpg'
    frame.to_csv(tmp_path/'source'/api.SOURCE,index=False)
    p=tmp_path/'source'/api.SOURCE_RECEIPT;r=json.loads(p.read_text())
    r['source_new_photo_candidates_SHA256']=api.digest(tmp_path/'source'/api.SOURCE)
    p.write_text(json.dumps(r))
    with pytest.raises(ValueError,match='filename'):
        api.load_source(tmp_path/'source')


def test_unclassified_not_treated_as_chromatic(tmp_path):
    api = load_module()
    frame=source(tmp_path/'source',api)
    measured=pd.DataFrame({'measurement_id':frame.measurement_id,
        'morph':['white','mixed_uncertain','mixed_uncertain'],
        'measurement_status':['classified_four_state_morph','not_evaluable_roi_or_flip_gate','image_acquisition_failed']})
    joined=api.join_labels(frame,measured)
    assert joined.classified.sum()==1 and joined.white.sum()==1
    assert int((joined.classified & ~joined.white).sum())==0
