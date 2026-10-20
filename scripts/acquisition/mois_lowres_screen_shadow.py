#!/usr/bin/env python3
"""Compare an in-memory 320px detector with the unchanged full ROI-v4 pipeline.

SHADOW ONLY: every acquired image receives the full frozen baseline, including
screen negatives and screen errors. Image download is shared, not repeated.
A successful screen comparison is NOT independent botanical validation.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import math
import time
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pandas as pd

HUES = ('white', 'yellow_orange', 'red_pink', 'blue_purple')
THRESHOLDS = (0.05, 0.10, 0.25)
IMAGE_SIZE = 320
SCREEN_COLUMNS = ('measurement_id', 'image_sha256', 'morph', 'screen_status',
                  'screen_max_confidence', 'screen_box_count', 'screen_seconds',
                  'baseline_roi_seconds', 'baseline_runtime_error', 'screen_error',
                  'production_filter_enabled')
SOURCE_BLOB = 'b1d7091db44682559bbe4b48a9f3b8db648501a6'


def require(ok, message):
    if not ok:
        raise ValueError(message)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


class ShadowEstimator:
    """Return the baseline result unchanged; screen cannot exclude a photograph."""
    def __init__(self, baseline, screen, records):
        self.baseline, self.screen, self.records = baseline, screen, records

    def measure(self, image):
        row = {'baseline_runtime_error': False, 'production_filter_enabled': False}
        start = time.perf_counter()
        try:
            result = self.baseline.measure(image)
        except Exception:
            row['baseline_runtime_error'] = True
            raise
        finally:
            row['baseline_roi_seconds'] = time.perf_counter() - start
            # Baseline runs first. The extra detector does not supply masks,
            # labels, predictor state, or decisions to the original estimator.
            start = time.perf_counter()
            try:
                confidence, count = self.screen(image)
                confidence, count = float(confidence), int(count)
                require(math.isfinite(confidence) and 0 <= confidence <= 1 and count >= 0,
                        'Invalid screen confidence/count')
                row.update(screen_status='ok', screen_max_confidence=confidence,
                           screen_box_count=count, screen_error='')
            except Exception as exc:
                row.update(screen_status='error', screen_max_confidence=None,
                           screen_box_count=None, screen_error=f'{type(exc).__name__}:{str(exc)[:250]}')
            row['screen_seconds'] = time.perf_counter() - start
            self.records.append(row)
        return result


def make_screen(weight, contract):
    from PIL import Image, ImageOps
    from ultralytics import YOLO
    detector = YOLO(str(weight))  # separate predictor; never reuse baseline state
    inference = contract['detector']['inference']

    def screen(image):
        preview = ImageOps.exif_transpose(image).convert('RGB')
        preview.thumbnail((IMAGE_SIZE, IMAGE_SIZE), Image.Resampling.BILINEAR)
        prediction = detector.predict(source=preview, imgsz=IMAGE_SIZE,
            conf=min(THRESHOLDS), iou=float(inference['nms_iou']),
            max_det=int(inference['maximum_detections']), augment=False,
            device='cpu', verbose=False, save=False, save_txt=False, save_crop=False)[0]
        if prediction.boxes is None or len(prediction.boxes) == 0:
            return 0.0, 0
        scores = prediction.boxes.conf.detach().cpu().numpy()
        return float(scores.max()), len(scores)
    return screen


def summarize(frame):
    required = {'measurement_id','morph','screen_status','screen_max_confidence',
                'baseline_roi_seconds','screen_seconds'}
    require(required <= set(frame), 'Missing shadow comparison columns')
    ids = frame.measurement_id.astype(str)
    require(frame.measurement_id.notna().all() and ids.str.len().gt(0).all()
            and ids.is_unique, 'Invalid shadow measurement identity')
    require(frame.morph.isin((*HUES,'mixed_uncertain')).all(), 'Unknown baseline hue')
    require(frame.screen_status.isin(('ok','error')).all(), 'Invalid screen status')
    for key in ('baseline_roi_seconds','screen_seconds'):
        values = pd.to_numeric(frame[key], errors='coerce')
        require(np.isfinite(values).all() and values.ge(0).all(), 'Invalid runtime timing')
    scores = pd.to_numeric(frame.screen_max_confidence, errors='coerce')
    ok = frame.screen_status.eq('ok')
    require(scores.loc[ok].between(0,1).all(), 'Invalid valid-screen confidence')
    classified = frame.morph.isin(HUES)
    report = {'schema':'mois_320px_screen_shadow_vs_frozen_ROIv4_v1',
        'photos_with_baseline_and_shadow':len(frame), 'screen_errors':int((~ok).sum()),
        'baseline_classifiable':int(classified.sum()),
        'baseline_unknown':int((~classified).sum()),
        'production_filter_enabled':False,
        'independent_botanical_flowering_validation_completed':False,
        'reference_is_automated_ROI_NOT_verified_flowering_or_petal_truth':True,
        'baseline_ROI_seconds':float(frame.baseline_roi_seconds.sum()),
        'extra_screen_seconds':float(frame.screen_seconds.sum()),
        'timing_scope':'ROI calls only, includes cold first inference; excludes downloads, initialization and palette classification',
        'thresholds':{}}
    for threshold in THRESHOLDS:
        negative = ok & scores.lt(threshold)  # screen errors must be kept
        n_classified = int(classified.sum())
        lost = int((negative & classified).sum())
        by_hue = {}
        for hue in HUES:
            same = frame.morph.eq(hue)
            n, missed = int(same.sum()), int((same & negative).sum())
            by_hue[hue] = {'baseline_n':n,'lost':missed,
                          'retention':(n-missed)/n if n else None}
        report['thresholds'][format(threshold,'.2f')] = {
            'hypothetical_skipped_images':int(negative.sum()),
            'baseline_classifiable_lost':lost,
            'baseline_unknown_in_skipped':int((negative & ~classified).sum()),
            'baseline_classifiable_retention':(n_classified-lost)/n_classified if n_classified else None,
            'projected_net_roi_seconds_saved_NOT_realized':float(
                frame.loc[negative,'baseline_roi_seconds'].sum()-frame.screen_seconds.sum()),
            'by_hue':by_hue}
    return report


def evaluate_shadow(worker, images, source, detector, sam, out, *, sidecar):
    """Call the actual frozen source's main, adding only an observation wrapper."""
    source, sidecar = Path(source), Path(sidecar)
    raw = source.read_bytes()
    require(hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest() == SOURCE_BLOB,
            'Frozen baseline source Git blob changed')
    spec = importlib.util.spec_from_file_location('_mois_original_frozen_measurement', source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    original_class = module.FrozenFlowerColourEstimator
    records = []

    def factory(weight, segmenter_dir, contract, *, torch_threads):
        baseline = original_class(weight, segmenter_dir, contract, torch_threads=torch_threads)
        try:
            screen = make_screen(weight, contract)
        except Exception as exc:
            message = f'{type(exc).__name__}:{exc}'
            def screen(image, error=message):
                raise RuntimeError('Screen initialization failed: '+error)
        return ShadowEstimator(baseline, screen, records)

    module.FrozenFlowerColourEstimator = factory
    module.parse_args = lambda: SimpleNamespace(worker_manifest=Path(worker), images_dir=Path(images),
        detector_weight=Path(detector), efficient_sam_dir=Path(sam), output_csv=Path(out),
        measurement_contract=module.DEFAULT_MEASUREMENT, execution_contract=module.DEFAULT_EXECUTION,
        roi_contract=module.DEFAULT_ROI, locked_result=module.DEFAULT_LOCKED)
    require(module.main() == 0, 'Original frozen baseline did not finish')
    baseline = pd.read_csv(out, dtype={'measurement_id':str})
    blind = pd.read_csv(worker, dtype=str, keep_default_na=False)
    require(len(records)==len(blind)==len(baseline) and blind.measurement_id.is_unique,
            'Shadow screen lacks an acquired original image')
    comparison = pd.DataFrame(records)
    comparison['measurement_id'] = blind.measurement_id.to_numpy()
    comparison = comparison.merge(baseline[['measurement_id','image_sha256','morph']],
                                  on='measurement_id', validate='one_to_one')
    comparison = comparison.loc[:,SCREEN_COLUMNS]
    write_sidecar(comparison, sidecar)
    return baseline  # unchanged original full-resolution ROI/palette outcomes


def write_sidecar(comparison, sidecar):
    sidecar = Path(sidecar)
    sidecar.mkdir(parents=True,exist_ok=True)
    file = sidecar/'screen_shadow.csv'
    require(not file.exists(), 'Refuse to overwrite shadow measurements')
    comparison.to_csv(file,index=False,lineterminator='\n')
    (sidecar/'screen_shadow_receipt.json').write_text(json.dumps({
        'schema':'mois_original_image_shared_download_shadow_screen_v1',
        'rows':len(comparison), 'csv_sha256':digest(file),
        'full_baseline_executed_for_all_acquired_images':True,
        'production_filter_enabled':False, 'screen_image_size':IMAGE_SIZE,
        'screen_thresholds_fixed_before_pixels':list(THRESHOLDS),
        'baseline_measurement_git_blob':SOURCE_BLOB},indent=2,sort_keys=True)+'\n')


def work(source, runtime, sam, group, output):
    from scripts.acquisition.mois_full_index_real_original_photo_measure import measure_group, plan, input_rows
    from scripts.acquisition.run_fcp_fresh_image_partition import run_batch, image_loader, PermanentOriginalPhotoError
    from scripts.acquisition.mois_prefetch import PrefetchLoader
    info = plan(source)
    units = range(group*4,min((group+1)*4,info['work_units']))
    urls=[]
    for unit in units:
        folder=Path(output)/f'shard_{unit:03d}'
        if (folder/'receipt.json').exists():
            r=json.loads((folder/'sealed/screen_shadow_receipt.json').read_text())
            require(digest(folder/'sealed/screen_shadow.csv')==r['csv_sha256'],
                    'Previously completed shadow data missing or altered')
        else:
            rows,_=input_rows(source,unit)
            urls.extend(rows.photo_url_large.tolist())

    def runner(*args, **kwargs):
        sidecar=Path(args[4])
        def evaluate(worker,images,source,detector,sam,out):
            return evaluate_shadow(worker,images,source,detector,sam,out,sidecar=sidecar)
        result=run_batch(*args,evaluator=evaluate,**kwargs)
        if not (sidecar/'screen_shadow.csv').exists():
            require(result['n_acquisition_failed']==result['n_source_photo_ids'],
                    'Missing shadow comparison for acquired images')
            write_sidecar(pd.DataFrame(columns=SCREEN_COLUMNS),sidecar)
        return result

    with PrefetchLoader(urls,image_loader,Path(output)/'prefetch',
                       permanent_errors=(PermanentOriginalPhotoError,)) as prefetch:
        return measure_group(source,group,runtime,sam,output,runner=runner,loader=prefetch)


def collect(source, workers, output):
    from scripts.acquisition.mois_full_index_real_original_photo_measure import collect as baseline_collect
    receipt=baseline_collect(source,workers,output)
    output=Path(output)
    expected=pd.read_parquet(output/'measured_photos.parquet')
    acquired=expected.loc[expected.measurement_status.ne('image_acquisition_failed')]
    frames=[]
    for file in sorted(Path(workers).rglob('screen_shadow.csv')):
        r=json.loads(file.with_name('screen_shadow_receipt.json').read_text())
        require(digest(file)==r['csv_sha256'], 'Shadow source checksum changed')
        require(r['full_baseline_executed_for_all_acquired_images'] is True
                and r['production_filter_enabled'] is False,
                'Screen must never filter before independent validation')
        frame=pd.read_csv(file,dtype={'measurement_id':str})
        require(len(frame)==r['rows'], 'Shadow row count changed')
        frames.append(frame)
    require(bool(frames),'No shadow worker receipts')
    shadow=pd.concat(frames,ignore_index=True)
    require(shadow.measurement_id.is_unique and set(shadow.measurement_id)==set(acquired.measurement_id),
            'Shadow coverage missing/duplicated original acquired measurement identity')
    joined=acquired[['measurement_id','morph','image_sha256']].merge(shadow,
        on='measurement_id',validate='one_to_one',suffixes=('_original',''))
    require(joined.morph.eq(joined.morph_original).all()
            and joined.image_sha256.eq(joined.image_sha256_original).all(),
            'Original source image or frozen ROI outcome changed in shadow comparison')
    report=summarize(shadow)
    report.update(original_photo_attempts=receipt['attempted'],
                  original_acquisition_failures=receipt['acquisition_failed'],
                  source_plan_sha256=receipt['source_plan_sha256'],
                  frozen_original_ROI_method=receipt['frozen_original_ROI_method'])
    shadow.to_csv(output/'screen_shadow.csv',index=False)
    report['comparison_csv_sha256']=digest(output/'screen_shadow.csv')
    (output/'screen_shadow_summary.json').write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n')
    (output/'source.json').write_bytes((Path(source)/'manifest.json').read_bytes())
    print('MOIS_ACTUAL_SHADOW_SCREEN_COMPLETE_NO_PHOTOS_DROPPED',json.dumps(report,sort_keys=True),flush=True)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    subs=p.add_subparsers(dest='command',required=True)
    w=subs.add_parser('work')
    for name in ('source','runtime','sam','output'):w.add_argument('--'+name,type=Path,required=True)
    w.add_argument('--group',type=int,required=True)
    c=subs.add_parser('collect')
    for name in ('source','workers','output'):c.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    if a.command=='work':work(a.source,a.runtime,a.sam,a.group,a.output)
    else:collect(a.source,a.workers,a.output)

if __name__=='__main__':main()
