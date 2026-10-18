#!/usr/bin/env python3
"""Measure frozen H9 candidates once; recover valid pixel checkpoints on resume.

Uses the unchanged, source-pinned FCP loader and ROI-v4 evaluator. H9's
observation-level flowering filter is not proof of per-image open petals.
The output is a labelled acquisition stratum, not a prevalence estimate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import time
from datetime import date
from pathlib import Path
from urllib.parse import urlsplit

import numpy as np
import pandas as pd

SOURCE = 'mois_fcp_H9_original_flowering_photo_refills_first256_unknown.csv'
SOURCE_RECEIPT = 'mois_fcp_H9_flowering_refills_photo_source_provenance.json'
SOURCE_SCHEMA = 'mois_fcp_H9_flowering_metadata_rescue_original_unclassifiable256_pilot_v1'
RESULT = 'mois_h9_measured_photo_labels.csv'
RECEIPT = 'mois_h9_measured_photo_receipt.json'
HUES = ('white','yellow_orange','red_pink','blue_purple')
PALETTE = ('white','yellow','orange','red','pink','magenta','purple','blue','bronze')
LICENSES = {'cc0','cc-by','cc-by-sa','cc-by-nc','cc-by-nc-sa'}
WORKER = ['measurement_id','image_filename','photo_license']


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError(message)


def digest_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def digest(path: Path) -> str:
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def write_json(path: Path, data: dict) -> None:
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(data,indent=2,sort_keys=True,allow_nan=False)+'\n')
    tmp.replace(path)


def load_source(root: Path) -> tuple[pd.DataFrame, dict]:
    rec=json.loads((root/SOURCE_RECEIPT).read_text())
    require(rec['schema']==SOURCE_SCHEMA and rec['new_photo_pixels_actually_measured'] is False,
            'Not a frozen H9 source-only candidate receipt')
    require(digest(root/SOURCE)==rec['source_new_photo_candidates_SHA256'],
            'Source digest mismatch')
    rows=pd.read_csv(root/SOURCE,dtype=str,keep_default_na=False)
    needed=set(WORKER)|{'photo_id','photo_url_large','inat_taxon_id','species','observer_id',
        'observation_id','latitude','longitude','observed_on','source_sampling_arm','source_annotation'}
    require(needed.issubset(rows.columns) and not {'morph','classified','white'}.intersection(rows),
            'Missing metadata or leaked colour in source')
    require(0<len(rows)<=256 and len(rows)==rec['rescued_flowering_annotated_new_photo_SOURCE_ONLY_IDs'],
            'Empty or inconsistent source size')
    for key in ('measurement_id','photo_id','observation_id','inat_taxon_id'):
        require(rows[key].is_unique and rows[key].ne('').all(), 'Duplicated/empty '+key)
    for row in rows.to_dict('records'):
        require(re.fullmatch(r'MOIS-[A-F0-9]{24}',row['measurement_id']) is not None,
                'Invalid measurement identifier')
        require(row['image_filename']==row['measurement_id']+'.jpg', 'Unsafe image filename')
        require(row['photo_id'].isdecimal() and row['inat_taxon_id'].isdecimal(), 'Invalid provider identifier')
        url=urlsplit(row['photo_url_large'])
        require(url.scheme=='https' and url.netloc=='inaturalist-open-data.s3.amazonaws.com'
            and not url.query and not url.fragment
            and re.fullmatch(r'/photos/'+re.escape(row['photo_id'])+r'/large\.(jpg|jpeg|png)',url.path),
            'Unsafe original photo URL')
        require(row['photo_license'] in LICENSES, 'Unsupported source photo licence')
        require(row['source_sampling_arm']=='fcp_h9_flowering_annotated_rescue_after_unknown'
            and row['source_annotation']=='flowers_and_fruits:flowers_(12:13)', 'Wrong sampling stratum')
        lat,lon=float(row['latitude']),float(row['longitude'])
        require(np.isfinite([lat,lon]).all() and -90<=lat<=90 and -180<=lon<=180, 'Invalid coordinates')
        require(date(1800,1,1)<=date.fromisoformat(row['observed_on'])<=date.today(), 'Invalid observation date')
    return rows,rec


class OneBatch:
    def __init__(self, rows: pd.DataFrame):
        self.rows=rows

    def get_batch(self, part: int, batch: int, size: int):
        require((part,batch,size)==(0,0,len(self.rows)), 'Expected one whole photo batch')
        return self.rows[WORKER].copy(),self.rows[WORKER+['photo_url_large']].copy(),1


def join_labels(source: pd.DataFrame, measured: pd.DataFrame) -> pd.DataFrame:
    required={'measurement_id','morph','measurement_status'}
    require(required.issubset(measured), 'Missing pixel outcome columns')
    require(measured.measurement_id.is_unique and set(source.measurement_id)==set(measured.measurement_id),
            'Pixel identifiers differ from frozen source')
    known=measured.morph.isin(HUES)
    require(measured.morph.isin((*HUES,'mixed_uncertain')).all()
        and known.eq(measured.measurement_status.eq('classified_four_state_morph')).all(),
        'Unknown colour was converted to a classified hue')
    extra=[c for c in measured if c not in source or c=='measurement_id']
    joined=source.merge(measured[extra],on='measurement_id',validate='one_to_one',sort=False)
    joined['classified']=joined.morph.isin(HUES)
    joined['white']=joined.morph.eq('white')
    return joined


def run(source: Path, runtime: Path, sam: Path, out: Path, *, executor=None) -> dict:
    start=time.monotonic()
    rows,source_rec=load_source(source)
    index=OneBatch(rows)
    expected_worker_sha=digest_bytes(rows[WORKER].to_csv(index=False).encode())
    measure=runtime/'scripts/analysis/measure_random_photo_first_blind_partition.py'
    detector=runtime/'data/atlas/qualification/roi_v4_training/jrc_yolo11n_last_v4.pt'
    expected_method=None
    if executor is None:
        from scripts.acquisition.run_fcp_fresh_image_partition import run_batch, verify_frozen_method
        expected_method=verify_frozen_method(measure,detector,sam)
        executor=run_batch
    out.mkdir(parents=True,exist_ok=True)
    binding={'source_csv_sha256':digest(source/SOURCE),
             'source_receipt_sha256':digest(source/SOURCE_RECEIPT)}
    bind=out/'source_binding.json'
    if bind.exists():
        require(json.loads(bind.read_text())==binding,'Cannot resume a different frozen photo source')
    else:
        require(not any(out.iterdir()),'Unbound nonempty output directory')
        write_json(bind,binding)
    work=out/'work'; terminal=work/'terminal_p000_b0000000.csv'; seal=terminal.with_suffix('.json')
    reuse=terminal.exists() or seal.exists()
    pixel_start=time.monotonic()
    if reuse:
        require(terminal.is_file() and seal.is_file(), 'Incomplete pixel checkpoint; do not overwrite')
        rec=json.loads(seal.read_text())
    else:
        rec=executor(Path('.'),0,0,len(rows),work,source=measure,detector=detector,sam=sam,indexed_source=index)
    require(rec['schema']=='fcp_fresh_blind_terminal_partition_v1' and rec['status']=='COMPLETE_WORK_UNIT',
            'Pixel work unit did not complete')
    require(rec['n_source_photo_ids']==rec['n_terminal_photo_ids']==len(rows)
        and rec['number_batches_for_this_partition']==1
        and rec['compute_batch_size_not_photo_cap']==len(rows)
        and rec['image_bytes_persisted'] is False and rec['reused_historical_measurement_csv'] is False,
            'Actual FCP work unit count/provenance mismatch')
    require(rec['selected_source_worker_sha256']==expected_worker_sha, 'Pixel source digest mismatch')
    require(rec['terminal_csv_sha256']==digest(terminal), 'Pixel terminal digest mismatch')
    if expected_method is not None:
        require(rec['method']==expected_method, 'Frozen ROI method changed')
    measured=pd.read_csv(terminal,dtype=str,keep_default_na=False)
    joined=join_labels(rows,measured)
    failed=int(joined.measurement_status.eq('image_acquisition_failed').sum())
    classified=int(joined.classified.sum())
    require(classified==rec['n_classifiable'] and failed==rec['n_acquisition_failed'], 'Outcome count mismatch')
    # FCP also validates the palette in production; check persisted seals on resume.
    fractions=[f'flower_fraction_{p}' for p in PALETTE]
    require(set(fractions).issubset(measured), 'Missing fixed-palette measurements')
    values=measured.loc[measured.morph.isin(HUES),fractions].apply(pd.to_numeric,errors='raise').to_numpy(float)
    require(np.isfinite(values).all() and (values>=0).all()
        and np.allclose(values.sum(axis=1),1,atol=1e-7), 'Invalid fixed-palette composition')
    pixel_seconds=time.monotonic()-pixel_start
    tmp=out/(RESULT+'.tmp'); joined.to_csv(tmp,index=False,lineterminator='\n'); tmp.replace(out/RESULT)
    summary={'schema':'mois_h9_one_process_measurement_v1',**binding,
        'attempted_photos':len(rows),'classified_photos':classified,'white_count':int(joined.white.sum()),
        'download_failures':failed,'hue_unknown_after_download':len(rows)-classified-failed,
        'classification_yield_per_attempt':classified/len(rows),
        'hue_counts':{h:int(joined.morph.eq(h).sum()) for h in HUES},
        'new_classifier_evaluator_invocations':0 if reuse or failed==len(rows) else 1,
        'reused_sealed_work_unit':reuse,'model_changed':False,
        'method_attestation':rec['method'],'outcome_csv_sha256':digest(out/RESULT),
        'pixel_stage_seconds_this_run':pixel_seconds,'elapsed_seconds_this_run':time.monotonic()-start,
        'annotation_stratum_not_unbiased_prevalence':True,
        'timing_is_not_a_paired_speedup_benchmark':True}
    write_json(out/RECEIPT,summary)
    print('MOIS_H9_REAL_MEASUREMENT '+json.dumps(summary,sort_keys=True),flush=True)
    return summary


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--runtime',type=Path,required=True)
    p.add_argument('--sam',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    args=p.parse_args()
    run(args.source,args.runtime,args.sam,args.out)
