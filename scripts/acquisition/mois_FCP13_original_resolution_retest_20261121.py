#!/usr/bin/env python3
"""MOIS: Compare frozen FCP ROI-v4 on true /original image bytes vs 500px large
variant for 13 source-verified failures. Preserve strict botanical and SHA gates.

The same detector and thresholds are used in both arms. This test is NOT an
expert perianth/flower confirmation and never updates the 37,104-species atlas.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import re
from pathlib import Path
from urllib.parse import urlsplit,urlunsplit
from zipfile import ZipFile

import numpy as np
import pandas as pd
from scripts.acquisition.mois_full_index_real_original_photo_measure import (
    SCHEMA as FCP_SOURCE_SCHEMA,
)
from scripts.analysis.mois_species_first_palette9 import PALETTE,class_labels

ZIP_SHA='92315ef6679ced282fd822dcb08336a541add546afa949fb657cda76f6fd90ea'
CSV_SHA='4a45d8bdd6e1900c3597d19f2c39239f7e81cda967e4d6cb9f589442fefacb50'
SOURCE='source/shards/shard-000.csv'
ORIGINAL_WORKER='pixel/shard_000/sealed/terminal_p000_b0000000.csv'
ORIGINAL_WORKER_RECEIPT='pixel/shard_000/receipt.json'
BASELINE_COUNT=19
FAIL_COUNT=13


def check(ok,msg):
    if not bool(ok):raise ValueError('MOIS_FCP13_FULL_ORIGINAL_VARIANT: '+msg)


def sha(raw):
    if isinstance(raw,(str,Path)):
        raw=Path(raw).read_bytes()
    return hashlib.sha256(raw).hexdigest()


def read_audited_source(previous,fail_csv):
    check(sha(previous)==ZIP_SHA,'Signed 19 original FCP worker ZIP SHA drift')
    check(sha(fail_csv)==CSV_SHA,'Frozen 13 source failure CSV SHA drift')
    with ZipFile(previous) as z:
        manifest=json.loads(z.read('source/manifest.json'))
        old_receipt=json.loads(z.read(ORIGINAL_WORKER_RECEIPT))
        old_bytes=z.read(SOURCE)
        prior_term_bytes=z.read(ORIGINAL_WORKER)
        source_csv_sha=sha(old_bytes)
        check(source_csv_sha==manifest['photo_shard_sha256']['0']
              and sha(prior_term_bytes)==old_receipt['terminal_sha256'],
              'Source shard or FCP previous pixel result receipt violated')
        check(old_receipt['attempted']==BASELINE_COUNT
              and old_receipt['classifiable']==6
              and old_receipt['unclassifiable_acquired']==FAIL_COUNT
              and old_receipt['acquisition_failed']==0,
              'Historical 19 original photos/13 technical failures not same cohort')
    source=pd.read_csv(io.BytesIO(old_bytes),dtype=str,keep_default_na=False)
    term=pd.read_csv(io.BytesIO(prior_term_bytes),dtype=str,keep_default_na=False)
    fail=pd.read_csv(fail_csv,dtype=str,keep_default_na=False)
    check(len(source)==len(term)==BASELINE_COUNT and len(fail)==FAIL_COUNT
          and source.photo_id.is_unique and fail.photo_id.is_unique
          and source.observation_id.is_unique and fail.observation_uuid.is_unique,
          'Original photo and observation uniqueness changed')
    d=fail.merge(source,on='photo_id',how='inner',validate='one_to_one')
    check(len(d)==FAIL_COUNT
          and d.observation_uuid.eq(d.observation_id).all()
          and d.species.eq(d.candidate_species).all()
          and d.source_full_original_SHA256.eq(d.original_image_sha256).all()
          and d.photo_license_x.eq(d.photo_license_y).all()
          and d.fcp_large_variant_URL.eq(d.photo_url_large).all(),
          'Original 13 IDs, UUIDs, licence, image URL or full SHA do not match source')
    term['photo_id']=term.measurement_id.map(
        source.set_index('measurement_id').photo_id)
    old=term.loc[term.photo_id.isin(fail.photo_id)].copy()
    check(len(old)==FAIL_COUNT and old.measurement_status.isin(
        ['not_evaluable_roi_or_flip_gate','not_evaluable_ambiguous_palette_composition']).all(),
        'Some original technical failures were incorrectly promoted to colours')
    check(old.image_sha256.map(lambda x:bool(re.fullmatch('[a-f0-9]{64}',x))).all()
          and not old.image_sha256.eq(d.set_index('photo_id').loc[
              old.photo_id,'original_image_sha256'].to_numpy()).any(),
          'Original large and actual original image SHA unexpectedly identical')
    return d,old,old_receipt


def original_url(url,photo_id):
    u=urlsplit(url)
    check(u.scheme=='https' and u.netloc=='inaturalist-open-data.s3.amazonaws.com'
          and not u.query and not u.fragment,
          'Previous photo URL host / transport changed')
    m=re.fullmatch(r'/photos/(\d+)/large\.(jpg|jpeg|png)',u.path)
    check(m is not None and m.group(1)==str(photo_id),
          'Invalid large variant URL / photo ID binding')
    return urlunsplit((u.scheme,u.netloc,f'/photos/{photo_id}/original.{m.group(2)}','',''))


def write_source(args):
    root=Path(args.out)
    check(not root.exists(),'Never overwrite signed source image shard')
    d,old,prior_receipt=read_audited_source(args.previous,args.failed)
    new=d[['measurement_id','image_filename','photo_id','observation_id',
            'inat_taxon_id','candidate_species','photo_url_large','photo_license_y',
            'shard','original_image_sha256','thumbnail_visual_triage']].copy()
    new.rename(columns={'photo_license_y':'photo_license'},inplace=True)
    new['previous_large_url']=new.photo_url_large
    new['photo_url_large']=[original_url(u,p) for u,p in zip(new.photo_url_large,new.photo_id)]
    new['original_resolution_variant']=True
    new['prior_large_measurement_status']=d.measurement_status.to_numpy()
    new['full_original_byte_SHA256_expected']=d.source_full_original_SHA256.to_numpy()
    check(new.measurement_id.is_unique and new.photo_id.is_unique
          and new.observation_id.is_unique and new.photo_url_large.str.contains('/original\\.').all(),
          'New full-res original image input identity corrupt')
    root.mkdir(parents=True)
    folder=root/'shards'
    folder.mkdir()
    file=folder/'shard-000.csv'
    new.to_csv(file,index=False,lineterminator='\n')
    manifest={'schema':FCP_SOURCE_SCHEMA,
        'source_arm':'MOIS_FCP13_SAME_IMAGES_TRUE_ORIGINAL_VARIANT_FROZEN_METHOD',
        'colour_measured':False,
        'photo_shard_sha256':{'0':sha(file)},
        'work_units':1,'worker_ids':[0],'new_photo_attempts_planned':FAIL_COUNT,
        'prior_true_FCP19_artifact_SHA256':ZIP_SHA,
        'prior_13_technical_unknown_source_SHA256':CSV_SHA,
        'photo_ID_and_observation_UUID_unchanged_from_previous':True,
        'original_full_photo_byte_SHA256_expected_present':True,
        'original_pixels_fetched_by_this_source':False,
        'FCP_model_and_strict_gates_relaxed':False,
        'true_botanical_petals_reviewed':False,
        'primary_37104_species_catalogue_append_authorized':False}
    (root/'manifest.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n')
    (root/'source_receipt.json').write_text(json.dumps({
         'n_original_fullres_source_photos':len(new),
         'n_original_identical_photos_with_frozen_large_variant_failures':FAIL_COUNT,
         'old_FCP_method_attestation':prior_receipt['method'],
         'source_shard_SHA256':sha(file),'manifest_SHA256':sha(root/'manifest.json')
    },indent=2,sort_keys=True)+'\n')
    print('MOIS_FCP13_FULLRES_SOURCE_VERIFIED',json.dumps({
         'n_source_photos':len(new),'shard_SHA256':sha(file)}),flush=True)


def audit(args):
    root=Path(args.out)
    check(not root.exists(),'Never overwrite actual original pixel retest results')
    d,old,prev_receipt=read_audited_source(args.previous,args.failed)
    folder=Path(args.source)
    manifest=json.loads((folder/'manifest.json').read_text())
    shard=folder/'shards/shard-000.csv'
    check(manifest['schema']==FCP_SOURCE_SCHEMA and manifest['new_photo_attempts_planned']==FAIL_COUNT
          and not manifest['FCP_model_and_strict_gates_relaxed']
          and sha(shard)==manifest['photo_shard_sha256']['0'],
          'Frozen source manifest not authentic or thresholds relaxed')
    expected=pd.read_csv(shard,dtype=str,keep_default_na=False)
    work=Path(args.worker)/'shard_000'
    receipt=json.loads((work/'receipt.json').read_text())
    outcome=work/'photos.csv'
    terminal=work/'sealed/terminal_p000_b0000000.csv'
    check(receipt['schema']=='mois_full_original93m_real_RGB_photo_unit_v1'
          and receipt['attempted']==FAIL_COUNT
          and receipt['source_sha256']==sha(shard)
          and receipt['photos_sha256']==sha(outcome)
          and receipt['terminal_sha256']==sha(terminal)
          and receipt['method']==prev_receipt['method'],
          'Source/shard/terminal or unchanged original FCP method invalid')
    measured=pd.read_csv(outcome,dtype=str,keep_default_na=False)
    tech=pd.read_csv(terminal,dtype=str,keep_default_na=False)
    check(len(measured)==len(tech)==FAIL_COUNT
          and measured.measurement_id.is_unique
          and set(expected.measurement_id)==set(measured.measurement_id)
          and set(tech.measurement_id)==set(expected.measurement_id),
          'Original source / terminal result ID difference')
    # For duplicate image/photo/observation fields, all existing identities must match.
    joined=expected.merge(measured,on='measurement_id',how='inner',suffixes=('','_worker'),validate='one_to_one')
    for col in ['photo_id','observation_id','inat_taxon_id','candidate_species','photo_license','photo_url_large']:
        check(joined[col].astype(str).equals(joined[col+'_worker'].astype(str)),
              'Worker mutated photo identity, licensing or URL: '+col)
    joined=joined.merge(tech[['measurement_id','roi_status','failure_reasons',
               'flower_effective_pixels','horizontal_flip_mask_iou',
               'horizontal_flip_colour_delta_e']],on='measurement_id',validate='one_to_one')
    joined=joined.merge(old[['measurement_id','measurement_status','morph',
               'image_sha256','flower_effective_pixels','horizontal_flip_mask_iou',
               'horizontal_flip_colour_delta_e','failure_reasons']],
               on='measurement_id',suffixes=('_fullres','_large'),validate='one_to_one')
    check(len(joined)==FAIL_COUNT and joined.photo_id.is_unique and joined.observation_id.is_unique,
          'Unexpected result join expansion')
    sha_ok=joined.image_sha256_fullres.eq(joined.full_original_byte_SHA256_expected)
    acquired=~joined.measurement_status_fullres.eq('image_acquisition_failed')
    check(sha_ok[acquired].all(),
          'Actual frozen FCP loader photo bytes differ from original full-resolution source SHA')
    check((~acquired|sha_ok).all(), 'Byte SHA mismatch of successful FCP images')
    frac=['flower_fraction_'+p for p in PALETTE]
    check(set(frac).issubset(joined.columns),'Original full-res numerical FCP nine components absent')
    vals=joined[frac].apply(pd.to_numeric,errors='coerce').to_numpy(float)
    has_roi=joined.morph_fullres.isin(['white','yellow_orange','red_pink','blue_purple']).to_numpy(bool)
    good=(acquired.to_numpy(bool)&has_roi&np.isfinite(vals).all(axis=1)
        & (vals>=0).all(axis=1)&(np.abs(vals.sum(axis=1)-1)<.001))
    check(int(np.count_nonzero(good))==receipt['classifiable'],
          'Original fullres classified ROI count not numeric-consistent')
    provisional=np.full(FAIL_COUNT,'',dtype=object)
    if good.any():provisional[good]=class_labels(vals[good])
    joined['original_variant_byte_SHA_matches_prior_photo_source']=sha_ok
    joined['fullres_numeric9_FCP_valid']=good
    joined['fullres_provisional_palette_label']=provisional
    joined['botanical_true_perianth_color_confirmed']=False
    joined['allow_to_change_species_global_photo9']=False
    joined['FCP_method_thresholds_unchanged']=True
    root.mkdir(parents=True)
    output=root/'FCP13_large_vs_full_original_with_paired_failure_reasons.csv'
    joined.to_csv(output,index=False,lineterminator='\n')
    types=joined.groupby('prior_large_measurement_status').agg(
        n_original_source_photos=('photo_id','size'),
        n_frozen_fullres_ROI_numeric=('fullres_numeric9_FCP_valid','sum'),
        n_exact_full_original_SHA=('original_variant_byte_SHA_matches_prior_photo_source','sum')
        ).reset_index().to_dict('records')
    report={'source_13_provisional_thumbnail_candidates':FAIL_COUNT,
        'prior_signed_FCP19_zip_SHA256':ZIP_SHA,'prior13_failure_csv_SHA256':CSV_SHA,
        'n_same_photo_identity':len(joined),'n_full_original_FCP_worker_image_downloaded':int(acquired.sum()),
        'n_original_FCP_worker_SHA_match':int(sha_ok[acquired].sum()),
        'n_original_variant_FCP_9_numeric_ROI':int(good.sum()),
        'n_original_variant_strict_label':int(np.count_nonzero(provisional!='')),
        'n_remaining_unclassifiable':int((acquired.to_numpy()&~good).sum()),
        'n_actual_original_photo_reacquisition_failed':int((~acquired).sum()),
        'n_were_13_large_variant_technical_failures':len(joined),
        'n_frozen_model_rescues':int(good.sum()),
        'n_photo_fullres_Sha_mismatch':int((acquired&~sha_ok).sum()),
        'baseline_same_FCP_method':True,
        'n_botanically_confirmed_petal_colour':0,
        'n_species_global_photo9_newly_promoted':0,
        'same_photo_no_new_independent_observations':True,
        'original_minus_large_is_not_an_independent_study':True,
        'type_summaries':types,
        'method_attestation':receipt['method'],
        'output_CSV_sha256':sha(output)}
    (root/'receipt.json').write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n')
    print('MOIS_FCP13_TRUE_FULL_ORIGINAL_BYTES_REAL_RETEST',
          json.dumps({k:v for k,v in report.items() if k!='method_attestation'},sort_keys=True),flush=True)

def work(args):
    """Use FCP's documented loader injection, preserving every model/check.

    Frozen image_loader only allows /large.* paths. A pinned same-image
    /original.* loader is authorized here by verifying provider host and photo
    identity, 12 MiB cap, decoded JPEG and signed full image raw-byte SHA.
    """
    from urllib.request import Request,urlopen
    from urllib.error import HTTPError,URLError
    from io import BytesIO
    from PIL import Image,UnidentifiedImageError
    import time
    from scripts.acquisition.run_fcp_fresh_image_partition import (
        run_batch,verify_frozen_method,_pace_provider_request,
        MAX_IMAGE_BYTES,PermanentOriginalPhotoError,ProviderTransientAccessError,
    )
    from scripts.acquisition.mois_full_index_real_original_photo_measure import measure_group
    src=Path(args.source)
    manifest=json.loads((src/'manifest.json').read_text())
    shard=src/'shards/shard-000.csv'
    check(manifest['schema']==FCP_SOURCE_SCHEMA
          and manifest['new_photo_attempts_planned']==FAIL_COUNT
          and manifest['photo_shard_sha256']['0']==sha(shard),
          'Pinned original FCP source is missing or mutated')
    items=pd.read_csv(shard,dtype=str,keep_default_na=False)
    expected=dict(zip(items.photo_url_large,items.full_original_byte_SHA256_expected))
    check(len(expected)==FAIL_COUNT and len(items)==FAIL_COUNT,
          'Original photo URL and source full image SHA count mismatch')
    original_runtime=Path(args.runtime)
    sam=Path(args.sam)
    measured_method=verify_frozen_method(
        original_runtime/'scripts/analysis/measure_random_photo_first_blind_partition.py',
        original_runtime/'data/atlas/qualification/roi_v4_training/jrc_yolo11n_last_v4.pt',sam)
    check(measured_method==json.loads((src/'source_receipt.json').read_text())['old_FCP_method_attestation'],
          'Frozen image ROI model/method differs from historical 19-photo comparison')

    def pinned_true_original_loader(url):
        check(url in expected,'FCP requested a URL not pinned to a source photo')
        u=urlsplit(url)
        check(u.scheme=='https' and u.netloc=='inaturalist-open-data.s3.amazonaws.com'
              and not u.query and not u.fragment
              and re.fullmatch(r'/photos/[0-9]+/original\.(?:jpg|jpeg|png)',u.path),
              'Non original-resolution provider URL in fullres FCP batch')
        req=Request(url,headers={'Accept':'image/*',
                  'User-Agent':'MOIS-botanically-unverified-original13-FCP/1.0'})
        for attempt in range(3):
            _pace_provider_request()
            try:
                with urlopen(req,timeout=45) as reply:
                    raw=reply.read(MAX_IMAGE_BYTES+1)
                if not raw or len(raw)>MAX_IMAGE_BYTES:
                    raise PermanentOriginalPhotoError('Original image outside frozen 12MiB bounded-source cap')
                try:
                    with Image.open(BytesIO(raw)) as img:
                        img.verify()
                except (UnidentifiedImageError,OSError,Image.DecompressionBombError) as exc:
                    raise PermanentOriginalPhotoError('Full original image failed decode') from exc
                # The exact bytes of this original /original photo were already
                # downloaded/sha-certified by the first 46-photo artifact.
                check(sha(raw)==expected[url],
                      'Re-fetched actual original JPEG byte SHA disagrees with original source receipt')
                return raw
            except HTTPError as exc:
                if exc.code in (404,410):
                    raise PermanentOriginalPhotoError(
                        f'Exact S3 original image HTTP {exc.code}') from exc
                if exc.code not in (408,425,429,500,502,503,504) or attempt==2:
                    raise ProviderTransientAccessError(
                        f'HTTP {exc.code} for prior signed original image; stop work unit') from exc
            except (URLError,TimeoutError,ConnectionError) as exc:
                if attempt==2:
                    raise ProviderTransientAccessError(
                        'Source provider transient failure for previously signed original') from exc
            time.sleep(2**attempt)
        raise ProviderTransientAccessError('Exhausted signed original image retry budget')

    results=measure_group(src,0,original_runtime,sam,Path(args.out),
                          runner=run_batch,loader=pinned_true_original_loader)
    check(len(results)==1 and results[0]['attempted']==FAIL_COUNT,
          'FCP true-original source/worker output incomplete')
    print('MOIS_FCP13_TRUE_ORIGINAL_CUSTOM_LOADER_SOURCE_PINNED_WORKER',
          json.dumps({'n_photos':FAIL_COUNT,
                      'n_classifiable':results[0]['classifiable'],
                      'n_acquisition_failed':results[0]['acquisition_failed'],
                      'byte_SHA_expected_for_each_original':True,
                      'frozen_FCP_method_reused':True},sort_keys=True),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    sub=p.add_subparsers(dest='cmd',required=True)
    for cmd in ('source','work','audit'):
        t=sub.add_parser(cmd)
        arguments=(('previous','failed','out') if cmd=='source'
                   else ('source','runtime','sam','out') if cmd=='work'
                   else ('previous','failed','source','worker','out'))
        for name in arguments:
            t.add_argument('--'+name,required=True)
    args=p.parse_args()
    if args.cmd=='source':write_source(args)
    elif args.cmd=='work':work(args)
    else:audit(args)