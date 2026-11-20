#!/usr/bin/env python3
"""MOIS 6 positive-control FCP ROI9 photos: paired /large versus /original.

Source: immutable 19-photo signed FCP ROI-v4 output, 6 previously numerically
valid images. Same image ID, frozen FCP method/ROI gates; fully SHA-certified
originals. These are pixel diagnostics, not confirmed true perianth colours.
"""
from __future__ import annotations
import argparse,io,json
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd
from scripts.acquisition.mois_full_index_real_original_photo_measure import (
    SCHEMA as FCP_SOURCE_SCHEMA
)
from scripts.acquisition.mois_FCP13_original_resolution_retest_20261121 import (
    ZIP_SHA, SOURCE, ORIGINAL_WORKER, ORIGINAL_WORKER_RECEIPT, sha, check,
    original_url
)
from scripts.analysis.mois_species_first_palette9 import PALETTE,class_labels

NUM_SOURCE=19
NUM_CONTROLS=6


def source_original_six(path):
    check(sha(path)==ZIP_SHA,'Previous FCP19 signed source ZIP drift')
    with ZipFile(path) as z:
        manifest=json.loads(z.read('source/manifest.json'))
        rec=json.loads(z.read(ORIGINAL_WORKER_RECEIPT))
        b=z.read(SOURCE)
        c=z.read(ORIGINAL_WORKER)
    check(sha(b)==manifest['photo_shard_sha256']['0']
          and sha(c)==rec['terminal_sha256']
          and rec['attempted']==NUM_SOURCE and rec['classifiable']==NUM_CONTROLS,
          'Original FCP19 source/terminal authenticity drift')
    source=pd.read_csv(io.BytesIO(b),dtype=str,keep_default_na=False)
    terminal=pd.read_csv(io.BytesIO(c),dtype=str,keep_default_na=False)
    valid=terminal.loc[terminal.measurement_status.eq('classified_four_state_morph')].copy()
    check(len(valid)==NUM_CONTROLS and valid.measurement_id.is_unique,
          'Original valid six FCP IDs drift')
    d=source.merge(valid,on='measurement_id',how='inner',suffixes=('','_large'),validate='one_to_one')
    check(len(d)==NUM_CONTROLS and d.photo_id.is_unique and d.observation_id.is_unique
          and d.original_image_sha256.str.len().eq(64).all(),
          'Original six photo/observation identity drift')
    return d,rec


def prep(args):
    root=Path(args.out)
    check(not root.exists(),'Source six cannot overwrite earlier genuine source')
    d,rec=source_original_six(args.previous)
    source_cols=['measurement_id','image_filename','photo_id','observation_id',
        'inat_taxon_id','candidate_species','photo_url_large','photo_license',
        'shard','original_image_sha256','thumbnail_visual_triage']
    out=d[source_cols].copy()
    out['previous_large_url']=out.photo_url_large.copy()
    out['photo_url_large']=[original_url(u,p) for u,p in zip(out.photo_url_large,out.photo_id)]
    out['full_original_byte_SHA256_expected']=d.original_image_sha256.to_numpy()
    out['original_resolution_variant']=True
    out['prior_large_measurement_status']=d.measurement_status.to_numpy()
    check(out.measurement_id.is_unique and out.observation_id.is_unique
          and out.photo_url_large.str.contains('/original\\.').all(),
          'True original six source integrity failed')
    root.mkdir(parents=True)
    (root/'shards').mkdir()
    shard=root/'shards/shard-000.csv'
    out.to_csv(shard,index=False,lineterminator='\n')
    manifest={'schema':FCP_SOURCE_SCHEMA,
       'source_arm':'MOIS_FCP6_VALID_LARGE_TO_TRUE_ORIGINAL_PAIRED_CONTROLS',
       'colour_measured':False,'work_units':1,'worker_ids':[0],
       'photo_shard_sha256':{'0':sha(shard)},
       'new_photo_attempts_planned':NUM_CONTROLS,
       'prior_TRUE_FCP19_ZIP_SHA256':ZIP_SHA,
       'original_pixels_fetched_by_this_source':False,
       'original_full_image_byte_SHA_pinned':True,
       'prior_large_numeric_nineROI_valid_all6':True,
       'FCP_model_weights_thresholds_unchanged':True,
       'true_botanical_petal_color_verified':False,
       'global_source_37104_photo9_changed':False}
    (root/'manifest.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n')
    (root/'source_receipt.json').write_text(json.dumps({
       'n_original_fullres_source_photos':NUM_CONTROLS,
       'old_FCP_method_attestation':rec['method'],
       'source_shard_SHA256':sha(shard),
       'manifest_SHA256':sha(root/'manifest.json')},indent=2,sort_keys=True)+'\n')
    print('MOIS_FCP6_PRIOR_PASS_CONTROLS_PINNED_SOURCE',
          json.dumps({'n_controls':NUM_CONTROLS,'shard_SHA256':sha(shard)}),flush=True)


def work(args):
    # Reuse the source-verifying custom full-original loader of the 13-failure
    # audit. The expected source batch cardinality is the only isolated change.
    from scripts.acquisition import mois_FCP13_original_resolution_retest_20261121 as same_photo
    check(same_photo.FAIL_COUNT==13,'Unexpected change to full original comparison loader')
    same_photo.FAIL_COUNT=NUM_CONTROLS
    same_photo.work(args)


def audit(args):
    root=Path(args.out)
    check(not root.exists(),'Refuse overwrite of positive-control FCP results')
    old,old_receipt=source_original_six(args.previous)
    source=Path(args.source)
    manifest=json.loads((source/'manifest.json').read_text())
    shard=source/'shards/shard-000.csv'
    check(manifest['schema']==FCP_SOURCE_SCHEMA and manifest['new_photo_attempts_planned']==NUM_CONTROLS
          and sha(shard)==manifest['photo_shard_sha256']['0'],
          'Positive source was not pinned')
    expect=pd.read_csv(shard,dtype=str,keep_default_na=False)
    folder=Path(args.worker)/'shard_000'
    rec=json.loads((folder/'receipt.json').read_text())
    measured=folder/'photos.csv'
    tech=folder/'sealed/terminal_p000_b0000000.csv'
    check(rec['attempted']==NUM_CONTROLS
          and rec['photos_sha256']==sha(measured)
          and rec['terminal_sha256']==sha(tech)
          and rec['method']==old_receipt['method']
          and rec['source_sha256']==sha(shard),
          'FCP full original six controls not exact model/source replica')
    worker=pd.read_csv(measured,dtype=str,keep_default_na=False)
    output=expect.merge(worker,on='measurement_id',suffixes=('','_original_worker'),validate='one_to_one')
    check(len(output)==NUM_CONTROLS and output.photo_id.is_unique,
          'Expected six cross-image pairs lost')
    for key in ['photo_id','observation_id','photo_license','photo_url_large']:
        check(output[key].eq(output[key+'_original_worker']).all(),
              'Worker changed frozen per-photo ID or URL '+key)
    output=output.merge(old[['measurement_id','measurement_status','morph','image_sha256',
        *['flower_fraction_'+p for p in PALETTE]]],on='measurement_id',
        how='inner',suffixes=('_original','_large'),validate='one_to_one')
    check(len(output)==NUM_CONTROLS,'Results not paired with same large variant')
    acquired=~output.measurement_status_original.eq('image_acquisition_failed')
    matching_sha=output.image_sha256_original.eq(output.full_original_byte_SHA256_expected)
    check(matching_sha[acquired].all(),
          'Original full JPEG byte SHA does not match previous 46-image provenance')
    arrays=['flower_fraction_'+x for x in PALETTE]
    v=output[[k+'_original' for k in arrays]].apply(pd.to_numeric,errors='coerce').to_numpy(float)
    w=output[[k+'_large' for k in arrays]].apply(pd.to_numeric,errors='coerce').to_numpy(float)
    labels_old=class_labels(w)
    eligible=output.morph_original.isin(['white','yellow_orange','red_pink','blue_purple']).to_numpy(bool)
    numeric=acquired.to_numpy(bool)&eligible&np.isfinite(v).all(axis=1)&(v>=0).all(axis=1)&(abs(v.sum(axis=1)-1)<.001)
    check(int(np.count_nonzero(numeric))==rec['classifiable'],
          'Original fullres numerical-valid control count inconsistent')
    labels_new=np.full(NUM_CONTROLS,'',dtype=object)
    if numeric.any(): labels_new[numeric]=class_labels(v[numeric])
    output['large_photo9_strict_label']=labels_old
    output['original_photo9_strict_label']=labels_new
    output['FCP_fullres_nine_numeric_valid']=numeric
    output['actual_full_original_sha_matches_prior46']=matching_sha
    output['photo9_categorical_label_stable_nonempty']=(labels_old==labels_new)&(labels_new!='')
    output['numeric9_L1_original_minus_large']=np.where(numeric,abs(w-v).sum(axis=1),np.nan)
    output['actual_petal_organ_independently_verified']=False
    output['species_photo_colour_would_be_promoted']=False
    root.mkdir(parents=True)
    fp=root/'same_six_control_original_large_RGB9_photo_paired.csv'
    output.to_csv(fp,index=False,lineterminator='\n')
    report={'same_six_photo_IDs_and_observation_UUIDs':NUM_CONTROLS,
      'original_fullres_successful_downloads':int(acquired.sum()),
      'original_image_byte_SHA_matches_prior46':int(matching_sha[acquired].sum()),
      'prior_large_image_9colour_numeric':NUM_CONTROLS,
      'original_image_9colour_numeric':int(numeric.sum()),
      'prior_large_single_9hue_label':int(np.count_nonzero(labels_old!='')),
      'original_single_9hue_label':int(np.count_nonzero(labels_new!='')),
      'n_stable_nonempty_paired_photo9_strict_labels':int(output.photo9_categorical_label_stable_nonempty.sum()),
      'n_original_failed_acquisition':int((~acquired).sum()),
      'same_original_frozen_FCP_method_attestation':True,
      'true_petal_organ_verified':False,'new_species_colour_labels_promoted':0,
      'n_photo_ID_and_observation_UUID_disjoint_from_13_failures':NUM_CONTROLS,
      'original_source19_ZIP_SHA256':ZIP_SHA,
      'output_CSV_SHA256':sha(fp)}
    (root/'receipt.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print('MOIS_FCP6_PRIOR_VALID_FULL_ORIGINAL_CONTROL_REAL_RESULT',
          json.dumps(report,sort_keys=True),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser()
    sub=p.add_subparsers(dest='cmd',required=True)
    for cmd in ('source','work','audit'):
        v=sub.add_parser(cmd)
        req=(('previous','out') if cmd=='source' else
             ('source','runtime','sam','out') if cmd=='work' else
             ('previous','source','worker','out'))
        for field in req: v.add_argument('--'+field,required=True)
    a=p.parse_args()
    if a.cmd=='source':prep(a)
    elif a.cmd=='work':work(a)
    else:audit(a)