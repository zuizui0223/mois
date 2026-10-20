#!/usr/bin/env python3
"""Prepare authentic RGB processing from ONLY one sealed NEW source photo tier.

No 93m-index rescan, source photo download or colour inference. Every photo ID
remains in the full target register. Compute buckets are exhaustive scheduling
partitions, not botanical/taxon/sampling caps, and all attached photos of one
original observation share a compute bucket. A source-photo must be absent
from the prior actual-image receipt before it can be queued again.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

from scripts.acquisition.mois_onepass_full93m_target_tiers import (
    COLUMNS, EXPECTED_PHOTOS, EXPECTED_SHARDS, EXPECTED_TAXA,
    RESULT_SCHEMA, TIER_ENDS, obs_draw, require, sha256,
    normalized_extension,
)
WORKER_COLUMNS = (
    "measurement_id", "photo_id", "inat_taxon_id", "species", "observation_id",
    "observer_id", "observed_on", "latitude", "longitude", "positional_accuracy",
    "photo_license", "photo_url_large", "image_filename", "source_arm",
    "selection_probability", "source_snapshot_date", "shard",
)
OUT_SCHEMA = "mois_full_source_observation_sampling_cohort_v1"
FCP_COUNT_SHA = "ef26c16252dd785f2d08bde3b18938c7624d73a54386fa3d3974dceb09ed97a9"


def observation_compute_bucket(uuid, n=64):
    require(1 <= n <= 1024, "Compute partition count is invalid")
    raw=hashlib.sha256(("mois-original-compute-bucket-v1|"+str(uuid)).encode()).digest()[:8]
    return int.from_bytes(raw,"big")%n

def read_prior_ids(paths):
    all_ids=set()
    for p in paths:
        with Path(p).open(newline="",encoding="utf-8") as f:
            reader=csv.DictReader(f)
            require("photo_id" in (reader.fieldnames or []), "Prior measured-photo receipt lacks photo_id")
            for r in reader:
                pid=str(r["photo_id"])
                require(pid.isdecimal(), "Previously measured original photo ID malformed")
                require(pid not in all_ids,"Duplicate photo ID in previous-measured source inputs")
                all_ids.add(pid)
    return all_ids


def read_taxa(path,*,check_sha=True):
    # The source count CSV is a separate independent FCP signed 42,111 taxa table.
    if check_sha:
        require(sha256(path)==FCP_COUNT_SHA, "Original source taxon name/count table SHA changed")
    with Path(path).open(newline="",encoding="utf-8") as f:
        reader=csv.DictReader(f)
        require({'inat_taxon_id','species'} <= set(reader.fieldnames or []),
                'Original FCP source taxon names missing')
        names={}
        for row in reader:
            tax=str(row['inat_taxon_id'])
            require(tax.isdecimal() and tax not in names and row['species'],
                    'Original source taxon ID duplicated/unresolved')
            names[tax]=row['species']
    return names


def queue_tier(registry,group_receipts,targets,counts,prior,out,*,
               tier,compute_bucket,compute_bucket_count=64,batch_size=256,
               check_taxa_sha=True,expected_taxa=EXPECTED_TAXA,
               expected_shards=EXPECTED_SHARDS):
    """Build ROI-v4 source queue from one tier and one exhaustive compute bucket.

    This is an acquisition/processing batch, not a scientific sample. Do not
    estimate global white composition from a partly completed compute bucket.
    """
    out=Path(out); targets=Path(targets)
    require(not out.exists() or not any(out.iterdir()),
            "Existing original RGB source roster is immutable")
    require(0<=tier<len(TIER_ENDS), "Photo-target tier index is invalid")
    require(0<=compute_bucket<compute_bucket_count<=1024,
            "Invalid exhaustive observation compute bucket")
    require(1<=batch_size<=256, "Original FCP ROI batch can contain at most 256 photos")
    full=json.loads(Path(registry).read_text())
    require(full.get('schema')=='mois_complete_full93m_original_photo_target_registry_v1'
            and full.get('exact_original_photo_IDs_in_all_disjoint_tiers')==EXPECTED_PHOTOS
            and full.get('original_named_taxa_from_full_source')==EXPECTED_TAXA
            and full.get('native_vs_introduced_not_hard_filtered') is True
            and full.get('source_only_NO_photo_pixels_or_hue_measured') is True,
            "No complete original-photo 93m precolour target registry")
    previous=read_prior_ids(prior)
    names=read_taxa(counts,check_sha=check_taxa_sha)
    require(len(names)==expected_taxa, "Incomplete original 42,111 taxon source names")
    audited_parts={}
    source_group_total=0
    for path in group_receipts:
        group=json.loads(Path(path).read_text())
        require(group.get('schema')==RESULT_SCHEMA
                and group.get('source_snapshot_SHA256')==full['full_source_manifest_sha256'],
                'Original photo tier group does not correspond to signed source index')
        for part in group['parts']:
            n=int(part['part'])
            require(n not in audited_parts and 0<=n<expected_shards,
                    'Overlapping or impossible 93m original photo target shard')
            audited_parts[n]=part
        source_group_total+=int(group['counts_by_disjoint_tier'][str(tier)])
    require(len(audited_parts)==expected_shards and set(audited_parts)==set(range(expected_shards)),
            'Complete 128 original-photo target shards must all be present')
    require(source_group_total==full['new_target_counts_by_tier'][str(tier)],
            'Source tier count differs between full frame and grouped signed registry')
    selected=[]; unsupported=[]; skipped=0; tier_photo_total=0; bucket_photo_total=0
    seen=set()
    for n in range(expected_shards):
        path=targets/f'tier-{tier:02d}'/f'part-{n:03d}.csv.gz'
        require(sha256(path)==audited_parts[n]['ranked_target_csv_sha256'][str(tier)],
                'Original precolour photo tier source bytes changed: '+str(n))
        with gzip.open(path,'rt',encoding='utf-8',newline='') as f:
            r=csv.DictReader(f)
            require(tuple(r.fieldnames or [])==COLUMNS,
                    'Original target tier columns changed')
            part_count=0
            for row in r:
                part_count+=1
                pid=row['photo_id'];uuid=row['observation_uuid']
                require(pid.isdecimal() and row['inat_taxon_id'] in names and uuid,
                        'Original photo ID/taxon/observation missing')
                require(row['tier']==str(tier) and float(row['observation_inclusion_score'])==obs_draw(uuid),
                        'Original deterministic observation hash selection drift')
                tier_photo_total+=1
                if observation_compute_bucket(uuid,compute_bucket_count)!=compute_bucket:
                    continue
                require(pid not in seen, 'Duplicate actual source photo ID within the compute bucket')
                seen.add(pid);bucket_photo_total+=1
                if pid in previous:
                    skipped+=1
                    continue
                ext=normalized_extension(row['image_extension'])
                if ext is None:
                    unsupported.append({
                        'photo_id':pid,'observation_uuid':uuid,
                        'source_image_extension':row['image_extension'],
                        'reason':'unsupported_by_frozen_FCP_ROI_image_loader',
                    })
                    continue
                mid='MOIS-'+hashlib.sha256(('mois-original-global-RGB-v1|'+pid).encode()).hexdigest()[:24].upper()
                selected.append({
                    'measurement_id':mid,'photo_id':pid,
                    'inat_taxon_id':row['inat_taxon_id'],
                    'species':names[row['inat_taxon_id']],
                    'observation_id':uuid,
                    'observer_id':row['observer_id'],
                    'observed_on':row['observed_on'],
                    'latitude':row['latitude'],
                    'longitude':row['longitude'],
                    'positional_accuracy':row['positional_accuracy'],
                    'photo_license':row['photo_license'],
                    'photo_url_large':f'https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.{ext}',
                    'image_filename':mid+'.jpg',
                    'source_arm':f'global93m_onepass_target_tier_{tier:02d}',
                    'selection_probability':str(TIER_ENDS[tier]),
                    'source_snapshot_date':'2026-09-27',
                    'shard':'',
                })
            require(part_count==audited_parts[n]['counts_by_disjoint_tier'][str(tier)],
                    'Signed tier photo IDs were added or lost in source')
    require(tier_photo_total==source_group_total
            and len(selected)+len(unsupported)+skipped==bucket_photo_total,
            'Onepass source tier compute bucket reconciliation failed')
    selected.sort(key=lambda d: hashlib.sha256(
        ('mois-prioritise-observation|'+d['photo_id']).encode()).digest())
    units=math.ceil(len(selected)/batch_size)
    out.mkdir(parents=True,exist_ok=True)
    (out/'shards').mkdir(exist_ok=False)
    unit_sha={}
    for unit in range(units):
        subset=selected[unit*batch_size:(unit+1)*batch_size]
        for row in subset:row['shard']=str(unit)
        file=out/'shards'/f'shard-{unit:03d}.csv'
        with file.open('w',encoding='utf-8',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=WORKER_COLUMNS,lineterminator='\n')
            writer.writeheader();writer.writerows(subset)
        unit_sha[str(unit)]=sha256(file)
    unsupported_path=out/'unsupported_selected_original_photo_formats.jsonl'
    with unsupported_path.open('w',encoding='utf-8') as f:
        for r in sorted(unsupported,key=lambda x:int(x['photo_id'])):
            f.write(json.dumps(r,sort_keys=True)+'\n')
    report={
        'schema':OUT_SCHEMA,
        'full_frame_photo_ids':EXPECTED_PHOTOS,
        'full_source_manifest_sha256':full['full_source_manifest_sha256'],
        'source_target_registry_sha256':sha256(registry),
        'full_source_original_tier':tier,
        'source_tier_probability_lower_bound':0 if tier==0 else TIER_ENDS[tier-1],
        'selection_probability':TIER_ENDS[tier],
        'selection_unit':'observation_uuid',
        'compute_bucket':compute_bucket,
        'compute_bucket_count':compute_bucket_count,
        'compute_bucket_is_exhaustive_processing_partition_not_statistical_inclusion_cap':True,
        'source_tier_total_photo_candidates':tier_photo_total,
        'source_tier_compute_bucket_total_photos':bucket_photo_total,
        'sampled_photo_ids_before_prior_dedup':bucket_photo_total,
        'previously_measured_photo_ids_skipped':skipped,
        'selected_new_photos_unsupported_by_frozen_loader':len(unsupported),
        'selected_new_photo_formats_unsupported_sha256':sha256(unsupported_path),
        'new_photo_attempts_planned':len(selected),
        'work_units':units,
        'worker_ids':list(range(units)),
        'group_ids':list(range(math.ceil(units/4))),
        'units_per_group':4,
        'photo_shard_sha256':unit_sha,
        'original_model_batch_size':batch_size,
        'sampled_before_any_hue_measurement':True,
        'source_original_images_measured':False,
        'colour_measured':False,
        'technical_unknowns_not_assumed_chromatic':True,
        'all_unselected_source_photos_eligible_at_higher_probability':True,
        'no_per_taxon_or_observation_photo_count_cap':True,
        'native_photo_status_not_inferred_or_hard_filtered':True,
        'ecological_analysis_run':False,
    }
    (out/'manifest.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print('MOIS_SOURCE_TIER_REUSE_ONLY_NEW_ORIGINAL_PHOTO_IDS',json.dumps({
        k:report[k] for k in ('full_source_original_tier','compute_bucket',
                             'source_tier_total_photo_candidates',
                             'previously_measured_photo_ids_skipped',
                             'new_photo_attempts_planned')},sort_keys=True),flush=True)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('registry','targets','counts','output'):
        p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--group-receipt',type=Path,action='append',required=True)
    p.add_argument('--previous',type=Path,action='append',default=[])
    p.add_argument('--tier',type=int,required=True)
    p.add_argument('--compute-bucket',type=int,required=True)
    p.add_argument('--compute-bucket-count',type=int,default=64)
    a=p.parse_args()
    queue_tier(a.registry,a.group_receipt,a.targets,a.counts,a.previous,a.output,
               tier=a.tier,compute_bucket=a.compute_bucket,
               compute_bucket_count=a.compute_bucket_count)


if __name__=='__main__':main()