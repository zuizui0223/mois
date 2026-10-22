#!/usr/bin/env python3
"""Build a SHA-pinned real-photo ROI queue from already routed tier06 buckets.

Reads only the 33,464-original-photo signed tier06/group00 routed source
archive, not the 93m index or all 128 source taxon/observation partitions.
Every original photo in the fixed first 64 of 512 observation-hash buckets is
scheduled unless previously attempted, separately benchmarked or unsupported.
The 448 other buckets are preserved in the signed source release. This is a
bounded *processing* cohort, not a complete world-flora probability estimate.

Original fixed FCP ROI-v4, no lowres or nativeness hard exclusion. No pixels
are downloaded by this script. This metadata plan is immutable and source-
SHA verified before its first photo is ever fetched by the warm worker.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import tarfile

from scripts.acquisition.mois_onepass_full93m_target_tiers import (
    COLUMNS, EXPECTED_PHOTOS, EXPECTED_TAXA, TIER_ENDS,
    obs_draw, normalized_extension, require, sha256,
)
from scripts.acquisition.mois_queue_from_full93m_target_tiers import (
    FCP_COUNT_SHA, WORKER_COLUMNS, OUT_SCHEMA, read_taxa,
    observation_compute_bucket,
)
from scripts.acquisition.mois_singlepass_signed_tier_bucket_router import ROUTE_SCHEMA

SOURCE_TIER = 6
SOURCE_GROUP = 0
BUCKETS = 512
FIRST_BUCKETS = 64
ORIGINAL_BATCH_SIZE = 256
MAX_NEW_ATTEMPTS = 12_000
ROUTE_RECEIPT_SHA = "c493f29beb2782795bee19d4223c2c026cd59a5a4f479b39bca8594f8bccde6f"
ROUTED_TAR_SHA = "c1ff70564bf35c74c4d656511cfd51ad1682bd00b5bf607e9e739011230902f0"
TECHNICAL_RESERVATION_SHA = "7038130d6a3cc45100f11a3994da7d6735c554dcb8e6bc829c27f5bef82f2639"
PRIOR_RECEIPT_SHA = "91a9d72fdf19e298d9261818247cd73e07feb3b0493f3355c75cdafab2e23f62"
EXPECTED_PRIOR_IDS = 839_472
SOURCE_ARM = "global93m_tier06_group00_buckets000to063_warm_R1"


def sha_bytes(b):
    return hashlib.sha256(b).hexdigest()


def previous_ids(parquet, receipt, *, expected_count=EXPECTED_PRIOR_IDS,
                 expected_receipt_sha=PRIOR_RECEIPT_SHA):
    import pandas as pd
    if expected_receipt_sha is not None:
        require(sha256(receipt)==expected_receipt_sha,
                "Prior original pixel checkpoint receipt changed")
    report=json.loads(Path(receipt).read_text())
    require(report.get("schema")=="mois_sealed_actual_pixel_checkpoint_v1"
            and report.get("combined_unique_photo_attempts")==expected_count,
            "Unrecognized prior real original-photo checkpoint")
    require(sha256(parquet)==report["output_sha256"]["combined_actual_photo_checkpoint.parquet"],
            "Prior genuine-photo rows SHA mismatch")
    frame=pd.read_parquet(parquet,columns=["photo_id","observation_id"])
    ids=frame.photo_id.astype(str)
    require(len(ids)==expected_count and ids.is_unique and ids.str.fullmatch(r"\d+").all(),
            "Prior original photo identity malformed or duplicated")
    obs={str(x).strip() for x in frame.observation_id.dropna().tolist()
         if str(x).strip()}
    return set(ids), obs


def technical_reservations(path, *, expected_sha=TECHNICAL_RESERVATION_SHA):
    if expected_sha is not None:
        require(sha256(path)==expected_sha,"Prior tier06 genuine technical photo-ID receipt changed")
    photos=set()
    observations=set()
    with Path(path).open(newline="",encoding="utf-8") as f:
        r=csv.DictReader(f)
        require({"photo_id","observation_id","acquired","status"}<=set(r.fieldnames or ()),
                "Missing genuine tier06 technical source-ID reservation columns")
        for row in r:
            pid=str(row["photo_id"]).strip()
            obs=str(row["observation_id"]).strip()
            require(pid.isdecimal() and bool(obs) and pid not in photos,
                    "Duplicate or invalid technical source photo ID")
            photos.add(pid);observations.add(obs)
    require(bool(photos),"Missing technical original-photo source reservation")
    return photos, observations


def build_queue(route_receipt, routed_archive, counts_file,
                old_parquet, old_receipt, technical_csv, output, *,
                first_buckets=FIRST_BUCKETS,
                expected_route_sha=ROUTE_RECEIPT_SHA,
                expected_archive_sha=ROUTED_TAR_SHA,
                expected_technical_sha=TECHNICAL_RESERVATION_SHA,
                expected_prior_receipt_sha=PRIOR_RECEIPT_SHA,
                expected_prior_count=EXPECTED_PRIOR_IDS,
                expected_total=EXPECTED_PHOTOS,
                expected_taxa=EXPECTED_TAXA,
                check_taxa_sha=True, max_new=MAX_NEW_ATTEMPTS):
    output=Path(output)
    require(not output.exists(),"Original photo source cohort is immutable")
    require(1 <= first_buckets <= BUCKETS and max_new > 0,
            "Invalid pre-pixel computing bucket coverage")
    if expected_route_sha is not None:
        require(sha256(route_receipt)==expected_route_sha,"Routed source group receipt SHA changed")
    if expected_archive_sha is not None:
        require(sha256(routed_archive)==expected_archive_sha,"Routed source bucket archive SHA changed")
    info=json.loads(Path(route_receipt).read_text())
    require(info.get("schema")==ROUTE_SCHEMA
            and info.get("source_tier")==SOURCE_TIER
            and info.get("source_group")==SOURCE_GROUP
            and info.get("compute_bucket_count")==BUCKETS
            and info.get("full_original_source_photo_IDs")==expected_total
            and info.get("native_status_not_hard_filtered") is True
            and info.get("images_downloaded")==0
            and info.get("flower_pixels_inferred")==0,
            "Not the SHA-signed original tier06 full-frame metadata routing")
    names=read_taxa(counts_file,check_sha=check_taxa_sha)
    require(len(names)==expected_taxa,"Frozen complete taxon ID/name source missing")
    prior_photos,prior_obs=previous_ids(
        old_parquet,old_receipt,expected_count=expected_prior_count,
        expected_receipt_sha=expected_prior_receipt_sha)
    tech_photos,tech_obs=technical_reservations(
        technical_csv,expected_sha=expected_technical_sha)
    require(not tech_photos & prior_photos,
            "Separately benchmarked tier06 originals already occur in prior checkpoint")
    files={int(x["bucket"]):x for x in info["routed_buckets"]}
    require(len(files)==len(info["routed_buckets"])
            and all(0<=i<BUCKETS for i in files)
            and sum(int(x["original_photo_IDs"]) for x in files.values())
                == info["original_source_photo_IDs_routed"],
            "Routed photo-source bucket manifest inconsistent")
    expected={"routing_receipt.json"}
    expected.update("bucket_csv/"+x["source_file"] for x in files.values())
    selected=[]
    unsupported=[]
    exclusions=Counter()
    photo_ids=set()
    seen_observation_bucket={}
    original_total=0
    bucket_candidates=0
    with tarfile.open(routed_archive,"r:*") as tf:
        members={}
        for member in tf.getmembers():
            if member.isfile():
                key=member.name.removeprefix("./")
                require(key not in members and key in expected,
                        "Unknown/duplicate compressed routed source asset")
                members[key]=member
            else:
                require(member.isdir(),"Unexpected link/special file in routed original source tar")
        require(set(members)==expected,"Routed original source has missing/extra data")
        packaged=tf.extractfile(members["routing_receipt.json"]).read()
        require(sha_bytes(packaged)==sha256(route_receipt),
                "Routed source receipt differs from pinned in-archive copy")
        for bucket,meta in sorted(files.items()):
            name="bucket_csv/"+meta["source_file"]
            payload=tf.extractfile(members[name]).read()
            require(sha_bytes(payload)==meta["source_csv_gzip_sha256"],
                    "Routed original source bucket SHA changed")
            n=0
            with gzip.GzipFile(fileobj=io.BytesIO(payload)) as g:
                with io.TextIOWrapper(g,encoding="utf-8",newline="") as inp:
                    reader=csv.DictReader(inp)
                    require(tuple(reader.fieldnames or ())==COLUMNS,
                            "Original signed source column order changed")
                    for row in reader:
                        n+=1
                        pid=row["photo_id"];obs=row["observation_uuid"]
                        require(pid.isdecimal() and obs and pid not in photo_ids
                                and row["tier"]==str(SOURCE_TIER)
                                and row["inat_taxon_id"] in names,
                                "Duplicate or malformed original source photo ID")
                        photo_ids.add(pid)
                        require(observation_compute_bucket(obs,BUCKETS)==bucket
                                and float(row["observation_inclusion_score"])==obs_draw(obs),
                                "Original observation bucket or source hash differs")
                        require(obs not in seen_observation_bucket
                                or seen_observation_bucket[obs]==bucket,
                                "Original observation split between source compute buckets")
                        seen_observation_bucket[obs]=bucket
                        if bucket>=first_buckets:
                            continue
                        bucket_candidates+=1
                        if pid in prior_photos or obs in prior_obs:
                            exclusions["prior_actual_photo_or_observation"]+=1
                            continue
                        if pid in tech_photos or obs in tech_obs:
                            exclusions["prior_technical_photo_or_observation"]+=1
                            continue
                        ext=normalized_extension(row["image_extension"])
                        if ext is None:
                            exclusions["unsupported_loader_format"]+=1
                            unsupported.append({
                                "photo_id":pid,"observation_uuid":obs,
                                "source_image_extension":row["image_extension"],
                                "reason":"unsupported_original_image_extension_not_flower_status",
                            })
                            continue
                        mid="MOIS-"+hashlib.sha256(
                            ("mois-original-global-RGB-v1|"+pid).encode()).hexdigest()[:24].upper()
                        selected.append({
                            "measurement_id":mid,"photo_id":pid,
                            "inat_taxon_id":row["inat_taxon_id"],
                            "species":names[row["inat_taxon_id"]],
                            "observation_id":obs,"observer_id":row["observer_id"],
                            "observed_on":row["observed_on"],
                            "latitude":row["latitude"],"longitude":row["longitude"],
                            "positional_accuracy":row["positional_accuracy"],
                            "photo_license":row["photo_license"],
                            "photo_url_large":f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.{ext}",
                            "image_filename":mid+".jpg",
                            "source_arm":SOURCE_ARM,
                            # Joint *metadata design* propensity within group00:
                            # 0.0064 tier-width x 64/512 bucket subhash.
                            # Ineligible taxa outside group00 have probability 0.
                            "selection_probability":format(
                                (TIER_ENDS[SOURCE_TIER]-TIER_ENDS[SOURCE_TIER-1])
                                *first_buckets/BUCKETS,".17g"),
                            "source_snapshot_date":"2026-09-27","shard":"",
                        })
            require(n==int(meta["original_photo_IDs"]),"Signed compressed source bucket row count changed")
            original_total+=n
    require(original_total==info["original_source_photo_IDs_routed"],
            "All 33,464 original source photo IDs not accounted for")
    require(bucket_candidates==len(selected)+len(unsupported)+sum(exclusions.values())
            -exclusions["unsupported_loader_format"],
            "All selected original photo IDs must end in one source disposition")
    require(0<len(selected)<=max_new,"New original-photo work queue absent/oversized")
    # Prioritization is observation-based; retain every nonexcluded photo from
    # each accepted observation and do not cap per species or photo number.
    selected.sort(key=lambda x:(
        hashlib.sha256(("mois-tier06-group00-observation-order-v1|"
                         +x["observation_id"]).encode()).hexdigest(),int(x["photo_id"])))
    units=math.ceil(len(selected)/ORIGINAL_BATCH_SIZE)
    output.mkdir(parents=True,exist_ok=False)
    (output/"shards").mkdir()
    sha_by_unit={}
    for unit in range(units):
        members=selected[unit*ORIGINAL_BATCH_SIZE:(unit+1)*ORIGINAL_BATCH_SIZE]
        for row in members:row["shard"]=str(unit)
        path=output/"shards"/f"shard-{unit:03d}.csv"
        with path.open("w",newline="",encoding="utf-8") as f:
            w=csv.DictWriter(f,fieldnames=WORKER_COLUMNS,lineterminator="\n")
            w.writeheader();w.writerows(members)
        sha_by_unit[str(unit)]=sha256(path)
    sidecar=output/"unsupported_selected_original_photo_formats.jsonl"
    sidecar.write_text("".join(json.dumps(x,sort_keys=True)+"\n" for x in unsupported))
    manifest={
        "schema":OUT_SCHEMA,
        "full_frame_photo_ids":expected_total,
        "source_tier":SOURCE_TIER,"source_group":SOURCE_GROUP,
        "source_arm":SOURCE_ARM,
        "full_source_manifest_sha256":info["full_source_registry_sha256"],
        "source_target_registry_sha256":info["full_source_registry_sha256"],
        "source_routing_receipt_sha256":sha256(route_receipt),
        "source_routed_archive_sha256":sha256(routed_archive),
        "prior_original_pixel_checkpoint_sha256":sha256(old_receipt),
        "technical_original_photo_reservation_sha256":sha256(technical_csv),
        "frozen_taxon_names_sha256":sha256(counts_file),
        "full_source_original_tier":SOURCE_TIER,
        "source_tier_probability_lower_bound":TIER_ENDS[SOURCE_TIER-1],
        "selection_probability":(TIER_ENDS[SOURCE_TIER]-TIER_ENDS[SOURCE_TIER-1])
                                *first_buckets/BUCKETS,
        "source_group_restricted_taxon_partitions":True,
        "not_probability_complete_global_source":True,
        "not_a_world_flora_accepted_species_prevalence_sample":True,
        "selection_unit":"observation_uuid",
        "compute_bucket_count":BUCKETS,"admitted_bucket_count":first_buckets,
        "admitted_buckets":list(range(first_buckets)),
        "compute_buckets_are_processing_subsample_not_taxon_caps":True,
        "source_group_original_photo_IDs":original_total,
        "selected_bucket_photo_IDs_before_dedup":bucket_candidates,
        "prior_and_technical_excluded_source_photo_IDs":dict(exclusions),
        "new_photo_attempts_planned":len(selected),
        "selected_new_photos_unsupported_by_frozen_loader":len(unsupported),
        "selected_new_photo_formats_unsupported_sha256":sha256(sidecar),
        "work_units":units,"worker_ids":list(range(units)),
        "group_ids":list(range(math.ceil(units/4))),
        "units_per_group":4,"original_model_batch_size":ORIGINAL_BATCH_SIZE,
        "photo_shard_sha256":sha_by_unit,
        "sampled_before_any_hue_measurement":True,
        "source_original_images_measured":False,"colour_measured":False,
        "technical_unknowns_not_assumed_chromatic":True,
        "no_per_taxon_or_observation_photo_count_cap":True,
        "native_photo_status_not_inferred_or_hard_filtered":True,
        "ecological_analysis_run":False,
        "previous_source_photo_IDs_not_redownloaded":True,
    }
    (output/"manifest.json").write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
    print("MOIS_SIGNED_TIER06_ROUTED_WARM_PRECOLOUR_PHOTO_QUEUE",
          json.dumps({k:manifest[k] for k in (
              "source_group_original_photo_IDs",
              "selected_bucket_photo_IDs_before_dedup",
              "prior_and_technical_excluded_source_photo_IDs",
              "new_photo_attempts_planned","work_units","group_ids"
          )},sort_keys=True),flush=True)
    return manifest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("route-receipt","routed-archive","counts-file",
                "old-parquet","old-receipt","technical-csv","output"):
        p.add_argument("--"+key,type=Path,required=True)
    p.add_argument("--first-buckets",type=int,default=FIRST_BUCKETS)
    a=p.parse_args()
    build_queue(a.route_receipt,a.routed_archive,a.counts_file,
                a.old_parquet,a.old_receipt,a.technical_csv,a.output,
                first_buckets=a.first_buckets)


if __name__=="__main__":
    main()
