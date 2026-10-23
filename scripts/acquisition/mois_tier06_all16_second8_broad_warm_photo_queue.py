#!/usr/bin/env python3
"""Freeze a *second*, non-overlapping 8-of-512 broad original-photo source queue.

All 597,070 tier06 source photo IDs are already routed once into 16 taxonomic
source groups x 512 observation-hash buckets. Select only observation buckets
4,5,...,11: NO photo from the prior all16 first4-bucket original ROI cohort
can be selected again. Verify newest pinned six-arm 887600 baseline, old 4138
warm group00 photo outcomes, prior 4468 all16 photos, technical holdout and
immutable source group SHA before allowing any original photo pixels to open.

This is an actual new original-photo target queue, never botanical flowering
ground truth, a native-only hard filter or global white-flower prevalence.
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
    obs_draw, normalized_extension, sha256, require,
)
from scripts.acquisition.mois_queue_from_full93m_target_tiers import (
    WORKER_COLUMNS, OUT_SCHEMA, read_taxa, observation_compute_bucket,
)
from scripts.acquisition.mois_singlepass_signed_tier_bucket_router import ROUTE_SCHEMA

TIER=6
SOURCE_GROUPS=16
BUCKETS=512
FIRST_BUCKETS=8
SOURCE_ARM="global93m_tier06_all16_buckets004to011_broad_warm_R2"
MAX_NEW=12_000
# GitHub published asset SHA256 digests, source locked BEFORE new colour pixels.
PINNED={
    "00":("c493f29beb2782795bee19d4223c2c026cd59a5a4f479b39bca8594f8bccde6f","c1ff70564bf35c74c4d656511cfd51ad1682bd00b5bf607e9e739011230902f0"),
    "01":("131d67b63c1ec85744c1c4e53ec932deae18df9a7a8841df1d4cc82b83df3c76","e9317629d15d85dc2a0d05e9c08b568e3d68cea35e342482090208a4ddbf1219"),
    "02":("4f25f155a5d1dfb62f801d7778549063645a52013984f0d2cbab3b98c6e99d00","faa7647cdc6d556b95dff929a351d70e32b83fc8340cb47861291c65fe71a715"),
    "03":("623dde9412fc7e8dda80505c48b2688a0d46ca0bcff2b0661e05f41d8d6d3b5a","e0ad6e1bffaa50d240ca8851ba5b0a3ffcec9ad82d7e73a4ca3c2a9fd749c3ca"),
    "04":("5a93d4dfe5f491f0b2122be850d72d4716f3b61b09d40d4df5fe20b48667e776","cee8c23f51907472b915c2b4aad7367718c74abccc09eac72ce96c6a17713142"),
    "05":("a6ea71dc5dfc3a54fc9fe2382f390fdbb63c21ed52a43cb43057a15b7e63af6b","5b78e8a47b16290afc9339e7419ec18a905d5fcd87b025d828165a6a12f6440a"),
    "06":("c6dd61b050933698f1872d44e38de66c19c152d7ddd5f3871d98d112a20748a3","656e4bf63c2dffc2d7c7bb7c3541575dba0593ae4a890be2357f7855d518eead"),
    "07":("a4dd5b1faa814cf4e73a8028bbfd9db933858acdc222a5468565fcf27b162732","7201990b000cd14ea960eda990dcc8113c30b2a887bfff92b6e0735f556a9481"),
    "08":("a82e24d004c9abbb8dea30d1220268652839740a268226b9dfd486eb78297c44","4bcbda7f46cbebbe5dc74dacfb95123c0442958c050ba53e3b83ae0f85377225"),
    "09":("abc591e8271c2d3573b3d2aa551b582831fdf6dd6335dfa1c3d700e1af369646","24d21b46395bfac98a0aeaddb9fad44df62fccf9ded0f5f6f01657fc355a6480"),
    "10":("661866c2dac25df00edc1b2004493500d092c3f52b3b45ff3e526ac51a7f9a30","d659b73d4c047d6219ad7322377aab1fff894bb88bb10475d63b739000a1cd69"),
    "11":("648d8e046f6953b7d1f426b672432f79f1be75b07452814344ef4d263918884a","67e55ce8112f2f9555598db0d4d379dbfb41a68d35d6fa8571dfa7e2c981ef73"),
    "12":("d4f2cbb7be605cba7e3ab6dfdeb17fa3ce60127e8fae12b7864ae88af279f2b2","047d98f07f53dbd86480c17bfa378f028443c589ec923ba412604f28fc7b709b"),
    "13":("823a65d2fd291f7436ced14944aab3c3bdeb506a3dfa957e0de7a7444861fb3a","982453c63d68cf961dc72c87993ff4b7740d2ffbedf6fc099e05a35a91424562"),
    "14":("bec3571ceb2b4138c7554284b9fe41c683aa4f306d01e4371efc5be4db399a98","50f5bc99ebdc8278e1817363cb219d5da8ad854b0319085aa8ed595c2682d328"),
    "15":("288c5ebd3c128d6349261d613238247888037fca5feb87af6c2ec2809d3ade35","beaa7b3c5be829ab93452567bc6a66834349e9f30ec7b4bf7f9c688dc0e344a8"),
}
REFERENCE_MAIN_ATTEMPTED=897_840
REFERENCE_MAIN_CLASSIFIED=404_870
REFERENCE_TIER06_00_ATTEMPTED=4_138
REFERENCE_TIER06_00_CLASSIFIED=1_573
REFERENCE_TIER06_SOURCE_TOTAL=597_070
REFERENCE_PREVIOUS_BROAD_ATTEMPTED=4_468
REFERENCE_PREVIOUS_BROAD_CLASSIFIED=1_760
REFERENCE_TECH_SHA="7038130d6a3cc45100f11a3994da7d6735c554dcb8e6bc829c27f5bef82f2639"


def identity_from_parquet(parquet, receipt, *, expected_attempted, expected_classified, full=False):
    import pandas as pd
    d=json.loads(Path(receipt).read_text())
    if full:
        require(d.get("schema")=="mois_sealed_actual_pixel_checkpoint_v1"
                and d.get("combined_unique_photo_attempts")==expected_attempted
                and d.get("combined_four_hue_classified")==expected_classified
                and sha256(parquet)==d["output_sha256"][Path(parquet).name],
                "Invalid signed main original-photo baseline identity")
    else:
        require(d.get("schema")=="mois_global93m_observation_baseline_real_RGB_photo_results_v1"
                and d.get("attempted")==expected_attempted
                and d.get("classifiable")==expected_classified
                and sha256(parquet)==d["measured_photo_parquet_sha256"],
                "Invalid signed completed tier06 source image identity")
    rows=pd.read_parquet(parquet,columns=["photo_id","observation_id"])
    ids=rows.photo_id.astype(str)
    require(len(rows)==expected_attempted and ids.is_unique
            and ids.str.fullmatch(r"\d+").all(),"Original photos duplicated in previous sources")
    observations={s for s in rows.observation_id.fillna("").astype(str) if s.strip()}
    return set(ids),observations


def technical_reservations(path):
    require(sha256(path)==REFERENCE_TECH_SHA,"Original tier06 real-photo technical reservation SHA changed")
    photos=set();observations=set()
    with Path(path).open(newline="",encoding="utf-8") as stream:
        reader=csv.DictReader(stream)
        require({"photo_id","observation_id","acquired","status"}<=set(reader.fieldnames or ()),
                "Original technical source/observation reservation columns missing")
        for row in reader:
            pid=row["photo_id"].strip();obs=row["observation_id"].strip()
            require(pid.isdecimal() and obs and pid not in photos,
                    "Original technical photo/observation identity invalid")
            photos.add(pid);observations.add(obs)
    require(len(photos)>=16,"Incomplete original technical reservation")
    return photos,observations


def build(registry,routed,main_parquet,main_receipt,group00_parquet,group00_receipt,
          broad_previous_parquet,broad_previous_receipt,
          technical_csv,counts_file,output,*,first_buckets=FIRST_BUCKETS,
          expected_full=EXPECTED_PHOTOS,expected_taxa=EXPECTED_TAXA,
          expected_tier_total=REFERENCE_TIER06_SOURCE_TOTAL,max_new=MAX_NEW,
          pinned=PINNED,check_taxa_sha=True):
    output=Path(output);routed=Path(routed)
    require(not output.exists(),"Cannot change original before-colour source selection")
    require(first_buckets==8 and max_new>=1,
            "Second tier06 broad original-photo wave freezes exactly buckets 4 through 11")
    names=read_taxa(counts_file,check_sha=check_taxa_sha)
    require(len(names)==expected_taxa,"Incomplete signed original plant taxon-name source")
    reg=json.loads(Path(registry).read_text())
    require(reg.get("schema")=="mois_complete_full93m_original_photo_target_registry_v1"
            and reg.get("exact_original_photo_IDs_in_all_disjoint_tiers")==expected_full
            and int(reg.get("new_target_counts_by_tier",{}).get("6",-1))==expected_tier_total
            and reg.get("native_vs_introduced_not_hard_filtered") is True,
            "Missing complete source-photo ID register or original tier06 count")
    source_registry_sha=sha256(registry)
    old_main,old_main_obs=identity_from_parquet(
        main_parquet,main_receipt,expected_attempted=REFERENCE_MAIN_ATTEMPTED,
        expected_classified=REFERENCE_MAIN_CLASSIFIED,full=True)
    group00,group00_obs=identity_from_parquet(
        group00_parquet,group00_receipt,expected_attempted=REFERENCE_TIER06_00_ATTEMPTED,
        expected_classified=REFERENCE_TIER06_00_CLASSIFIED)
    require(not old_main & group00 and not old_main_obs & group00_obs,
            "Already published tier06 group00 observations overlap old main source")
    broad_previous,broad_previous_obs=identity_from_parquet(
        broad_previous_parquet,broad_previous_receipt,
        expected_attempted=REFERENCE_PREVIOUS_BROAD_ATTEMPTED,
        expected_classified=REFERENCE_PREVIOUS_BROAD_CLASSIFIED)
    require(not broad_previous & (old_main | group00)
            and not broad_previous_obs & (old_main_obs | group00_obs),
            "The previous all16 real-original broad photo/source arms overlap historic photo identities")
    tech,tech_obs=technical_reservations(technical_csv)
    require(not tech & (old_main|group00),
            "Technical tier06 photo IDs overlap already published source photos")
    old_photos=old_main|group00|broad_previous
    old_obs=old_main_obs|group00_obs|broad_previous_obs

    group_counts={};source_photos=0;candidate_total=0
    chosen=[];unsupported=[];skipped=Counter();ids=set();seen_obs_group={}
    for n in range(SOURCE_GROUPS):
        group=f"{n:02d}";root=routed/f"group{group}"
        receipt=root/"routing_receipt.json"
        archive=root/f"tier06-group{group}-512-source-buckets.tar"
        expected_receipt,expected_archive=pinned[group]
        require(sha256(receipt)==expected_receipt and sha256(archive)==expected_archive,
                "Immutable signed tier06 group source Release digest mismatch: "+group)
        d=json.loads(receipt.read_text())
        require(d.get("schema")==ROUTE_SCHEMA and d.get("source_tier")==TIER
                and d.get("source_group")==n and d.get("compute_bucket_count")==BUCKETS
                and d.get("full_original_source_photo_IDs")==expected_full
                and d.get("full_source_registry_sha256")==source_registry_sha
                and d.get("source_group_read_exactly_once") is True
                and d.get("images_downloaded")==0 and d.get("flower_pixels_inferred")==0
                and d.get("native_status_not_hard_filtered") is True,
                "Source-group routing receipt is not the frozen full original photo frame")
        buckets={int(b["bucket"]):b for b in d["routed_buckets"]}
        require(len(buckets)==len(d["routed_buckets"])
                and sum(int(b["original_photo_IDs"]) for b in buckets.values())
                    ==int(d["original_source_photo_IDs_routed"]),
                "Signed source group contains lost/duplicate original bucket identities")
        group_counts[group]=int(d["original_source_photo_IDs_routed"])
        source_photos+=group_counts[group]
        with tarfile.open(archive,"r:") as t:
            members={}
            for member in t.getmembers():
                if member.isfile():
                    key=member.name.removeprefix("./")
                    require(key not in members and
                            (key=="routing_receipt.json" or
                             key in {"bucket_csv/"+b["source_file"] for b in buckets.values()}),
                            "Unrecognized/duplicate original routed photo source tar member")
                    members[key]=member
                else:
                    require(member.isdir(),"Non-source file/symlink in immutable photo metadata tar")
            require(len(members)==1+len(buckets)
                    and hashlib.sha256(t.extractfile(members["routing_receipt.json"]).read()).hexdigest()
                       ==expected_receipt,
                    "Archive and released routed source manifest differ")
            for index in range(4,4+first_buckets):
                b=buckets.get(index)
                if b is None:
                    continue
                member=members["bucket_csv/"+b["source_file"]]
                data=t.extractfile(member).read()
                require(hashlib.sha256(data).hexdigest()==b["source_csv_gzip_sha256"],
                        "Original source bucket photo metadata SHA changed")
                rows=0
                with gzip.GzipFile(fileobj=io.BytesIO(data),mode="rb") as stream:
                    with io.TextIOWrapper(stream,encoding="utf-8",newline="") as f:
                        csvreader=csv.DictReader(f)
                        require(tuple(csvreader.fieldnames or ())==COLUMNS,
                                "Original plant photo metadata field order changed")
                        for row in csvreader:
                            rows+=1
                            pid=row["photo_id"];obs=row["observation_uuid"]
                            require(pid.isdecimal() and obs and pid not in ids
                                    and row["inat_taxon_id"] in names
                                    and row["tier"]==str(TIER),
                                    "Original photo ID/taxon/observation lost or duplicated")
                            ids.add(pid)
                            require(observation_compute_bucket(obs,BUCKETS)==index
                                    and float(row["observation_inclusion_score"])==obs_draw(obs),
                                    "Observation original hash changed or wrong source bucket")
                            require(obs not in seen_obs_group or seen_obs_group[obs]==n,
                                    "Same original observation split among source taxon groups")
                            seen_obs_group[obs]=n
                            candidate_total+=1
                            if pid in old_photos or obs in old_obs:
                                skipped["previously_finished_photo_or_observation"]+=1
                                continue
                            if pid in tech or obs in tech_obs:
                                skipped["technical_source_photo_or_observation"]+=1
                                continue
                            ext=normalized_extension(row["image_extension"])
                            if ext is None:
                                skipped["unsupported_original_image_extension"]+=1
                                unsupported.append({
                                    "photo_id":pid,"observation_uuid":obs,
                                    "source_image_extension":row["image_extension"],
                                    "reason":"frozen_loader_unsupported_metadata_not_flower_absence",
                                })
                                continue
                            mid="MOIS-"+hashlib.sha256(
                                ("mois-original-global-RGB-v1|"+pid).encode()).hexdigest()[:24].upper()
                            chosen.append({
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
                                "selection_probability":format(
                                    (TIER_ENDS[TIER]-TIER_ENDS[TIER-1])
                                    *first_buckets/BUCKETS,".17g"),
                                "source_snapshot_date":"2026-09-27",
                                "shard":"",
                            })
                require(rows==int(b["original_photo_IDs"]),
                        "Original bucket source row count was tampered")
    require(source_photos==expected_tier_total,
            "All 16 routed source-photo groups do not sum to complete signed tier06 frame")
    require(candidate_total==len(chosen)+sum(skipped.values())
            and 0<len(chosen)<=max_new,
            "Pre-colour source-photo IDs lost/duplicated or new batch exceeded limit")

    chosen.sort(key=lambda row:(
        hashlib.sha256(("mois-broad-tier06-observation-shard-v1|"+row["observation_id"]).encode()).hexdigest(),
        int(row["photo_id"])))
    units=math.ceil(len(chosen)/256)
    output.mkdir(parents=True)
    (output/"shards").mkdir()
    shard_sha={}
    for unit in range(units):
        members=chosen[unit*256:(unit+1)*256]
        for row in members:row["shard"]=str(unit)
        f=output/"shards"/f"shard-{unit:03d}.csv"
        with f.open("w",encoding="utf-8",newline="") as stream:
            writer=csv.DictWriter(stream,fieldnames=WORKER_COLUMNS,lineterminator="\n")
            writer.writeheader();writer.writerows(members)
        shard_sha[str(unit)]=sha256(f)
    sidecar=output/"unsupported_selected_original_photo_formats.jsonl"
    sidecar.write_text("".join(json.dumps(x,sort_keys=True)+"\n" for x in unsupported))
    receipt={
        "schema":OUT_SCHEMA,
        "full_frame_photo_ids":expected_full,
        "source_tier":TIER,"source_group":"all16",
        "source_arm":SOURCE_ARM,
        "full_source_original_tier":TIER,
        "full_source_manifest_sha256":source_registry_sha,
        "source_target_registry_sha256":source_registry_sha,
        "complete_tier06_source_photo_IDs":source_photos,
        "all16_immutable_signed_source_group_counts":group_counts,
        "source_original_images_measured":False,
        "source_tier_probability_lower_bound":TIER_ENDS[TIER-1],
        "selection_probability":(TIER_ENDS[TIER]-TIER_ENDS[TIER-1])
                                  *first_buckets/BUCKETS,
        "selection_unit":"observation_uuid",
        "compute_bucket_count":BUCKETS,"admitted_bucket_count":first_buckets,
        "admitted_bucket_start":4,
        "admitted_buckets":list(range(4,4+first_buckets)),
        "previous_broad_source_parquet_sha256":sha256(broad_previous_parquet),
        "previous_broad_original_photo_attempts":len(broad_previous),
        "previous_broad_observation_overlap_with_new_bucket_range_forbidden":True,
        "metadata_original_photos_in_admitted_buckets":candidate_total,
        "previously_measured_photo_ids_skipped":skipped["previously_finished_photo_or_observation"],
        "previous_photo_or_observation_exclusions":dict(skipped),
        "selected_new_photos_unsupported_by_frozen_loader":len(unsupported),
        "selected_new_photo_formats_unsupported_sha256":sha256(sidecar),
        "new_photo_attempts_planned":len(chosen),
        "work_units":units,"worker_ids":list(range(units)),
        "group_ids":list(range(math.ceil(units/4))),
        "units_per_group":4,"original_model_batch_size":256,
        "photo_shard_sha256":shard_sha,
        "sampled_before_any_hue_measurement":True,
        "colour_measured":False,
        "no_per_taxon_or_observation_photo_count_cap":True,
        "native_photo_status_not_inferred_or_hard_filtered":True,
        "technical_unknowns_not_assumed_chromatic":True,
        "prior_photo_and_observation_exclusion_verified":True,
        "not_probability_complete_global_source":True,
        "not_a_world_flora_accepted_species_prevalence_sample":True,
        "ecological_analysis_run":False,
    }
    (output/"manifest.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print("MOIS_BROAD_ALL16_TIER06_FROZEN_ORIGINAL_RGB_SOURCE_READY",
          json.dumps({k:receipt[k] for k in (
              "complete_tier06_source_photo_IDs",
              "metadata_original_photos_in_admitted_buckets",
              "previous_photo_or_observation_exclusions",
              "new_photo_attempts_planned","work_units","group_ids")
          },sort_keys=True),flush=True)
    return receipt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("registry","routed","main-parquet","main-receipt",
                "group00-parquet","group00-receipt",
                "broad-previous-parquet","broad-previous-receipt",
                "technical-csv","counts-file","output"):
        p.add_argument("--"+key,type=Path,required=True)
    p.add_argument("--first-buckets",type=int,default=FIRST_BUCKETS)
    a=p.parse_args()
    build(a.registry,a.routed,a.main_parquet,a.main_receipt,
          a.group00_parquet,a.group00_receipt,
          a.broad_previous_parquet,a.broad_previous_receipt,
          a.technical_csv,a.counts_file,a.output,first_buckets=a.first_buckets)


if __name__=="__main__":
    main()
