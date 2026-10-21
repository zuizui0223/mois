#!/usr/bin/env python3
"""Route each signed source tier group into observation-stable compute buckets ONCE.

Without this router, opening 128 signed source photo parts separately for every
one of 512 compute buckets scales as a repeated scan of the full tier.
Here the source group archives are read once, bucket files are compressed on
disk with bounded open handles, and later workers need only their bucket files.

This is source metadata only: no photos downloaded, no early native filtering,
no flower-stage decisions, no ROI/hue estimation or changed selection probability.
"""
from __future__ import annotations

import argparse
from collections import Counter, OrderedDict
import csv
import gzip
import json
import os
from pathlib import Path
import re
import shutil
import tarfile
import tempfile

from scripts.acquisition.mois_onepass_full93m_target_tiers import (
    COLUMNS, EXPECTED_PHOTOS, EXPECTED_TAXA, RESULT_SCHEMA, TIER_ENDS,
    obs_draw, require, sha256,
)
from scripts.acquisition.mois_queue_from_full93m_target_tiers import observation_compute_bucket


ROUTE_SCHEMA = "mois_93m_singlepass_signed_source_group_bucket_routing_v1"
SOURCE_GROUP = re.compile(r"^tier-(\d{2})-group-(\d{2})\.tar$")


def route_signed_group(registry, receipt, archive, output, *,
                       bucket_count=512, max_open=32,
                       expected_total=EXPECTED_PHOTOS,
                       expected_taxa=EXPECTED_TAXA,
                       parts_per_group=8):
    registry, receipt, archive, output = map(Path, (registry, receipt, archive, output))
    require(not output.exists(), "Original signed bucket output may not be overwritten")
    require(1 <= bucket_count <= 1024 and 1 <= max_open <= 64,
            "Invalid bounded compute bucket routing dimensions")
    name = SOURCE_GROUP.fullmatch(archive.name)
    require(name is not None, "Unsupported source tier group archive name")
    tier, group_no = map(int, name.groups())
    require(0 <= tier < len(TIER_ENDS),
            "Source bucket router received an unknown observation hash tier")
    full = json.loads(registry.read_text(encoding="utf-8"))
    group = json.loads(receipt.read_text(encoding="utf-8"))
    require(full.get("schema") == "mois_complete_full93m_original_photo_target_registry_v1"
            and full.get("exact_original_photo_IDs_in_all_disjoint_tiers") == expected_total
            and full.get("original_named_taxa_from_full_source") == expected_taxa
            and full.get("source_only_NO_photo_pixels_or_hue_measured") is True
            and full.get("native_vs_introduced_not_hard_filtered") is True,
            "Not the immutable complete original photo ID opportunity frame")
    require(group.get("schema") == RESULT_SCHEMA
            and group.get("source_snapshot_SHA256") == full["full_source_manifest_sha256"]
            and group.get("all_source_photo_IDs_retained_in_one_and_only_one_tier") is True
            and group.get("source_original_photo_IDs_uniqueness_proved_by_upstream_primary_key") is True
            and group.get("same_observation_uses_same_selection_score") is True
            and group.get("images_downloaded") == 0
            and group.get("colour_pixels_measured") == 0,
            "One-pass source group identity or uniqueness proof missing")
    require(archive.name in full["published_tier_archive_SHA256"]
            and sha256(archive) == full["published_tier_archive_SHA256"][archive.name],
            "Frozen source tier archive SHA changed")
    parts = {int(x["part"]):x for x in group["parts"]}
    expected_parts = set(range(group_no*parts_per_group,(group_no+1)*parts_per_group))
    require(set(parts) == expected_parts,
            "Unrecognized or missing fixed original photo source partitions")

    output.parent.mkdir(parents=True,exist_ok=True)
    directory=Path(tempfile.mkdtemp(prefix="mois-signed-route-",dir=output.parent))
    open_files=OrderedDict()
    counts=Counter()
    first_opened=set()
    total=0
    try:
        # No extraction of raw photos: this tar contains only source CSV metadata.
        source_dir=directory/"source"
        source_dir.mkdir()
        expected_names={
            f"tier-{tier:02d}/part-{number:03d}.csv.gz"
            for number in expected_parts
        }
        with tarfile.open(archive) as handle:
            entries=handle.getmembers()
            files={m.name for m in entries if m.isfile()}
            require(files == expected_names
                    and all(m.isfile() or m.isdir() for m in entries)
                    and len(files)==parts_per_group,
                    "Source tar does not match exactly the original signed tier")
            handle.extractall(source_dir, filter="data")
        target=directory/"bucket_csv"
        target.mkdir()

        def write_bucket(number,row):
            if number in open_files:
                f,writer=open_files[number]
                open_files.move_to_end(number)
            else:
                if len(open_files) >= max_open:
                    _,(old_file,_) = open_files.popitem(last=False)
                    old_file.close()
                file=target/f"bucket-{number:04d}.csv.gz"
                first=(number not in first_opened)
                f=gzip.open(file,"wt" if first else "at",encoding="utf-8",newline="")
                writer=csv.DictWriter(f,fieldnames=COLUMNS,lineterminator="\n")
                if first:
                    writer.writeheader()
                    first_opened.add(number)
                open_files[number]=(f,writer)
            writer.writerow(row)
            counts[number]+=1

        for number,source_part in sorted(parts.items()):
            path=source_dir/f"tier-{tier:02d}/part-{number:03d}.csv.gz"
            require(sha256(path)==source_part["ranked_target_csv_sha256"][str(tier)],
                    "An individual original source-photo part SHA changed")
            rows=0
            with gzip.open(path,"rt",encoding="utf-8",newline="") as f:
                reader=csv.DictReader(f)
                require(tuple(reader.fieldnames or ())==COLUMNS,
                        "Source photo column provenance changed")
                for row in reader:
                    rows+=1
                    pid=row["photo_id"]
                    uuid=row["observation_uuid"]
                    require(pid.isdecimal() and bool(uuid) and row["tier"]==str(tier),
                            "Source photo identity or disjoint tier assignment missing")
                    require(float(row["observation_inclusion_score"])==obs_draw(uuid),
                            "Original observation inclusion hash differs from source registry")
                    compute=observation_compute_bucket(uuid,bucket_count)
                    write_bucket(compute,row)
            require(rows == int(source_part["counts_by_disjoint_tier"][str(tier)]),
                    "Signed source-photo partition count changed")
            total+=rows
        for f,_ in open_files.values():
            f.close()
        open_files.clear()
        require(total == int(group["counts_by_disjoint_tier"][str(tier)])
                and sum(counts.values()) == total,
                "Original source group photos were lost while routing")
        # Source-only input has already proved global photo ID uniqueness through
        # its SQLite primary key and signed part SHA; no full-93m Python set.
        bucket_receipts=[]
        for number in sorted(counts):
            path=target/f"bucket-{number:04d}.csv.gz"
            bucket_receipts.append({
                "bucket":number,
                "original_photo_IDs":int(counts[number]),
                "source_csv_gzip_sha256":sha256(path),
                "source_file":path.name,
            })
        report={
            "schema":ROUTE_SCHEMA,
            "source_tier":tier,
            "source_group":group_no,
            "source_tier_probability_interval":[
                0 if tier==0 else TIER_ENDS[tier-1],
                TIER_ENDS[tier],
            ],
            "full_original_source_photo_IDs":expected_total,
            "full_source_registry_sha256":sha256(registry),
            "group_source_receipt_sha256":sha256(receipt),
            "group_source_archive_sha256":sha256(archive),
            "compute_bucket_count":bucket_count,
            "maximum_open_compressed_source_files":max_open,
            "original_source_photo_IDs_routed":total,
            "routed_buckets":bucket_receipts,
            "observation_compute_hash_unchanged":True,
            "original_photo_ID_uniqueness_inherited_from_pinned_upstream_primary_key":True,
            "source_photo_and_observation_grain_preserved":True,
            "native_status_not_hard_filtered":True,
            "images_downloaded":0,
            "flower_pixels_inferred":0,
            "source_group_read_exactly_once":True,
            "no_repeated_source_tier_scan_per_bucket":True,
        }
        (directory/"routing_receipt.json").write_text(
            json.dumps(report,indent=2,sort_keys=True)+"\n")
        # Source files are already verified by signed receipt and need not be
        # duplicated in the *result* archive.
        shutil.rmtree(source_dir)
        os.replace(directory,output)
        print("MOIS_SIGNED_93M_ONEPASS_BUCKET_ROUTING",
              json.dumps({"tier":tier,"group":group_no,"photos":total,
                          "buckets":len(bucket_receipts)},sort_keys=True),flush=True)
        return report
    except BaseException:
        for f,_ in open_files.values():
            f.close()
        shutil.rmtree(directory,ignore_errors=True)
        raise


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("registry","receipt","archive","output"):
        p.add_argument("--"+k,type=Path,required=True)
    p.add_argument("--bucket-count",type=int,default=512)
    p.add_argument("--max-open",type=int,default=32)
    a=p.parse_args()
    route_signed_group(a.registry,a.receipt,a.archive,a.output,
                       bucket_count=a.bucket_count,max_open=a.max_open)


if __name__=="__main__":
    main()
