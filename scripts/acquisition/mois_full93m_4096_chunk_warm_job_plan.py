#!/usr/bin/env python3
"""Schedule the exact 93,209,351 signed iNaturalist source photos as bounded jobs.

The existing 128×14 source-gzip manifest is SHA signed. Preserve every
original photo ID, observation and species without any colour-based prefilter.
A *planned* photo chunk has no pixel data or child SHA until actually split.
Partitioned execution: <=4096 genuine source photo IDs per bounded chunk;
up to 8 sequential original chunks per persistent warm YOLO worker.

No S3/image API calls, no original hue inference, no invented child hashes.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import os
import shutil
import tempfile

SCHEMA="mois_full93m_signed_source_4096_photo_chunks_8chunk_warm_worker_plan_v1"
EXPECTED_SOURCE_SHA="6b607a4b7b0863d65a0eb514d4158434453e7ae4d909f4cf1e3372b424eb9279"
PHOTO_COUNT=93_209_351
PARTS=128
TIERS=14
COLUMNS=("part","tier","relative_source_gzip_path","source_gzip_sha256",
         "expected_photo_rows","photo_level_sampling_probability")
CHUNK_FIELDS=("tier","part","chunk_index","relative_source_gzip_path",
 "source_gzip_sha256","source_row_offset_inclusive","planned_original_photo_ids",
 "child_chunk_sha256_available_only_after_real_source_split")
BATCH_FIELDS=("tier","part","warm_batch_index","first_chunk_index",
              "child_chunks_in_batch","expected_original_photo_ids",
              "source_gzip_sha256","one_model_bootstrap_target",
              "actual_source_chunk_shas_required_before_image_work")

def require(x,msg):
    if not x: raise ValueError("MOIS_FULL93M_WARM_PLAN: "+msg)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(4 << 20),b""):h.update(b)
    return h.hexdigest()

def run(source,out,*,expected_sha=EXPECTED_SOURCE_SHA,
        photo_count=PHOTO_COUNT,parts=PARTS,tiers=TIERS,
        chunk_size=4096,chunks_per_warm=8):
    source,out=Path(source),Path(out)
    require(not out.exists(),"Never overwrite a previously signed workload plan")
    require(1<=chunk_size<=8192 and 1<=chunks_per_warm<=8,
            "Unsafe/unbounded source photo chunk width or warm batch size")
    require(sha(source)==expected_sha,"Complete 1792 source photo work-plan SHA mismatch")
    with source.open(newline="",encoding="utf-8") as inp:
        reader=csv.DictReader(inp)
        require(tuple(reader.fieldnames or ())==COLUMNS,"Wrong original photo plan grain")
        source_jobs=list(reader)
    require(len(source_jobs)==parts*tiers, "Source part×tier photo plan missing")
    source_jobs.sort(key=lambda r:(int(r["tier"]),int(r["part"])))
    require({(int(r["tier"]),int(r["part"])) for r in source_jobs}
            == {(t,p) for t in range(tiers) for p in range(parts)},
            "Some signed original source part×tier IDs are duplicate or missing")
    src_count=0
    out.parent.mkdir(parents=True,exist_ok=True)
    stage=Path(tempfile.mkdtemp(prefix="mois-full93m-chunk-all-",dir=out.parent))
    chunk_path=stage/"all_photo_source_4096_chunk_tasks.csv"
    warm_path=stage/"all_one_warm_model_8chunk_batches.csv"
    try:
        chunks=batches=0
        by_tier={str(t):{"photo_ids":0,"chunks":0,"batches":0}
                 for t in range(tiers)}
        with (chunk_path.open("w",newline="",encoding="utf-8") as fp,
              warm_path.open("w",newline="",encoding="utf-8") as fw):
            wc=csv.DictWriter(fp,fieldnames=CHUNK_FIELDS,lineterminator="\n")
            wb=csv.DictWriter(fw,fieldnames=BATCH_FIELDS,lineterminator="\n")
            wc.writeheader();wb.writeheader()
            for r in source_jobs:
                t,p=int(r["tier"]),int(r["part"])
                n=int(r["expected_photo_rows"])
                require(n>=0 and r["relative_source_gzip_path"]==
                        "tier-%02d/part-%03d.csv.gz"%(t,p),
                        "Photo source path/count not SHA-pinned")
                digest=r["source_gzip_sha256"]
                require(len(digest)==64 and all(c in "0123456789abcdef" for c in digest),
                        "Source original photo part SHA invalid")
                require(r["photo_level_sampling_probability"]=="1",
                        "Do not silently drop original source photographs")
                src_count+=n
                by_tier[str(t)]["photo_ids"]+=n
                num_chunks=math.ceil(n/chunk_size)
                for k in range(num_chunks):
                    width=min(chunk_size,n-k*chunk_size)
                    require(width>0,"Zero-sized actual photo chunk")
                    wc.writerow({
                        "tier":t,"part":p,"chunk_index":k,
                        "relative_source_gzip_path":r["relative_source_gzip_path"],
                        "source_gzip_sha256":digest,
                        "source_row_offset_inclusive":k*chunk_size,
                        "planned_original_photo_ids":width,
                        "child_chunk_sha256_available_only_after_real_source_split":"true",
                    })
                    chunks+=1;by_tier[str(t)]["chunks"]+=1
                for batch in range(math.ceil(num_chunks/chunks_per_warm)):
                    first=batch*chunks_per_warm
                    last=min(num_chunks,first+chunks_per_warm)
                    width=min(n,chunk_size*last)-chunk_size*first
                    require(0<width<=chunks_per_warm*chunk_size,
                            "Malformed warm batch source photo count")
                    wb.writerow({
                        "tier":t,"part":p,"warm_batch_index":batch,
                        "first_chunk_index":first,
                        "child_chunks_in_batch":last-first,
                        "expected_original_photo_ids":width,
                        "source_gzip_sha256":digest,
                        "one_model_bootstrap_target":"true",
                        "actual_source_chunk_shas_required_before_image_work":"true",
                    })
                    batches+=1;by_tier[str(t)]["batches"]+=1
        require(src_count==photo_count,"Whole original 93m photo ID counts changed")
        result={
            "schema":SCHEMA,
            "original_signed_1792_job_source_sha256":expected_sha,
            "original_source_photo_IDs_entire_universe":src_count,
            "original_source_part_tier_jobs":parts*tiers,
            "max_original_photos_per_child_chunk":chunk_size,
            "max_consecutive_chunks_per_warm_worker":chunks_per_warm,
            "planned_child_photo_chunk_tasks":chunks,
            "planned_one_warm_model_worker_batches":batches,
            "by_tier":by_tier,
            "chunk_task_csv_sha256":sha(chunk_path),
            "warm_batch_csv_sha256":sha(warm_path),
            "each_source_part_id_retained_once_per_tier":True,
            "original_photo_IDs_preselected_by_hue_or_flowering":False,
            "output_is_plan_only_no_child_sha_or_image_bytes_yet":True,
            "images_downloaded":0,
            "flower_colour_classes_inferred":0,
            "actual_photo_scans_completed_by_this_plan":0,
        }
        (stage/"all93m_chunk_workplan_receipt.json").write_text(
            json.dumps(result,indent=2,sort_keys=True)+"\n")
        os.replace(stage,out)
        return result
    except BaseException:
        shutil.rmtree(stage,ignore_errors=True)
        raise

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original-1792-jobs",type=Path,required=True)
    p.add_argument("--source-sha256",default=EXPECTED_SOURCE_SHA)
    p.add_argument("--chunk-size",type=int,default=4096)
    p.add_argument("--chunks-per-warm",type=int,default=8)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    r=run(a.original_1792_jobs,a.output,expected_sha=a.source_sha256,
          chunk_size=a.chunk_size,chunks_per_warm=a.chunks_per_warm)
    print("MOIS_FULL93M_EXHAUSTIVE_WARM_PHOTO_WORKER_PLAN",json.dumps(r,sort_keys=True))

if __name__=="__main__":
    main()
