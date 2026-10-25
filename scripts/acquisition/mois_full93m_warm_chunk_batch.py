#!/usr/bin/env python3
"""Process consecutive 93m signed photo chunks under ONE warm YOLO model.

A chunk is 4096 original photo IDs (at most), not a source part or global
sample. This avoids a cold YOLO restart for each of ~24k prospective chunks.
No independent model/ROI output becomes an original 4-hue class. For large
jobs use durable persistent workers near the Open Data source and explicit
cluster-wide provider traffic control; this script has a per-process cap.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path

from scripts.acquisition.mois_full93m_signed_chunk_pixel_scan import (
    run as process_one, sha, CHUNK_SCHEMA
)
from scripts.acquisition.mois_full93m_medium_pixel_scan import (
    Pacer, WarmFlowerBoxDetector, download_medium
)

SCHEMA = "mois_original_full93m_signed_multi_chunk_one_warm_model_v1"
def require(x,why):
    if not x:
        raise ValueError("MOIS_FULL93M_BATCH_WARM: "+why)

def run(chunks,output,*,first_chunk,chunk_count,
        detector,download,max_real_photo_ids=32768,
        workers=8,batch_size=512):
    chunks,output=Path(chunks),Path(output)
    require(0<=first_chunk and 1<=chunk_count<=8
            and 1<=max_real_photo_ids<=32768
            and 1<=workers<=12 and 1<=batch_size<=1024,
            "Bounded worker photo/chunk resource ceiling missing")
    source=json.loads((chunks/"chunk_receipt.json").read_text())
    require(source.get("schema")==CHUNK_SCHEMA and
            source.get("source_fully_partitioned_once_into_nonoverlapping_photo_chains") is True,
            "Signed original part chunk source unavailable")
    manifest=chunks/"chunk_manifest.csv"
    require(sha(manifest)==source["manifest_sha256"],"Source chunk manifest SHA changed")
    with manifest.open(newline="",encoding="utf-8") as fh:
        rows=list(csv.DictReader(fh))
    require(len(rows)==source["chunks"]
            and first_chunk+chunk_count<=len(rows),
            "Requested worker chunks exceed genuine source partition")
    selected=rows[first_chunk:first_chunk+chunk_count]
    total=sum(int(x["expected_photo_ids"]) for x in selected)
    require(0<total<=max_real_photo_ids,
            "Chunk batch would exceed explicitly bounded actual original photo count")
    require(all(int(row["chunk_index"])==first_chunk+i
                and int(row["part"])==source["original_part"]
                and int(row["tier"])==source["original_tier"]
                and row["original_source_part_sha256"]==
                    source["original_source_part_sha256"]
                for i,row in enumerate(selected)),
            "Chunk rows do not belong to one SHA-verified original photo parent")
    output.mkdir(parents=True,exist_ok=True)
    completed=[]
    for row in selected:
        idx=int(row["chunk_index"])
        folder=output/f"chunk-{idx:05d}"
        proof=process_one(chunks,folder,idx,detector=detector,
                          download=download,workers=workers,
                          batch_size=batch_size)
        require(proof["source_chunk_all_photo_ids_scanned"] is True
                and proof["actual_new_or_resumed_medium_photo_attempts"]==
                    int(row["expected_photo_ids"])
                and proof["entire_parent_original_photo_part_scanned"] is False,
                "Only all original photo IDs in each selected child may be counted")
        completed.append({
            "chunk_index":idx,
            "original_photo_ids":int(row["expected_photo_ids"]),
            "actual_medium_pixels_decoded":proof["medium_decoded_photo_count"],
            "source_child_pixel_receipt_sha256":sha(folder/"complete_receipt.json"),
            "pixel_status_counts":proof["screen_counts"],
        })
    require(sum(x["original_photo_ids"] for x in completed)==total,
            "Original image chunks no longer sum to the scheduled source total")
    digest = sha(manifest)
    summary={
        "schema":SCHEMA,
        "original_source_photo_part_sha256":source["original_source_part_sha256"],
        "original_source_photo_part_total_metadata_IDs":source["actual_original_photo_IDs"],
        "original_photo_part_complete":False,
        "source_93m_all_pixels_complete":False,
        "chunk_manifest_sha256":digest,
        "source_chunk_start":first_chunk,
        "source_chunks_in_batch":chunk_count,
        "original_source_part_total_chunks":source["chunks"],
        "photo_IDs_fully_pixel_attempted_in_this_batch":total,
        "new_model_bootstraps_in_this_batch":1 if isinstance(detector,WarmFlowerBoxDetector) else 0,
        "original_source_four_hue_ROI_outcomes_unchanged":True,
        "completed_chunk_receipts":completed,
    }
    file=output/(f"warm_batch_{first_chunk:05d}_through_{first_chunk+chunk_count-1:05d}.json")
    if file.exists():
        require(json.loads(file.read_text())==summary,
                "Repeat batch receipt changed")
    else:
        file.write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    return summary

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--chunks",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--first-chunk",type=int,required=True)
    p.add_argument("--chunk-count",type=int,required=True)
    p.add_argument("--max-photo-ids",type=int,default=32768)
    p.add_argument("--weights",type=Path,required=True)
    p.add_argument("--weights-sha256",required=True)
    p.add_argument("--device",default="cpu")
    p.add_argument("--requests-per-second",type=float,default=2.5)
    p.add_argument("--workers",type=int,default=8)
    p.add_argument("--batch-size",type=int,default=512)
    args=p.parse_args()
    require(0<args.requests_per_second<=3.0,
            "Source provider pace must be explicitly approved before increasing")
    detector=WarmFlowerBoxDetector(args.weights,args.weights_sha256,
                                   0,"item",device=args.device,conf=.02)
    pace=Pacer(args.requests_per_second)
    result=run(args.chunks,args.output,first_chunk=args.first_chunk,
               chunk_count=args.chunk_count,detector=detector,
               download=lambda url:download_medium(url,pace),
               max_real_photo_ids=args.max_photo_ids,
               workers=args.workers,batch_size=args.batch_size)
    print("MOIS_REAL_93M_SOURCE_ONE_WARM_MULTICHUNK_PIXELS",
          json.dumps(result,sort_keys=True))

if __name__=="__main__":
    main()
