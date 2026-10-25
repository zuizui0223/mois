#!/usr/bin/env python3
"""Read actual medium pixels from ONE bounded chunk of a signed MOIS 93m source part.

Exact source-part -> signed 4096-photo chunk -> downloaded 500px pixels.
Keeps parent photo part INCOMPLETE until all its child chunks are completed.
One warm YOLO bootstrap can be shared by repeated run(...) calls in a worker.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import json
from pathlib import Path

from scripts.acquisition.mois_full93m_source_part_chunks import (
    sha, SCHEMA as CHUNK_SCHEMA
)
from scripts.acquisition.mois_full93m_medium_pixel_scan import (
    Pacer, WarmFlowerBoxDetector, download_medium, scan_source
)

SCHEMA = "mois_signed_full93m_parent_photo_part_bounded_actual_pixel_chunk_v1"

def require(cond, msg):
    if not cond:
        raise ValueError("MOIS_93M_REAL_PIXEL_CHUNK: " + msg)

def run(chunks, output, chunk_index, *, detector,
        download, workers=8, batch_size=512, max_packets=None):
    chunks, output = Path(chunks), Path(output)
    report = json.loads((chunks / "chunk_receipt.json").read_text())
    require(report.get("schema") == CHUNK_SCHEMA
            and report.get("source_fully_partitioned_once_into_nonoverlapping_photo_chains") is True
            and report.get("floral_hue_or_bloom_inferred") is False
            and report.get("image_pixels_downloaded") == 0,
            "Only source-verified complete original photo chunk roster is accepted")
    manifest = chunks / "chunk_manifest.csv"
    require(sha(manifest) == report["manifest_sha256"],
            "Bounded source chunk manifest SHA differs")
    with manifest.open(newline="", encoding="utf-8") as f:
        table = list(csv.DictReader(f))
    require(len(table) == report["chunks"]
            and set(int(x["chunk_index"]) for x in table) == set(range(report["chunks"])),
            "Missing or duplicate source photo chunks")
    require(isinstance(chunk_index,int) and 0 <= chunk_index < len(table),
            "Requested non-existent source photo chunk")
    row = table[chunk_index]
    tier = int(row["tier"]); part = int(row["part"])
    parent = report["original_source_part_sha256"]
    require(part == report["original_part"]
            and tier == report["original_tier"]
            and row["original_source_part_sha256"] == parent
            and row["chunk_index"] == str(chunk_index)
            and 1 <= int(row["expected_photo_ids"]) <= report["chunk_size_maximum"]
            and row["chunk_relative_path"] ==
                "chunks/chunk-%05d.csv.gz" % chunk_index,
            "Source chunk not tied to original signed photo/observation part")
    source = chunks / row["chunk_relative_path"]
    require(source.is_file() and sha(source) == row["chunk_gzip_sha256"],
            "Source photo-ID chunk file not SHA verified")
    with gzip.open(source, "rt", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        photo_ids = [x["photo_id"] for x in reader]
    require(len(photo_ids) == int(row["expected_photo_ids"])
            and len(set(photo_ids)) == len(photo_ids)
            and photo_ids[0] == row["first_photo_id"]
            and photo_ids[-1] == row["last_photo_id"],
            "Source chunk photo IDs/count/first-last changed")
    result = scan_source(
        source, output,
        expected_source_sha256=row["chunk_gzip_sha256"],
        expected_count=int(row["expected_photo_ids"]),
        part=part, tier=tier,
        source_role="full93m_chunk",
        parent_source_sha256=parent,
        chunk_index=chunk_index,
        detector=detector, download=download,
        workers=workers, batch_size=batch_size, max_packets=max_packets,
    )
    complete = bool(result["all_original_input_photos_attempted"])
    require(result["source_role"] == "full93m_chunk"
            and result["source_parent_sha256"] == parent
            and result["source_chunk_index"] == chunk_index
            and result["entire_original_source_part_completed"] is False
            and result["actual_unique_source_photo_ids_attempted"] <=
                int(row["expected_photo_ids"]),
            "A partial photo chunk must never be recorded as a whole source part")
    summary = {
        "schema": SCHEMA,
        "parent_source_sha256": parent,
        "original_source_tier": tier,
        "original_source_part": part,
        "source_chunk_index": chunk_index,
        "source_chunk_manifest_sha256": report["manifest_sha256"],
        "original_source_part_total_photo_ids": report["actual_original_photo_IDs"],
        "source_chunk_expected_photo_ids": int(row["expected_photo_ids"]),
        "actual_new_or_resumed_medium_photo_attempts": result["actual_unique_source_photo_ids_attempted"],
        "source_chunk_all_photo_ids_scanned": complete,
        "entire_parent_original_photo_part_scanned": False,
        "source_93m_all_images_scanned": False,
        "original_four_hue_ROI_classification_modified": False,
        "medium_decoded_photo_count": result["medium_image_pixels_decoded"],
        "screen_counts": result["screen_counts"],
    }
    progress = output / ("chunk_progress_%05d.json" % result["actual_unique_source_photo_ids_attempted"])
    if progress.exists():
        require(json.loads(progress.read_text()) == summary,
                "Previously sealed photo-chunk progress receipt changed")
    else:
        progress.write_text(json.dumps(summary, indent=2,sort_keys=True) + "\n")
    return summary

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--chunks",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--chunk-index",type=int,required=True)
    p.add_argument("--weights",type=Path,required=True)
    p.add_argument("--weights-sha256",required=True)
    p.add_argument("--device",default="cpu")
    p.add_argument("--workers",type=int,default=8)
    p.add_argument("--batch-size",type=int,default=512)
    p.add_argument("--requests-per-second",type=float,default=2.5)
    p.add_argument("--max-packets",type=int)
    a=p.parse_args()
    require(0<a.requests_per_second<=3 and 1<=a.workers<=12,
            "Pilot and provider-rate safety gate")
    detector=WarmFlowerBoxDetector(a.weights,a.weights_sha256,0,"item",
                                   device=a.device,conf=.02)
    pacer=Pacer(a.requests_per_second)
    res=run(a.chunks,a.output,a.chunk_index,
            detector=detector,
            download=lambda url:download_medium(url,pacer),
            workers=a.workers,batch_size=a.batch_size,
            max_packets=a.max_packets)
    print("MOIS_REAL_SIGNED_PARENT_93M_SOURCE_CHUNK_PIXEL_PROGRESS",
          json.dumps(res,sort_keys=True))

if __name__=="__main__":
    main()
