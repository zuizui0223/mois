#!/usr/bin/env python3
"""Split ONE SHA-verified MOIS original 93m source photo shard into resumable chunks.

Streaming and outcome-blind. No photo bytes, model inference, photo/species caps
or geographic filters. This solves the 594,457-row largest tier/part job timeout:
a signed 4,096-original-photo chunk is an independent bounded image task.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import tempfile

from scripts.acquisition.mois_onepass_full93m_target_tiers import COLUMNS

SCHEMA = "mois_93m_exact_full_source_part_to_bounded_photo_chunks_v1"
FIELDS = (
    "tier", "part", "chunk_index", "original_source_part_sha256",
    "chunk_relative_path", "chunk_gzip_sha256",
    "expected_photo_ids", "first_photo_id", "last_photo_id",
)
def require(value, reason):
    if not value:
        raise ValueError("MOIS_93M_SOURCE_CHUNKS: " + reason)

def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for piece in iter(lambda: stream.read(4 << 20), b""):
            h.update(piece)
    return h.hexdigest()

def run(source, output, *, expected_source_sha256,
        expected_photo_count, tier, part, chunk_size=4096):
    source, output = Path(source), Path(output)
    require(not output.exists(), "Refuse overwrite of immutable source chunk registry")
    require(0 <= tier < 14 and 0 <= part < 128 and 1 <= chunk_size <= 8192,
            "Invalid tier, part or bounded chunk size")
    require(type(expected_photo_count) is int and expected_photo_count > 0,
            "Expected source photo-ID count missing")
    require(sha(source) == expected_source_sha256,
            "Original published source part SHA256 differs")
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix="mois-93m-photo-chunks-", dir=output.parent))
    chunks = staging / "chunks"
    chunks.mkdir()
    count = 0
    manifest = []
    prev_ids = set()
    try:
        with gzip.open(source, "rt", encoding="utf-8", newline="") as inp:
            reader = csv.DictReader(inp)
            require(tuple(reader.fieldnames or ()) == tuple(COLUMNS),
                    "Not a frozen MOIS source tier photo/observation CSV")
            iterator = iter(reader)
            while True:
                rows = []
                for _ in range(chunk_size):
                    row = next(iterator, None)
                    if row is None:
                        break
                    pid = row["photo_id"]
                    require(pid.isdecimal() and row["observation_uuid"]
                            and row["inat_taxon_id"].isdecimal()
                            and row["tier"] == str(tier),
                            "Original photo, observation, taxon or tier invalid")
                    require(pid not in prev_ids,
                            "Repeated original photo ID in a signed source part")
                    prev_ids.add(pid)
                    rows.append(row)
                if not rows:
                    break
                index = len(manifest)
                name = f"chunk-{index:05d}.csv.gz"
                dest = chunks / name
                with dest.open("wb") as raw:
                    with gzip.GzipFile(filename="", fileobj=raw,
                                       mode="wb", compresslevel=2, mtime=0) as zipped:
                        with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as text:
                            writer = csv.DictWriter(text, fieldnames=COLUMNS,
                                                    lineterminator="\n")
                            writer.writeheader()
                            writer.writerows(rows)
                count += len(rows)
                manifest.append({
                    "tier": tier, "part": part, "chunk_index": index,
                    "original_source_part_sha256": expected_source_sha256,
                    "chunk_relative_path": "chunks/" + name,
                    "chunk_gzip_sha256": sha(dest),
                    "expected_photo_ids": len(rows),
                    "first_photo_id": rows[0]["photo_id"],
                    "last_photo_id": rows[-1]["photo_id"],
                })
        require(count == expected_photo_count,
                "Signed source part photo count does not reconcile across chunks")
        require(len(prev_ids) == count, "Original photo IDs are not unique")
        with (staging / "chunk_manifest.csv").open("w", encoding="utf-8",
                                                   newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=FIELDS, lineterminator="\n")
            writer.writeheader()
            writer.writerows(manifest)
        receipt = {
            "schema": SCHEMA,
            "original_source_part_sha256": expected_source_sha256,
            "original_tier": tier, "original_part": part,
            "expected_original_photo_IDs": expected_photo_count,
            "actual_original_photo_IDs": count,
            "chunk_size_maximum": chunk_size,
            "chunks": len(manifest),
            "first_photo_id": manifest[0]["first_photo_id"],
            "last_photo_id": manifest[-1]["last_photo_id"],
            "manifest_sha256": sha(staging / "chunk_manifest.csv"),
            "source_fully_partitioned_once_into_nonoverlapping_photo_chains": True,
            "image_pixels_downloaded": 0,
            "floral_hue_or_bloom_inferred": False,
            "native_or_taxon_filt": False,
        }
        (staging / "chunk_receipt.json").write_text(
            json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        os.replace(staging, output)
        return receipt
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--source-sha256", required=True)
    parser.add_argument("--expected-count", type=int, required=True)
    parser.add_argument("--tier", type=int, required=True)
    parser.add_argument("--part", type=int, required=True)
    parser.add_argument("--chunk-size", type=int, default=4096)
    parser.add_argument("--output", type=Path, required=True)
    a = parser.parse_args()
    r = run(a.source, a.output, expected_source_sha256=a.source_sha256,
            expected_photo_count=a.expected_count, tier=a.tier,
            part=a.part, chunk_size=a.chunk_size)
    print("MOIS_93M_RESUMABLE_SOURCE_PHOTO_CHUNKS", json.dumps(r, sort_keys=True))

if __name__ == "__main__":
    main()
