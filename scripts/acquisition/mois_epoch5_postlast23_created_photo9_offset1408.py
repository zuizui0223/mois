#!/usr/bin/env python3
"""MOIS epoch5 offset1408: disjoint 128 zero-numeric taxa, after preceding source-zero terminal (formal atlas remains at offset896).

No selected candidate is a numeric ROI without terminal original RGB, independent
audit and formal append. The full 1,231,346 photo-ID/observation UUID ledger is
checked. Only 102 recent raw image SHAs are known, not 1.23 million.
"""
import argparse
import hashlib
import json
from pathlib import Path
from scripts.acquisition.mois_20261115_created_epoch2_disjoint_original_photo9_after2048 import source as core_source
from scripts.acquisition import mois_epoch4_post2048_last23_original_photo9 as predecessor

m = predecessor.m
LAST_DELTA_PATH = "old/epoch5_offset896/new_real_original_nine_photo9_delta.parquet"
LAST_DELTA_SHA = "545d034d3a396db769a9ca26640405bfd9d9b0349074306259cfb906f0a52b11"
LAST_ATLAS_SHA = "c8c9084ff6ba6d587c41ffd168b942fa21a2fc9932520092598b420dc89bf2b1"
m.SCHEMA = "mois_incremental_20261116_epoch5_post_offset896_1231346_original_photo9_v1"
m.SOURCE_ARM = "MOIS_EPOCH5_OFFSET1408_FULL_1231346_ID_UUID_RECENT102_SHA"
m.MIN_CREATED = "2026-11-15T21:27:57+00:00"
m.CREATED_D1 = "2026-11-15"
m.MAX_SOURCE_TAXA = 128
m.ATLAS_SHA = LAST_ATLAS_SHA
m.OLD_PHOTOS = 1231346
m.OLD_NUMERIC_PHOTOS = 467566
m.OLD_NUMERIC_TAXA = 39912
m.OLD_CLEAR_TAXA = 39173
m.DELTA_FILES = m.DELTA_FILES + (LAST_DELTA_PATH,)
m.DELTA_N = m.DELTA_N + (1,)

def history_1408():
    checks = (
        (predecessor.prev.EPOCH3_DELTA_PATH, predecessor.prev.EPOCH3_DELTA_SHA),
        (predecessor.prev.PREV_DELTA_PATH, predecessor.prev.PREV_DELTA_SHA),
        (predecessor.prev.NEW_DELTA_PATH, predecessor.prev.NEW_DELTA_SHA),
        (predecessor.LAST_DELTA, predecessor.LAST_DELTA_SHA),
        (LAST_DELTA_PATH, LAST_DELTA_SHA),
    )
    for path, sha in checks:
        m.require(Path(path).is_file() and m.sha256(path) == sha,
                  "signed historical independent delta absent/altered: " + path)
    ids, observations, recent_raw = predecessor.prev._full_history()
    m.require(
        len(ids) == 1231346 and len(recent_raw) == 102 and
        "764119551" in ids and
        "48da2604-b864-4216-a73c-8d290724a61d" in observations and
        "0746e7b90e71677e9df3145a7b0c53918db67cee8f3c321d51a97747cd82b1c1" in recent_raw and
        "763811321" in ids and
        "183a6f88-4367-417b-a0b3-404d0a0ccd0b" in observations and
        "fb25d663cd95983e6e44078efe8a29504a2f3ad1d1371a9bedb550913ebffe56" in recent_raw,
        "formal post-offset896 history/102 recent raw SHA incomplete")
    return ids, observations, recent_raw

def select_taxa_1408(atlas, n=128, offset=1408):
    z = atlas.loc[atlas.after_strict9.eq(0)].copy()
    m.require(len(atlas) == 42111 and len(z) == 2199 and
              atlas.inat_taxon_id.is_unique and
              z.inat_taxon_id.str.fullmatch(r"\d+").all() and
              n == 128 and offset == 1408,
              "frozen species frame/next offset invalid")
    z["priority"] = z.inat_taxon_id.map(lambda s: hashlib.sha256(
        ("MOIS-20261116-EPOCH5-POST-21:27:57Z-DISJOINT-2199|" + s).encode()).hexdigest())
    ordered = z.sort_values("priority", kind="stable").inat_taxon_id.tolist()
    selected = ordered[1408:1536]
    m.require(len(ordered) == len(set(ordered)) == 2199 and
              len(selected) == len(set(selected)) == 128 and
              not set(selected).intersection(ordered[:1408]),
              "next cohort reused previous 1,408 taxon IDs")
    return selected

m.history = history_1408
m.select_taxa = select_taxa_1408
m.source = core_source

def main():
    parser = argparse.ArgumentParser()
    stages = parser.add_subparsers(dest="stage", required=True)
    for stage, fields in (("source", ("atlas", "output")),
                          ("audit", ("source", "worker", "output")),
                          ("append", ("atlas", "source", "worker", "audit", "output"))):
        sub = stages.add_parser(stage)
        for flag in fields:
            sub.add_argument("--" + flag, required=True)
        if stage == "source":
            sub.add_argument("--offset", type=int, default=1408)
    args = parser.parse_args()
    if args.stage == "source":
        m.source(args.atlas, args.output, quota=128, offset=args.offset)
        receipt = Path(args.output) / "receipt.json"
        r = json.loads(receipt.read_bytes())
        m.require((r["hash_disjoint_taxon_offset"], r["searched_zero_numeric_taxa"],
                   r["source_API_error_count"], r["old_original_photo_attempt_ID_count"],
                   r["old_original_photo_obs_ID_overlap"]) == (1408, 128, 0, 1231346, 0)
                  and r["old_all_image_raw_SHA_disjointness_proven"] is False
                  and r["original_pixels_fetched"] == r["strict9_numeric_photos"] == 0
                  and r["created_utc_lower_bound"] == m.MIN_CREATED,
                  "signed source not eligible for downstream worker")
        r["recent_known_original_image_raw_SHA_count"] = 102
        r.pop("recent_97_nonempty_raw_image_SHA_count", None)
        receipt.write_text(json.dumps(r, sort_keys=True, indent=2) + "\n")
        print("MOIS_EPOCH5_OFFSET1408_SOURCE_REBASED_SIGNED", json.dumps(r, sort_keys=True), flush=True)
    elif args.stage == "audit":
        m.audit(args.source, args.worker, args.output)
    else:
        m.append(args.atlas, args.source, args.worker, args.audit, args.output)

if __name__ == "__main__":
    main()
