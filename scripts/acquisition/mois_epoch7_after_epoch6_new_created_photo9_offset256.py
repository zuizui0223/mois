#!/usr/bin/env python3
"""MOIS epoch7: new-created interval after signed complete epoch6.

This is a new temporal source interval, not a replay of any epoch6 source SHA.
The source selects 128 taxon IDs from the 2,197 still-zero-strict9 taxa in a
frozen epoch7 hash order; previously queried taxa are eligible only for
observations created after the epoch6 final source receipt. All older photo
IDs/observation UUIDs are excluded. Raw image SHA disjointness against the
full historical ~1.23 million images is NOT certified.
"""
import argparse
import hashlib
import json
from pathlib import Path

from scripts.acquisition import (
    mois_epoch6_post_complete2199_new_created_original_photo9_offset2176 as predecessor
)
from scripts.acquisition import mois_20261115_created_epoch2_disjoint_original_photo9_after2048 as core

m = predecessor.m
m.SCHEMA = "mois_incremental_20261117_epoch7_offset256_new_time_window_1231352_original9_v1"
m.SOURCE_ARM = "MOIS_EPOCH7_OFFSET256_AFTER_TERMINAL_EPOCH6_2199_NEW_CREATED_INTERVAL"
# Epoch6 final source artifact was frozen at 2026-11-17T11:47:05Z.
m.MIN_CREATED = "2026-11-17T11:47:06+00:00"
m.CREATED_D1 = "2026-11-17"
m.MAX_SOURCE_TAXA = 128
m.ATLAS_SHA = "8589c7ec5cfafcf22141bdfcf8c66e70b1cc7a0b0c5b5c7b8e6c0adef4d99b30"
m.OLD_PHOTOS = 1231352
m.OLD_NUMERIC_PHOTOS = 467568
m.OLD_NUMERIC_TAXA = 39914
m.OLD_CLEAR_TAXA = 39174

def select_taxa_epoch7(atlas, n=128, offset=256):
    m.require(len(atlas) == 42111 and atlas.inat_taxon_id.is_unique,
              "nominal atlas must contain exactly 42,111 distinct taxa")
    z = atlas.loc[atlas.after_strict9.eq(0)].copy()
    m.require(len(z) == 2197 and z.inat_taxon_id.str.fullmatch(r"\d+").all(),
              "epoch7 must select only the 2,197 presently unmeasured taxa")
    z["priority"] = z.inat_taxon_id.map(
        lambda tid: hashlib.sha256(
            ("MOIS-20261117-EPOCH7-AFTER-EPOCH6-TERMINAL-TEMPORAL-NEW-CREATED|" + tid).encode()
        ).hexdigest())
    ordered = z.sort_values("priority", kind="stable").inat_taxon_id.tolist()
    selected = ordered[offset:offset+n]
    m.require(offset == 256 and n == 128
              and len(selected) == len(set(selected)) == n
              and not set(selected).intersection(ordered[:offset]),
              "epoch7 third disjoint taxon cohort changed")
    return selected

# Existing signed history includes the 1,231,352 photo-ID/UUID ledger, the
# independently audited epoch6 offset1920 technical unknown, and 108 recent
# actual raw-image hashes. The zero-source offset2176 added no photo IDs.
m.history = predecessor.history_after_offset1920
m.select_taxa = select_taxa_epoch7
m.source = core.source

def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="stage", required=True)
    for stage, keys in (
        ("source", ("atlas", "output")),
        ("audit", ("source", "worker", "output")),
        ("append", ("atlas", "source", "worker", "audit", "output")),
    ):
        q = sub.add_parser(stage)
        for key in keys:
            q.add_argument("--" + key, required=True)
        if stage == "source":
            q.add_argument("--offset", type=int, default=128)
    args = parser.parse_args()
    if args.stage == "source":
        m.source(args.atlas, args.output, quota=128, offset=args.offset)
        receipt_path = Path(args.output) / "receipt.json"
        r = json.loads(receipt_path.read_bytes())
        m.require(
            (r["hash_disjoint_taxon_offset"], r["searched_zero_numeric_taxa"],
             r["source_API_error_count"], r["old_original_photo_attempt_ID_count"],
             r["old_original_photo_obs_ID_overlap"]) == (256, 128, 0, 1231352, 0)
            and r["old_all_image_raw_SHA_disjointness_proven"] is False
            and r["created_utc_lower_bound"] == m.MIN_CREATED
            and r["original_pixels_fetched"] == r["strict9_numeric_photos"] == 0,
            "new temporal source-only receipt gate failed")
        r["recent_known_original_image_raw_SHA_count"] = 108
        r["prior_signed_epoch6_taxon_frame_complete"] = 2199
        r["currently_zero_numeric_taxa"] = 2197
        r["previous_epoch6_terminal_run"] = 44076781021
        r["previous_epoch7_offset0_run"] = 44084420377
        r["previous_epoch7_offset0_source_ZIP_SHA256"] = "8eaee465a4d4b31c318a0abd26303cebcedaa23345fa32d57c16cbdd84ac28d2"
        r["previous_epoch7_offset0_source_candidate_count"] = 0
        r["previous_epoch7_offset128_run"] = 44088047232
        r["previous_epoch7_offset128_source_ZIP_SHA256"] = "a1fdd2b04e03a7d8455c7bd09c0959016a75427bd1c6d2a6c2c35dc63d5ddf46"
        r["previous_epoch7_offset128_source_candidate_count"] = 0
        r["previous_epoch6_terminal_source_ZIP_SHA256"] = (
            "3ce764b6d837348d02009173513539772d71de569088832a2f768d09c22fb87e"
        )
        r.pop("recent_97_nonempty_raw_image_SHA_count", None)
        receipt_path.write_text(json.dumps(r, sort_keys=True, indent=2) + "\n")
        print("MOIS_EPOCH7_OFFSET256_NEW_TEMPORAL_SOURCE_SIGNED",
              json.dumps(r, sort_keys=True), flush=True)
    elif args.stage == "audit":
        m.audit(args.source, args.worker, args.output)
    else:
        m.append(args.atlas, args.source, args.worker, args.audit, args.output)

if __name__ == "__main__":
    main()
