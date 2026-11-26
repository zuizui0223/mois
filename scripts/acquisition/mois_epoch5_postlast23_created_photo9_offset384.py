#!/usr/bin/env python3
"""MOIS epoch5 offset384: fresh observation-created window, disjoint from audited offsets0,128,256.

Re-query a *new time window*, in a new deterministically hash-partitioned taxon
order. Source candidates are never counted as measured ROI. Every candidate
photo ID and independent observation UUID must be absent from the signed
1,231,345-image history; compare image raw SHA only with the 101 recent known
originals, never claim full 1.23M raw-byte disjointness.
"""
import argparse
import hashlib
import json
from pathlib import Path

# Capture the unmodified core function BEFORE predecessor import monkey-patches its module.
from scripts.acquisition.mois_20261115_created_epoch2_disjoint_original_photo9_after2048 import source as core_source
from scripts.acquisition import mois_epoch4_post2048_last23_original_photo9 as predecessor

m = predecessor.m
LAST_DELTA_SHA = predecessor.LAST_DELTA_SHA
LAST_ATLAS_SHA = predecessor.LAST_ATLAS_SHA
history_last23 = predecessor.history_last23
m.SCHEMA = "mois_incremental_20261116_epoch5_post_terminal2199_window_cc_original_photo9_v1"
m.SOURCE_ARM = "MOIS_EPOCH5_2199_ZERO_NUMERIC_NEW_CREATED_WINDOW_DISJOINT_ORIGINAL"
m.MIN_CREATED = "2026-11-15T21:27:57+00:00"
m.CREATED_D1 = "2026-11-15"
m.MAX_SOURCE_TAXA = 128

def select_epoch5_taxa(atlas, n=128, offset=384):
    z = atlas.loc[atlas.after_strict9.eq(0)].copy()
    m.require(len(atlas) == 42111 and len(z) == 2199 and
              atlas.inat_taxon_id.is_unique and
              z.inat_taxon_id.str.fullmatch(r"\d+").all(),
              "epoch5 signed 42,111 / 2,199 frozen taxon frame invalid")
    m.require(n == 128 and offset == 384, "epoch5 fourth 128 must use offset384 exactly")
    z["priority"] = z.inat_taxon_id.map(lambda s: hashlib.sha256(
        ("MOIS-20261116-EPOCH5-POST-21:27:57Z-DISJOINT-2199|" + s).encode()).hexdigest())
    ordered = z.sort_values("priority", kind="stable").inat_taxon_id.tolist()
    selected = ordered[384:512]
    m.require(not set(selected).intersection(ordered[:384]), "epoch5 offset384 overlaps signed offset0/128/256 taxa")
    m.require(len(ordered) == len(set(ordered)) == 2199 and
              len(selected) == len(set(selected)) == 128, "epoch5 duplicate taxonomy")
    return selected

m.select_taxa = select_epoch5_taxa
m.history = history_last23
m.source = core_source

def main():
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="stage", required=True)
    for stage, fields in (
        ("source", ("atlas", "output")), 
        ("audit", ("source", "worker", "output")),
        ("append", ("atlas", "source", "worker", "audit", "output"))
    ):
        q = sub.add_parser(stage)
        for f in fields:
            q.add_argument("--" + f, required=True)
        if stage == "source":
            q.add_argument("--offset", type=int, default=256)
    a = p.parse_args()
    if a.stage == "source":
        m.source(a.atlas, a.output, quota=128, offset=a.offset)
        receipt_path = Path(a.output) / "receipt.json"
        r = json.loads(receipt_path.read_bytes())
        m.require(r["hash_disjoint_taxon_offset"] == 384 and
                  r["searched_zero_numeric_taxa"] == 128 and
                  r["source_API_error_count"] == 0 and
                  r["old_original_photo_attempt_ID_count"] == 1231345 and
                  r["old_original_photo_obs_ID_overlap"] == 0 and
                  r["old_all_image_raw_SHA_disjointness_proven"] is False and
                  r["original_pixels_fetched"] == r["strict9_numeric_photos"] == 0 and
                  r["created_utc_lower_bound"] == m.MIN_CREATED,
                  "epoch5 signed source expected terminal eligibility failed")
        r["recent_known_original_image_raw_SHA_count"] = 101
        r.pop("recent_97_nonempty_raw_image_SHA_count", None)
        receipt_path.write_text(json.dumps(r, sort_keys=True, indent=2) + "\n")
        print("MOIS_EPOCH5_OFFSET256_SOURCE_128_AFTER_SIGNED_OFFSET256_SUCCESS", json.dumps(r, sort_keys=True), flush=True)
    elif a.stage == "audit":
        m.audit(a.source, a.worker, a.output)
    else:
        m.append(a.atlas, a.source, a.worker, a.audit, a.output)

if __name__ == "__main__":
    main()
