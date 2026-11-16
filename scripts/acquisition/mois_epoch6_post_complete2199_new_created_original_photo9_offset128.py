#!/usr/bin/env python3
"""MOIS epoch6 offset128: overlap-safe post-sweep new iNat original-photo opportunities.

No new ROI9 is counted without original RGB worker terminal, independent audit and
append. The historical 1,231,346 photo-ID/UUID exclusions and 102 known raw
image SHAs remain mandatory; all historical raw SHAs are not available.
A conservative overlapping created-after window prevents missed late arrivals
for taxa checked early in epoch5. The epoch6 taxon hash order is fresh.
"""
import argparse
import hashlib
import json
from pathlib import Path
from scripts.acquisition import mois_epoch5_postlast23_created_photo9_offset2176 as previous

m = previous.m
m.SCHEMA = "mois_incremental_20261117_epoch6_offset128_after_terminal2199_1231346_original_photo9_v1"
m.SOURCE_ARM = "MOIS_EPOCH6_OFFSET128_POST_TERMINAL2199_OVERLAP_SAFE_NEW_ORIGINAL"
m.MIN_CREATED = "2026-11-15T21:27:57+00:00"
m.CREATED_D1 = "2026-11-15"
m.MAX_SOURCE_TAXA = 128
m.history = previous.history_2176

def select_taxa_epoch6(atlas, n=128, offset=128):
    z = atlas.loc[atlas.after_strict9.eq(0)].copy()
    m.require(
        len(atlas) == 42111 and len(z) == 2199 and atlas.inat_taxon_id.is_unique
        and z.inat_taxon_id.str.fullmatch(r"\d+").all()
        and n == 128 and offset == 128,
        "epoch6 full atlas / initial 128 taxon frame invalid",
    )
    z["priority"] = z.inat_taxon_id.map(lambda taxon: hashlib.sha256(
        ("MOIS-20261117-EPOCH6-POST-2199-DISJOINT-NEW-CREATED|" + taxon).encode()
    ).hexdigest())
    ordered = z.sort_values("priority", kind="stable").inat_taxon_id.tolist()
    chosen = ordered[128:256]
    m.require(len(ordered) == len(set(ordered)) == 2199
              and len(chosen) == len(set(chosen)) == 128 and not set(chosen).intersection(ordered[:128]),
              "epoch6 selected taxon IDs are not unique")
    return chosen

m.select_taxa = select_taxa_epoch6
m.source = previous.core_source

def main():
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="stage", required=True)
    for stage, fields in (
        ("source", ("atlas", "output")),
        ("audit", ("source", "worker", "output")),
        ("append", ("atlas", "source", "worker", "audit", "output")),
    ):
        q = sub.add_parser(stage)
        for field in fields:
            q.add_argument("--" + field, required=True)
        if stage == "source":
            q.add_argument("--offset", type=int, default=128)
    a = p.parse_args()
    if a.stage == "source":
        m.source(a.atlas, a.output, quota=128, offset=a.offset)
        receipt_file = Path(a.output) / "receipt.json"
        receipt = json.loads(receipt_file.read_bytes())
        m.require(
            (receipt["hash_disjoint_taxon_offset"],
             receipt["searched_zero_numeric_taxa"],
             receipt["source_API_error_count"],
             receipt["old_original_photo_attempt_ID_count"],
             receipt["old_original_photo_obs_ID_overlap"]) == (128, 128, 0, 1231346, 0)
            and receipt["old_all_image_raw_SHA_disjointness_proven"] is False
            and receipt["created_utc_lower_bound"] == m.MIN_CREATED
            and receipt["original_pixels_fetched"] == receipt["strict9_numeric_photos"] == 0,
            "epoch6 source stage did not satisfy provenance gate",
        )
        receipt["recent_known_original_image_raw_SHA_count"] = 102
        receipt.pop("recent_97_nonempty_raw_image_SHA_count", None)
        receipt_file.write_text(json.dumps(receipt, sort_keys=True, indent=2) + "\n")
        print("MOIS_EPOCH6_OFFSET128_SIGNED_SOURCE", json.dumps(receipt, sort_keys=True), flush=True)
    elif a.stage == "audit":
        m.audit(a.source, a.worker, a.output)
    else:
        m.append(a.atlas, a.source, a.worker, a.audit, a.output)

if __name__ == "__main__":
    main()
