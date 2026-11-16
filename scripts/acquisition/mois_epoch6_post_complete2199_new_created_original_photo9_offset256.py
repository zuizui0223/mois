#!/usr/bin/env python3
"""MOIS epoch6 offset256: hash-disjoint 128 zero-ROI taxa after signed offset128.

Use the terminal 1,231,347-photo atlas and full historical photo/observation ledger.
New image source candidates are metadata until original RGB worker, independent
strict9 audit and formal append succeed. Historical ~1.23M raw-image SHA
disjointness is not established.
"""
import argparse
import hashlib
import json
from pathlib import Path
from scripts.acquisition import mois_epoch6_post_complete2199_new_created_original_photo9_offset128 as predecessor
from scripts.acquisition import mois_20261115_created_epoch2_disjoint_original_photo9_after2048 as core

m = predecessor.m
NEW_DELTA = "old/epoch6_offset128/new_real_original_nine_photo9_delta.parquet"
NEW_DELTA_SHA = "117d262b70866a367f11684a5340a285e30c898810a29e67b69c59578bc54d4a"
m.SCHEMA = "mois_incremental_20261117_epoch6_offset256_after_terminal128_1231347_original_photo9_v1"
m.SOURCE_ARM = "MOIS_EPOCH6_OFFSET256_POST_SIGNED128_FULL_1231347_DISJOINT"
m.MIN_CREATED = "2026-11-15T21:27:57+00:00"
m.CREATED_D1 = "2026-11-15"
m.MAX_SOURCE_TAXA = 128
m.ATLAS_SHA = "fe2799372abcd0223328cafa8dc37a3a2113230bdd0e7832d281eee2882a2a24"
m.OLD_PHOTOS = 1231347
m.OLD_NUMERIC_PHOTOS = 467566
m.OLD_NUMERIC_TAXA = 39912
m.OLD_CLEAR_TAXA = 39173
m.DELTA_FILES = m.DELTA_FILES + (NEW_DELTA,)
m.DELTA_N = m.DELTA_N + (1,)

def history_after_offset128():
    # Validate the old historical SHA chain and *new* independent audit bytes.
    for p, sha in (
        (predecessor.previous.LAST_DELTA_PATH,
         predecessor.previous.LAST_DELTA_SHA),
    ):
        m.require(Path(p).is_file() and m.sha256(p) == sha,
                  "earlier signed independent audit changed: " + p)
    m.require(Path(NEW_DELTA).is_file() and m.sha256(NEW_DELTA) == NEW_DELTA_SHA,
              "offset128 signed audit delta absent or altered")
    ids, observations, raw = core.history()
    m.require(len(ids) == 1231347 and len(raw) == 103 and
              "764358617" in ids and
              "8137069d-6363-4326-a9ab-e1547ae8006f" in observations and
              "94915a0f4042e92445e9ab890ab66db84acc2e83f603a01377cf17deb79b6444" in raw and
              "764119551" in ids and
              "0746e7b90e71677e9df3145a7b0c53918db67cee8f3c321d51a97747cd82b1c1" in raw,
              "offset128 photo/observation/103-known-raw-SHA not in full historical ledger")
    return ids, observations, raw

def select_taxa_epoch6(atlas, n=128, offset=256):
    z = atlas.loc[atlas.after_strict9.eq(0)].copy()
    m.require(len(atlas) == 42111 and len(z) == 2199 and
              atlas.inat_taxon_id.is_unique and
              z.inat_taxon_id.str.fullmatch(r"\d+").all() and
              n == 128 and offset == 256,
              "epoch6 full atlas or offset256 cohort invalid")
    z["priority"] = z.inat_taxon_id.map(lambda taxon: hashlib.sha256(
        ("MOIS-20261117-EPOCH6-POST-2199-DISJOINT-NEW-CREATED|" + taxon).encode()
    ).hexdigest())
    ordered = z.sort_values("priority", kind="stable").inat_taxon_id.tolist()
    selected = ordered[256:384]
    m.require(len(ordered) == len(set(ordered)) == 2199 and
              len(selected) == len(set(selected)) == 128 and
              not set(selected).intersection(ordered[:256]),
              "epoch6 offset256 reuses taxon from signed offsets0/128")
    return selected

m.history = history_after_offset128
m.select_taxa = select_taxa_epoch6
m.source = core.source

def main():
    parser = argparse.ArgumentParser()
    stages = parser.add_subparsers(dest="stage", required=True)
    for stage, fields in (
        ("source", ("atlas", "output")),
        ("audit", ("source", "worker", "output")),
        ("append", ("atlas", "source", "worker", "audit", "output")),
    ):
        p = stages.add_parser(stage)
        for key in fields:
            p.add_argument("--" + key, required=True)
        if stage == "source":
            p.add_argument("--offset", type=int, default=256)
    args = parser.parse_args()
    if args.stage == "source":
        m.source(args.atlas, args.output, quota=128, offset=args.offset)
        receipt_path = Path(args.output) / "receipt.json"
        r = json.loads(receipt_path.read_bytes())
        m.require(
            (r["hash_disjoint_taxon_offset"], r["searched_zero_numeric_taxa"],
             r["source_API_error_count"], r["old_original_photo_attempt_ID_count"],
             r["old_original_photo_obs_ID_overlap"]) == (256, 128, 0, 1231347, 0)
            and r["old_all_image_raw_SHA_disjointness_proven"] is False
            and r["created_utc_lower_bound"] == m.MIN_CREATED
            and r["original_pixels_fetched"] == r["strict9_numeric_photos"] == 0,
            "epoch6 offset256 source provenance gate failed",
        )
        r["recent_known_original_image_raw_SHA_count"] = 103
        r.pop("recent_97_nonempty_raw_image_SHA_count", None)
        receipt_path.write_text(json.dumps(r, sort_keys=True, indent=2) + "\n")
        print("MOIS_EPOCH6_OFFSET256_SOURCE_SIGNED", json.dumps(r, sort_keys=True), flush=True)
    elif args.stage == "audit":
        m.audit(args.source, args.worker, args.output)
    else:
        m.append(args.atlas, args.source, args.worker, args.audit, args.output)

if __name__ == "__main__":
    main()
