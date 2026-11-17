#!/usr/bin/env python3
"""MOIS epoch6 offset768: hash-disjoint 128 zero-ROI taxa after signed offset640.

Use the terminal 1,231,348-photo atlas and full historical photo/observation ledger.
New image source candidates are metadata until original RGB worker, independent
strict9 audit and formal append succeed. Historical ~1.23M raw-image SHA
disjointness is not established.
"""
import argparse
import hashlib
import json
from pathlib import Path
from scripts.acquisition import mois_epoch6_post_complete2199_new_created_original_photo9_offset640 as predecessor
import pandas as pd
from scripts.acquisition import mois_20261115_created_epoch2_disjoint_original_photo9_after2048 as core

m = predecessor.m
NEW_DELTA = "old/epoch6_offset640/new_real_original_nine_photo9_delta.parquet"
NEW_DELTA_SHA = "9b27a4535fd4f6f5db6e0f2329449730301e83a3cc6b47f1380b1e371f002a59"
m.SCHEMA = "mois_incremental_20261117_epoch6_offset768_after_terminal640_1231348_original_photo9_v1"
m.SOURCE_ARM = "MOIS_EPOCH6_OFFSET768_POST_SIGNED128_FULL_1231348_DISJOINT"
m.MIN_CREATED = "2026-11-15T21:27:57+00:00"
m.CREATED_D1 = "2026-11-15"
m.MAX_SOURCE_TAXA = 128
m.ATLAS_SHA = "ee28c331a58c1f16fcd91de31b189366cf16ba76750737741f6fb5b36306b923"
m.OLD_PHOTOS = 1231348
m.OLD_NUMERIC_PHOTOS = 467566
m.OLD_NUMERIC_TAXA = 39912
m.OLD_CLEAR_TAXA = 39173

def history_after_offset640():
    # Load earlier independently audited 1,231,347 originals exactly once.
    # Add the signed offset640 original RGB measurement delta, not source-only metadata.
    m.require(Path(NEW_DELTA).is_file() and m.sha256(NEW_DELTA) == NEW_DELTA_SHA,
              "offset640 independent original-pixel audit delta absent or altered")
    current = m.OLD_PHOTOS
    try:
        m.OLD_PHOTOS = 1231347
        ids, observations, raw = predecessor.history_after_offset512()
    finally:
        m.OLD_PHOTOS = current
    m.require(len(ids) == 1231347 and len(raw) == 103,
              "offset640 predecessor photo or raw-SHA ledger incorrect")
    d = pd.read_parquet(NEW_DELTA)
    m.require(len(d) == 1 and {"photo_id", "observation_id", "image_sha256", "strict9"}.issubset(d.columns),
              "offset640 signed pixel audit delta missing identity or strict9 fields")
    row = d.iloc[0]
    photo, obs, img = str(row["photo_id"]), str(row["observation_id"]), str(row["image_sha256"])
    m.require(photo == "764361971" and obs == "0f629749-ebc5-4cb3-9f6f-ed991530a9bb" and
              img == "562aaa43736460cceb74e25705f6ccdd0c27a518f70dfaed3dcc645914db3ede" and
              not bool(row["strict9"]) and
              photo not in ids and obs not in observations and img not in raw,
              "offset640 new original photo, observation, or SHA not disjoint")
    ids.add(photo)
    observations.add(obs)
    raw.add(img)
    m.require(len(ids) == 1231348 and len(raw) == 104,
              "offset768 full historical photo ID/UUID/raw image-SHA counts incorrect")
    return ids, observations, raw

def select_taxa_epoch6(atlas, n=128, offset=768):
    z = atlas.loc[atlas.after_strict9.eq(0)].copy()
    m.require(len(atlas) == 42111 and len(z) == 2199 and
              atlas.inat_taxon_id.is_unique and
              z.inat_taxon_id.str.fullmatch(r"\d+").all() and
              n == 128 and offset == 768,
              "epoch6 full atlas or offset768 cohort invalid")
    z["priority"] = z.inat_taxon_id.map(lambda taxon: hashlib.sha256(
        ("MOIS-20261117-EPOCH6-POST-2199-DISJOINT-NEW-CREATED|" + taxon).encode()
    ).hexdigest())
    ordered = z.sort_values("priority", kind="stable").inat_taxon_id.tolist()
    selected = ordered[768:896]
    m.require(len(ordered) == len(set(ordered)) == 2199 and
              len(selected) == len(set(selected)) == 128 and
              not set(selected).intersection(ordered[:768]),
              "epoch6 offset768 reuses taxon from signed offsets0/128")
    return selected

m.history = history_after_offset640
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
            p.add_argument("--offset", type=int, default=768)
    args = parser.parse_args()
    if args.stage == "source":
        m.source(args.atlas, args.output, quota=128, offset=args.offset)
        receipt_path = Path(args.output) / "receipt.json"
        r = json.loads(receipt_path.read_bytes())
        m.require(
            (r["hash_disjoint_taxon_offset"], r["searched_zero_numeric_taxa"],
             r["source_API_error_count"], r["old_original_photo_attempt_ID_count"],
             r["old_original_photo_obs_ID_overlap"]) == (768, 128, 0, 1231348, 0)
            and r["old_all_image_raw_SHA_disjointness_proven"] is False
            and r["created_utc_lower_bound"] == m.MIN_CREATED
            and r["original_pixels_fetched"] == r["strict9_numeric_photos"] == 0,
            "epoch6 offset768 source provenance gate failed",
        )
        r["recent_known_original_image_raw_SHA_count"] = 104
        r.pop("recent_97_nonempty_raw_image_SHA_count", None)
        receipt_path.write_text(json.dumps(r, sort_keys=True, indent=2) + "\n")
        print("MOIS_EPOCH6_OFFSET768_SOURCE_SIGNED", json.dumps(r, sort_keys=True), flush=True)
    elif args.stage == "audit":
        m.audit(args.source, args.worker, args.output)
    else:
        m.append(args.atlas, args.source, args.worker, args.audit, args.output)

if __name__ == "__main__":
    main()
