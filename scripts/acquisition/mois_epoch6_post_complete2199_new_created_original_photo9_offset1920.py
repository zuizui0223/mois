#!/usr/bin/env python3
"""MOIS epoch6 offset1920, 15 signed disjoint predecessor source cohorts.

The frozen 2,199-taxon ordering predates first-ROI gains at offsets 768 and 1792.
Restore both newly covered taxa into the original priority frame before slicing
the distinct 1920:2048 cohort. Never equate original photo RGB9 with species-wide
petal colour or a GBIF occurrence coordinate. Full historical raw SHA comparison
is unavailable, even after photo-ID and observation UUID census exclusions.
"""
import argparse
import hashlib
import json
from pathlib import Path
import pandas as pd
from scripts.acquisition import mois_epoch6_post_complete2199_new_created_original_photo9_offset1792 as predecessor
from scripts.acquisition import mois_20261115_created_epoch2_disjoint_original_photo9_after2048 as core

m = predecessor.m
NEW_DELTA = "old/epoch6_offset1792/new_real_original_nine_photo9_delta.parquet"
NEW_DELTA_SHA = "dd3f86884ec5c3890371f0d53b0f52eb83bc6bd8d23c5d341c2aad6452e58453"
PREVIOUS_NEW_NUMERIC = ("1249987", "560402")
m.SCHEMA = "mois_incremental_20261117_epoch6_offset1920_after_terminal1792_1231351_original_photo9_v1"
m.SOURCE_ARM = "MOIS_EPOCH6_OFFSET1920_POST_SIGNED1792_FULL_1231351_FROZEN2199_DISJOINT"
m.MIN_CREATED = "2026-11-15T21:27:57+00:00"
m.CREATED_D1 = "2026-11-15"
m.MAX_SOURCE_TAXA = 128
m.ATLAS_SHA = "6aa39e76a039e0c3161a5d115f5dd27fd0e91e438ae947df7f284fbeccf2796f"
m.OLD_PHOTOS = 1231351
m.OLD_NUMERIC_PHOTOS = 467568
m.OLD_NUMERIC_TAXA = 39914
m.OLD_CLEAR_TAXA = 39174

def history_after_offset1792():
    m.require(Path(NEW_DELTA).is_file() and m.sha256(NEW_DELTA) == NEW_DELTA_SHA,
              "offset1792 independently signed real-pixel delta absent or changed")
    current = m.OLD_PHOTOS
    try:
        m.OLD_PHOTOS = 1231350
        ids, observations, raw = predecessor.history_after_offset1152()
    finally:
        m.OLD_PHOTOS = current
    m.require(len(ids) == 1231350 and len(raw) == 106,
              "predecessor photo-ID and recent image raw-SHA census drift")
    d = pd.read_parquet(NEW_DELTA)
    m.require(len(d) == 1 and {"photo_id","observation_id","image_sha256",
                                "strict9","inat_taxon_id"}.issubset(d.columns),
              "offset1792 independently audited pixel delta identity absent")
    r = d.iloc[0]
    photo, obs, img = str(r.photo_id), str(r.observation_id), str(r.image_sha256)
    m.require(photo == "764438797"
              and obs == "949fcd16-9591-4f94-b0d0-291cab53125c"
              and img == "b06119dcc68a0231a63c9fd18120af2394914e42f096c964c48e112f783965db"
              and str(r.inat_taxon_id) == "560402" and bool(r.strict9)
              and photo not in ids and obs not in observations and img not in raw,
              "offset1792 original image ID/UUID/raw SHA or ROI9 status invalid")
    fractions = [r["strict9_" + s] for s in m.PALETTE]
    m.require(len(fractions) == 9
              and all(pd.notna(x) and float(x) >= 0 for x in fractions)
              and abs(sum(float(x) for x in fractions)-1.0) <= 0.001,
              "offset1792 genuine RGB9 finite/nonnegative/sum1 invariant failed")
    ids.add(photo); observations.add(obs); raw.add(img)
    m.require(len(ids) == 1231351 and len(raw) == 107,
              "complete independent original photo ID or available raw SHA count mismatch")
    return ids, observations, raw

def select_taxa_epoch6(atlas, n=128, offset=1920):
    missing = atlas.loc[atlas.after_strict9.eq(0)].copy()
    two_new = atlas.loc[atlas.inat_taxon_id.isin(PREVIOUS_NEW_NUMERIC)].copy()
    m.require(len(atlas) == 42111 and atlas.inat_taxon_id.is_unique
              and len(missing) == 2197 and len(two_new) == 2
              and two_new.after_strict9.eq(1).all()
              and not set(PREVIOUS_NEW_NUMERIC).intersection(missing.inat_taxon_id)
              and n == 128 and offset == 1920,
              "restore exactly two numeric first-coverage taxa before frozen slice")
    frozen = pd.concat([missing, two_new], ignore_index=True)
    m.require(len(frozen) == 2199 and frozen.inat_taxon_id.is_unique
              and frozen.inat_taxon_id.str.fullmatch(r"\d+").all(),
              "frozen original 2,199-taxon frame inconsistent")
    frozen["priority"] = frozen.inat_taxon_id.map(
        lambda tid: hashlib.sha256(
            ("MOIS-20261117-EPOCH6-POST-2199-DISJOINT-NEW-CREATED|" + tid).encode()
        ).hexdigest()
    )
    ordered = frozen.sort_values("priority", kind="stable").inat_taxon_id.tolist()
    m.require(PREVIOUS_NEW_NUMERIC[0] in ordered[768:896]
              and PREVIOUS_NEW_NUMERIC[1] in ordered[1792:1920],
              "previously measured taxa have shifted from their signed source cohorts")
    selected = ordered[offset:offset+n]
    m.require(len(selected) == len(set(selected)) == n
              and not set(selected).intersection(ordered[:offset])
              and not set(PREVIOUS_NEW_NUMERIC).intersection(selected),
              "offset1920 cohort repeats already searched taxon")
    return selected

m.history = history_after_offset1792
m.select_taxa = select_taxa_epoch6
m.source = core.source

def main():
    p = argparse.ArgumentParser()
    stages = p.add_subparsers(dest="stage", required=True)
    for stage, keys in (
        ("source", ("atlas", "output")),
        ("audit", ("source", "worker", "output")),
        ("append", ("atlas", "source", "worker", "audit", "output")),
    ):
        q = stages.add_parser(stage)
        for key in keys:
            q.add_argument("--" + key, required=True)
        if stage == "source":
            q.add_argument("--offset", type=int, default=1920)
    args = p.parse_args()
    if args.stage == "source":
        m.source(args.atlas, args.output, quota=128, offset=args.offset)
        rp = Path(args.output) / "receipt.json"
        r = json.loads(rp.read_bytes())
        m.require(
            (r["hash_disjoint_taxon_offset"], r["searched_zero_numeric_taxa"],
             r["source_API_error_count"], r["old_original_photo_attempt_ID_count"],
             r["old_original_photo_obs_ID_overlap"]) == (1920, 128, 0, 1231351, 0)
            and r["old_all_image_raw_SHA_disjointness_proven"] is False
            and r["created_utc_lower_bound"] == m.MIN_CREATED
            and r["original_pixels_fetched"] == r["strict9_numeric_photos"] == 0,
            "offset1920 frozen source-only gate failed",
        )
        r["recent_known_original_image_raw_SHA_count"] = 107
        r["frozen_original_zero_numeric_taxon_frame"] = 2199
        r["currently_zero_numeric_taxa"] = 2197
        r["new_numeric_taxa_removed_from_sampling_only"] = list(PREVIOUS_NEW_NUMERIC)
        r.pop("recent_97_nonempty_raw_image_SHA_count", None)
        rp.write_text(json.dumps(r, sort_keys=True, indent=2) + "\n")
        print("MOIS_EPOCH6_OFFSET1920_SOURCE_SIGNED",
              json.dumps(r, sort_keys=True), flush=True)
    elif args.stage == "audit":
        m.audit(args.source, args.worker, args.output)
    else:
        m.append(args.atlas, args.source, args.worker, args.audit, args.output)

if __name__ == "__main__":
    main()
