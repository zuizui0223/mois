#!/usr/bin/env python3
"""Epoch3 after first genuine offset1792 ROI: immutable cohort, rebased evidence.

Frozen epoch2 atlas determines the taxon-hash order so a first newly measured
species does not shift subsequent 128-taxon cohorts. The newer signed epoch3
atlas determines which photo IDs/UUIDs and nine-colour numerical traits exist.
No photo pixel, original resolution, botanical petal, or all-history raw SHA
certainty is inferred from a source roster.
"""
import argparse
from pathlib import Path
import pandas as pd
from scripts.acquisition import mois_epoch3_new_created_disjoint_photo9 as previous
m = previous.m

# A validated terminal 4/4-success stage and exact byte hash, not a planned job.
NEW_ATLAS_SHA = "7a2cd8cab556f4ec3f7496d0bc62e77dc9ec82f834ccf0cbd1ce25a315afe5ee"
FROZEN_ATLAS_SHA = "6d08bfc203f317f1fc4a8d3d5f22015a7d0b393cba8ce9fcb18dd0ea74fe1909"
NEW_DELTA_SHA = "537aa3b8419411ec82cecfc6a9768cfb313d65364a4b85afc1edb2953211cfb9"
NEW_DELTA_PATH = "old/epoch3_offset1792/new_real_original_nine_photo9_delta.parquet"
m.SCHEMA = "mois_incremental_20261115_epoch3_rebased_after_offset1792_original_photo9_v2"
m.SOURCE_ARM = "MOIS_EPOCH3_POST1792_SIGNED_UNUSED_ORIGINAL"
m.ATLAS_SHA = NEW_ATLAS_SHA
m.OLD_PHOTOS = 1231342
m.OLD_NUMERIC_PHOTOS = 467566
m.OLD_NUMERIC_TAXA = 39912
m.OLD_CLEAR_TAXA = 39173
m.MAX_SOURCE_TAXA = 128
m.DELTA_FILES = m.DELTA_FILES + (NEW_DELTA_PATH,)
m.DELTA_N = m.DELTA_N + (1,)

_previous_history = m.history

def history():
    m.require(Path(NEW_DELTA_PATH).is_file() and m.sha256(NEW_DELTA_PATH) == NEW_DELTA_SHA,
              "epoch3 offset1792 independently audited byte delta missing/altered")
    ids, observations, recent_raw = _previous_history()
    m.require(len(ids) == 1231342 and
              "763768595" in ids and
              "821d0884-737a-4644-8f8e-6e36c278b9fb" in observations and
              "581ab6b386ab9a2667121a6250ee3d409f5bd5cac914c0777bdab3e0527d1d84" in recent_raw and
              len(recent_raw) == 98,
              "post1792 photo/UUID/raw-SHA identity ledger incomplete")
    return ids, observations, recent_raw
m.history = history

def select_taxa(atlas, n=128, offset=0):
    frozen = Path("frozen/epoch2_atlas.csv.gz")
    m.require(frozen.is_file() and m.sha256(frozen) == FROZEN_ATLAS_SHA,
              "frozen pre-epoch3 taxon cohort atlas absent/altered")
    old = pd.read_csv(frozen, dtype={"inat_taxon_id": str})
    m.require(len(old) == 42111 and old.inat_taxon_id.is_unique and
              int(old.after_strict9.eq(0).sum()) == 2200,
              "frozen epoch2 target state changed")
    cohort = previous.select_epoch3_taxa(old, n, offset)
    now = atlas.set_index("inat_taxon_id").loc[cohort]
    m.require(len(cohort) == n and now.after_strict9.eq(0).all(),
              "cohort includes already measured taxa after epoch3 offset1792")
    return cohort
m.select_taxa = select_taxa

def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="stage", required=True)
    for stage, fields in (("source", ("atlas", "output")),
                          ("audit", ("source", "worker", "output")),
                          ("append", ("atlas", "source", "worker", "audit", "output"))):
        p = sub.add_parser(stage)
        for field in fields:
            p.add_argument("--" + field, required=True)
        if stage == "source":
            p.add_argument("--offset", type=int, default=1920)
    args = parser.parse_args()
    if args.stage == "source":
        m.source(args.atlas, args.output, quota=128, offset=args.offset)
    elif args.stage == "audit":
        m.audit(args.source, args.worker, args.output)
    else:
        m.append(args.atlas, args.source, args.worker, args.audit, args.output)

if __name__ == "__main__":
    main()
