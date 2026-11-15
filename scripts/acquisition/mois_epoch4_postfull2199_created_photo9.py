#!/usr/bin/env python3
"""Epoch4: only genuinely unused original photos for the 2,199 still-zero ROI9 taxa.

This is a new hash-disjoint *query* epoch after the completed epoch3 source pass.
Taxa may have been queried in epoch3; photo IDs and independent observation UUIDs
MUST never have been attempted. Do not claim all-history raw-byte SHA coverage.
"""
import argparse
import hashlib
from scripts.acquisition import mois_epoch3_post1792_rebased_photo9 as previous

m = previous.m
m.SCHEMA = "mois_incremental_20261115_epoch4_postfull2199_cc_original_photo9_v1"
m.SOURCE_ARM = "MOIS_EPOCH4_2199_ZERO_NUMERIC_REQUERY_UNUSED_ORIGINAL"
m.MIN_CREATED = "2026-11-14T18:43:28+00:00"
m.CREATED_D1 = "2026-11-14"
m.MAX_SOURCE_TAXA = 128

def select_epoch4_taxa(atlas, n=128, offset=0):
    z = atlas.loc[atlas.after_strict9.eq(0)].copy()
    m.require(len(atlas)==42111 and len(z)==2199 and
              atlas.inat_taxon_id.is_unique and
              z.inat_taxon_id.str.fullmatch(r"\d+").all(),
              "epoch4 exact signed 2,199 zero-numeric taxon frame changed")
    m.require(n==128 and offset>=0 and offset%128==0 and
              offset+n<=len(z),
              "epoch4 128-taxon disjoint source cohort not available")
    z["priority"] = z.inat_taxon_id.map(lambda s: hashlib.sha256(
        ("MOIS-20261115-EPOCH4-HASH-DISJOINT-2199|"+s).encode()).hexdigest())
    ids = z.sort_values("priority",kind="stable").iloc[offset:offset+n].inat_taxon_id.tolist()
    m.require(len(ids)==n and len(set(ids))==n,
              "epoch4 query taxon ID duplication")
    return ids

m.select_taxa = select_epoch4_taxa

def main():
    parser = argparse.ArgumentParser()
    stages = parser.add_subparsers(dest="stage", required=True)
    for stage, required in (
        ("source", ("atlas", "output")),
        ("audit", ("source", "worker", "output")),
        ("append", ("atlas", "source", "worker", "audit", "output"))
    ):
        p=stages.add_parser(stage)
        for name in required:
            p.add_argument("--"+name, required=True)
        if stage=="source":
            p.add_argument("--offset",type=int,default=0)
    args=parser.parse_args()
    if args.stage=="source":
        m.source(args.atlas,args.output,quota=128,offset=args.offset)
    elif args.stage=="audit":
        m.audit(args.source,args.worker,args.output)
    else:
        m.append(args.atlas,args.source,args.worker,args.audit,args.output)

if __name__ == "__main__":
    main()
