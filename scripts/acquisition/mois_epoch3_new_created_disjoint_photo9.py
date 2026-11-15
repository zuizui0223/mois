#!/usr/bin/env python3
"""One-shot epoch3: bounded distinct taxa and newly created observation photo9.

Chronological overlap with epoch2 begins at its earliest queried moment, so
observations appearing after an early group was scanned are not lost. All
previous 1,231,341 photo IDs, observation UUIDs and 97 available raw-image
digests remain excluded. Only source candidates pass to original-pixel FCP
after a separately signed source artifact; unknown nine-vectors stay missing.
This cannot certify raw-byte uniqueness against all historical 1.23m photos.
"""
import argparse
import hashlib
from scripts.acquisition import mois_20261115_created_epoch2_disjoint_original_photo9_after2048 as m

m.SCHEMA = "mois_incremental_20261115_epoch3_post_earliest_epoch2_scan_cc_original_photo9_v1"
m.SOURCE_ARM = "MOIS_EPOCH3_NEW_CREATED_HISTORY_EXCLUDED_ORIGINALS"
m.MIN_CREATED = "2026-11-14T18:43:28+00:00"
m.CREATED_D1 = "2026-11-14"
m.MAX_SOURCE_TAXA = 128

def select_epoch3_taxa(a, n=128, offset=0):
    z = a.loc[a.after_strict9.eq(0)].copy()
    m.require(len(z)==2200 and z.inat_taxon_id.str.fullmatch(r"\d+").all(),
              "epoch3 target 2,200 prior zero-numeric taxa changed")
    z["priority"] = z.inat_taxon_id.map(
        lambda s: hashlib.sha256(("MOIS-20261115-EPOCH3-DISJOINT-CREATED|"+s).encode()).hexdigest())
    taxa=z.sort_values("priority", kind="stable").iloc[offset:offset+n].inat_taxon_id.tolist()
    m.require(len(taxa)==n and len(set(taxa))==n, "epoch3 taxon cohort not unique")
    return taxa

m.select_taxa=select_epoch3_taxa

def main():
    p=argparse.ArgumentParser()
    sub=p.add_subparsers(dest="stage",required=True)
    for stage, flags in [("source",["atlas","output"]),("audit",["source","worker","output"]),
                         ("append",["atlas","source","worker","audit","output"])]:
        q=sub.add_parser(stage)
        for f in flags:q.add_argument("--"+f,required=True)
        if stage=="source":q.add_argument("--offset",type=int,default=0)
    a=p.parse_args()
    if a.stage=="source":m.source(a.atlas,a.output,quota=128,offset=a.offset)
    elif a.stage=="audit":m.audit(a.source,a.worker,a.output)
    else:m.append(a.atlas,a.source,a.worker,a.audit,a.output)
if __name__=="__main__":
    main()
