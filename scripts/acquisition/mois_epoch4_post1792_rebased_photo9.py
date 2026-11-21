#!/usr/bin/env python3
"""Epoch4 source after independently audited offset1792 technical-unknown original.

Continue the frozen zero-ROI taxon-hash cohort without remeasuring historical
photos. Enforce full original-photo ID and independent observation UUID history
including the two epoch4 technical unknowns (offset896 and offset1792).
Only recent image raw SHA is comparable: all ~1.23m historical image bytes
are NOT claimed to have been checked.
"""
import argparse
from pathlib import Path
from scripts.acquisition import mois_epoch4_post896_rebased_photo9 as prev
from scripts.acquisition import mois_epoch3_post1792_rebased_photo9 as base

m=prev.m
NEW_ATLAS_SHA="55e643c0e64a93c68a64e6b25ada1e6f71caf2d93741b0973577e915dd26931b"
NEW_DELTA_SHA="c15c9accbf55306b46a24bb1e915ade0fbb29324f4d2bb633e0269c6c9c7ff9c"
NEW_DELTA_PATH="old/epoch4_offset1792/new_real_original_nine_photo9_delta.parquet"
PREV_DELTA_SHA="1daa87b330c4664940293b11ca4e9d2a5bd9f7a3ee21e4aa44124206bc4a1cd8"
PREV_DELTA_PATH="old/epoch4_offset896/new_real_original_nine_photo9_delta.parquet"
EPOCH3_DELTA_SHA="537aa3b8419411ec82cecfc6a9768cfb313d65364a4b85afc1edb2953211cfb9"
EPOCH3_DELTA_PATH="old/epoch3_offset1792/new_real_original_nine_photo9_delta.parquet"

m.SCHEMA="mois_incremental_20261116_epoch4_rebased_after_offset1792_original_photo9_v3"
m.SOURCE_ARM="MOIS_EPOCH4_POST1792_SIGNED_UNUSED_ORIGINAL_ID_DISJOINT"
m.ATLAS_SHA=NEW_ATLAS_SHA
m.OLD_PHOTOS=1231344
m.OLD_NUMERIC_PHOTOS=467566
m.OLD_NUMERIC_TAXA=39912
m.OLD_CLEAR_TAXA=39173
m.MAX_SOURCE_TAXA=128
m.DELTA_FILES=m.DELTA_FILES+(NEW_DELTA_PATH,)
m.DELTA_N=m.DELTA_N+(1,)

# Skip the earlier *fixed-count* rebase guards while retaining original
# full historical ID/UUID/delta bytes reader. All added delta SHA are checked.
_full_history=base._previous_history

def history():
    for path,sha in (
        (EPOCH3_DELTA_PATH,EPOCH3_DELTA_SHA),
        (PREV_DELTA_PATH,PREV_DELTA_SHA),
        (NEW_DELTA_PATH,NEW_DELTA_SHA),
    ):
        m.require(Path(path).is_file() and m.sha256(path)==sha,
                  "signed historical delta missing or modified: "+path)
    ids,observations,recent_raw=_full_history()
    m.require(
        len(ids)==1231344 and len(recent_raw)==100 and
        "763768595" in ids and
        "821d0884-737a-4644-8f8e-6e36c278b9fb" in observations and
        "581ab6b386ab9a2667121a6250ee3d409f5bd5cac914c0777bdab3e0527d1d84" in recent_raw and
        "763837356" in ids and
        "993e4dfc-a8c4-4162-9328-d7facffee028" in observations and
        "aa82fdc048f13a6d4f6729cd2bd841d382ffb618ba33db07031647c0a3f51faf" in recent_raw and
        "763563677" in ids and
        "37c68c6f-d877-48aa-b266-586ab6e969e6" in observations and
        "6825d888f96dc177cde4f87396d1bfe317c3547cf4f987a9bfd30c4584c297d0" in recent_raw,
        "post1792 signed photo/observation/raw identity ledger incomplete"
    )
    return ids,observations,recent_raw
m.history=history
m.select_taxa=prev.prev.select_epoch4_taxa

def main():
    parser=argparse.ArgumentParser()
    stages=parser.add_subparsers(dest="stage",required=True)
    for stage,fields in (("source",("atlas","output")),("audit",("source","worker","output")),
                         ("append",("atlas","source","worker","audit","output"))):
        p=stages.add_parser(stage)
        for name in fields: p.add_argument("--"+name,required=True)
        if stage=="source":p.add_argument("--offset",type=int,default=1920)
    a=parser.parse_args()
    if a.stage=="source":m.source(a.atlas,a.output,quota=128,offset=a.offset)
    elif a.stage=="audit":m.audit(a.source,a.worker,a.output)
    else:m.append(a.atlas,a.source,a.worker,a.audit,a.output)

if __name__=="__main__":
    main()
