#!/usr/bin/env python3
"""Post-offset896: new attempted original retained as technical unknown, source history rebased.

The signed 42,111-taxon atlas and 1,231,343 independent photo-ID attempts,
including Quercus kerrii (technical unknown), must be present before source.
The epoch4 2,199-taxon hash order stays fixed; source selections are disjoint.
No assertion of all-history raw-image SHA uniqueness or true petal colour.
"""
import argparse
from pathlib import Path
from scripts.acquisition import mois_epoch4_postfull2199_created_photo9 as prev
from scripts.acquisition import mois_epoch3_post1792_rebased_photo9 as base

m=prev.m
NEW_ATLAS_SHA="6c3aa46ba089649774b81d80a0fb8fbc49c33d71684a9058022eda5751d29309"
NEW_DELTA_SHA="1daa87b330c4664940293b11ca4e9d2a5bd9f7a3ee21e4aa44124206bc4a1cd8"
NEW_DELTA_PATH="old/epoch4_offset896/new_real_original_nine_photo9_delta.parquet"
PRIOR_EPOCH3_DELTA="old/epoch3_offset1792/new_real_original_nine_photo9_delta.parquet"
PRIOR_EPOCH3_SHA="537aa3b8419411ec82cecfc6a9768cfb313d65364a4b85afc1edb2953211cfb9"

m.SCHEMA="mois_incremental_20261116_epoch4_rebased_after_offset896_original_photo9_v2"
m.SOURCE_ARM="MOIS_EPOCH4_POST896_SIGNED_UNMEASURED_ORIGINAL_ID_DISJOINT"
m.ATLAS_SHA=NEW_ATLAS_SHA
m.OLD_PHOTOS=1231343
m.OLD_NUMERIC_PHOTOS=467566
m.OLD_NUMERIC_TAXA=39912
m.OLD_CLEAR_TAXA=39173
m.MAX_SOURCE_TAXA=128
m.DELTA_FILES=m.DELTA_FILES+(NEW_DELTA_PATH,)
m.DELTA_N=m.DELTA_N+(1,)

# The original complete-history reader, before epoch3 installed a 1,231,342 check:
_original_history=base._previous_history

def history():
    m.require(Path(NEW_DELTA_PATH).is_file() and m.sha256(NEW_DELTA_PATH)==NEW_DELTA_SHA,
              "offset896 independent audit delta missing/modified")
    m.require(Path(PRIOR_EPOCH3_DELTA).is_file() and
              m.sha256(PRIOR_EPOCH3_DELTA)==PRIOR_EPOCH3_SHA,
              "epoch3 offset1792 full-history independent delta missing/modified")
    ids,observations,recent_raw=_original_history()
    m.require(len(ids)==1231343 and len(recent_raw)==99 and
              "763837356" in ids and
              "993e4dfc-a8c4-4162-9328-d7facffee028" in observations and
              "aa82fdc048f13a6d4f6729cd2bd841d382ffb618ba33db07031647c0a3f51faf" in recent_raw and
              "763768595" in ids and
              "821d0884-737a-4644-8f8e-6e36c278b9fb" in observations,
              "post896 all historical photo/UUID/new raw SHA chain incomplete")
    return ids,observations,recent_raw

m.history=history
m.select_taxa=prev.select_epoch4_taxa

def main():
    parser=argparse.ArgumentParser()
    stages=parser.add_subparsers(dest="stage",required=True)
    for stage,req in (("source",("atlas","output")),("audit",("source","worker","output")),("append",("atlas","source","worker","audit","output"))):
        p=stages.add_parser(stage)
        for name in req: p.add_argument("--"+name,required=True)
        if stage=="source":p.add_argument("--offset",type=int,default=1024)
    args=parser.parse_args()
    if args.stage=="source":
        m.source(args.atlas,args.output,quota=128,offset=args.offset)
    elif args.stage=="audit":
        m.audit(args.source,args.worker,args.output)
    else:
        m.append(args.atlas,args.source,args.worker,args.audit,args.output)

if __name__=="__main__":
    main()
