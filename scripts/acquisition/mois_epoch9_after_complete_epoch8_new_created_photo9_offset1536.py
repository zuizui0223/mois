#!/usr/bin/env python3
"""MOIS epoch9 offset1536: 128 never-before-queried taxa, post-epoch8 created observations.

Use only authenticated original pixels from disjoint photo ID / observation UUID
history; no raw-SHA guarantee for the complete 1.23M old images.
"""
import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from scripts.acquisition import mois_epoch9_after_complete_epoch8_new_created_photo9_offset1408 as prior

m=prior.m
m.SCHEMA="mois_incremental_20261119_epoch9_offset1536_post_terminal_offset1408_1231357_original9_v1"
m.SOURCE_ARM="MOIS_EPOCH9_OFFSET1536_AFTER_TERMINAL_OFFSET1408_NEW_CREATED"
m.MIN_CREATED="2026-11-18T19:07:09+00:00"
m.CREATED_D1="2026-11-18"
m.MAX_SOURCE_TAXA=128
m.ATLAS_SHA="e9091728b95f556a95c1e65253a9f407688f2ed7131e4cf606c2bfe8e81b9393"
m.OLD_PHOTOS=1231357
m.OLD_NUMERIC_PHOTOS=467568
m.OLD_NUMERIC_TAXA=39914
m.OLD_CLEAR_TAXA=39174

def history_after_epoch9_offset1408():
    old_photo_count, old_atlas_sha = m.OLD_PHOTOS, m.ATLAS_SHA
    try:
        m.OLD_PHOTOS=1231356
        m.ATLAS_SHA="bc646ee29c02fa02111e073f096b2c7a63594d7294d7e658c2c7d6e8424bc154"
        ids, observations, raw = prior.history_after_epoch8_offset1920()
    finally:
        m.OLD_PHOTOS, m.ATLAS_SHA=old_photo_count, old_atlas_sha
    m.require(len(ids)==1231356 and len(raw)==112,
              "prior full photo/observation-ID and recent raw-image-SHA history drift")
    path=Path("old/epoch9_offset1408/new_real_original_nine_photo9_delta.parquet")
    expected_sha="85fa012a826f73fc08b5de147f9897d15583a79bab62069af236500c29ee3de1"
    m.require(path.is_file() and m.sha256(path)==expected_sha,
              "signed independently audited offset1408 delta missing or changed")
    delta=pd.read_parquet(path)
    m.require(len(delta)==1 and
              {"photo_id","observation_id","image_sha256","strict9","inat_taxon_id"}.issubset(delta.columns),
              "invalid prior independent audited delta schema")
    row=delta.iloc[0]
    pid,oid,h,tid=(str(row.photo_id),str(row.observation_id),
                   str(row.image_sha256),str(row.inat_taxon_id))
    m.require(pid=="765143736" and
              oid=="852cbc3b-5d7d-4239-82c8-e95ef2d15bd8" and
              h=="c9c99ffe4846e5c3428a36700c827bfd83531c7530b58d8107b09c5dfdb97e64" and
              tid=="209990" and not bool(row.strict9) and
              pid not in ids and oid not in observations and h not in raw,
              "prior source photo ID, UUID, raw-SHA, taxon or strict9 status drift")
    m.require(all(pd.isna(row["strict9_"+k]) for k in m.PALETTE),
              "technical unknown must have nine missing numerical palette components")
    ids.add(pid); observations.add(oid); raw.add(h)
    m.require(len(ids)==1231357 and len(raw)==113,
              "new signed historical photo/raw-SHA cardinality mismatch")
    return ids,observations,raw

def select_taxa(atlas,n=128,offset=1536):
    m.require(len(atlas)==42111 and atlas.inat_taxon_id.is_unique,
              "frozen atlas cardinality/uniqueness mismatch")
    missing=atlas.loc[atlas.after_strict9.eq(0)].copy()
    m.require(len(missing)==2197 and missing.inat_taxon_id.str.fullmatch(r"\d+").all(),
              "frozen unmeasured taxon set changed")
    missing["priority"]=missing.inat_taxon_id.map(
        lambda t: hashlib.sha256(
          ("MOIS-20261118-EPOCH8-AFTER-EPOCH7-2197-TEMPORAL-NEW-CREATED|"+t).encode()
        ).hexdigest())
    ordered=missing.sort_values("priority",kind="stable").inat_taxon_id.tolist()
    selection=ordered[offset:offset+n]
    m.require((n,offset,len(selection),len(set(selection)))==(128,1536,128,128) and
              not set(selection).intersection(ordered[:offset]),
              "priorly queried taxa or missing/disordered next128 selection")
    return selection

m.history=history_after_epoch9_offset1408
m.select_taxa=select_taxa

def main():
    p=argparse.ArgumentParser()
    sub=p.add_subparsers(dest="stage",required=True)
    for stage,keys in (("source",("atlas","output")),
                       ("audit",("source","worker","output")),
                       ("append",("atlas","source","worker","audit","output"))):
        q=sub.add_parser(stage)
        for k in keys:q.add_argument("--"+k,required=True)
        if stage=="source":q.add_argument("--offset",type=int,default=1536)
    args=p.parse_args()
    if args.stage=="source":
        m.source(args.atlas,args.output,quota=128,offset=args.offset)
        receipt=Path(args.output)/"receipt.json"
        r=json.loads(receipt.read_bytes())
        m.require(
            (r["hash_disjoint_taxon_offset"],r["searched_zero_numeric_taxa"],
             r["source_API_error_count"],r["old_original_photo_attempt_ID_count"],
             r["old_original_photo_obs_ID_overlap"])==(1536,128,0,1231357,0)
            and r["old_all_image_raw_SHA_disjointness_proven"] is False
            and r["created_utc_lower_bound"]==m.MIN_CREATED
            and r["original_pixels_fetched"]==r["strict9_numeric_photos"]==0,
            "new source-only receipt does not satisfy disjoint/source-only contract")
        r["recent_known_original_image_raw_SHA_count"]=113
        r["preceding_epoch9_offset1408_run"]=44465833745
        r["preceding_epoch9_offset1408_source_artifact"]=14501094720
        r["preceding_epoch9_offset1408_worker_artifact"]=14502466602
        r["preceding_epoch9_offset1408_independent_audit_artifact"]=14502092599
        r["preceding_epoch9_offset1408_formal_append_artifact"]=14502741240
        r["preceding_epoch9_offset1408_source_ZIP_SHA256"]="e358fb493073c2e35a59bdc829a053ac3cc4dfacbd2a57a46a4c8bb86e7fb31e"
        r["preceding_epoch9_offset1408_worker_ZIP_SHA256"]="a5ebe8b0fc4461897a7d81544355b50594cc34bb3a09b19f0b6ca95416d1b5ec"
        r["preceding_epoch9_offset1408_audit_ZIP_SHA256"]="7592113b8cc81c389ae72c14ac2a04117dc310def2a801fa82e2b646ff034b9d"
        r["preceding_epoch9_offset1408_append_ZIP_SHA256"]="7cb7ca7c0773d5f45313803e88aaa5ea6c5cf0fde1e639464dbafad700828740"
        r["preceding_epoch9_offset1408_audited_delta_SHA256"]="85fa012a826f73fc08b5de147f9897d15583a79bab62069af236500c29ee3de1"
        r["preceding_epoch9_offset1408_original_image_raw_SHA256"]="c9c99ffe4846e5c3428a36700c827bfd83531c7530b58d8107b09c5dfdb97e64"
        r["preceding_epoch9_offset1408_original_attempts_increment"]=1
        r["preceding_epoch9_offset1408_technical_unknown_increment"]=1
        r["preceding_epoch9_offset1408_strict9_increment"]=0
        receipt.write_text(json.dumps(r,sort_keys=True,indent=2)+"\n")
        print("EPOCH9_OFFSET1536_SOURCE_FROZEN",json.dumps(r,sort_keys=True),flush=True)
    elif args.stage=="audit":
        m.audit(args.source,args.worker,args.output)
    else:
        m.append(args.atlas,args.source,args.worker,args.audit,args.output)

if __name__=="__main__":
    main()
