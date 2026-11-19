#!/usr/bin/env python3
"""MOIS epoch9 offset1664: next 128 distinct zero-ROI taxa, after signed offset1536.

Uses the same certified 1,231,357 photo-ID / independent observation-UUID history
and only the 113 known recent raw SHA values; never claims global raw SHA uniqueness.
"""
import argparse
import hashlib
import json
from pathlib import Path
from scripts.acquisition import mois_epoch9_after_complete_epoch8_new_created_photo9_offset1536 as prev

m = prev.m
m.SCHEMA = "mois_incremental_20261119_epoch9_offset1664_post_terminal_offset1536_1231357_original9_v1"
m.SOURCE_ARM = "MOIS_EPOCH9_OFFSET1664_AFTER_TERMINAL_OFFSET1536_NEW_CREATED"
m.MIN_CREATED = "2026-11-18T19:07:09+00:00"
m.CREATED_D1 = "2026-11-18"
m.MAX_SOURCE_TAXA = 128
m.ATLAS_SHA = "e9091728b95f556a95c1e65253a9f407688f2ed7131e4cf606c2bfe8e81b9393"
m.OLD_PHOTOS = 1231357
m.OLD_NUMERIC_PHOTOS = 467568
m.OLD_NUMERIC_TAXA = 39914
m.OLD_CLEAR_TAXA = 39174

def select_taxa(atlas, n=128, offset=1664):
    m.require(len(atlas)==42111 and atlas.inat_taxon_id.is_unique,
              "atlas cardinality / unique iNaturalist taxon IDs drift")
    missing=atlas.loc[atlas.after_strict9.eq(0)].copy()
    m.require(len(missing)==2197 and missing.inat_taxon_id.str.fullmatch(r"\d+").all(),
              "unmeasured taxon selection population changed")
    missing["priority"]=missing.inat_taxon_id.map(
        lambda tid: hashlib.sha256(
            ("MOIS-20261118-EPOCH8-AFTER-EPOCH7-2197-TEMPORAL-NEW-CREATED|"+tid).encode()
        ).hexdigest())
    ordered=missing.sort_values("priority",kind="stable").inat_taxon_id.tolist()
    picked=ordered[offset:offset+n]
    m.require((n,offset,len(picked),len(set(picked)))==(128,1664,128,128)
              and not set(picked).intersection(ordered[:offset]),
              "epoch9 next cohort is nonexclusive or incomplete")
    return picked

m.select_taxa=select_taxa

def main():
    p=argparse.ArgumentParser()
    sub=p.add_subparsers(dest="stage",required=True)
    for stage,keys in (("source",("atlas","output")),
                       ("audit",("source","worker","output")),
                       ("append",("atlas","source","worker","audit","output"))):
        q=sub.add_parser(stage)
        for k in keys:q.add_argument("--"+k,required=True)
        if stage=="source":q.add_argument("--offset",type=int,default=1664)
    args=p.parse_args()
    if args.stage=="source":
        m.source(args.atlas,args.output,quota=128,offset=args.offset)
        receipt=Path(args.output)/"receipt.json"
        r=json.loads(receipt.read_bytes())
        m.require(
            (r["hash_disjoint_taxon_offset"],r["searched_zero_numeric_taxa"],
             r["source_API_error_count"],r["old_original_photo_attempt_ID_count"],
             r["old_original_photo_obs_ID_overlap"])==(1664,128,0,1231357,0)
            and r["old_all_image_raw_SHA_disjointness_proven"] is False
            and r["created_utc_lower_bound"]==m.MIN_CREATED
            and r["original_pixels_fetched"]==r["strict9_numeric_photos"]==0,
            "source-only receipt drift")
        r["recent_known_original_image_raw_SHA_count"]=113
        r["preceding_epoch9_offset1536_run"]=44467964373
        r["preceding_epoch9_offset1536_source_artifact"]=14502659671
        r["preceding_epoch9_offset1536_source_ZIP_SHA256"]="f14347c1431100f25f179e89bfaff119fd9973d71b9d0b43cde6ebbdc2ad97f2"
        r["preceding_epoch9_offset1536_source_receipt_SHA256"]="6ef01785f23b6fad72630a0188044638052307eb2d6bba37484c5fc814faa110"
        r["preceding_epoch9_offset1536_source_status_SHA256"]="4a8914d379ebb7dc7431945b18db8855abf45ca129e387d6e8d97c98dd6db939"
        r["preceding_epoch9_offset1536_new_source_candidates"]=0
        receipt.write_text(json.dumps(r,sort_keys=True,indent=2)+"\n")
        print("EPOCH9_OFFSET1664_SOURCE_FROZEN",json.dumps(r,sort_keys=True),flush=True)
    elif args.stage=="audit":
        m.audit(args.source,args.worker,args.output)
    else:
        m.append(args.atlas,args.source,args.worker,args.audit,args.output)

if __name__=="__main__":
    main()
