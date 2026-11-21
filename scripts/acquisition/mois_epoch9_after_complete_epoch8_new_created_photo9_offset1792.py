#!/usr/bin/env python3
"""MOIS epoch9 offset1792: next 128 distinct zero-ROI taxa, after signed offset1664.

Uses the same certified 1,231,357 photo-ID / independent observation-UUID history
and only the 113 known recent raw SHA values; never claims global raw SHA uniqueness.
"""
import argparse
import hashlib
import json
from pathlib import Path
from scripts.acquisition import mois_epoch9_after_complete_epoch8_new_created_photo9_offset1664 as prev

m = prev.m
m.SCHEMA = "mois_incremental_20261119_epoch9_offset1792_post_terminal_offset1664_1231357_original9_v1"
m.SOURCE_ARM = "MOIS_EPOCH9_OFFSET1792_AFTER_TERMINAL_OFFSET1664_NEW_CREATED"
m.MIN_CREATED = "2026-11-18T19:07:09+00:00"
m.CREATED_D1 = "2026-11-18"
m.MAX_SOURCE_TAXA = 128
m.ATLAS_SHA = "e9091728b95f556a95c1e65253a9f407688f2ed7131e4cf606c2bfe8e81b9393"
m.OLD_PHOTOS = 1231357
m.OLD_NUMERIC_PHOTOS = 467568
m.OLD_NUMERIC_TAXA = 39914
m.OLD_CLEAR_TAXA = 39174

def select_taxa(atlas, n=128, offset=1792):
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
    m.require((n,offset,len(picked),len(set(picked)))==(128,1792,128,128)
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
        if stage=="source":q.add_argument("--offset",type=int,default=1792)
    args=p.parse_args()
    if args.stage=="source":
        m.source(args.atlas,args.output,quota=128,offset=args.offset)
        receipt=Path(args.output)/"receipt.json"
        r=json.loads(receipt.read_bytes())
        m.require(
            (r["hash_disjoint_taxon_offset"],r["searched_zero_numeric_taxa"],
             r["source_API_error_count"],r["old_original_photo_attempt_ID_count"],
             r["old_original_photo_obs_ID_overlap"])==(1792,128,0,1231357,0)
            and r["old_all_image_raw_SHA_disjointness_proven"] is False
            and r["created_utc_lower_bound"]==m.MIN_CREATED
            and r["original_pixels_fetched"]==r["strict9_numeric_photos"]==0,
            "source-only receipt drift")
        r["recent_known_original_image_raw_SHA_count"]=113
        r["preceding_epoch9_offset1664_run"]=44481265003
        r["preceding_epoch9_offset1664_source_artifact"]=14509181621
        r["preceding_epoch9_offset1664_source_ZIP_SHA256"]="30a8efc63fab6295ccc7d93a2c95779b10f55f6f6a9514610220323cf7ec9e63"
        r["preceding_epoch9_offset1664_source_receipt_SHA256"]="3446c3ede4c2185e1701fcf55bd137c9ab5d564b63ff5da722be88017be3ac7e"
        r["preceding_epoch9_offset1664_source_status_SHA256"]="be7de4701df8f17144f3c000052754f076a6b757690aca76a402c35c2f0e27b5"
        r["preceding_epoch9_offset1664_new_source_candidates"]=0
        receipt.write_text(json.dumps(r,sort_keys=True,indent=2)+"\n")
        print("EPOCH9_OFFSET1792_SOURCE_FROZEN",json.dumps(r,sort_keys=True),flush=True)
    elif args.stage=="audit":
        m.audit(args.source,args.worker,args.output)
    else:
        m.append(args.atlas,args.source,args.worker,args.audit,args.output)

if __name__=="__main__":
    main()
