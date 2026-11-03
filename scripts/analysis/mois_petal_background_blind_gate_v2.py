#!/usr/bin/env python3
"""Metadata-minimized MOIS 128 botanical-review screen and private readout v2.

This supplements v1 without changing historic original ROI9 or creating
botanical truth. Real human responses remain local and are never sent to CI.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path

from scripts.analysis import mois_petal_background_contrast_review_v1 as base

SCHEMA="mois_existing128_petal_background_metadata_minimized_review_v2"
BLIND_COLUMNS=("photo_id","source_roster_sha256","reviewer_slot","reviewer_id",
    "flower_visible","petal_visible","background_visible",
    "photo_file_sha256","image_width","image_height","image_file_name",
    "fg_points_json","bg_points_json","fg_r","fg_g","fg_b",
    "bg_r","bg_g","bg_b","deltaE76","deltaLstar",
    "sRGB_luminance_ratio","notes")
FORBIDDEN=("observation_id","inat_taxon_id","source_design","latitude_zone",
           "local_season","postoutcome_source_gap","screen_status",
           "technical_diagnostic_class","review_priority")
ORIGINAL_FIELDS_PATTERN=r"const fields=\[[^\]]+\];"
OLD_STORAGE="MOIS_128_PETAL_BG_MANUAL_V1|"
NEW_STORAGE="MOIS_128_PETAL_BG_MANUAL_V2|"


def preflight(rows):
    base.require(len(rows)>0,"Original source review roster is empty")
    ids=[x["photo_id"] for x in rows]
    base.require(len(set(ids))==len(rows) and all(x.isdecimal() for x in ids),
                 "Invalid/duplicated original photo IDs")
    def counted(field):
        return dict(sorted(Counter(str(r[field]) for r in rows).items()))
    return {"schema":SCHEMA,"photo_ids":len(rows),
            "regional_counts":counted("latitude_zone"),
            "local_season_counts":counted("local_season"),
            "source_design_counts":counted("source_design"),
            "technical_screen_counts":counted("screen_status"),
            "botanically_verified_count":0,
            "original_strict9_increment":0,
            "island_mainland_status_in_roster":False,
            "source_purposeful_ROI9_failure_sample":True,
            "ecological_white_or_dull_flower_rate_estimated":False,
            "nothing_read_from_new_photo_bytes":True}


def build_blind(rows, template, output, *, source_sha=base.ROSTER_SHA):
    """Only photo IDs enter reviewer HTML; locked coordinator roster stays separate."""
    output=Path(output)
    base.require(not output.exists(),"Do not overwrite existing reviewer packet")
    base.require(len(set(x["photo_id"] for x in rows))==len(rows),
                 "Duplicate original source photo IDs")
    base.require(len({x["observation_id"] for x in rows})==len(rows) and
                 len({x["inat_taxon_id"] for x in rows})==len(rows),
                 "Repeat original observation/taxon in source")
    base.require(len(source_sha)==64 and all(x in "0123456789abcdef" for x in source_sha),
                 "Missing authentic source roster hash")
    markup=Path(template).read_text(encoding="utf-8")
    base.require(markup.count("__SOURCE_ROWS__")==1 and
                 markup.count("__SOURCE_SHA__")==1 and
                 "__REVIEWER_SLOT__" in markup,"Source review HTML template changed")
    markup,n=re.subn(ORIGINAL_FIELDS_PATTERN,'const fields=["photo_id"];',markup)
    base.require(n==1 and markup.count(OLD_STORAGE)==1,
                 "Cannot minimize original pre-review metadata fields")
    markup=markup.replace(OLD_STORAGE,NEW_STORAGE)
    for word in FORBIDDEN:
        base.require(word not in markup,"Human reviewer template exposes "+word)
    output.mkdir(parents=True)
    digests={}
    for slot in "AB":
        ordered=sorted(rows,key=lambda r:hashlib.sha256(
            ("MOIS_128_CONTRAST_REVIEW_BLIND_"+slot+"|"+r["photo_id"]).encode()).hexdigest())
        visible=[{"photo_id":r["photo_id"]} for r in ordered]
        json_rows=json.dumps(visible,ensure_ascii=False,separators=(",",":"))
        json_rows=json_rows.replace("<","\\u003c").replace(">","\\u003e").replace("&","\\u0026")
        html=markup.replace("__SOURCE_ROWS__",json_rows)
        html=html.replace("__SOURCE_SHA__",source_sha).replace("__REVIEWER_SLOT__",slot)
        for word in FORBIDDEN:
            base.require(word not in html,"Reviewer output contains identity metadata "+word)
        file=output/("mois_blind_contrast_reviewer_"+slot+".html")
        file.write_text(html,encoding="utf-8")
        digests[file.name]=base.sha(file)
    receipt={"schema":SCHEMA,"source_roster_sha256":source_sha,
             "reviewer_photo_ids":len(rows),"independent_reviewer_packets":2,
             "included_source_fields":["photo_id"],
             "no_region_or_model_metadata_embedded_in_reviewer_html":True,
             "not_fully_blind_to_recognizable_photo_or_original_photo_id":True,
             "original_strict9_unknown_unchanged":True,
             "botanical_verified_photos":0,"contrast_measurements_from_real_photos":0,
             "no_source_photo_reacquired_or_inferred":True,
             "source_data_remain_coordinator_only":True,
             "reviewer_exports_contain_no_geography_taxonomy_or_detector_state":True,
             "original_photo_rows_are_ROI9_failure_enriched_not_floristic_sample":True,
             "outputs_sha256":digests}
    (output/"blind128_pre_review_receipt.json").write_text(
        json.dumps(receipt,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    (output/"blind128_source_support_counts.json").write_text(
        json.dumps(preflight(rows),sort_keys=True,indent=2)+"\n",encoding="utf-8")
    return receipt


def blank_form(rows,slot,reviewer_id):
    """Logistical blank form for synthetic tests, never botanical labels."""
    base.require(slot in ("A","B"),"Invalid reviewer slot")
    return [{key:({"photo_id":r["photo_id"],"source_roster_sha256":base.ROSTER_SHA,
                  "reviewer_slot":slot,"reviewer_id":reviewer_id}.get(key,""))
             for key in BLIND_COLUMNS} for r in rows]


def expand_blind_review(original,records,slot,*,source_sha=base.ROSTER_SHA):
    base.require(slot in ("A","B"),"Invalid reviewer slot")
    ref={r["photo_id"]:r for r in original}
    base.require(len(ref)==len(original),"Original source photo duplicates")
    base.require(len(records)==len(ref),"Reviewer photo count changed")
    seen=set()
    expanded=[]
    for row in records:
        base.require(set(row)==set(BLIND_COLUMNS),"Unexpected blind reviewer CSV columns")
        id=row["photo_id"]
        base.require(id in ref and id not in seen,"Blind reviewer source photo identity drift")
        seen.add(id)
        base.require(row["reviewer_slot"]==slot,"Reviewer slot changed")
        base.require(row["source_roster_sha256"]==source_sha,"Original source roster SHA changed")
        expanded.append({**ref[id],**row})
    base.require(seen==set(ref),"Reviewer source roster differs from frozen photos")
    base.require(all(set(r)==set(base.REVIEW_COLUMNS) for r in expanded),
                 "Expanded private review has invalid source fields")
    return expanded


def reconcile_blind(original,first,second,*,source_sha=base.ROSTER_SHA):
    a=expand_blind_review(original,first,"A",source_sha=source_sha)
    b=expand_blind_review(original,second,"B",source_sha=source_sha)
    return base.reconcile(original,a,b,expected_roster_sha=source_sha)


def aggregate_only(original, result):
    """No private photo/observation IDs or individual reviewer decisions exported."""
    n=len(original)
    success=int(result["confirmed_optical_photo_ids"])
    nonpetal=int(result["confirmed_nonpetal_cases"])
    unresolved=int(result["unresolved_or_disagreement"])
    different=int(result["different_image_bytes"])
    base.require(success+nonpetal+unresolved+different==n,
                 "Private review statuses fail count reconciliation")
    return {"schema":SCHEMA+"_aggregate",
            "source_photo_attempts":n,"verified_photo_n":success,
            "reviewer_agreed_no_petaloid_background_n":nonpetal,
            "unresolved_or_review_disagreement_n":unresolved,
            "identical_original_photo_bytes_disagreed_n":different,
            "unreviewed_or_unresolved_photo_n":unresolved+different,
            "review_coverage_fraction":(n-unresolved)/n,
            "original_strict9_success_increment":0,
            "regional_descriptive_estimates_eligible":False,
            "source_selection_conditioned_on_original_ROI9_failure":True,
            "not_native_global_flora_frequency":True,
            "optical_values_are_uncalibrated_visible_image_sRGB_not_pollinator_vision":True,
            "raw_individual_review_data_public":False,
            "requires_private_photo_file_authentication":True}


def _read_blind_csv(path):
    with Path(path).open(newline="",encoding="utf-8-sig") as stream:
        reader=csv.DictReader(stream)
        base.require(tuple(reader.fieldnames or ())==BLIND_COLUMNS,
                     "Unexpected blind CSV exported columns")
        return list(reader)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest="task",required=True)
    p=sub.add_parser("build")
    p.add_argument("--roster",type=Path,required=True)
    p.add_argument("--template",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    q=sub.add_parser("validate-private")
    q.add_argument("--roster",type=Path,required=True)
    q.add_argument("--reviewer-a",type=Path,required=True)
    q.add_argument("--reviewer-b",type=Path,required=True)
    q.add_argument("--output",type=Path,required=True)
    opts=parser.parse_args()
    source=base._load_roster(opts.roster,base.ROSTER_SHA,base.ROSTER_N)
    if opts.task=="build":
        result=build_blind(source,opts.template,opts.output)
        print("MOIS128_METADATA_MINIMIZED_REVIEW_READY",json.dumps(result,sort_keys=True))
    else:
        result=reconcile_blind(source,_read_blind_csv(opts.reviewer_a),
                              _read_blind_csv(opts.reviewer_b))
        public=aggregate_only(source,result)
        base.require(not opts.output.exists(),"Never overwrite genuine private human review")
        opts.output.mkdir(parents=True)
        # Private per-photo table stays on coordinator computer. Never upload to CI.
        if result["confirmed_records"]:
            file=opts.output/"PRIVATE_verified_photo_contrasts.csv"
            with file.open("w",newline="",encoding="utf-8") as stream:
                writer=csv.DictWriter(stream,fieldnames=list(result["confirmed_records"][0]))
                writer.writeheader();writer.writerows(result["confirmed_records"])
        (opts.output/"PRIVATE_only_aggregate.json").write_text(
            json.dumps(public,sort_keys=True,indent=2)+"\n",encoding="utf-8")
        print("MOIS128_PRIVATE_REVIEW_GATED",json.dumps(public,sort_keys=True))


if __name__=="__main__":
    main()
