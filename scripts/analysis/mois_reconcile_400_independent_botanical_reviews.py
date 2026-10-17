#!/usr/bin/env python3
"""MOIS: reconcile 400 source-photo botanical reviews without a colour outcome.

The two forms are blind to observation date, season, geography, humidity and
automatic photo hue; join frozen source metadata ONLY AFTER both independent
reviewers have submitted decisions. Third-person adjudication is optional but
only allowed for genuine disagreement. Review scores never replace a source
photo with a more convenient image and never manufacture flower colour.

DO NOT use this as genetic flower-colour polymorphism evidence.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

from scripts.analysis.mois_prepare_400_botanical_review_handoff import (
    FIELDS, STAGE, TISSUE, TOKEN_SEED, SOURCE_SHA, token, SOURCE_MAP,
)

ATTESTED_HANDoFF = "400_ORIGINAL_PHOTO_REVIEW_FORMS_BOTH_UNFILLED_NOT_COLOUR"
EXPECTED_SOURCE = "ORIGINAL_PHOTO_SOURCE_ONLY_NOT_A_BLOOM_OR_COLOUR"
ADJ_FIELDS = ("review_token","adjudicator_id","adjudicated_open_flower_stage",
              "adjudicated_roi_organ","adjudication_notes")
RESULT_FIELDS = ("inat_taxon_id","species","photo_id","observation_id","observer_id",
                 "observed_on","season_year","source_site_lat_cell_0p25",
                 "source_site_lon_cell_0p25","source_local_season",
                 "review_token","reviewer_A_id","reviewer_B_id",
                 "reviewer_A_stage","reviewer_A_roi_organ",
                 "reviewer_B_stage","reviewer_B_roi_organ",
                 "adjudicator_id","botanical_photo_stage_status",
                 "eligible_independently_verified_open_petaloid_photo_for_later_roi",
                 "original_RGB_photo_colour_measured_in_MOIS",
                 "botanically_verified_white_flower_hue")
SUMMARY_FIELDS = ("inat_taxon_id","species","season_year","site_lat_cell_0p25",
                 "site_lon_cell_0p25","source_local_season",
                 "source_original_photo_opportunities","reviewed_by_both_independent_botanists",
                 "reviewer_disagreement_needing_adjudication",
                 "unreviewed_or_one_reviewer_only",
                 "reviewer_agreement_nonopen_or_nonpetaloid",
                 "independently_verified_open_petaloid_eligible_for_later_ROI",
                 "genuine_original_image_ROI_pixel_colour_attempts_completed",
                 "valid_open_white_flower_colour_outcomes",
                 "valid_open_chromatic_flower_colour_outcomes")
STATUSES = (
    "awaiting_two_independent_reviews",
    "awaiting_second_independent_review",
    "reviewer_disagreement_pending",
    "botanical_open_petaloid_agreed",
    "botanical_open_petaloid_adjudicated",
    "other_open_flower_organ_not_primary",
    "nonopen_bud_fruit_or_no_visible_flower",
    "noncomparable_bract_seed_leaf_background",
    "uncertain_flower_or_tissue",
)


def sha(file:Path)->str:
    return hashlib.sha256(file.read_bytes()).hexdigest()


def read_rows(file:Path,required_header:tuple[str,...])->list[dict]:
    with file.open(newline="",encoding="utf-8") as h:
        reader=csv.DictReader(h)
        if tuple(reader.fieldnames or ()) != required_header:
            raise ValueError("Unexpected review or provenance columns (possible season/colour leakage)")
        return list(reader)


def review_state(r:dict,who:str)->tuple[str,str,str]:
    reviewer_id=r["reviewer_id"].strip()
    stage=r["open_flower_stage"].strip()
    tissue=r["roi_organ"].strip()
    values=(reviewer_id,stage,tissue)
    if not any(values):
        return ("","","")
    if not all(values):
        raise ValueError("Partial reviewer reply; reviewer ID, stage and ROI organ required together: "+who)
    if stage not in STAGE or tissue not in TISSUE:
        raise ValueError("Unregistered botanical stage or ROI organ label")
    if len(reviewer_id)>128 or "\n" in reviewer_id:
        raise ValueError("Unsafe botanical reviewer identity")
    return values


def classify(stage:str,tissue:str,*,adjudicated:bool)->str:
    if not stage or not tissue:
        raise ValueError("Botanical stage must be fully specified before classification")
    if stage=="uncertain" or tissue=="uncertain":
        return "uncertain_flower_or_tissue"
    if stage!="open_flower":
        return "nonopen_bud_fruit_or_no_visible_flower"
    if tissue=="petal_tepal":
        return ("botanical_open_petaloid_adjudicated" if adjudicated
                else "botanical_open_petaloid_agreed")
    if tissue=="other_flower_organ":
        return "other_open_flower_organ_not_primary"
    return "noncomparable_bract_seed_leaf_background"


def reconcile(
    photo_source:Path, photo_selection_receipt:Path,
    form_a:Path, form_b:Path, review_handoff_receipt:Path,
    out:Path, *, adjudication:Path|None=None, _synthetic_test:bool=False,
)->dict:
    if out.exists() and any(out.iterdir()):
        raise ValueError("Never overwrite previously sealed botanical review results")
    source_doc=json.loads(photo_selection_receipt.read_text(encoding="utf-8"))
    handoff_doc=json.loads(review_handoff_receipt.read_text(encoding="utf-8"))
    if (source_doc.get("status")!=EXPECTED_SOURCE
        or source_doc.get("selected_photo_manifest_sha256")!=sha(photo_source)
        or handoff_doc.get("status")!=ATTESTED_HANDoFF
        or handoff_doc.get("source_photo_manifest_SHA256")!=sha(photo_source)
        or handoff_doc.get("reviewer_forms_season_date_climate_colour_blind") is not True
        or handoff_doc.get("n_actual_independent_botanical_reviews_completed")!=0
        or handoff_doc.get("new_photo_pixels_measured")!=0):
        raise ValueError("Original source, season-blind review forms or source SHA are not attested")
    if not _synthetic_test and (
        sha(photo_source)!=SOURCE_SHA or
        handoff_doc.get("n_original_photos")!=400 or
        source_doc.get("selected_original_photo_ids")!=400 or
        source_doc.get("selected_distinct_source_observers")!=400 or
        handoff_doc.get("n_original_named_taxa")!=5):
        raise ValueError("Real source can only admit the exact pre-colour 400 photo IDs / five taxa")

    # Historical source roster's exact field order and original photo-ID
    # namespace, observer identifiers and dates are never changed by reviews.
    with photo_source.open(newline="",encoding="utf-8") as h:
        inp=csv.DictReader(h)
        required={
            "inat_taxon_id","species","photo_id","observation_id","observer_id",
            "observed_on","photo_license","photo_url_large","attribution",
            "site_latcell","site_loncell","season_year","selected_season",
        }
        if not required.issubset(inp.fieldnames or []):
            raise ValueError("Original photo-only seasonal source schema incomplete")
        source=list(inp)
    if len(source)!=(4 if _synthetic_test else 400):
        raise ValueError("Not all original source photos are present")
    if len({r["photo_id"] for r in source})!=len(source) or (
       len({r["observation_id"] for r in source})!=len(source) or
       len({r["observer_id"] for r in source})!=len(source)):
        raise ValueError("One original photo/observation/observer reused across seasons")
    token_source={}
    species_season=Counter()
    for r in source:
        pid=r["photo_id"]
        if not pid.isascii() or not pid.isdecimal():
            raise ValueError("Original photo ID invalid")
        t=token(pid)
        if t in token_source:
            raise ValueError("Blind review-token collision from original photo IDs")
        if r["selected_season"] not in ("winter","spring","summer","autumn") or (
            not r["season_year"].isdigit() or not r["species"] or not r["inat_taxon_id"].isdigit()):
            raise ValueError("Frozen source season/year/taxon identity invalid")
        if not re.fullmatch(r"https://inaturalist-open-data[.]s3[.]amazonaws[.]com/photos/"
                            +re.escape(pid)+r"/large[.](?:jpg|jpeg|png)",r["photo_url_large"],re.I):
            raise ValueError("Original licensed photo URL no longer matches original photo ID")
        token_source[t]=r
        species_season[(r["inat_taxon_id"],r["selected_season"])]+=1
    if not _synthetic_test and (
        len(species_season)!=10 or set(species_season.values())!={40} or
        len({r["inat_taxon_id"] for r in source})!=5):
        raise ValueError("Exact source 5 taxa × 2 seasons ×40 original opportunities changed")

    a=read_rows(form_a,FIELDS)
    b=read_rows(form_b,FIELDS)
    forms={}
    reviewer_labels={}
    for who,entries in (("A",a),("B",b)):
        if len(entries)!=len(source):
            raise ValueError("Reviewer "+who+" did not retain every source photo including unreviewed")
        maps={}
        ids=set()
        for row in entries:
            t=row["review_token"]
            if t in maps or t not in token_source:
                raise ValueError("Duplicate or non-original blind photo token in reviewer "+who)
            orig=token_source[t]
            if any((
                 row["species"]!=orig["species"],
                 row["original_photo_url"]!=orig["photo_url_large"],
                 row["photo_license"]!=orig["photo_license"],
                 row["attribution"]!=orig["attribution"],
            )):
                raise ValueError("Original photo metadata changed between signed reviewer sheets")
            u,stage,organ=review_state(row,who)
            if u:
                ids.add(u)
            maps[t]=(u,stage,organ)
        if set(maps)!=set(token_source):
            raise ValueError("Reviewer "+who+" form omits original photos")
        if len(ids)>1:
            raise ValueError("Reviewer "+who+" form has multiple identities; not one independent reviewer")
        forms[who]=maps
        reviewer_labels[who]=next(iter(ids),"")
    if reviewer_labels["A"] and reviewer_labels["B"] and reviewer_labels["A"]==reviewer_labels["B"]:
        raise ValueError("Reviewers A and B must be independent humans")

    adjudication_map={}
    if adjudication is not None:
        for row in read_rows(adjudication,ADJ_FIELDS):
            t=row["review_token"]
            if t in adjudication_map or t not in token_source:
                raise ValueError("Duplicate/alien photo adjudication token")
            who,stage,organ=row["adjudicator_id"].strip(),row["adjudicated_open_flower_stage"].strip(),row["adjudicated_roi_organ"].strip()
            if not all((who,stage,organ)) or stage not in STAGE or organ not in TISSUE:
                raise ValueError("Adjudicated botanical state is incomplete or invalid")
            if who in (reviewer_labels["A"],reviewer_labels["B"]):
                raise ValueError("Third-party adjudicator cannot be either original reviewer")
            adjudication_map[t]=(who,stage,organ)
    output=[]
    status_counts=Counter()
    by_group=defaultdict(Counter)
    for t,orig in sorted(token_source.items(),key=lambda pair:int(pair[1]["photo_id"])):
        A=forms["A"][t]
        B=forms["B"][t]
        if not A[0] and not B[0]:
            status="awaiting_two_independent_reviews"
        elif not A[0] or not B[0]:
            status="awaiting_second_independent_review"
        elif A[1:]==B[1:]:
            status=classify(A[1],A[2],adjudicated=False)
        else:
            status="reviewer_disagreement_pending"
        adjudicator=""
        if t in adjudication_map:
            if status!="reviewer_disagreement_pending":
                raise ValueError("Third reviewer cannot alter agreeing, absent or partial reviews")
            adjudicator,stage,organ=adjudication_map[t]
            status=classify(stage,organ,adjudicated=True)
        eligible=status in ("botanical_open_petaloid_agreed","botanical_open_petaloid_adjudicated")
        status_counts[status]+=1
        key=(orig["inat_taxon_id"],orig["species"],orig["season_year"],
             orig["site_latcell"],orig["site_loncell"],orig["selected_season"])
        counters=by_group[key]
        counters["n"]+=1
        counters["both"]+=int(bool(A[0] and B[0]))
        counters["pending"]+=int(status=="reviewer_disagreement_pending")
        counters["partial"]+=int(status in ("awaiting_two_independent_reviews","awaiting_second_independent_review"))
        counters["noneligible"]+=int(status in (
             "other_open_flower_organ_not_primary",
             "nonopen_bud_fruit_or_no_visible_flower",
             "noncomparable_bract_seed_leaf_background",
             "uncertain_flower_or_tissue",
        ))
        counters["eligible"]+=int(eligible)
        output.append(dict(zip(RESULT_FIELDS,(
            orig["inat_taxon_id"],orig["species"],orig["photo_id"],
            orig["observation_id"],orig["observer_id"],orig["observed_on"],
            orig["season_year"],orig["site_latcell"],orig["site_loncell"],
            orig["selected_season"],t,A[0],B[0],A[1],A[2],B[1],B[2],
            adjudicator,status,eligible,False,False
        ))))
    if len(output)!=len(source) or sum(status_counts.values())!=len(source):
        raise ValueError("Botanical photo-state denominator/identity not complete")

    out.mkdir(parents=True,exist_ok=True)
    results_csv=out/"mois_400_source_photo_botanical_qualification_COLOUR_UNMEASURED.csv"
    with results_csv.open("x",encoding="utf-8",newline="") as h:
        w=csv.DictWriter(h,fieldnames=RESULT_FIELDS,lineterminator="\n")
        w.writeheader();w.writerows(output)
    summary=out/"mois_400_taxon_siteyear_season_reviewer_sample_flow.csv"
    with summary.open("x",encoding="utf-8",newline="") as h:
        w=csv.DictWriter(h,fieldnames=SUMMARY_FIELDS,lineterminator="\n")
        w.writeheader()
        for key,c in sorted(by_group.items()):
            w.writerow(dict(zip(SUMMARY_FIELDS,(
                 *key,c["n"],c["both"],c["pending"],c["partial"],
                 c["noneligible"],c["eligible"],0,0,0
            ))))
    full=(len(output) if _synthetic_test else 400)
    no_pending=not any(status_counts[s] for s in (
        "awaiting_two_independent_reviews","awaiting_second_independent_review",
        "reviewer_disagreement_pending"))
    status=("BOTANICAL_QUALIFICATION_COMPLETE_BUT_COLOUR_NOT_MEASURED"
            if no_pending else "INCOMPLETE_OR_NOT_STARTED_BOTANICAL_REVIEW_NOT_WHITE_RATE")
    report={
        "schema":"mois_source_only_400_photo_two_blind_botanists_adjudication_v1",
        "status":status,
        "source_original_photo_manifest_sha256":sha(photo_source),
        "photo_selection_receipt_sha256":sha(photo_selection_receipt),
        "original_blind_review_handoff_receipt_sha256":sha(review_handoff_receipt),
        "reviewer_A_form_SHA256":sha(form_a),
        "reviewer_B_form_SHA256":sha(form_b),
        "third_party_adjudication_form_SHA256":sha(adjudication) if adjudication else None,
        "selected_original_photo_opportunities":full,
        "independent_source_observations":len(source),
        "independent_source_photographers":len({r["observer_id"] for r in source}),
        "source_siteyear_season_groups":len(by_group),
        "n_botanist_A_reviewed_photos":sum(bool(forms["A"][t][0]) for t in token_source),
        "n_botanist_B_reviewed_photos":sum(bool(forms["B"][t][0]) for t in token_source),
        "n_reviewed_by_both_independently":sum(c["both"] for c in by_group.values()),
        "botanical_dispositions":dict(sorted(status_counts.items())),
        "n_confirmed_open_petaloid_photos_FOR_LATER_ROI":sum(status_counts[s] for s in (
            "botanical_open_petaloid_agreed","botanical_open_petaloid_adjudicated")),
        "new_colour_ROI_measurements_from_these_400":0,
        "independently_confirmed_white_flower_COLOUR_COUNT":0,
        "independently_confirmed_chromatic_flower_COLOUR_COUNT":0,
        "not_yet_measured_colour_is_not_nonwhite":True,
        "photographs_with_unknown_botanical_stage_are_not_nonflowering":True,
        "reviewer_inputs_exclude_month_season_geo_climate_and_pixel_hue":True,
        "photo_colour_white_frequency_estimated":False,
        "climate_seasonal_flower_colour_effect_estimated":False,
        "sample_flow_csv_sha256":sha(summary),
        "photo_qualification_csv_sha256":sha(results_csv),
    }
    (out/"mois_400_two_botanist_photo_review_reconciliation_receipt.json").write_text(
         json.dumps(report,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    return report


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source400",type=Path,required=True)
    p.add_argument("--source400-receipt",type=Path,required=True)
    p.add_argument("--reviewer-A",type=Path,required=True)
    p.add_argument("--reviewer-B",type=Path,required=True)
    p.add_argument("--review-handoff-receipt",type=Path,required=True)
    p.add_argument("--third-adjudication",type=Path)
    p.add_argument("--outdir",type=Path,required=True)
    args=p.parse_args(argv)
    print(json.dumps(reconcile(
        args.source400,args.source400_receipt,args.reviewer_A,args.reviewer_B,
        args.review_handoff_receipt,args.outdir,adjudication=args.third_adjudication
    ),sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
