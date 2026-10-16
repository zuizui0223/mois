#!/usr/bin/env python3
"""Prepare two independent season/colour/climate-blinded review forms for 400 MOIS photos.

Photos are ONLY source opportunities and remain UNREVIEWED. This script neither
retrieves photo pixels nor assigns blooming status or visible petal hue.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

SOURCE_SHA="b961b50ee6efe9ec08a5575a8eb1e34b1edf7b7fcb426226fb11ee901cfd3329"
TOKEN_SEED="mois-20261017-400-opening-petaloid-two-botanists-blind-v1"
FIELDS=("review_token","species","original_photo_url","photo_license","attribution",
        "reviewer_id","open_flower_stage","roi_organ","notes")
SOURCE_MAP=("inat_taxon_id","species","photo_id","observation_id","observer_id",
            "observed_on","photo_license","photo_url_large","source_rank",
            "site_latcell","site_loncell","season_year","selected_season","review_token")
STAGE=("open_flower","buds_only","fruit_only","no_flower_visible","uncertain")
TISSUE=("petal_tepal","other_flower_organ","spikelet_glume_bract","fruit_seed",
        "leaf_stem","background","uncertain")


def sha(b:bytes)->str:
    return hashlib.sha256(b).hexdigest()


def token(pid:str)->str:
    return sha((TOKEN_SEED+"|"+pid).encode())[:20]


def sheet(path:Path, rows:list[dict], fields=FIELDS)->str:
    with path.open("x",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields,lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    return sha(path.read_bytes())


def make(source:Path,out:Path,*,fixture:bool=False)->dict:
    if out.exists() and any(out.iterdir()):
        raise ValueError("Do not overwrite source-identified human reviewers' handoff")
    raw=source.read_bytes()
    if not fixture and sha(raw)!=SOURCE_SHA:
        raise ValueError("Source sample not the exact independently observer-disjoint 400 IDs")
    with source.open(newline="",encoding="utf-8") as f:
        original=list(csv.DictReader(f))
    if not fixture and len(original)!=400:
        raise ValueError("Review stage cannot replace/subsample original 400 photo IDs")
    ids=[r["photo_id"] for r in original]
    if (len(ids)!=len(set(ids)) or
        len({r["observer_id"] for r in original})!=len(original) or
        len({r["observation_id"] for r in original})!=len(original)):
        raise ValueError("Photo, observation or observer reused in the blind primary wave")
    if not all(r["photo_url_large"].startswith("https://") and
               r["photo_license"] and r["species"] and
               r["photo_id"] and r["observed_on"] for r in original):
        raise ValueError("Missing authenticated source image/identity/licence fields")
    if len({r["species"] for r in original})!=(2 if fixture else 5):
        raise ValueError("Taxon composition not the preserved blinded five-taxon source")
    review_order=sorted(original,key=lambda r:sha(
        (TOKEN_SEED+"|revieworder|"+r["photo_id"]).encode()))
    records=[]
    coordinator=[]
    for r in review_order:
        t=token(r["photo_id"])
        records.append({"review_token":t,"species":r["species"],
                        "original_photo_url":r["photo_url_large"],
                        "photo_license":r["photo_license"],"attribution":r["attribution"],
                        "reviewer_id":"","open_flower_stage":"","roi_organ":"","notes":""})
        coordinator.append({**{k:r[k] for k in SOURCE_MAP if k!="review_token"},
                            "review_token":t})
    if len({r["review_token"] for r in records})!=len(records):
        raise ValueError("Collision in blinded nonidentifying review tokens")
    out.mkdir(parents=True,exist_ok=True)
    summaries={}
    for who in ("A","B"):
        order=sorted(records,key=lambda r:sha(
            (f"mois-20261017-independent-review-{who}|"+r["review_token"]).encode()))
        name=f"mois_400_botanical_open_flower_reviewer_{who}_UNFILLED.csv"
        summaries[who]=sha_name=sheet(out/name,order)
    map_name="mois_400_review_token_source_mapping_COORDINATOR_ONLY.csv"
    map_sha=sheet(out/map_name,coordinator,SOURCE_MAP)
    metadata={
        "schema":"mois_five_taxa_independent_open_flower_handoff_v1",
        "status":"400_ORIGINAL_PHOTO_REVIEW_FORMS_BOTH_UNFILLED_NOT_COLOUR",
        "source_photo_manifest_SHA256":sha(raw),
        "n_original_photos":len(original),
        "distinct_observers":len(original),
        "n_original_named_taxa":len({r["species"] for r in original}),
        "reviewer_form_A_SHA256":summaries["A"],
        "reviewer_form_B_SHA256":summaries["B"],
        "coordinator_only_photo_to_token_mapping_SHA256":map_sha,
        "reviewer_forms_season_date_climate_colour_blind":True,
        "reviewer_form_columns":list(FIELDS),
        "allowed_review_stage_codes":list(STAGE),
        "allowed_roi_organ_codes":list(TISSUE),
        "all_review_response_fields_blank":True,
        "n_actual_independent_botanical_reviews_completed":0,
        "n_true_open_flower_photos_confirmed":0,
        "new_photo_pixels_measured":0,
        "seasonal_white_fraction_estimated":False,
        "human_responses_must_be_independent_before_adjudication":True,
        "source_only_not_a_globally_representative_flower_sample":True,
    }
    (out/"MOIS_400_INDEPENDENT_REVIEW_HANDoFF_RECEIPT.json").write_text(
        json.dumps(metadata,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    return metadata


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--source-original-400",type=Path,required=True)
    parser.add_argument("--outdir",type=Path,required=True)
    a=parser.parse_args()
    print("MOIS_400_SOURCE_ONLY_BOTANIST_BLINDED_FORMS",
          json.dumps(make(a.source_original_400,a.outdir),sort_keys=True))
