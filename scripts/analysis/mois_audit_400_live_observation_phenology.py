#!/usr/bin/env python3
"""MOIS real 400-photo original iNaturalist observation phenology availability.

- Original selected 400 photo/observation/observer IDs are SOURCE frozen and
  outcome-blind, distinct from old MOIS 100,543 images and FCP pilot 64.
- The live iNaturalist API is used only after selection for a current
  OBSERVATION-level ancillary phenology/missingness audit. It does not tell
  whether the particular selected original photo depicts open petaloid tissue.
- It does not alter the original 400 selection, download pixels, classify
  floral colour, join weather, or estimate seasonal flower-white frequencies.
- Current API observations can differ from the dated 2026-09-27 source.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import time
import uuid
from collections import Counter,defaultdict
from datetime import date,datetime,timezone
from pathlib import Path

import requests

from scripts.analysis.mois_audit_live_inat_observation_phenology_2025 import (
    API, VALUE, AGENT, plant_status, request,
)
from scripts.analysis.mois_observer_disjoint_photo_wave_20261017 import PHOTO_COLUMNS

SELECTED_SHA="b961b50ee6efe9ec08a5575a8eb1e34b1edf7b7fcb426226fb11ee901cfd3329"
ORIGINAL_SNAPSHOT="2026-09-27"
SOURCE_STATUS="ORIGINAL_PHOTO_SOURCE_ONLY_NOT_A_BLOOM_OR_COLOUR"
ROWCOLS=(
    "original_source_photo_id","original_source_observation_UUID",
    "original_source_observed_on","original_taxon_id","original_species",
    "original_site_0p25_lat_cell","original_site_0p25_lon_cell",
    "original_source_season_year","original_local_season",
    "live_API_observation_id","live_API_UUID_matched_original",
    "live_API_observed_on_exactly_matches_frozen_source",
    "same_original_photo_ID_still_listed_on_live_observation",
    "current_plant_phenology_value_ids","current_observation_phenology_NOT_photo",
    "public_API_metadata_NOT_frozen_on_original_date",
    "original_picture_OPEN_petal_tepal_independently_reviewed",
)
GROUPCOLS=(
    "inat_taxon_id","species","source_site_lat0p25","source_site_lon0p25",
    "original_source_season_year","original_local_season",
    "original_photo_observation_opportunities","current_UUID_resolved",
    "current_photo_and_original_observed_on_both_match",
    "current_observation_flowering_NOT_photo","current_observation_fruiting_NOT_photo",
    "current_observation_flower_budding_NOT_photo",
    "current_observation_no_evidence_of_flowering_NOT_photo",
    "current_observation_phenology_missing_or_unknown",
    "current_observation_other_or_conflicting_annotation",
    "original_picture_OPEN_petal_tepal_independently_reviewed",
)
CATEGORIES={
 "current_observation_annotation_flowering_NOT_photo":"flowering",
 "current_observation_annotation_fruiting_NOT_photo":"fruiting",
 "current_observation_annotation_flower_budding_NOT_photo":"flower_budding",
 "current_observation_annotation_no_evidence_of_flowering_NOT_photo":"no_evidence",
}


def sha(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for part in iter(lambda:f.read(1024*1024),b""):h.update(part)
    return h.hexdigest()


def original_400(source:Path,receipt:Path,*,_synthetic_test:bool=False)->list[dict]:
    doc=json.loads(receipt.read_text(encoding="utf-8"))
    if (doc.get("status")!=SOURCE_STATUS or
        doc.get("selected_photo_manifest_sha256")!=sha(source) or
        doc.get("did_not_infer_flower_opening_or_white_fraction") is not True or
        doc.get("did_not_open_or_reclassify_photo_pixels") is not True):
        raise ValueError("Only original SHA-certified 400 photo-ID source opportunities allowed")
    if not _synthetic_test and (
        sha(source)!=SELECTED_SHA or
        doc.get("selected_original_photo_ids")!=400 or
        doc.get("selected_distinct_source_observers")!=400 or
        doc.get("selected_taxa")!=5):
        raise ValueError("Not exact five-taxon 400 original photo source roster")
    with source.open(encoding="utf-8",newline="") as f:
        q=csv.DictReader(f)
        if tuple(q.fieldnames or ())!=PHOTO_COLUMNS:
            raise ValueError("Original 400 photo source and blind metadata header changed")
        rows=list(q)
    expected=4 if _synthetic_test else 400
    if len(rows)!=expected or any(
        len({r[k] for r in rows})!=len(rows)
        for k in ("photo_id","observation_id","observer_id")
    ):
        raise ValueError("Original source photo/observation/observer count is incomplete or reused")
    groups=Counter()
    seen_taxa={}
    for r in rows:
        photo,observer,obs=r["photo_id"],r["observer_id"],r["observation_id"]
        if not photo.isascii() or not photo.isdecimal() or not observer or not obs:
            raise ValueError("Invalid original source identity")
        try:
            if str(uuid.UUID(obs))!=obs:
                raise ValueError("UUID not normalised in original source")
            y=date.fromisoformat(r["observed_on"]).year
            year=int(r["season_year"])
            ilat,ilon=int(r["site_latcell"]),int(r["site_loncell"])
        except (ValueError,TypeError) as e:
            raise ValueError("Original observed_on UUID/source site-year invalid") from e
        if (r["source_snapshot_date"]!=ORIGINAL_SNAPSHOT or
            not r["species"] or not r["inat_taxon_id"].isdecimal() or
            not 0<=ilat<=719 or not 0<=ilon<=1439 or
            not 2020<=y<=2026 or r["selected_season"] not in
               ("winter","spring","summer","autumn")):
            raise ValueError("Original site-year or source observation identity drift")
        # December belongs to the next meteorological season-year.
        if year!=y+int(r["observed_on"][5:7]=="12"):
            raise ValueError("Source season-year disagrees with immutable observed-on date")
        if r["inat_taxon_id"] in seen_taxa and seen_taxa[r["inat_taxon_id"]]!=r["species"]:
            raise ValueError("One original taxon id has multiple source names")
        seen_taxa[r["inat_taxon_id"]]=r["species"]
        if not re.fullmatch(
            r"https://inaturalist-open-data[.]s3[.]amazonaws[.]com/photos/"
            +re.escape(photo)+r"/large[.](?:jpg|jpeg|png)",r["photo_url_large"],re.I):
            raise ValueError("Original photo URL no longer matches frozen source photo ID")
        groups[(r["inat_taxon_id"],r["site_latcell"],r["site_loncell"],
                r["season_year"],r["selected_season"])]+=1
    if not _synthetic_test and (
        len(groups)!=10 or set(groups.values())!={40} or len(seen_taxa)!=5):
        raise ValueError("The five-taxon × two-season frozen source sample changed")
    return rows


def collect(
    source:Path, receipt:Path, out:Path,*, batch_size:int=8,
    pause:float=.65, _synthetic_test:bool=False, _request=None,
)->dict:
    if out.exists() and any(out.iterdir()):
        raise ValueError("Cannot overwrite previously observed public API annotation evidence")
    if not 1<=batch_size<=8 or (pause<.25 and not _synthetic_test):
        raise ValueError("Live provider API throttling/size bounds must be respected")
    rows=original_400(source,receipt,_synthetic_test=_synthetic_test)
    uuids=sorted({r["observation_id"] for r in rows})
    fetch=_request
    with requests.Session() as session:
        if fetch is None:
            fetch=lambda url,params=None:request(session,url,params)
        v2={}
        for offset in range(0,len(uuids),batch_size):
            wanted=uuids[offset:offset+batch_size]
            result=fetch(API+"/v2/observations/"+",".join(wanted),
                         {"fields":"uuid,id,observed_on"})
            items=result.get("results")
            if not isinstance(items,list):
                raise ValueError("Unexpected provider UUID-keyed observation result schema")
            for obj in items:
                uid=obj.get("uuid")
                if uid not in wanted or uid in v2:
                    raise ValueError("Provider UUID API returned alien or repeated original observation")
                num=obj.get("id")
                if not str(num).isdigit():
                    raise ValueError("Provider current original observation numeric ID invalid")
                v2[uid]=int(num)
            if pause:time.sleep(pause)
        numbers=sorted(set(v2.values()))
        v1={}
        for offset in range(0,len(numbers),batch_size):
            wanted=numbers[offset:offset+batch_size]
            result=fetch(API+"/v1/observations/"+",".join(map(str,wanted)))
            items=result.get("results")
            if not isinstance(items,list):
                raise ValueError("Unexpected current V1 public plant annotation result")
            for obj in items:
                num=obj.get("id")
                if num not in wanted or num in v1:
                    raise ValueError("Public annotation API returned alien/repeated original observation")
                v1[num]=obj
            if pause:time.sleep(pause)
    results=[]
    categories=Counter()
    group=defaultdict(Counter)
    matched=0
    for r in sorted(rows,key=lambda r:int(r["photo_id"])):
        key=(r["inat_taxon_id"],r["species"],r["site_latcell"],r["site_loncell"],
             r["season_year"],r["selected_season"])
        entry=group[key]
        entry["n"]+=1
        uid=r["observation_id"]
        current=v1.get(v2.get(uid,-1))
        if current is None:
            cat="no_current_public_observation_unknown"
            oid=""
            day_match=False
            original_photo_match=False
            values=[]
        else:
            if current.get("uuid")!=uid:
                raise ValueError("Current public iNat observation UUID changed from original source")
            oid=current["id"]
            day_match=current.get("observed_on")==r["observed_on"]
            original_photo_match=any(str(p.get("id"))==r["photo_id"]
                for p in (current.get("photos") or []) if isinstance(p,dict))
            cat,values=plant_status(current.get("annotations"))
            # A changed live source day/photo breaks the historical context,
            # but never removes the original source opportunity from totals.
            if not day_match:
                cat="live_observation_day_changed_since_source_annotation_NOT_USABLE"
            elif not original_photo_match:
                cat="original_photo_not_currently_on_live_observation_UNKNOWN"
        categories[(r["selected_season"],cat)]+=1
        entry["live"]+=int(current is not None)
        entry["both_match"]+=int(day_match and original_photo_match)
        if cat in CATEGORIES:
            entry[CATEGORIES[cat]]+=1
        elif cat in ("no_current_public_observation_unknown",
                   "no_current_observation_phenology_field_unknown",
                   "no_current_plant_phenology_annotation_unknown",
                   "original_photo_not_currently_on_live_observation_UNKNOWN",
                   "live_observation_day_changed_since_source_annotation_NOT_USABLE"):
            entry["unknown"]+=1
        else:
            entry["other"]+=1
        matched+=int(day_match and original_photo_match)
        results.append(dict(zip(ROWCOLS,(
            r["photo_id"],uid,r["observed_on"],r["inat_taxon_id"],r["species"],
            r["site_latcell"],r["site_loncell"],r["season_year"],
            r["selected_season"],oid,current is not None,day_match,
            original_photo_match,";".join(map(str,values)),cat,True,False
        ))))
    if len(results)!=len(rows) or sum(v["n"] for v in group.values())!=len(rows):
        raise ValueError("Current source annotation audit dropped original source photo denominator")
    out.mkdir(parents=True,exist_ok=True)
    all_csv=out/"current_live_iNaturalist_400_photo_original_observation_phenology_NOT_PHOTO.csv"
    with all_csv.open("x",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=ROWCOLS,lineterminator="\n")
        writer.writeheader();writer.writerows(results)
    summary=out/"current_live_iNaturalist_400_taxon_siteyear_season_phenology_availability.csv"
    with summary.open("x",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=GROUPCOLS,lineterminator="\n")
        writer.writeheader()
        for k,c in sorted(group.items()):
            writer.writerow(dict(zip(GROUPCOLS,(
                *k,c["n"],c["live"],c["both_match"],c["flowering"],c["fruiting"],
                c["flower_budding"],c["no_evidence"],c["unknown"],c["other"],0
            ))))
    report={
        "schema":"mois_original_400_observation_only_iNat_live_plant_phenology_v1",
        "status":"CURRENT_PUBLIC_OBSERVATION_PHENOLOGY_EFFORT_ONLY_NOT_BOTANICAL_OPEN_FLOWERS",
        "source_20260927_original_400_SHA256":sha(source),
        "source_selection_receipt_sha256":sha(receipt),
        "current_API_accessed_UTC":datetime.now(timezone.utc).isoformat(),
        "n_original_outcome_blind_selected_photo_opportunities":len(rows),
        "n_original_independent_observations":len(rows),
        "n_original_source_independent_photographers":len(rows),
        "n_source_taxon_siteyear_seasons":len(group),
        "n_live_current_v2_uuid_ids_resolved":len(v2),
        "n_live_current_v1_full_annotation_records":len(v1),
        "n_live_original_observed_on_and_photo_id_both_match":matched,
        "phenology_label_code_12_13_14_15_21":VALUE,
        "current_observation_status_by_source_local_season":{
            s:{name:n for (season,name),n in categories.items() if season==s}
            for s in sorted({r["selected_season"] for r in rows})},
        "source_original_photo_ids_never_replaced_or_filtered":True,
        "api_source_is_LIVE_NOT_20260927_frozen_annotations":True,
        "observation_level_annotation_does_not_confirm_open_petal_in_original_photo":True,
        "missing_live_annotation_is_unknown_not_nonflowering":True,
        "n_independently_botanically_verified_open_photo_petals":0,
        "no_colour_ROI_photo_measurements_from_these_400":True,
        "no_new_white_flower_frequency_or_weather_effect_estimated":True,
        "source_400_original_photo_current_observation_annotations_sha256":sha(all_csv),
        "taxon_season_phenology_availability_sha256":sha(summary),
    }
    (out/"mois_400_original_phenology_live_observation_ancillary_receipt.json").write_text(
        json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original-400-photo-manifest",type=Path,required=True)
    p.add_argument("--original-400-source-receipt",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    a=p.parse_args()
    r=collect(a.original_400_photo_manifest,a.original_400_source_receipt,a.outdir)
    print("MOIS_400_SOURCE_LIVE_PUBLIC_I_NAT_OBSERVATION_PHENOLOGY_QC",json.dumps({
        "original_source_photos":r["n_original_outcome_blind_selected_photo_opportunities"],
        "source_group_count":r["n_source_taxon_siteyear_seasons"],
        "live_resolved":r["n_live_current_v2_uuid_ids_resolved"],
        "original_photo_date_both_match":r["n_live_original_observed_on_and_photo_id_both_match"],
        "by_season":r["current_observation_status_by_source_local_season"],
        "no_validated_photo_open_flowers":True,
    },sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
