#!/usr/bin/env python3
"""MOIS first16 priority species: true original CC photographs from exact
iNaturalist Flowers (12|13) annotations, excluding all 48 earlier observations.

Original byte/SHA/UUID provenance, botanical photo contact sheets and review
templates are the deliverable. A flower annotation is *not* a verified petal
ROI, a pigment measurement, a species-wide fixed hue or an LGM refugium.
"""
from __future__ import annotations

import argparse
import io
import json
import re
import time
from collections import Counter
from pathlib import Path
from urllib.error import HTTPError,URLError
from urllib.parse import urlencode
from urllib.request import Request,urlopen
from zipfile import ZipFile

import pandas as pd

from scripts.acquisition.mois_20261119_targeted_16_native_photo_contact_proof import (
    API, TARGET_OBS, TARGET_SPECIES_MAX, read_source_queue,
    source_species_map, fetch_original, thumbnail, contact_sheet,
    sha_file, sha_bytes, history, require, strict_new_photo
)

SCHEMA="mois_native_photo9_first16_explicit_flowering_annotated_original_20261120_v1"
PRIOR_ZIP_SHA="a6c96f5e2db118bef282688a237d720362d651a16a611ca01907ad12c5cb0998"
PRIOR_CSV="actual_CC_original_photo_UUID_SHA_metadata.csv"
FLOWER_ATTRIBUTE=12
FLOWER_VALUE=13
MAX_PAGES=3
PAGE_SIZE=24


def read_first48(path):
    require(sha_file(path)==PRIOR_ZIP_SHA,
            "Previous true original 48-photo ZIP digest drifted")
    with ZipFile(path) as z:
        receipt=json.loads(z.read("receipt.json"))
        require(receipt["n_exact_species_targeted"]==16
                and receipt["n_authentic_CC_original_image_SHA_verified"]==48
                and receipt["n_species_with_3_independent_photos"]==16
                and receipt["source_api_gate_errors"]==0,
                "First16 prior proof counts changed")
        for name,digest in receipt["outputs_sha256"].items():
            require(sha_bytes(z.read(name))==digest,
                    "Previous original 48 image metadata/source receipt changed "+name)
        for name,digest in receipt["contact_sha256"].items():
            require(sha_bytes(z.read(name))==digest,
                    "Previously reviewed contact sheet SHA changed")
        d=pd.read_csv(io.BytesIO(z.read(PRIOR_CSV)),
                      dtype={"photo_id":str,"observation_id":str,"inat_taxon_id":str})
    require(len(d)==48 and d.photo_id.nunique()==48
            and d.observation_id.nunique()==48
            and d.candidate_species.nunique()==16,
            "First48 original image IDs or source taxonomy not unique")
    require(d.original_image_sha256.map(
                lambda s: bool(re.fullmatch("[a-f0-9]{64}",str(s)))).all(),
            "First48 original image SHA missing")
    d["flowers_annotated"]=[flower_annotation(json.loads(s))
                    for s in d.plant_phenology_annotations_json]
    require(int(d.flowers_annotated.sum())==5,
            "Previous first48 floral annotation baseline unexpectedly changed")
    return d


def flower_annotation(items):
    """Never accept an API query parameter alone as observation evidence."""
    return any(isinstance(a,dict)
               and str(a.get("controlled_attribute_id"))==str(FLOWER_ATTRIBUTE)
               and str(a.get("controlled_value_id"))==str(FLOWER_VALUE)
               for a in items)


def request_flower_page(taxon_id,page):
    params={
        "taxon_id":str(taxon_id), "term_id":FLOWER_ATTRIBUTE,
        "term_value_id":FLOWER_VALUE, "photos":"true",
        "photo_licensed":"true", "quality_grade":"research",
        "per_page":PAGE_SIZE, "page":page, "order_by":"observed_on",
        "order":"desc", "captive":"false"
    }
    url=API+"?"+urlencode(params)
    req=Request(url,headers={
        "Accept":"application/json",
        "User-Agent":"MOIS-research-real-flowering-photo-proof/1.0"
    })
    last=None
    for attempt in range(3):
        try:
            with urlopen(req,timeout=35) as response:
                require(response.status==200,"iNaturalist API non-200")
                raw=response.read(12_000_001)
            require(len(raw)<=12_000_000,
                    "iNaturalist metadata response exceeds 12MB")
            data=json.loads(raw)
            require(isinstance(data,dict) and isinstance(data.get("results"),list),
                    "Unexpected iNaturalist observations API response")
            return url,data
        except (HTTPError,URLError,TimeoutError,ValueError,json.JSONDecodeError) as exc:
            last=str(exc)
            if attempt<2:time.sleep(attempt+1)
    raise ValueError("iNaturalist annotated page failed: "+str(last)[:180])


def eligible_from_page(data,tid,used_ids,used_obs,seen,counts):
    output=[]
    for obs in data["results"]:
        counts["api_observations_seen"]+=1
        taxon=obs.get("taxon") or {}
        if not isinstance(taxon,dict) or str(taxon.get("id"))!=str(tid):
            counts["wrong_taxon_or_descendant"]+=1
            continue
        if obs.get("captive") is True:
            counts["captive"]+=1
            continue
        annotations=obs.get("annotations") or []
        if not flower_annotation(annotations):
            counts["API_flower_annotation_not_present_in_record"]+=1
            continue
        counts["exact_taxon_flowering_annotated"]+=1
        uuid=str(obs.get("uuid") or "")
        if not re.fullmatch(r"[0-9a-fA-F-]{36}",uuid):
            counts["bad_observation_UUID"]+=1
            continue
        if uuid in used_obs or uuid in seen:
            counts["historical_or_previous_wave_UUID"]+=1
            continue
        photo=None
        for media in obs.get("photos") or []:
            checked=strict_new_photo(media,tid)
            if checked is not None and checked["photo_id"] not in used_ids:
                photo=checked;break
        if photo is None:
            counts["no_unused_licensed_original_photo"]+=1
            continue
        if photo["photo_id"] in used_ids:
            counts["photo_id_already_used"]+=1
            continue
        photo.update({
            "observation_id":uuid,
            "observation_numeric_id":str(obs.get("id") or ""),
            "inat_taxon_id":str(tid),
            "observed_on":str(obs.get("observed_on") or ""),
            "observer_login":str((obs.get("user") or {}).get("login") or ""),
            "plant_phenology_annotations_json":json.dumps(annotations),
            "iNaturalist_flowers_annotation_12_13_in_response":True,
            "iNaturalist_research_grade_reported":str(obs.get("quality_grade") or ""),
            "source_api_flowers_annotation_query":True,
            "true_botanical_petal_ROI_verified":False,
            "new_image_9hue_ROI_measured":False,
            "species_colour_confirmed":False,
            "original_raw_bytes_retained":False
        })
        output.append(photo)
        seen.add(uuid)
        counts["unused_licensed_flowers_observations"]+=1
    return output


def run(args):
    root=Path(args.out)
    require(not root.exists(),"Refuse overwrite of real photo source proof")
    queue=read_source_queue(args.targets)
    target=source_species_map(args.source_env,queue)
    old_ids,old_obs,old_raw=history()
    prior=read_first48(args.prior48)
    used_ids=old_ids|set(prior.photo_id.astype(str))
    used_obs=old_obs|set(prior.observation_id.astype(str))
    used_sha=old_raw|set(prior.original_image_sha256.astype(str))
    require(len(old_ids)==1231337 and len(prior)==48
            and not(set(prior.photo_id)&old_ids)
            and not(set(prior.observation_id)&old_obs),
            "Previous first48 cannot be disjoint from historical original ID ledger")

    root.mkdir(parents=True)
    statuses=[];photorecords=[];errors=[];review=[]
    seen_photos=set(); seen_obs=set(); seen_sha=set()
    for item in target.itertuples(index=False):
        tid=str(item.inat_taxon_id)
        state={"rank":int(item.rank),"candidate_species":item.candidate_species,
               "inat_taxon_id":tid,"original_proxy_hue9":item.source_photo9_label,
               "old_48_photos_excluded":48,
               "n_flowering_annotated_original_photos":0,
               "n_independent_observation_UUID":0,
               "n_API_pages_queried":0,
               "record_state":"NO_NEW_ANNOTATED_FLOWER_PHOTO"}
        counts=Counter()
        picked=[];visited=set()
        for page in range(1,MAX_PAGES+1):
            try:
                api_url,payload=request_flower_page(tid,page)
                state["n_API_pages_queried"]+=1
                state["API_flowering_total_reported"]=int(
                    payload.get("total_results") or 0)
                candidates=eligible_from_page(payload,tid,used_ids,
                                               used_obs,visited,counts)
                for photo in candidates:
                    if len(picked)>=TARGET_OBS:break
                    if photo["photo_id"] in seen_photos or photo["observation_id"] in seen_obs:
                        continue
                    try:
                        raw,picture,original_url=fetch_original(photo)
                        digest=sha_bytes(raw)
                        if digest in used_sha or digest in seen_sha:
                            counts["byte_SHA_duplicate"]+=1
                            continue
                        photo.update({
                           "original_image_sha256":digest,
                           "original_image_bytes":len(raw),
                           "original_image_dimensions":[picture.width,picture.height],
                           "original_url":original_url,
                           "candidate_species":item.candidate_species,
                           "rank":int(item.rank),
                           "thumbnail":thumbnail(picture)
                        })
                        picked.append(photo)
                        seen_photos.add(photo["photo_id"])
                        seen_obs.add(photo["observation_id"])
                        seen_sha.add(digest)
                    except Exception as exc:
                        errors.append({"species":item.candidate_species,
                           "photo_id":photo["photo_id"],
                           "error":str(exc)[:180]})
                if len(picked)>=TARGET_OBS:break
                if len(payload["results"])<PAGE_SIZE:break
            except Exception as exc:
                errors.append({"species":item.candidate_species,"page":page,
                                "error":str(exc)[:220]})
                counts["API_page_error"]+=1
                break
            finally:
                time.sleep(args.pace)
        picked=picked[:TARGET_OBS]
        state.update({k:int(v) for k,v in counts.items()})
        state["n_flowering_annotated_original_photos"]=len(picked)
        state["n_independent_observation_UUID"]=len(set(p["observation_id"] for p in picked))
        if picked:
            folder=root/f"species_{int(item.rank):03d}_{tid}"
            folder.mkdir()
            sheet=contact_sheet(item.candidate_species,picked,folder)
            state["contact_sheet"]=str(sheet.relative_to(root))
            state["contact_SHA256"]=sha_file(sheet)
            state["record_state"]="FLOWERING_ANNOTATED_SOURCE_READY_PETAL_REVIEW_REQUIRED"
        elif counts.get("API_page_error"):
            state["record_state"]="API_FAILURE_OR_RETRY_REQUIRED"
        for p in picked:
            rec={k:v for k,v in p.items() if k!="thumbnail"}
            photorecords.append(rec)
            review.append({
                "rank":int(item.rank),"species":item.candidate_species,
                "photo_id":p["photo_id"],"observation_id":p["observation_id"],
                "iNat_flowering_annotation":True,
                "flower_organ_visible_REVIEW_REQUIRED":"",
                "corolla_or_perianth_visible_REVIEW_REQUIRED":"",
                "floral_ROI_measurement_eligible_REVIEW_REQUIRED":"",
                "reviewed_hue9_REVIEW_REQUIRED":"",
                "reviewer":"","review_notes":"",
                "botanical_decision_not_yet_made":True
            })
        statuses.append(state)
        print("MOIS_FLOWER_ANNOTATION_SOURCE_ROW",
              json.dumps({"species":item.candidate_species,
                          "n_images":len(picked),
                          "state":state["record_state"],
                          "api_pages":state["n_API_pages_queried"]}),flush=True)

    photos=pd.DataFrame(photorecords)
    if photos.empty:
        photos=pd.DataFrame(columns=[
            "candidate_species","inat_taxon_id","photo_id","observation_id",
            "original_image_sha256","iNaturalist_flowers_annotation_12_13_in_response"
        ])
    status=pd.DataFrame(statuses)
    pfile=root/"flowering_annotation_original_photo_SHA_UUID_sources.csv"
    sfile=root/"source_species_flowering_query_status.csv"
    rfile=root/"BOTANICAL_PETAL_REVIEW_NOT_YET_FILLED.csv"
    efile=root/"original_image_or_API_failures.json"
    pfile.write_text(photos.to_csv(index=False))
    sfile.write_text(status.to_csv(index=False))
    rfile.write_text(pd.DataFrame(review).to_csv(index=False))
    efile.write_text(json.dumps(errors,indent=2)+"\n")

    html=["<!doctype html><html lang='en'><meta charset='utf-8'>",
          "<title>MOIS Flowers Annotation Source QA</title>",
          "<body><h1>Source-annotated FLOWERS — NOT validated perianths</h1>"]
    for state in statuses:
        if state.get("contact_sheet"):
            html.append("<p>"+state["candidate_species"]+" — "
              +str(state["n_flowering_annotated_original_photos"])
              +" independent observations</p><img width='980' src='"
              +state["contact_sheet"]+"'>")
    html.append("</body></html>")
    index=root/"contact_sheet_index.html"
    index.write_text("\n".join(html)+"\n")
    require(photos.photo_id.is_unique and photos.observation_id.is_unique,
            "New flowering wave repeats photo or observation")
    if len(photos):
        require(photos.original_image_sha256.is_unique,
                "New flowering wave contains byte-identical original photos")
        require(photos.iNaturalist_flowers_annotation_12_13_in_response.eq(True).all(),
                "Source flowering annotation missing from some accepted record")
    receipt={
        "schema":SCHEMA,
        "source_original_PRIOR48_ZIP_SHA256":PRIOR_ZIP_SHA,
        "source_queue_199_genera_github_blob_sha1":
           "efc85eb9cc9fa842d630d69ab1962fa7a5029ab1",
        "source_37104_accepted_species_ENV_ZIP_SHA256":
           "d8b97528188d9f8d63cb802d5b3f6e8f7f251937fd30af42aeae47e43ebad1fc",
        "first16_prior_true_original_photos":len(prior),
        "first16_prior_flowering_annotations":int(prior.flowers_annotated.sum()),
        "base_source_photo_ids_excluded":len(old_ids),
        "base_source_observation_UUIDs_excluded":len(old_obs),
        "new_annotated_original_photo_SHA_checked":len(photos),
        "new_species_with_any_annotated_image":int(
             status.n_flowering_annotated_original_photos.gt(0).sum()),
        "new_species_with_three_annotations":int(
             status.n_flowering_annotated_original_photos.ge(3).sum()),
        "n_target_species":TARGET_SPECIES_MAX,
        "n_annotation_API_page_errors":int(status.API_page_error.fillna(0).sum())
            if "API_page_error" in status else 0,
        "n_original_image_fetch_errors":int(sum("photo_id" in e for e in errors)),
        "former_48_and_history_disjoint_by_photo_and_observation_UUID":True,
        "later_MOIS_2026_waves_not_fully_excluded":True,
        "source_annotated_flowering_is_not_visible_petal":True,
        "petal_organ_manual_review_has_not_been_completed":True,
        "new_palette9_colours_measured":0,
        "strict_species_colour_confirmed":0,
        "outputs_SHA256":{
            name.name:sha_file(name) for name in
            (pfile,sfile,rfile,efile,index)
        },
        "contact_SHA256":{
            x["contact_sheet"]:x["contact_SHA256"]
            for x in statuses if x.get("contact_sheet")
        }
    }
    (root/"receipt.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print("MOIS_NATIVE_FLOWERS_ANNOTATED_FIRST16_AUTHENTIC_SOURCE",
          json.dumps(receipt,sort_keys=True),flush=True)
    require(receipt["new_annotated_original_photo_SHA_checked"]>0,
            "No new verified flowers-annotated original pixels; do not report success")

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--targets",required=True)
    parser.add_argument("--source-env",required=True)
    parser.add_argument("--prior48",required=True)
    parser.add_argument("--out",required=True)
    parser.add_argument("--pace",type=float,default=0.8)
    args=parser.parse_args()
    require(0.6<=args.pace<=3.0,"Unbounded API pace")
    run(args)
