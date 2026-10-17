#!/usr/bin/env python3
"""Restore date *candidates* for original historical MOIS photos from public iNat.

Immutable 267425 original ROI photo IDs/hues are NEVER reassigned or replaced.
A blind, fixed 200-photo SOURCE sample tests the current public observation API
before committing to a 100543-photo historical date join. API observation-level
plant phenology tags are NOT independently verified photographed open petals.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import time
from datetime import date, datetime, timezone
from pathlib import Path
import pandas as pd
import requests
from scripts.analysis.mois_seasonal_267k_climate_normals import (
    MOIS_ROWS, MOIS_HUES, PHOTO_DATA, PHOTO_RECEIPT, HUES, read_signed, require,
)

PILOT_N = 200
BATCH = 50
DELAY_SECONDS = 1.25
SOURCE_N_HISTORICAL = 100543
SOURCE_PHOTO_SHA256 = "fbb54b58dc91e40e2e166081212588acc18adb19af4bdbcd1c7fe1778bee28d9"
API_URL = "https://api.inaturalist.org/v1/observations"


def deterministic_blind_sample(d:pd.DataFrame, n:int=PILOT_N) -> pd.DataFrame:
    """Hash original photo ID only, without colour, date, site or observed taxa."""
    require(n > 0 and n <= PILOT_N, "Only a bounded blind historic date-recovery pilot is supported")
    require(d.photo_id.is_unique and d.observation_id.is_unique,
            "Cannot pilot ambiguous/duplicate source photo or observation identity")
    require(not d.photo_id.isna().any() and not d.observation_id.isna().any(),
            "Original photo and observation IDs must not be missing")
    keys=d.photo_id.map(lambda x:hashlib.sha256(
        ("mois_historical_dated_source_pilot_v1:"+str(x)).encode()).hexdigest())
    chosen=d.assign(blind_photo_hash=keys).sort_values(
        ["blind_photo_hash","photo_id"],kind="stable").head(n)
    require(len(chosen)==n and chosen.observation_id.is_unique, "Blind source sample failed")
    return chosen.drop(columns="blind_photo_hash").reset_index(drop=True)


def normalize_public_observation(source:dict, observation:dict|None) -> dict:
    pid=str(source["photo_id"]);oid=str(source["observation_id"])
    result={"original_photo_id":pid,"original_observation_id":oid,
            "source_taxon_id":str(source["inat_taxon_id"]),
            "current_public_taxon_id":None,"source_photo_reappeared_in_api":False,
            "taxon_id_now_matches_source":None,"retrieved_observed_on":None,
            "has_current_observation_level_flowering_annotation":False,
            "current_observation_phenology_labels":None,
            "provenance_status":"not_found_in_current_public_api"}
    if observation is None:
        return result
    if str(observation.get("id"))!=oid:
        raise ValueError("API returned the wrong observation for original identity")
    nowtax=(observation.get("taxon") or {}).get("id")
    result["current_public_taxon_id"]=None if nowtax is None else str(nowtax)
    result["taxon_id_now_matches_source"]=(str(nowtax)==str(source["inat_taxon_id"]))
    photo_set={str(photo.get("id")) for photo in (observation.get("photos") or [])}
    if pid not in photo_set:
        result["provenance_status"]="observation_found_but_original_photo_not_in_current_photo_set"
        return result
    result["source_photo_reappeared_in_api"]=True
    current=observation.get("observed_on")
    try:
        dt=date.fromisoformat(str(current)) if current is not None else None
    except (ValueError,TypeError):
        dt=None
    if dt is None or dt>datetime.now(timezone.utc).date() or dt.year < 1800:
        result["provenance_status"]="original_photo_verified_but_current_observation_date_unavailable_or_invalid"
        return result
    result["retrieved_observed_on"]=dt.isoformat()
    result["provenance_status"]="source_photo_and_observation_verified_with_current_public_date"
    labels=sorted(set(str((ann.get("controlled_value") or {}).get("label",""))
                       for ann in (observation.get("annotations") or [])
                       if (ann.get("controlled_value") or {}).get("label")))
    result["current_observation_phenology_labels"]=";".join(labels) if labels else None
    result["has_current_observation_level_flowering_annotation"]=(
        "Flowers" in labels or "Flowers and Fruits" in labels)
    return result


def public_observations_batch(ids:list[str], *, session:requests.Session) -> tuple[dict,dict]:
    require(0<len(ids)<=BATCH and len(set(ids))==len(ids),
            "Duplicate or over-limit public observation API batch")
    headers={"User-Agent":"MOIS research date-provenance pilot (https://github.com/zuizui0223/mois)"}
    for retry in range(3):
        try:
            r=session.get(API_URL,params={"id":",".join(ids),"per_page":len(ids)},
                          headers=headers,timeout=30)
            if r.status_code in (429,500,502,503,504):
                if retry==2:
                    r.raise_for_status()
                time.sleep(min(12,2**retry*3))
                continue
            r.raise_for_status()
            data=r.json()
            break
        except (requests.RequestException,ValueError):
            if retry==2:
                raise
            time.sleep(2**retry)
    require(isinstance(data,dict) and isinstance(data.get("results"),list),
            "iNaturalist API payload schema changed")
    require(len(data["results"])<=len(ids), "API returned too many observations")
    require(data.get("total_results",0)<=len(ids),
            "API results paginated or query filter ignored; cannot silently lose IDs")
    lookup={}
    for observation in data["results"]:
        key=str(observation["id"])
        require(key in ids and key not in lookup,
                "Public API observation ID mismatch or duplicate")
        lookup[key]=observation
    canonical=json.dumps(data,ensure_ascii=False,sort_keys=True,
                         separators=(",",":")).encode()
    return lookup,{"request_count":1,"returned_observations":len(lookup),
                   "original_public_API_response_SHA256":hashlib.sha256(canonical).hexdigest()}


def run(original_zip:Path,output:Path,*,sample_n:int=PILOT_N) -> dict:
    require(not output.exists() or not any(output.iterdir()),
            "Never rewrite a completed source original sample or its receipts")
    d,receipt,source_sha=read_signed(original_zip,PHOTO_DATA,PHOTO_RECEIPT,
        "source_photo_level_csv_sha256",
        ["photo_id","observation_id","inat_taxon_id","old_mois_or_new_fcp",
         "observed_on","morph"])
    require(source_sha==SOURCE_PHOTO_SHA256,"Source original 267425 photo SHA drift")
    require(len(d)==MOIS_ROWS and d.morph.isin(HUES).sum()==MOIS_HUES,
            "Original signed MOIS measured pixel and hue census drift")
    historic=d[d.old_mois_or_new_fcp=="historical"].copy()
    require(len(historic)==SOURCE_N_HISTORICAL and
            historic.observed_on.isna().all() and
            historic.photo_id.is_unique and historic.observation_id.is_unique,
            "Original historical 100543-photo missing-date cohort changed")
    chosen=deterministic_blind_sample(historic,sample_n)
    receipt_ids=hashlib.sha256(("\n".join(chosen.photo_id)+"\n").encode()).hexdigest()
    lookup={};batch_ledger=[]
    with requests.Session() as client:
        ids=chosen.observation_id.tolist()
        for start in range(0,len(ids),BATCH):
            if start:
                time.sleep(DELAY_SECONDS)  # fewer than 60 requests/minute, total <=4 requests
            chunk=ids[start:start+BATCH]
            returned,info=public_observations_batch(chunk,session=client)
            overlap=set(lookup)&set(returned)
            require(not overlap,"Duplicate original observation results in separate batches")
            lookup.update(returned)
            batch_ledger.append({"batch_number":len(batch_ledger)+1,
                                 "requested_original_observations":len(chunk),**info})
    result_rows=[normalize_public_observation(r,lookup.get(str(r["observation_id"])))
                 for r in chosen.to_dict(orient="records")]
    output.mkdir(parents=True,exist_ok=True)
    out=output/"mois_historical_original_photo_observation_date_blind200_current_iNat.csv"
    pd.DataFrame(result_rows).to_csv(out,index=False,lineterminator="\n")
    success=sum(x["provenance_status"]=="source_photo_and_observation_verified_with_current_public_date"
                for x in result_rows)
    same_tax=sum(x["taxon_id_now_matches_source"] is True for x in result_rows)
    current_flowers=sum(x["has_current_observation_level_flowering_annotation"] for x in result_rows)
    summary={
        "schema":"mois_historical_original_photo_date_public_iNat_pilot_20261017_v1",
        "scientific_owner":"zuizui0223/mois",
        "original_full_source_photo_SHA256":source_sha,
        "historical_photo_count_without_source_calendar_date":SOURCE_N_HISTORICAL,
        "blind_sample_photo_count":len(chosen),"blind_source_photo_ids_SHA256":receipt_ids,
        "current_public_API_url":API_URL,
        "current_public_API_fetch_date_UTC":datetime.now(timezone.utc).date().isoformat(),
        "api_batches":batch_ledger,"n_actual_public_API_calls":len(batch_ledger),
        "n_returned_observation_records":len(lookup),
        "n_original_photo_ID_matched_and_current_observation_date_present":success,
        "n_original_current_taxon_ID_matches_source":same_tax,
        "n_current_observation_level_flower_tags_without_photo_specific_verification":current_flowers,
        "date_recovery_photo_csv_sha256":hashlib.sha256(out.read_bytes()).hexdigest(),
        "no_original_source_photo_classification_modified":True,
        "no_source_taxonomic_assignment_modified":True,
        "no_missing_calendar_date_imputation":True,
        "no_retrieved_observation_plant_phenology_assumed_to_be_photo_specific_open_flowers":True,
        "historical_photo_date_always_current_observation_metadata_NOT_original_archived_observed_on":True,
        "global_seasonal_white_effects_not_yet_refitted":True,
        "rate_limit":"4 API calls max, 50 IDs each, >=1.25 seconds between calls",
        "source_select_blind_to_colour_and_environment":True
    }
    (output/"mois_historical_blind200_observation_date_QC_receipt.json").write_text(
        json.dumps(summary,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("MOIS_ORIGINAL_HISTORICAL_PUBLIC_DATE_PILOT",json.dumps({
        "sample":len(chosen),"api_observations":len(lookup),
        "exact_photo_and_observation_date":success,
        "original_source_taxon_still_matches":same_tax,
        "current_plant_flowering_tags_not_independently_reviewed":current_flowers},
        sort_keys=True),flush=True)
    return summary


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-archive",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--sample-n",type=int,default=PILOT_N)
    a=p.parse_args()
    run(a.source_archive,a.output,sample_n=a.sample_n)
