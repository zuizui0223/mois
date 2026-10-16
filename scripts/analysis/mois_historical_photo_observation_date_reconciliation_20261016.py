#!/usr/bin/env python3
"""Historical 100,543 photo original-source DATE recovery; no new image downloads.

Uses the exact original FCP photo-site ZIP, SHA256-pinned. Source sampling:
cryptographic hash of ONLY original photo IDs, 7,200 of 100,348 original
PUBLIC VALID site photo opportunities. Colour is NOT loaded until the
sampled independent observation IDs and API date responses are frozen.

iNaturalist officially supports GET /v1/observations/{id1,id2,...} up to 200.
This pilot sends <=60 serial 120-ID requests, >=1.4 sec between requests,
and 2 limited retries for transient access. Uses current (October 2026)
observed_on, not a proven snapshot of the archived September 2026 source.
Joins dates ONLY if the originally selected photo ID remains attached to
the returned observation. Logs no individual photo ID, username or location.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import time
import urllib.error
import urllib.request
import zipfile
from collections import Counter, defaultdict
from datetime import date, datetime, timezone
from pathlib import Path

SCHEMA="mois_real_historical_100543_photo_original_obs_date_api_hash_sample_v1"
STATUS="CURRENT_INAT_OBS_DATE_JOINED_HISTORICAL_FCP_CLASSIFIED_PHOTOS_NOT_NEW_IMAGE_PIXEL_MEASUREMENTS"
SOURCE_SHA="2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0"
ZIP_MEMBER="fcp42111_real_site_abiotic/all_original_photo_site_abiotic_coverage.csv.gz"
BREADTH="fcp42111_real_site_abiotic/breadth_flower_colour_climate_soil.csv.gz"
REGIONAL="fcp42111_real_site_abiotic/taxon_cell_flower_colour_climate_soil.csv.gz"
SEED="MOIS_HISTORICAL_PHOTO_OBS_DATE_OUTCOME_BLIND_20261016_V1"
N_SOURCE=100543
N_VALID_PUBLIC=100348
N_SAMPLE=7200
BATCH=120
REQUEST_BUDGET=60
GAP=1.4
MAX_BODY=9*1024*1024
API="https://api.inaturalist.org/v1/observations/"
PHOTO_FOUR={"white","yellow_orange","red_pink","blue_purple"}
SOURCE_SNAPSHOT=date(2026,9,27)
MONTH_NAMES=("JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC")
PRIVATE_FIELDS=[
    "inat_taxon_id","observation_id","photo_id","latitude","longitude",
    "observed_on_current_api","observed_month","observed_year",
    "api_validation_status","api_photos_matched_to_original","api_current_taxon_same",
    "api_observer_id","current_flowers_annotation","measurement_status","morph",
    "photo_inclusion_probability","photo_design_weight","geographic_domain",
]


def digest(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(4*1024*1024),b""):
            h.update(block)
    return h.hexdigest()


def _stratum(lat:float)->str:
    return ("TROPICAL" if abs(lat)<23.5 else
            "NORTH_EXTRATROPICAL" if lat>=23.5 else
            "SOUTH_EXTRATROPICAL")


def _local_season(lat:float,month:int)->str:
    north={
        12:"WINTER",1:"WINTER",2:"WINTER",
        3:"SPRING",4:"SPRING",5:"SPRING",
        6:"SUMMER",7:"SUMMER",8:"SUMMER",
        9:"AUTUMN",10:"AUTUMN",11:"AUTUMN"}
    if abs(lat)<23.5:return f"CALENDAR_Q{(month-1)//3+1}"
    if lat>=23.5:return north[month]
    return {"WINTER":"SUMMER","SPRING":"AUTUMN",
            "SUMMER":"WINTER","AUTUMN":"SPRING"}[north[month]]


def _historical_rows(archive:Path, *, synthetic=False, max_photos=N_SAMPLE)->tuple[list[dict],int]:
    if not synthetic and digest(archive)!=SOURCE_SHA:
        raise ValueError("Original FCP 100543 source photo-site ZIP SHA256 mismatch")
    import pandas as pd
    with zipfile.ZipFile(archive) as z:
        if z.testzip() is not None or not {ZIP_MEMBER,BREADTH,REGIONAL}.issubset(z.namelist()):
            raise ValueError("Original photo-site/image-colour ZIP corrupted")
        cols=["inat_taxon_id","observation_id","photo_id","latitude","longitude","site_geo_status"]
        source=pd.read_csv(z.open(ZIP_MEMBER),compression="gzip",usecols=cols)
        if not synthetic:
            if (len(source)!=N_SOURCE or source.photo_id.nunique()!=N_SOURCE
                or source.observation_id.nunique()!=N_SOURCE):
                raise ValueError("Original 100543 independent source photo/observation census changed")
        good=source.loc[source.site_geo_status.eq("VALID_PUBLIC_ORIGINAL_PHOTO_POINT")].copy()
        if not synthetic and len(good)!=N_VALID_PUBLIC:
            raise ValueError("Original unobscured public photo location denominator changed")
        if len(good)<max_photos or len(good)==0 or good.photo_id.duplicated().any():
            raise ValueError("Original public site source insufficient for chosen independent photo ID count")
        good["hash_sort"]=[
            hashlib.sha256((SEED+"/"+str(int(x))).encode("ascii")).hexdigest()
            for x in good.photo_id]
        sampled=good.sort_values(["hash_sort","photo_id"],kind="stable").head(max_photos)
        assert len(sampled)==max_photos and sampled.observation_id.is_unique
        # IMPORTANT: photo hue remains untouched until AFTER source
        # admission and outcome-blind observation-date ID selection.
        chosen=set(sampled.photo_id.astype("int64").tolist())
        labels={}
        for member in (BREADTH,REGIONAL):
            t=pd.read_csv(z.open(member),compression="gzip",
               usecols=["photo_id","observation_id","inat_taxon_id",
                        "morph","measurement_status"])
            for row in t.itertuples(index=False):
                photo=int(row.photo_id)
                if photo not in chosen:continue
                entry=(int(row.observation_id),int(row.inat_taxon_id),
                       str(row.measurement_status),str(row.morph))
                if photo in labels and labels[photo]!=entry:
                    raise ValueError("Same photo has incompatible historical ROI/four-colour labels")
                labels[photo]=entry
        if chosen!=set(labels):
            raise ValueError("Original FCP photo coverage rows not all mapped to the actual pixel colour archive")
        result=[]
        for row in sampled.itertuples(index=False):
            pid=int(row.photo_id)
            obs=int(row.observation_id)
            tax=int(row.inat_taxon_id)
            x=labels[pid]
            if x[:2]!=(obs,tax) or not (math.isfinite(row.latitude) and math.isfinite(row.longitude)):
                raise ValueError("Source original photo/taxon/site lineage discrepancy")
            if (x[2]=="classified_four_state_morph")!=(x[3] in PHOTO_FOUR):
                raise ValueError("FCP original real image morph and ROI quality status drifted")
            result.append({
                "inat_taxon_id":tax,"observation_id":obs,"photo_id":pid,
                "latitude":float(row.latitude),"longitude":float(row.longitude),
                "measurement_status":x[2],"morph":x[3]
            })
    return result,len(good)


def _valid_date(value)->str|None:
    if not isinstance(value,str) or not re.fullmatch(r"\d{4}-\d\d-\d\d",value):
        return None
    try:d=date.fromisoformat(value)
    except ValueError:return None
    return value if d<=SOURCE_SNAPSHOT else None


def _annotation(result:dict)->str:
    a=result.get("annotations")
    if not isinstance(a,list):
        return "API_ANNOTATIONS_NOT_INCLUDED"
    stages=[]
    for x in a:
        if not isinstance(x,dict) or x.get("controlled_attribute_id")!=12:continue
        val=x.get("controlled_value_id")
        score=x.get("vote_score")
        if val in (13,14,15,21) and (
            score is None or isinstance(score,int) and score>0):
            stages.append(val)
    return ("FLOWERS" if 13 in stages else "BUDS" if 15 in stages
            else "FRUITS" if 14 in stages else
            "NO_FLOWERS_FRUITS" if 21 in stages else "NO_ANNOTATION")


def _response(ids:list[int],opener=urllib.request.urlopen)->tuple[dict, str]:
    if not 1<=len(ids)<=BATCH or any(not isinstance(v,int) or v<=0 for v in ids):
        raise ValueError("Only <=120 selected original observation IDs per provider request")
    url=API+",".join(str(v) for v in ids)
    req=urllib.request.Request(url,headers={
        "User-Agent":"Mois-2026-archived-photo-calendar-pilot/1.0 (+https://github.com/zuizui0223/mois)",
        "Accept":"application/json"})
    with opener(req,timeout=75) as resp:
        if getattr(resp,"status",200)!=200:
            raise ValueError("iNaturalist API server did not serve HTTP 200")
        b=resp.read(MAX_BODY+1)
    if len(b)>MAX_BODY:raise ValueError("Original observation response exceeded protected request byte limit")
    obj=json.loads(b)
    if (not isinstance(obj,dict) or not isinstance(obj.get("results"),list)
        or not isinstance(obj.get("total_results"),int)
        or obj["total_results"]<0 or obj["total_results"]>len(ids)):
        raise ValueError("Unexpected iNaturalist /observations/{ids} source API schema")
    matched={}
    selected=set(ids)
    for row in obj["results"]:
        if not isinstance(row,dict) or type(row.get("id")) is not int or row["id"] not in selected:
            raise ValueError("Extraneous returned public observation is not in original selected source IDs")
        if row["id"] in matched:raise ValueError("Duplicate original observation in one provider response")
        matched[row["id"]]=row
    return matched,hashlib.sha256(b).hexdigest()


def _summarize(rows:list[dict], *, min_public=25)->list[dict]:
    grouped=defaultdict(list)
    for r in rows:
        if r["api_validation_status"]!="MATCHED_ORIGINAL_PHOTO_WITH_VALID_DATE":continue
        group=(r["geographic_domain"],r["observed_month"])
        grouped[group].append(r)
    summaries=[]
    for (zone,month),subset in sorted(grouped.items()):
        n=len(subset)
        if n<min_public:
            summaries.append({"zone":zone,"month":month,
                              "public_suppressed_below":min_public})
            continue
        classify=[r for r in subset if r["measurement_status"]=="classified_four_state_morph"]
        whites=sum(r["morph"]=="white" for r in classify)
        by_tax=defaultdict(lambda:[0,0])
        for row in classify:
            tax=row["inat_taxon_id"]
            by_tax[tax][0]+=1
            by_tax[tax][1]+=int(row["morph"]=="white")
        observers={str(r["api_observer_id"]) for r in subset if r["api_observer_id"]}
        summaries.append({
            "zone":zone,"month":month,"original_photo_sample_obs_with_verified_date":n,
            "n_classifiable_original_photo_hues":len(classify),
            "n_white_among_classifiable":whites,
            "roi_classifiability_fraction":len(classify)/n,
            "photo_white_fraction_given_classifiable":(
                whites/len(classify) if classify else None),
            "species_equal_white_fraction_given_classifiable":(
                sum(w/count for count,w in by_tax.values())/len(by_tax)
                if by_tax else None),
            "n_distinct_original_taxa_with_classifiable_photo":len(by_tax),
            "n_distinct_identified_current_API_observers":len(observers),
            "positive_Flowers_annotation_source_obs":sum(
                r["current_flowers_annotation"]=="FLOWERS" for r in subset),
            "annotations_verified_as_of_archive_snapshot":False
        })
    return summaries


def run(archive:Path,output_private:Path,output_public:Path,*,synthetic=False,
        max_photos=N_SAMPLE,batch_size=BATCH,pause=GAP,
        getter=_response,sleep=time.sleep)->dict:
    if any(p.exists() for p in (output_private,output_public)):
        raise ValueError("Refuse overwrite of original source observation dates or outcomes")
    if not 1<=max_photos<=N_SAMPLE or not 1<=batch_size<=BATCH or \
       math.ceil(max_photos/batch_size)>REQUEST_BUDGET or pause<GAP:
        raise ValueError("Exceeded original-provider 60-query/safe-pacing colour-blind date pilot")
    rows,population=_historical_rows(archive,synthetic=synthetic,max_photos=max_photos)
    obs={r["observation_id"]:r for r in rows}
    receipts=[];lookups={}
    batches=[rows[j:j+batch_size] for j in range(0,len(rows),batch_size)]
    for i,batch in enumerate(batches):
        if i:sleep(pause)
        ids=[r["observation_id"] for r in batch]
        # Transient GET request failures are retried at limited rate;
        # persistent failures MUST NOT create a falsely complete source.
        for retry in range(3):
            try:
                found,api_sha=getter(ids)
                break
            except (urllib.error.URLError,TimeoutError) as e:
                if retry>=2:raise
                sleep(pause*(retry+1))
        if not set(found).issubset(set(ids)):
            raise ValueError("Unexpected IDs were returned from original source batch")
        if len(api_sha)!=64:
            raise ValueError("Raw HTTP response source SHA256 required for every observation batch")
        lookups.update(found)
        receipts.append({"batch_index":i,"selected_ids":len(ids),
                         "returned_observations":len(found),
                         "api_response_sha256":api_sha})
    status=Counter()
    observer_ids=set()
    by_year=Counter()
    for r in rows:
        entry=lookups.get(r["observation_id"])
        r.update({"observed_on_current_api":"",
                  "observed_month":"","observed_year":"",
                  "api_validation_status":"NOT_RETURNED_CURRENT_PUBLIC_API",
                  "api_photos_matched_to_original":False,
                  "api_current_taxon_same":None,
                  "api_observer_id":"",
                  "current_flowers_annotation":"NOT_VERIFIED",
                  "photo_inclusion_probability":max_photos/population,
                  "photo_design_weight":population/max_photos,
                  "geographic_domain":_stratum(r["latitude"])})
        if entry is None:
            status[r["api_validation_status"]]+=1
            continue
        current_photos=entry.get("photos")
        if not isinstance(current_photos,list) or not any(
            isinstance(p,dict) and p.get("id")==r["photo_id"]
            for p in current_photos):
            r["api_validation_status"]="ORIGINAL_SOURCE_PHOTO_NOT_ATTACHED_ON_CURRENT_API"
            status[r["api_validation_status"]]+=1
            continue
        r["api_photos_matched_to_original"]=True
        dt=_valid_date(entry.get("observed_on"))
        if dt is None:
            r["api_validation_status"]="OBSERVED_DATE_MISSING_INVALID_OR_AFTER_ORIGINAL_SNAPSHOT"
            status[r["api_validation_status"]]+=1
            continue
        r["observed_on_current_api"]=dt
        r["observed_month"]=int(dt[5:7])
        r["observed_year"]=int(dt[:4])
        r["api_validation_status"]="MATCHED_ORIGINAL_PHOTO_WITH_VALID_DATE"
        tax=entry.get("taxon")
        r["api_current_taxon_same"]=(
            isinstance(tax,dict) and tax.get("id")==r["inat_taxon_id"])
        user=entry.get("user")
        if isinstance(user,dict) and type(user.get("id")) is int:
            r["api_observer_id"]=user["id"]
            observer_ids.add(user["id"])
        r["current_flowers_annotation"]=_annotation(entry)
        by_year[r["observed_year"]]+=1
        status[r["api_validation_status"]]+=1
    if sum(status.values())!=max_photos:
        raise ValueError("Original color-blind photo source date coverage ledger incomplete")
    for p in (output_private,output_public):
        p.parent.mkdir(parents=True,exist_ok=True)
    tmp=output_private.with_name(output_private.name+".INCOMPLETE")
    with tmp.open("w",newline="",encoding="utf8") as f:
        writer=csv.DictWriter(f,fieldnames=PRIVATE_FIELDS)
        writer.writeheader()
        for r in rows:
            writer.writerow({k:r[k] for k in PRIVATE_FIELDS})
    report={
        "schema":SCHEMA,"status":STATUS,
        "source_original_100543_photo_zip_sha256":(
            SOURCE_SHA if not synthetic else digest(archive)),
        "source_eligible_original_public_site_photo_opportunities":population,
        "source_original_archived_photo_opportunities":(
            N_SOURCE if not synthetic else None),
        "outcome_blind_original_photo_hash_sampling_salt":SEED,
        "n_random_hash_sampled_original_observations":max_photos,
        "sample_original_photographic_inclusion_probability":max_photos/population,
        "n_serial_public_API_observation_batch_queries":len(batches),
        "max_original_ids_per_API_request":batch_size,
        "provider_minimum_inter_request_seconds":pause,
        "source_original_API_response_sha256_by_batch":receipts,
        "current_api_original_photo_and_date_status":dict(sorted(status.items())),
        "n_exact_original_photo_present_with_valid_current_API_date":status[
            "MATCHED_ORIGINAL_PHOTO_WITH_VALID_DATE"],
        "n_qualified_current_API_observers":len(observer_ids),
        "matched_source_observed_year_counts":dict(sorted(by_year.items())),
        "source_observed_month_by_hemisphere_photo_white":_summarize(rows),
        "original_selected_photo_private_date_join_sha256":digest(tmp),
        "original_API_observed_on_not_proven_frozen_at_archive_date":True,
        "original_four_hue_source_pixels_already_measured_not_newly_reprocessed":True,
        "source_not_a_probability_sample_of_all_world_flowering_plants":True,
        "source_historically_selected_100543_photos_not_all_93m_photo_opportunities":True,
        "original_FCP_42111_all_photo_ID_census_now_complete":False,
        "new_13_environment_climate_photo_regression_estimated":False,
        "original_iNaturalist_image_pixels_downloaded":False,
        "photographic_hue_and_date_not_used_for_sample_admission":True,
        "current_API_timestamp_UTC":datetime.now(timezone.utc).isoformat(),
        "synthetic_test_only":bool(synthetic),
        "note":"Old 100543 measured photo opportunity cohort, 7200 unique source original observation IDs outcome-blind sampled. Current iNaturalist API observed_on/photo-attachment verification, subject to platform changes and missingness. Month white fractions conditional on image ROI; not world flowering flora or causal climate slopes."
    }
    # Do not commit or upload private table with photos/sites/observer IDs.
    stage=output_public.with_name(output_public.name+".INCOMPLETE")
    stage.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n",encoding="utf8")
    os.replace(tmp,output_private)
    os.replace(stage,output_public)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original-historical-photo-source-zip",type=Path,required=True)
    p.add_argument("--out-private-current-date-photo-join",type=Path,required=True)
    p.add_argument("--out-public-date-month-receipt",type=Path,required=True)
    a=p.parse_args()
    outcome=run(a.original_historical_photo_source_zip,
                a.out_private_current_date_photo_join,
                a.out_public_date_month_receipt)
    print("MOIS_ACTUAL_HISTORICAL_ORIGINAL_PHOTOS_WITH_CURRENT_API_OBS_MONTHS",
          json.dumps({k:outcome[k] for k in (
              "source_eligible_original_public_site_photo_opportunities",
              "n_random_hash_sampled_original_observations",
              "n_exact_original_photo_present_with_valid_current_API_date",
              "current_api_original_photo_and_date_status",
              "source_observed_month_by_hemisphere_photo_white"
          )},sort_keys=True))


if __name__=="__main__":
    main()
