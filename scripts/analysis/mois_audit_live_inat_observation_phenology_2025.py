#!/usr/bin/env python3
"""MOIS supplemental iNaturalist PUBLIC observation phenology annotation census.

Reads only the frozen 2025 photo source. Fetches current public records read-only
to quantify availability of existing crowd annotations. Not historically frozen
in 2026-09-27 source and NOT photographic petal/tepal anthesis validation.

No photo/image pixel retrieval, no source-photo replacement, no colour
labels supplied to the annotator lookup or used for API selection.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,time,uuid,zipfile
from collections import Counter
from pathlib import Path

import requests

API="https://api.inaturalist.org"
AGENT={"User-Agent":"mois-original-dated-photo-phenology-availability-audit/1.0"}
SOURCE_CSV="postblind/first_original_source_matched_season_64_photo_roi_QC.csv"
SOURCE_RECEIPT="postblind/first_original_source_matched_season_64_roi_receipt.json"
EXPECTED_SHA="1b16cdabd2936567dee87b8344019f40f315600c6510dfc1f25a03835024b6a0"
FIELDS=("photo_id","original_observation_uuid","original_observation_day",
        "original_photo_season","live_2026_10_16_api_status",
        "live_numeric_observation_id","live_observed_on_matches_original",
        "original_photo_still_listed_with_live_observation",
        "plant_phenology_attribute_present_in_live_API",
        "live_observation_plant_phenology_value_ids",
        "plant_phenology_observation_annotation_status_NOT_PHOTO",
        "photo_open_flower_petal_tepal_independently_verified")
TERM=12
VALUE={13:"flowering",14:"fruiting",15:"flower_budding"}


def digest(raw:bytes)->str:
    return hashlib.sha256(raw).hexdigest()


def load_source(archive:Path)->list[dict]:
    with zipfile.ZipFile(archive) as z:
        raw=z.read(SOURCE_CSV)
        receipt=json.loads(z.read(SOURCE_RECEIPT))
    if (receipt.get("postblind_join_sha256")!=digest(raw)
        or digest(raw)!=EXPECTED_SHA or
        receipt.get("status")!="FIRST_ORIGINAL_SEASONAL_PIXELS_MEASURED_BOTANICAL_ANTHESIS_NOT_VALIDATED"
        or receipt.get("original_photos_measured")!=64):
        raise ValueError("Not the exact frozen FCP original 2025 summer/autumn 64 photo source")
    rows=list(csv.DictReader(io.StringIO(raw.decode("utf-8"))))
    if len(rows)!=64 or len({q["photo_id"] for q in rows})!=64 or len({q["observation_id"] for q in rows})!=64:
        raise ValueError("Original photo IDs or observation UUIDs not unique")
    for r in rows:
        if str(uuid.UUID(r["observation_id"]))!=r["observation_id"] or r["observed_on"][:7]!=r["observed_year_month"]:
            raise ValueError("Original signed observation identity/date malformed")
        if r["local_season"] not in ("summer","autumn"):
            raise ValueError("Non-original paired summer/autumn source month")
    return rows


def plant_status(annotations)->tuple[str,list[int]]:
    if not isinstance(annotations,list):
        return "no_current_observation_phenology_field_unknown",[]
    values=set()
    for a in annotations:
        if not isinstance(a,dict):
            raise ValueError("Unknown live API annotation record type")
        term=a.get("controlled_attribute_id",a.get("term_id"))
        val=a.get("controlled_value_id",a.get("term_value_id"))
        if str(term)==str(TERM):
            try:
                values.add(int(val))
            except (ValueError,TypeError) as e:
                raise ValueError("Unknown original provider phenology value code") from e
    if not values:
        return "no_current_plant_phenology_annotation_unknown",[]
    if len(values)>1:
        return "conflicting_current_observation_phenology_NOT_photo",sorted(values)
    v=next(iter(values))
    return (("current_observation_annotation_"+VALUE[v]+"_NOT_photo"
             if v in VALUE else "unrecognized_current_observation_phenology_NOT_photo"),[v])


def request(session:requests.Session,url:str,params:dict|None=None)->dict:
    error=None
    for attempt in range(3):
        try:
            response=session.get(url,params=params,timeout=(8,25),headers=AGENT)
            response.raise_for_status()
            if len(response.content)>2_500_000:
                raise ValueError("Unbounded live observation API response")
            return response.json()
        except (requests.RequestException,ValueError) as e:
            error=e
            if attempt<2:time.sleep(2*(attempt+1))
    raise RuntimeError("Public iNaturalist observation API temporarily unavailable; do not label missing as unannotated") from error


def run(source_zip:Path,out:Path,*,batch=8,delay=.6)->dict:
    if out.exists() and any(out.iterdir()):
        raise ValueError("Do not overwrite live independently sourced observations")
    if not 1<=batch<=8 or not .25<=delay<=3:
        raise ValueError("Keep source provider batch sizes and rate conservative")
    photos=load_source(source_zip)
    uuids={r["observation_id"]:r for r in photos}
    all_ids=sorted(uuids)
    live_id_by_uuid={}
    by_id={}
    with requests.Session() as s:
        for start in range(0,len(all_ids),batch):
            subset=all_ids[start:start+batch]
            path="/v2/observations/"+",".join(subset)
            data=request(s,API+path,{"fields":"uuid,id,observed_on"})
            results=data.get("results",[])
            if not isinstance(results,list):
                raise ValueError("Unexpected current iNaturalist v2 UUID response")
            for r in results:
                orig=r.get("uuid")
                if orig not in subset or orig in live_id_by_uuid:
                    raise ValueError("Current API observation UUID not from frozen photo source or duplicated")
                if not str(r.get("id","")).isdigit():
                    raise ValueError("Current API observation numeric ID is not verifiable")
                live_id_by_uuid[orig]=int(r["id"])
            time.sleep(delay)

        # Batch numeric-ID V1 is needed: V2 fields=all omits the annotations
        # key even when V1 contains actual crowd-supplied annotation records.
        numeric=sorted(set(live_id_by_uuid.values()))
        for start in range(0,len(numeric),batch):
            subset=numeric[start:start+batch]
            data=request(s,API+"/v1/observations/"+",".join(str(i) for i in subset))
            results=data.get("results",[])
            if not isinstance(results,list):
                raise ValueError("Unexpected live annotation API response")
            for r in results:
                oid=r.get("id")
                if oid not in subset or oid in by_id:
                    raise ValueError("Current phenology API returned wrong original observation identity")
                by_id[oid]=r
            time.sleep(delay)

    rows=[]
    by_season=Counter()
    statuses=Counter()
    for source in photos:
        ou=source["observation_id"]
        pid=source["photo_id"]
        current=by_id.get(live_id_by_uuid.get(ou,-1))
        if current is None:
            status="no_current_public_observation_unknown"
            anno_values=[]
            date_match=False
            photo_match=False
            obs_id=""
            present=False
        else:
            if current.get("uuid")!=ou:
                raise ValueError("Live V1 observation UUID differs from original 2026-09-27 source")
            obs_id=current.get("id")
            date_match=current.get("observed_on")==source["observed_on"]
            photo_match=any(str(p.get("id"))==pid for p in (current.get("photos") or [])
                            if isinstance(p,dict))
            present=isinstance(current.get("annotations"),list)
            status,anno_values=plant_status(current.get("annotations"))
            if not date_match:
                status="live_observation_date_changed_since_source_annotation_UNUSABLE"
            if not photo_match:
                status="original_source_photo_not_associated_with_live_observation_unknown"
        season=source["local_season"]
        by_season[season]+=1
        statuses[(season,status)]+=1
        rows.append(dict(zip(FIELDS,(
            pid,ou,source["observed_on"],season,
            "present_current_live_observation" if current else "unavailable_current_observation",
            obs_id,date_match,photo_match,present,
            ";".join(str(v) for v in anno_values),status,False))))
    if by_season!={"summer":32,"autumn":32}:
        raise ValueError("Current iNaturalist annotation API query changed sample membership")
    out.mkdir(parents=True,exist_ok=True)
    dest=out/"original_2025_64_live_iNaturalist_observation_plant_phenology_NOT_PHOTO.csv"
    with dest.open("x",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=FIELDS,lineterminator="\n")
        writer.writeheader()
        writer.writerows(sorted(rows,key=lambda r:int(r["photo_id"])))
    receipt={
        "schema":"mois_live_public_iNat_observation_annotation_supplement_to_frozen_2025_photos_v1",
        "status":"CURRENT_OBSERVATION_PHENOLOGY_METADATA_ONLY_NO_PHOTO_ANTHESIS_CONFIRMED",
        "frozen_original_photo_receipt_sha256":EXPECTED_SHA,
        "live_api_read_only":True,
        "live_observation_annotations_not_frozen_on_2026_09_27":True,
        "source_original_photos_not_changed_or_replaced":True,
        "n_original_source_photos":len(photos),
        "n_live_v2_uuid_resolved":len(live_id_by_uuid),
        "n_live_v1_observation_records":len(by_id),
        "by_season_source_photo_counts":dict(by_season),
        "by_season_current_annotation_status":{
            season:{status:n for (s,status),n in statuses.items() if s==season}
            for season in ("summer","autumn")},
        "flowering_term_id":12,
        "flowering_term_value_id":13,
        "fruiting_term_value_id":14,
        "flower_budding_term_value_id":15,
        "observation_annotation_does_not_certify_specific_photo_open_flower":True,
        "actual_photo_open_flower_petal_tepal_verified_count":0,
        "never_label_annotation_absence_nonflowering":True,
        "no_seasonal_white_flower_rate_or_climate_effect_inferred":True,
        "original_photo_observation_live_annotation_table_sha256":digest(dest.read_bytes()),
    }
    (out/"2025_Asclepias_live_public_observation_phenology_receipt.json").write_text(
        json.dumps(receipt,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    return receipt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fcp-original-2025-64-ROI-zip",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    args=p.parse_args()
    result=run(args.fcp_original_2025_64_ROI_zip,args.outdir)
    print("MOIS_CURRENT_LIVE_INAT_PLANT_PHENOLOGY_64_ORIGINAL_SOURCE_RESULTS",
          json.dumps({
            "n_live_v2_uuid_resolved":result["n_live_v2_uuid_resolved"],
            "n_live_v1_observation_records":result["n_live_v1_observation_records"],
            "statuses_by_season":result["by_season_current_annotation_status"],
          },sort_keys=True))


if __name__=="__main__":
    main()
