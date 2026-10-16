#!/usr/bin/env python3
"""Stream every ORIGINAL date-verified FCP photo into blinded 4096-photo ROI packets.

The full ORIGINAL source is the 2026-09-27 FCP 42,111 taxa/93,209,351 photo
opportunity census. This is DATE-FIRST COMPUTE PRIORITIZATION, not a 4096
photo sampling cap: every source photo with valid original observed_on from
an individually attested FCP date shard is dispatched exactly once.

The original iNaturalist four-table metadata has NO photo-specific blooming
annotation. Every queued photograph is BLOOM UNVERIFIED and needs original
ROI-v4 pixel inspection; no Flowers-observation tags, historical hues,
current-API dates or source-predicted seasonal white colours determine selection.

Each packet carries independent worker/acquisition/private-ID joins, with
photo, date, observer, species and geographic coordinates entirely hidden
from the image worker. The protected full queue is never published; no
provider image URL is opened by this code.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sqlite3
from pathlib import Path

from scripts.acquisition.mois_full_source_dated_flower_candidates_20261016 import (
    sha,SCHEMA as DATE_SCHEMA,STATUS as DATE_STATUS,
    source_observation_day,climatic_season,FULL_WORK_PARTITIONS,
    FRAME_SHA,
)

SCHEMA="mois_full_original_dated_photo_blind_roi_packet_source_v1"
STATUS="ALL_SOURCE_VALID_DATED_PHOTOS_DISPATCHED_BLOOM_UNVERIFIED_NO_PIXELS"
PACKET_SCHEMA="mois_full_original_dated_photo_blind_roi_packet_v1"
PACKET_STATUS="BLIND_ORIGINAL_DATE_ELIGIBLE_IMAGES_NO_PIXELS_OPENED"
MAX_PACKET=4096
DEFAULT_PACKET=4096
IMAGE_COMPUTE_BATCH=64
ID_SALT="MOIS_ALL_ORIGINAL_OBS_DATE_POSITIVE_PHOTO_ROI_V1"
URL_RE=re.compile(
    r"^https://inaturalist-open-data\.s3\.amazonaws\.com/photos/([0-9]+)/large\.(jpg|jpeg|png)$",re.I
)
WORKER=("measurement_id","image_filename","photo_license")
ACQUISITION=("measurement_id","image_filename","photo_url_large","photo_license")
PRIVATE=(
    "measurement_id","original_photo_id","original_observation_uuid",
    "inat_taxon_id","species","observer_id","observed_on","observed_month",
    "latitude","longitude","latitude_domain","local_season",
    "photo_license","source_snapshot_date","bloom_state_before_roi",
)


def blind_id(photo_id:str)->str:
    if not isinstance(photo_id,str) or not photo_id.isdecimal() or int(photo_id)<=0:
        raise ValueError("No original provider photo identity")
    digest=hashlib.sha256((ID_SALT+"\x1f"+photo_id).encode("ascii")).hexdigest()
    return "FCPR-"+digest.upper()[:24]


def _write_packet(parent:Path,index:int,source:list[dict],input_queue_sha:str)->dict:
    if not 1<=len(source)<=MAX_PACKET:
        raise ValueError("Zero or oversized source photo compute packets cannot be sealed")
    name=f"packet_{index:06d}"
    folder=parent/name
    folder.mkdir()
    worker=folder/"worker.csv"
    acquisition=folder/"acquisition.csv"
    private=folder/"private_source.csv"
    with worker.open("x",encoding="utf8",newline="") as w, \
         acquisition.open("x",encoding="utf8",newline="") as a, \
         private.open("x",encoding="utf8",newline="") as b:
        ww=csv.DictWriter(w,fieldnames=WORKER);ww.writeheader()
        aa=csv.DictWriter(a,fieldnames=ACQUISITION);aa.writeheader()
        bb=csv.DictWriter(b,fieldnames=PRIVATE);bb.writeheader()
        for row in source:
            photo=row["photo_id"];mid=blind_id(photo)
            original_name=mid+".jpg"
            basic={
                "measurement_id":mid,"image_filename":original_name,
                "photo_license":row["photo_license"],
            }
            ww.writerow(basic)
            aa.writerow({**basic,"photo_url_large":row["photo_url_large"]})
            bb.writerow({
                "measurement_id":mid,"original_photo_id":photo,
                "original_observation_uuid":row["observation_id"],
                "inat_taxon_id":row["inat_taxon_id"],
                "species":row["species"],
                "observer_id":row["observer_id"],
                "observed_on":row["observed_on"],
                "observed_month":row["observed_month"],
                "latitude":row["latitude"],"longitude":row["longitude"],
                "latitude_domain":row["latitude_domain"],
                "local_season":row["local_season"],
                "photo_license":row["photo_license"],
                "source_snapshot_date":row["source_snapshot_date"],
                "bloom_state_before_roi":row["original_photo_flowering_state"],
            })
    record={
        "schema":PACKET_SCHEMA,"status":PACKET_STATUS,
        "packet":index,"n_source_original_date_known_photos":len(source),
        "n_compute_batches_of_max_64_photos":(
            (len(source)+IMAGE_COMPUTE_BATCH-1)//IMAGE_COMPUTE_BATCH),
        "source_original_date_queue_sha256":input_queue_sha,
        "file_sha256":{
            "worker.csv":sha(worker),"acquisition.csv":sha(acquisition),
            "private_source.csv":sha(private)
        },
        "worker_only_fields":list(WORKER),
        "worker_has_observation_location_species_or_hue":False,
        "original_photo_selection_cap":None,
        "original_flower_roi_verified":False,
        "date_is_frozen_original_snapshot_not_live_api":True,
        "original_photo_pixels_opened":False,
    }
    (folder/"packet_receipt.json").write_text(
        json.dumps(record,indent=2,sort_keys=True)+"\n",encoding="utf8")
    return {"packet_index":index,"packet_receipt_sha256":sha(folder/"packet_receipt.json"),
            "n_source_dated_photos":len(source)}


def dispatch(date_queue:Path,date_receipt:Path,packet_dir:Path, *,
             packet_size:int=DEFAULT_PACKET,synthetic=False)->dict:
    if packet_size<1 or packet_size>MAX_PACKET:
        raise ValueError("4096 is a bounded compute packet maximum, not source sample cap")
    stage=packet_dir.with_name(packet_dir.name+".INCOMPLETE")
    if packet_dir.exists() or stage.exists():
        raise ValueError("No overwrite of already sealed or interrupted original photo ROI packets")
    source=json.loads(date_receipt.read_text(encoding="utf8"))
    if (source.get("schema")!=DATE_SCHEMA
        or source.get("status")!=DATE_STATUS
        or source.get("private_full_dated_source_photo_queue_sha256")!=sha(date_queue)
        or source.get("source_photo_cap_per_taxon","INVALID") is not None
        or source.get("photo_hue_classification_attempted") is not False
        or source.get("actual_flower_photographs_verified_from_pixels")!=0
        or source.get("metadata_date_known_is_not_blooming_photo_confirmation") is not True
        or source.get("synthetic_test_only") is not synthetic
        or not synthetic and (
            source.get("original_source_frame_sha256")!=FRAME_SHA
            or source.get("n_fcp_shards_expected_for_global_completion")
                !=FULL_WORK_PARTITIONS)):
        raise ValueError("Date-before-bloom source packet has no complete FCP source evidence")
    packet_dir.parent.mkdir(parents=True,exist_ok=True)
    stage.mkdir()
    sqlite=stage/"_source_identity_ledger.sqlite"
    count=0;seen_obs=0;packets=[];buffer=[]
    domains={}
    with sqlite3.connect(sqlite) as db:
        db.execute("PRAGMA journal_mode=DELETE")
        db.execute("""CREATE TABLE photo_ids(
           photo_id TEXT PRIMARY KEY, measurement_id TEXT UNIQUE NOT NULL)""")
        db.execute("""CREATE TABLE observations(
           uuid TEXT PRIMARY KEY,taxon INTEGER NOT NULL,day TEXT NOT NULL,
           lat REAL NOT NULL,lon REAL NOT NULL,observer TEXT NOT NULL)""")
        with date_queue.open(encoding="utf8") as f:
            for line in f:
                r=json.loads(line)
                if (r.get("original_photo_flowering_state")!="NOT_YET_PHOTO_VERIFIED"
                    or r.get("origin")!="ORIGINAL_DATED_FCP_SHARD_NO_PHOTO_PIXELS"
                    or not str(r.get("photo_id","")).isdecimal()
                    or not str(r.get("inat_taxon_id","")).isdecimal()
                    or not str(r.get("observation_id","")).strip()
                    or not str(r.get("observer_id","")).strip()):
                    raise ValueError("Original dated flower-unverified photo missing source identity")
                photo=str(r["photo_id"]);taxon=int(r["inat_taxon_id"])
                if not taxon or not int(photo):
                    raise ValueError("Nonpositive original date-known taxon or photo ID")
                day,state=source_observation_day(r.get("observed_on"))
                if day is None or state!="ORIGINAL_SOURCE_DATE_KNOWN" or \
                   r.get("observed_month")!=day.month or r.get("observed_year")!=day.year:
                    raise ValueError("Date-first ROI source must retain real original observed_on")
                lat=float(r["latitude"]);lon=float(r["longitude"])
                if not (-90<=lat<=90 and -180<=lon<=180):
                    raise ValueError("Original source photo has out-of-range geography")
                zone,season=climatic_season(lat,day.month)
                if zone!=r["latitude_domain"] or season!=r["local_season"]:
                    raise ValueError("Original hemisphere-local season was modified")
                url=str(r.get("photo_url_large",""))
                match=URL_RE.fullmatch(url)
                if not match or match.group(1)!=photo or not r.get("photo_license"):
                    raise ValueError("Original photo ID/URL/license does not match source archive")
                mid=blind_id(photo)
                try:
                    db.execute("INSERT INTO photo_ids VALUES (?,?)",(photo,mid))
                except sqlite3.IntegrityError as exc:
                    raise ValueError("Original provider photo ID or blinded ID duplicate") from exc
                observer=str(r["observer_id"])
                uuid=str(r["observation_id"])
                cur=db.execute("INSERT OR IGNORE INTO observations VALUES (?,?,?,?,?,?)",
                    (uuid,taxon,day.isoformat(),lat,lon,observer))
                if cur.rowcount:
                    seen_obs+=1
                else:
                    values=db.execute(
                        "SELECT taxon,day,lat,lon,observer FROM observations WHERE uuid=?",
                        (uuid,)).fetchone()
                    if values!=(taxon,day.isoformat(),lat,lon,observer):
                        raise ValueError("Original observation multi-photo site/date/taxon inconsistent")
                buffer.append(r)
                count+=1
                domains[zone]=domains.get(zone,0)+1
                if len(buffer)==packet_size:
                    packets.append(_write_packet(stage,len(packets),buffer,sha(date_queue)))
                    buffer=[]
                if count%10000==0:db.commit()
        if buffer:packets.append(_write_packet(stage,len(packets),buffer,sha(date_queue)))
        db.commit()
    sqlite.unlink()
    if count!=source.get("n_original_photo_observation_dates_usable") or \
       seen_obs>count or sum(x["n_source_dated_photos"] for x in packets)!=count:
        raise ValueError("Source photo date-eligible coverage was not conserved in blind dispatch")
    result={
        "schema":SCHEMA,"status":STATUS,
        "source_full_dated_shard_receipt_sha256":sha(date_receipt),
        "source_private_dated_photo_queue_sha256":sha(date_queue),
        "source_original_shard":source["shard_number"],
        "n_original_photo_opportunities_in_sealed_shard":source[
            "n_original_eligible_source_photographs_examined"],
        "n_date_known_original_photo_ids_dispatched_to_roi":count,
        "n_original_independent_observation_uuids_dispatched":seen_obs,
        "n_bounded_original_roi_packets":len(packets),
        "max_photo_count_per_packet":packet_size,
        "image_compute_batch_size_not_sample_cap":IMAGE_COMPUTE_BATCH,
        "n_original_photo_ids_omitted_despite_valid_source_date":0,
        "n_original_photos_with_undetermined_flowering_state":count,
        "n_original_photo_pixels_measured_at_dispatch":0,
        "photo_level_phenology_annotations_not_in_original_metadata":True,
        "photos_from_source_without_date_not_sent_to_image_evaluator":True,
        "all_photo_ids_with_source_date_included_without_sampling_cap":True,
        "per_original_packet":packets,
        "date_known_geographic_zone_candidate_counts":domains,
        "protected_source_photo_ids_and_coordinates_in_public_receipt":False,
        "synthetic_test_only":bool(synthetic)
    }
    (stage/"full_dated_roi_dispatch_receipt.json").write_text(
        json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf8")
    stage.replace(packet_dir)
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--original-private-dated-photo-queue",required=True,type=Path)
    ap.add_argument("--original-source-date-receipt",required=True,type=Path)
    ap.add_argument("--out-private-roi-packets",required=True,type=Path)
    ap.add_argument("--compute-packet-size",type=int,default=DEFAULT_PACKET)
    args=ap.parse_args()
    result=dispatch(args.original_private_dated_photo_queue,
                    args.original_source_date_receipt,
                    args.out_private_roi_packets,
                    packet_size=args.compute_packet_size)
    print(json.dumps({
       "status":result["status"],
       "n_original_date_known_photos_dispatched":result[
           "n_date_known_original_photo_ids_dispatched_to_roi"],
       "n_image_pixels_measured":0,
       "n_packets":result["n_bounded_original_roi_packets"]
    },sort_keys=True))


if __name__=="__main__":
    main()
