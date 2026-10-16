#!/usr/bin/env python3
"""Reassemble EVERY original date-valid photo ROI outcome in ONE FCP source shard.

Every source photographed opportunity with originally observed_on has one
sealed image ROI terminal, or this shard is explicitly incomplete. An
actual photo contributes to seasonal white/coloured only when original
ROI-v4 says classified_four_state_morph; no individual source photo is
pre-labelled blooming by its iNaturalist observation annotation.

The date-first process can save image work on date-absent photo IDs while
retaining their original source exclusion denominator. A single source
shard result is PARTIAL metadata and blooming-photo coverage of the full
42111/93209351 worldwide frame, never worldwide ecological inference.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sqlite3
from collections import Counter,defaultdict
from pathlib import Path

from scripts.acquisition.mois_full_source_dated_flower_candidates_20261016 import (
    sha,source_observation_day,climatic_season
)
from scripts.acquisition.mois_prepare_all_dated_original_roi_packets_20261016 import (
    SCHEMA as DISPATCH_SCHEMA,STATUS as DISPATCH_STATUS,
    PACKET_SCHEMA,PACKET_STATUS,blind_id,IMAGE_COMPUTE_BATCH,
    WORKER,PRIVATE
)
from scripts.acquisition.mois_execute_all_dated_original_roi_workunit_20261016 import (
    SCHEMA as WORK_SCHEMA,STATUS as WORK_STATUS,BIO,FCP_BACKEND_GIT_BLOB
)
from scripts.analysis.mois_original_full_fcp_dated_flower_ROI_20261016 import (
    MORPHS,MODEL_COMMIT,MEASUREMENT_BLOB,DETECTOR_SHA,
    SAM_ENCODER_SHA,SAM_DECODER_SHA
)

SCHEMA="mois_completed_one_original_date_known_photo_shard_ROI_bloom_v1"
STATUS="SHARD_COMPLETE_ORIGINAL_DATED_PHOTO_ROI_MEASURED_NOT_FULL_GLOBAL"
POSITIVE_FIELDS=(
    "original_photo_id","original_observation_uuid","original_taxon_id",
    "original_species_name","original_genus","original_observer",
    "original_observed_on","observed_month","latitude","longitude",
    "latitude_domain","local_season","original_ROI_v4_hue",
    "actual_photo_contains_classifiable_flower",
    "original_snapshot_20260927_date_proven"
)


def _validate_real_method(method:dict):
    if (method.get("source_commit")!=MODEL_COMMIT
        or method.get("measurement_git_blob")!=MEASUREMENT_BLOB
        or method.get("detector_sha256")!=DETECTOR_SHA
        or method.get("efficient_sam_encoder_sha256")!=SAM_ENCODER_SHA
        or method.get("efficient_sam_decoder_sha256")!=SAM_DECODER_SHA
        or method.get("location_blind") is not True
        or method.get("no_old_measurement_csv_inputs") is not True
        or method.get("persisted_image_pixels") is not False):
        raise ValueError("ROI method was not the original frozen FCP YOLO/EfficientSAM")


def reassemble(dispatch_root:Path,terminals_root:Path,
               out_private_positive:Path,out_public_receipt:Path,*,
               synthetic=False,min_public=10)->dict:
    if min_public<10:
        raise ValueError("Original source small photo/observer bins cannot be published")
    for path in (out_private_positive,out_public_receipt):
        if path.exists() or path.with_name(path.name+".INCOMPLETE").exists():
            raise ValueError("Cannot overwrite or replace original source ROI photo evidence")
    source_path=dispatch_root/"full_dated_roi_dispatch_receipt.json"
    dispatch=json.loads(source_path.read_text(encoding="utf8"))
    if (dispatch.get("schema")!=DISPATCH_SCHEMA
        or dispatch.get("status")!=DISPATCH_STATUS
        or dispatch.get("synthetic_test_only") is not synthetic
        or dispatch.get("n_original_photo_ids_omitted_despite_valid_source_date")!=0
        or dispatch.get("n_original_photo_pixels_measured_at_dispatch")!=0
        or dispatch.get("all_photo_ids_with_source_date_included_without_sampling_cap") is not True
        or dispatch.get("n_original_photo_opportunities_in_sealed_shard",0) <
            dispatch.get("n_date_known_original_photo_ids_dispatched_to_roi",0)):
        raise ValueError("Not the complete original source date-known photo dispatch")
    out_private_positive.parent.mkdir(parents=True,exist_ok=True)
    out_public_receipt.parent.mkdir(parents=True,exist_ok=True)
    stage=out_private_positive.with_name(out_private_positive.name+".INCOMPLETE")
    public_stage=out_public_receipt.with_name(out_public_receipt.name+".INCOMPLETE")
    private_month=out_private_positive.parent/"original_roi_month_effort_private.json"
    private_month_stage=private_month.with_name(private_month.name+".INCOMPLETE")
    if private_month.exists() or private_month_stage.exists():
        raise ValueError("Existing private original month ROI effort evidence cannot be overwritten")
    sqlite=out_private_positive.with_name(out_private_positive.name+".ids.sqlite.INCOMPLETE")
    original_n=0
    classified_n=0
    white_n=0
    terminal_n=0
    failure_n=0
    observation_ids=set()
    statuses=Counter()
    month_zone=defaultdict(Counter)
    packet_hashes=[]
    try:
        with sqlite3.connect(sqlite) as db,stage.open("w",encoding="utf8") as output:
            db.execute("PRAGMA journal_mode=DELETE")
            db.execute("PRAGMA synchronous=NORMAL")
            db.execute("CREATE TABLE ids (measurement TEXT PRIMARY KEY, photo_id TEXT UNIQUE NOT NULL)")
            db.execute("""CREATE TABLE observations(
                uuid TEXT PRIMARY KEY, taxon INTEGER, observer TEXT, day TEXT,
                lat REAL, lon REAL, n_classified INTEGER NOT NULL DEFAULT 0,
                n_white INTEGER NOT NULL DEFAULT 0)""")
            for packet in dispatch["per_original_packet"]:
                p=packet["packet_index"]
                folder=dispatch_root/f"packet_{p:06d}"
                receipt_path=folder/"packet_receipt.json"
                if sha(receipt_path)!=packet["packet_receipt_sha256"]:
                    raise ValueError("Original source photo packet receipt SHA changed")
                sealed=json.loads(receipt_path.read_text(encoding="utf8"))
                n=sealed["n_source_original_date_known_photos"]
                if (sealed.get("schema")!=PACKET_SCHEMA
                    or sealed.get("status")!=PACKET_STATUS
                    or n!=packet["n_source_dated_photos"]
                    or n<1
                    or sealed.get("source_original_date_queue_sha256")!=
                        dispatch["source_private_dated_photo_queue_sha256"]
                    or sealed.get("original_flower_roi_verified") is not False):
                    raise ValueError("Original dated source packet is incomplete/unrecognized")
                for name in ("worker.csv","acquisition.csv","private_source.csv"):
                    if sha(folder/name)!=sealed["file_sha256"][name]:
                        raise ValueError("Original date-first photo worker/private source was changed")
                with (folder/"worker.csv").open(newline="",encoding="utf8") as w,\
                     (folder/"private_source.csv").open(newline="",encoding="utf8") as f:
                    worker=csv.DictReader(w)
                    private=csv.DictReader(f)
                    if tuple(worker.fieldnames or ())!=WORKER or \
                       tuple(private.fieldnames or ())!=PRIVATE:
                        raise ValueError("Source worker no longer location-blind or private join changed")
                    assigned=list(worker)
                    identities=list(private)
                if len(assigned)!=n or len(identities)!=n:
                    raise ValueError("Original source packet membership count changed")
                by_measurement={}
                for w,original in zip(assigned,identities):
                    mid=w["measurement_id"]
                    if (mid!=original["measurement_id"]
                        or mid!=blind_id(original["original_photo_id"])
                        or original["bloom_state_before_roi"]!="NOT_YET_PHOTO_VERIFIED"):
                        raise ValueError("Original photo ID source/worker blinded mapping mismatch")
                    try:
                        db.execute("INSERT INTO ids VALUES (?,?)",(mid,original["original_photo_id"]))
                    except sqlite3.IntegrityError as exc:
                        raise ValueError("Original photo identity repeated across compute packets") from exc
                    by_measurement[mid]=original
                original_n+=n
                sub_receipts=[]
                for batch in range((n+IMAGE_COMPUTE_BATCH-1)//IMAGE_COMPUTE_BATCH):
                    unit=terminals_root/f"photo_roi_packet_{p:06d}_batch_{batch:03d}"
                    if not unit.is_dir():
                        raise ValueError("Original date-valid photo has unmeasured ROI image work unit")
                    work=unit/"original_dated_roi_work_unit_receipt.json"
                    terminal=unit/"original_fcp_date_valid_photo_ROI_v4_terminal.csv"
                    rec=json.loads(work.read_text(encoding="utf8"))
                    if (rec.get("schema")!=WORK_SCHEMA or rec.get("status")!=WORK_STATUS
                        or rec.get("source_packet_number")!=p
                        or rec.get("source_image_compute_batch")!=batch
                        or rec.get("source_packet_receipt_sha256")!=sha(receipt_path)
                        or rec.get("source_dispatch_receipt_sha256")!=sha(source_path)
                        or rec.get("original_terminal_csv_sha256")!=sha(terminal)
                        or rec.get("synthetic_test_only") is not synthetic
                        or rec.get("no_historical_image_measurements_reused") is not True
                        or rec.get("original_fcp_image_pixels_saved_after_evaluation") is not False
                        or rec.get("missing_photo_flower_hues_imputed") is not False):
                        raise ValueError("A source original photo ROI-v4 terminal was not sealed")
                    if not synthetic:
                        if rec.get("original_fcp_image_backend_git_blob")!=FCP_BACKEND_GIT_BLOB:
                            raise ValueError("Original source photo evaluator backend changed")
                        _validate_real_method(rec.get("original_source_frozen_roi_method",{}))
                    ids=[r["measurement_id"] for r in
                         assigned[batch*IMAGE_COMPUTE_BATCH:(batch+1)*IMAGE_COMPUTE_BATCH]]
                    unit_status=Counter()
                    with terminal.open(newline="",encoding="utf8") as f:
                        reader=csv.DictReader(f)
                        required={"measurement_id","morph","measurement_status",
                                  *(f"flower_fraction_{c}" for c in BIO)}
                        if not required.issubset(reader.fieldnames or ()) or \
                           any(k in (reader.fieldnames or []) for k in (
                               "original_photo_id","latitude","longitude",
                               "inat_taxon_id","observed_on","species","photo_url_large")):
                            raise ValueError("Original image terminal leaked source geographical/date identity")
                        seen=set()
                        for row in reader:
                            mid=row["measurement_id"]
                            if mid not in ids or mid in seen:
                                raise ValueError("Original photo ROI terminal missing/extra/repeated identifier")
                            seen.add(mid)
                            original=by_measurement[mid]
                            cl=row["measurement_status"]=="classified_four_state_morph"
                            hue=row["morph"]
                            if cl!=(hue in MORPHS) or (
                                not cl and hue!="mixed_uncertain"):
                                raise ValueError("Source original ROI-v4 flower hue gate changed")
                            unit_status["n"]+=1
                            unit_status["class"]+=int(cl)
                            unit_status["white"]+=int(cl and hue=="white")
                            unit_status[row["measurement_status"]]+=1
                            if cl:
                                fraction=[float(row[f"flower_fraction_{k}"]) for k in BIO]
                                if any(not math.isfinite(v) or v<0 for v in fraction) or \
                                   abs(sum(fraction)-1)>1e-7:
                                    raise ValueError("Actual original photo nine-colour ROI simplex invalid")
                            # Date+photo ROI is a DOUBLE positive source gate.
                            date_string=original["observed_on"]
                            dt,state=source_observation_day(date_string)
                            if dt is None or state!="ORIGINAL_SOURCE_DATE_KNOWN" or \
                               dt.month!=int(original["observed_month"]):
                                raise ValueError("Original archived date could not be restored to photo")
                            lat=float(original["latitude"]);lon=float(original["longitude"])
                            domain,season=climatic_season(lat,dt.month)
                            if domain!=original["latitude_domain"] or \
                               season!=original["local_season"]:
                                raise ValueError("Source photo original hemisphere or month changed")
                            source_zone=month_zone[(domain,dt.month)]
                            source_zone["all"]+=1
                            source_zone["class"]+=int(cl)
                            source_zone["white"]+=int(cl and hue=="white")
                            source_zone["permanent_photo_acquisition_failed"]+=int(
                                row["measurement_status"]=="image_acquisition_failed")
                            if cl:
                                # Count one source observation once across
                                # multiple original source photographs.
                                obs=original["original_observation_uuid"]
                                taxon=int(original["inat_taxon_id"])
                                observer=original["observer_id"]
                                cur=db.execute("""INSERT INTO observations
                                   (uuid,taxon,observer,day,lat,lon,n_classified,n_white)
                                   VALUES (?,?,?,?,?,?,1,?)
                                   ON CONFLICT(uuid) DO UPDATE SET
                                     n_classified=n_classified+1,
                                     n_white=n_white+excluded.n_white
                                   WHERE taxon=excluded.taxon
                                     AND observer=excluded.observer
                                     AND day=excluded.day
                                     AND lat=excluded.lat
                                     AND lon=excluded.lon""",
                                   (obs,taxon,observer,date_string,lat,lon,
                                    int(hue=="white")))
                                if cur.rowcount!=1:
                                    raise ValueError("Photos sharing source observation UUID disagree")
                                output.write(json.dumps({
                                    "original_photo_id":original["original_photo_id"],
                                    "original_observation_uuid":obs,
                                    "original_taxon_id":taxon,
                                    "original_species_name":original["species"],
                                    "original_genus":original["species"].split()[0],
                                    "original_observer":observer,
                                    "original_observed_on":date_string,
                                    "observed_month":dt.month,
                                    "latitude":lat,"longitude":lon,
                                    "latitude_domain":domain,"local_season":season,
                                    "original_ROI_v4_hue":hue,
                                    "actual_photo_contains_classifiable_flower":True,
                                    "original_snapshot_20260927_date_proven":True,
                                },sort_keys=True)+"\n")
                    if len(seen)!=len(ids) or (
                        unit_status["n"]!=rec["n_original_date_known_unverified_bloom_photos_submitted"]
                        or unit_status["class"]!=rec["n_source_roi_four_colour_classifiable_flower_photos"]
                        or unit_status["white"]!=rec["n_source_roi_white_classifiable_photos"]
                        or unit_status["image_acquisition_failed"]!=rec["n_source_image_acquisition_failed"]):
                        raise ValueError("Original ROI verified count/membership ledger not conserved")
                    terminal_n+=unit_status["n"]
                    classified_n+=unit_status["class"]
                    white_n+=unit_status["white"]
                    failure_n+=unit_status["image_acquisition_failed"]
                    statuses.update(unit_status)
                    sub_receipts.append(sha(work))
                packet_hashes.append({
                    "packet":p,"n_date_known_photos":n,
                    "n_bounded_roi_work_units":len(sub_receipts),
                    "sha256_roi_work_unit_receipts":sub_receipts
                })
                db.commit()
            n_observations=db.execute("SELECT COUNT(*) FROM observations").fetchone()[0]
        if (original_n!=dispatch["n_date_known_original_photo_ids_dispatched_to_roi"]
            or terminal_n!=original_n or classified_n>terminal_n
            or white_n>classified_n or n_observations>classified_n):
            raise ValueError("Original date-first photo dispatch and all measured ROI outcomes not conserved")
        if (sum(v["all"] for v in month_zone.values())!=terminal_n
            or sum(v["class"] for v in month_zone.values())!=classified_n
            or sum(v["white"] for v in month_zone.values())!=white_n
            or sum(v["permanent_photo_acquisition_failed"]
                   for v in month_zone.values())!=failure_n):
            raise ValueError("Original photographed month/photo-ROI status conservation failed")
        private_month_doc={
            "schema":"mois_private_original_source_photo_month_roi_effort_v1",
            "status":"SOURCE_COMPLETE_SHARD_ALL_DATE_KNOWN_PHOTO_MONTH_ROI_COUNTS",
            "source_original_roi_shard":dispatch["source_original_shard"],
            "source_original_dated_photo_dispatch_sha256":sha(source_path),
            "n_original_dated_photo_roi_terminals":terminal_n,
            "n_original_roi_four_hue_classifiable":classified_n,
            "n_original_roi_white":white_n,
            "n_permanently_unavailable_photo_pixels":failure_n,
            "latitude_domain_month_photo_effort":[
                {"latitude_domain":zone,"original_month":month,
                 "n_original_date_known_photo_roi_terminals":v["all"],
                 "n_original_four_hue_classifiable":v["class"],
                 "n_original_photo_roi_white":v["white"],
                 "n_source_image_permanent_acquisition_failures":
                     v["permanent_photo_acquisition_failed"],
                 "n_roi_other_unclassifiable":
                     v["all"]-v["class"]-
                     v["permanent_photo_acquisition_failed"]}
                for (zone,month),v in sorted(month_zone.items())
            ],
            "original_photo_identifiers_observers_and_precise_sites_not_in_this_file":True,
            "small_source_month_bins_private_only":True,
            "synthetic_test_only":bool(synthetic),
        }
        private_month_stage.write_text(
            json.dumps(private_month_doc,indent=2,sort_keys=True)+"\n",
            encoding="utf8")
        public_bins=[]
        for (zone,month),v in sorted(month_zone.items()):
            if v["all"]<min_public:
                public_bins.append({
                    "latitude_domain":zone,"original_calendar_month":month,
                    "public_suppressed_below":min_public,
                })
            else:
                public_bins.append({
                    "latitude_domain":zone,"original_calendar_month":month,
                    "n_original_date_known_photo_opportunities":v["all"],
                    "n_original_image_verified_four_hue_classifiable":v["class"],
                    "n_original_classified_white_photos":v["white"],
                    "photo_classifiable_fraction":v["class"]/v["all"],
                    "photo_white_fraction_given_classifiable":(
                        v["white"]/v["class"] if v["class"] else None)
                })
        result={
            "schema":SCHEMA,"status":STATUS,
            "source_original_dated_photo_dispatch_sha256":sha(source_path),
            "source_original_date_queue_sha256":dispatch[
                "source_private_dated_photo_queue_sha256"],
            "shard_number":dispatch["source_original_shard"],
            "original_photo_opportunities_in_source_shard":dispatch[
                "n_original_photo_opportunities_in_sealed_shard"],
            "n_original_date_known_photos_dispatched":original_n,
            "n_original_date_known_photo_ROI_terminals":terminal_n,
            "n_original_bloom_four_hue_classified_photo_pixels":classified_n,
            "n_original_bloom_white_four_hue_classified_photo_pixels":white_n,
            "n_original_distinct_observation_uuids_with_classifiable_flower":n_observations,
            "n_original_permanent_photo_acquisition_failures":failure_n,
            "n_original_image_roi_other_unclassifiable":terminal_n-classified_n-failure_n,
            "n_not_date_known_source_photo_opportunities_not_image_processed":
                dispatch["n_original_photo_opportunities_in_sealed_shard"]-original_n,
            "no_date_valid_original_image_omitted_from_roi_processing":True,
            "old_100543_photo_outcomes_ever_reused":False,
            "observation_level_iNaturalist_Flowers_annotations_used_as_photo_positives":False,
            "one_original_observation_uuid_counted_once_despite_multiple_photos":True,
            "source_photo_per_species_sample_limit":None,
            "n_finalized_original_photo_roi_packets":len(packet_hashes),
            "per_packet_roi_receipt_sha256":packet_hashes,
            "original_photo_observed_month_roi_quality_counts":public_bins,
            "private_all_month_roi_effort_schema":private_month_doc["schema"],
            "private_all_month_roi_effort_file_name":private_month.name,
            "private_all_month_roi_effort_sha256":sha(private_month_stage),
            "protected_positive_bloom_photo_month_ndjson_sha256":sha(stage),
            "this_shard_alone_is_NOT_complete_global_42111_source":True,
            "actual_all_128_original_photo_shards_complete":False,
            "synthetic_test_only":bool(synthetic),
        }
        public_stage.write_text(
            json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf8")
        private_month_stage.replace(private_month)
        stage.replace(out_private_positive)
        public_stage.replace(out_public_receipt)
        sqlite.unlink()
        return result
    except Exception:
        # Intentional quarantine: failures keep private incomplete files
        # for manual inspection and do not promote a fake complete ROI shard.
        raise


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--protected-date-first-packet-root",required=True,type=Path)
    p.add_argument("--protected-sealed-original-roi-terminal-root",required=True,type=Path)
    p.add_argument("--out-private-original-actual-bloom-photo-ndjson",required=True,type=Path)
    p.add_argument("--out-public-original-date-roi-shard-receipt",required=True,type=Path)
    a=p.parse_args()
    r=reassemble(a.protected_date_first_packet_root,
                 a.protected_sealed_original_roi_terminal_root,
                 a.out_private_original_actual_bloom_photo_ndjson,
                 a.out_public_original_date_roi_shard_receipt)
    print(json.dumps({
        "status":r["status"],
        "n_original_date_known_roi_terminal_photos":r[
            "n_original_date_known_photo_ROI_terminals"],
        "n_actual_photo_pixel_classified_flower_hues":r[
            "n_original_bloom_four_hue_classified_photo_pixels"],
        "n_original_actual_white_bloom_photos":r[
            "n_original_bloom_white_four_hue_classified_photo_pixels"]
    },sort_keys=True))


if __name__=="__main__":main()
