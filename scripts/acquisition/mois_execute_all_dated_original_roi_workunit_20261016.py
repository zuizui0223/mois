#!/usr/bin/env python3
"""Run one bounded DATE-FIRST original FCP photo ROI-v4 image batch.

All original 42,111 taxa are eligible with no per-species photo limit.
Packets are made ONLY from original observed_on-valid complete FCP photo
shards. Unlike FCP's separate 93M PHOTO colour-polymorphism study, Mois's
seasonal target does not need to download images whose ORIGINAL source
observed_on is missing. The metadata universe and exclusion ledger are
retained without sampling or presumed bloom status.

The actual unchanged original FCP ROI-v4, YOLO detector, EfficientSAM source
and model hashes must pass their original method fingerprint attestation.
Measurement backend gets only hashed measurement IDs, image filenames, and
image pixel bytes; no taxon/genus, source dates, climate or coordinates.
The provider must have coordinated global rate limiting and licensed access.

An unsuccessful/transient provider network request aborts the work batch;
only the source backend's explicit permanent original-image failure class
may be sealed as image acquisition failure. Files are atomically promoted
only after full 64-image compute-unit count and result validation.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
import tempfile
from pathlib import Path

import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from scripts.acquisition.mois_prepare_all_dated_original_roi_packets_20261016 import (
    WORKER,ACQUISITION,PACKET_SCHEMA,PACKET_STATUS,
    SCHEMA as DISPATCH_SCHEMA,STATUS as DISPATCH_STATUS,
    IMAGE_COMPUTE_BATCH,blind_id
)
from scripts.acquisition.mois_full_source_dated_flower_candidates_20261016 import sha
from scripts.acquisition.mois_execute_sample_roi_v4_20261016 import (
    backend_from_fcp,FCP_BACKEND_GIT_BLOB,VALID_ID
)

SCHEMA="mois_uncapped_original_dated_photo_roi_v4_work_unit_v1"
STATUS="COMPLETE_TRUE_ORIGINAL_DATED_PHOTO_PIXEL_ROI_WORK_UNIT_NOT_GLOBAL"
BIO=("white","yellow","orange","red","pink","magenta","purple","blue","bronze")


def read_batch(packet_dir:Path,dispatch_dir:Path,
               packet_number:int,batch:int)->tuple[dict,list[dict],list[dict]]:
    if packet_number<0 or batch<0:
        raise ValueError("ROI packet indices must be nonnegative")
    dispatch_file=dispatch_dir/"full_dated_roi_dispatch_receipt.json"
    source=json.loads(dispatch_file.read_text(encoding="utf8"))
    if (source.get("schema")!=DISPATCH_SCHEMA or source.get("status")!=DISPATCH_STATUS
        or source.get("n_original_photo_ids_omitted_despite_valid_source_date")!=0
        or source.get("all_photo_ids_with_source_date_included_without_sampling_cap") is not True
        or source.get("original_flower_roi_verified") is True
        or source.get("n_original_photo_pixels_measured_at_dispatch")!=0
        or source.get("n_original_photo_ids_omitted_despite_valid_source_date")!=0
        or packet_number>=source.get("n_bounded_original_roi_packets",-1)
        or source.get("image_compute_batch_size_not_sample_cap")!=IMAGE_COMPUTE_BATCH):
        raise ValueError("Not a full outcome-blind original date-source image packet dispatch")
    stamp=next((x for x in source["per_original_packet"]
               if x["packet_index"]==packet_number),None)
    if stamp is None or stamp["packet_receipt_sha256"]!=sha(packet_dir/"packet_receipt.json"):
        raise ValueError("Original source packet receipt SHA does not match dispatch")
    doc=json.loads((packet_dir/"packet_receipt.json").read_text(encoding="utf8"))
    if (doc.get("schema")!=PACKET_SCHEMA or doc.get("status")!=PACKET_STATUS
        or doc.get("packet")!=packet_number
        or doc.get("original_flower_roi_verified") is not False
        or doc.get("original_photo_selection_cap","INVALID") is not None
        or doc.get("original_photo_pixels_opened") is not False
        or doc.get("worker_has_observation_location_species_or_hue") is not False
        or doc.get("source_original_date_queue_sha256")!=
            source.get("source_private_dated_photo_queue_sha256")
        or doc.get("n_source_original_date_known_photos")!=stamp["n_source_dated_photos"]):
        raise ValueError("Original date-source image packet was already measured or reselected")
    source_packets=[]
    for name,fields in (("worker.csv",WORKER),("acquisition.csv",ACQUISITION)):
        if sha(packet_dir/name)!=doc["file_sha256"].get(name):
            raise ValueError("Original photo image packet bytes drifted")
        with (packet_dir/name).open("r",newline="",encoding="utf8") as f:
            reader=csv.DictReader(f)
            if tuple(reader.fieldnames or ())!=fields:
                raise ValueError("Original FCP worker location-blind columns do not match")
            rows=list(reader)
            if len(rows)!=doc["n_source_original_date_known_photos"]:
                raise ValueError("Original dated ROI packet incomplete")
            source_packets.append(rows)
    if sha(packet_dir/"private_source.csv")!=doc["file_sha256"].get("private_source.csv"):
        raise ValueError("Private original date/photo ID mapping corrupted")
    workers,acquisition=source_packets
    if len({v["measurement_id"] for v in workers})!=len(workers):
        raise ValueError("Blinded original photo IDs collide")
    for w,a in zip(workers,acquisition):
        if (not VALID_ID.fullmatch(w["measurement_id"])
            or w["image_filename"]!=w["measurement_id"]+".jpg"
            or not w["photo_license"]
            or any(w[x]!=a[x] for x in WORKER)
            or not a["photo_url_large"].startswith(
                "https://inaturalist-open-data.s3.amazonaws.com/photos/")):
            raise ValueError("Original photo blind image worker acquired source identity or shifted")
    pos=batch*IMAGE_COMPUTE_BATCH
    if pos>=len(workers):
        raise ValueError("Image batch number beyond original dated source packet")
    return doc,workers[pos:pos+IMAGE_COMPUTE_BATCH],acquisition[pos:pos+IMAGE_COMPUTE_BATCH]


def execute(packet_dir:Path,dispatch_dir:Path,
            output_root:Path,packet_number:int,batch:int,*,
            backend,source_script:Path,original_detector:Path,
            original_sam:Path,synthetic=False)->dict:
    original_receipt,workers,acq=read_batch(
        packet_dir,dispatch_dir,packet_number,batch)
    if not synthetic and not os.environ.get("FCP_PROVIDER_GLOBAL_RATE_FILE"):
        raise ValueError("Original image provider requires coordinated global pacing/access")
    fingerprint=backend.verify_frozen_method(source_script,original_detector,original_sam)
    if synthetic:
        if not fingerprint.get("synthetic_test_only"):
            raise ValueError("Mock ROI source cannot masquerade as certified real image method")
    elif (fingerprint.get("location_blind") is not True or
          fingerprint.get("no_old_measurement_csv_inputs") is not True):
        raise ValueError("Original frozen FCP ROI-v4 evaluator method provenance missing")
    target=output_root/f"photo_roi_packet_{packet_number:06d}_batch_{batch:03d}"
    stage=output_root/(target.name+".INCOMPLETE")
    if target.exists() or stage.exists():
        raise ValueError("No overwrite of a previous sealed or unresolved image batch")
    output_root.mkdir(parents=True,exist_ok=True)
    stage.mkdir()
    with tempfile.TemporaryDirectory(prefix="mois_real_date_positive_roi_") as tmp:
        tmp=Path(tmp);pixels=tmp/"images";pixels.mkdir()
        okay=[];failed=[]
        for ac in acq:
            try:
                body=backend.image_loader(ac["photo_url_large"])
                if not body:
                    raise ValueError("Empty original image data; source image request not authoritative")
                (pixels/ac["image_filename"]).write_bytes(body)
                okay.append({key:ac[key] for key in WORKER})
            except backend.PermanentOriginalPhotoError as exc:
                failed.append({
                    "measurement_id":ac["measurement_id"],
                    "morph":"mixed_uncertain",
                    "measurement_status":"image_acquisition_failed",
                    "failure_reasons":type(exc).__name__+":"+str(exc)[:160],
                    **{"flower_fraction_"+c:0. for c in BIO}
                })
        measured=pd.DataFrame(columns=[
            "measurement_id","morph","measurement_status",
            *["flower_fraction_"+c for c in BIO]])
        if okay:
            blind=tmp/"only_image_pixels_worker.csv"
            pd.DataFrame(okay,columns=list(WORKER)).to_csv(blind,index=False)
            measured=backend.evaluate_with_frozen_source(
                blind,pixels,source_script,original_detector,
                original_sam,tmp/"original_roi_hue.csv")
        all_results=pd.concat(
            [measured,pd.DataFrame(failed)],ignore_index=True,sort=False)
        expected=[row["measurement_id"] for row in workers]
        backend.validate_terminal(all_results,expected)
        all_results=all_results.sort_values("measurement_id",kind="stable")
        end=stage/"original_fcp_date_valid_photo_ROI_v4_terminal.csv"
        all_results.to_csv(end,index=False)
        passed=int(all_results.measurement_status.eq(
            "classified_four_state_morph").sum())
        white=int(all_results.morph.eq("white").sum())
        receipt={
            "schema":SCHEMA,"status":STATUS,
            "source_packet_number":packet_number,
            "source_image_compute_batch":batch,
            "n_original_date_known_unverified_bloom_photos_submitted":len(workers),
            "n_source_roi_four_colour_classifiable_flower_photos":passed,
            "n_source_roi_white_classifiable_photos":white,
            "n_source_image_acquisition_failed":len(failed),
            "n_classification_missing_or_ambiguous":len(workers)-passed-len(failed),
            "source_packet_receipt_sha256":sha(packet_dir/"packet_receipt.json"),
            "source_dispatch_receipt_sha256":sha(
                dispatch_dir/"full_dated_roi_dispatch_receipt.json"),
            "original_terminal_csv_sha256":sha(end),
            "original_source_frozen_roi_method":fingerprint,
            "original_fcp_image_backend_git_blob":(
                "SYNTHETIC_TEST_ONLY" if synthetic else FCP_BACKEND_GIT_BLOB),
            "original_image_ids_or_photo_geography_exposed_to_classifier":False,
            "all_source_date_known_photos_covered_by_this_unit":False,
            "missing_photo_flower_hues_imputed":False,
            "no_historical_image_measurements_reused":True,
            "original_fcp_image_pixels_saved_after_evaluation":False,
            "synthetic_test_only":bool(synthetic),
        }
        (stage/"original_dated_roi_work_unit_receipt.json").write_text(
            json.dumps(receipt,indent=2,sort_keys=True)+"\n",encoding="utf8")
    stage.replace(target)
    return receipt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--private-original-photo-roi-packet-dir",required=True,type=Path)
    p.add_argument("--private-original-date-dispatch-root",required=True,type=Path)
    p.add_argument("--private-original-roi-terminals-root",required=True,type=Path)
    p.add_argument("--packet-number",required=True,type=int)
    p.add_argument("--batch-number",required=True,type=int)
    p.add_argument("--fcp-source-repository",required=True,type=Path)
    p.add_argument("--original-fcp-roi-measure-script",required=True,type=Path)
    p.add_argument("--original-fcp-yolo-detector",required=True,type=Path)
    p.add_argument("--original-fcp-efficient-sam",required=True,type=Path)
    p.add_argument("--provider-access-reviewed",action="store_true")
    a=p.parse_args()
    if not a.provider_access_reviewed:
        p.error("Confirm source photo image provider rights, pacing and compute budget")
    backend=backend_from_fcp(a.fcp_source_repository)
    result=execute(a.private_original_photo_roi_packet_dir,
                   a.private_original_date_dispatch_root,
                   a.private_original_roi_terminals_root,
                   a.packet_number,a.batch_number,backend=backend,
                   source_script=a.original_fcp_roi_measure_script,
                   original_detector=a.original_fcp_yolo_detector,
                   original_sam=a.original_fcp_efficient_sam)
    print(json.dumps({k:result[k] for k in (
        "status","source_packet_number","source_image_compute_batch",
        "n_original_date_known_unverified_bloom_photos_submitted",
        "n_source_roi_four_colour_classifiable_flower_photos"
    )},sort_keys=True))


if __name__=="__main__":main()
