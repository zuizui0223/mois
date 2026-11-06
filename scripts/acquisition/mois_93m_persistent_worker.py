#!/usr/bin/env python3
"""Long-lived, restartable medium-image scanner for the entire MOIS photo census.

Replace the previous pattern of one cold GitHub Actions worker per ~1,000
photographs. A worker can scan MANY independent 4,096-photo signed chunks
with a single frozen YOLO instance, one shared network pool, and one reused
image-decode thread executor. No network requests occur on import or plan mode.

Census inclusion = ALL of the 93,209,351 source photo IDs. Chunks, worker
slots, season, species and original observations are PROCESSING partitions,
not biological/photo selection rules. Source negatives, errors and quality
unknowns are retained without reinterpreting them as petal colours.

A user-supplied *aggregate* source request budget is divided by worker slots.
This is an engineering ceiling, not a license to exceed provider policies.
Deploy cloud workers with one authoritative slot each and near the S3 region.
"""
from __future__ import annotations
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import time

from scripts.acquisition.mois_full93m_medium_pixel_scan import (
    HOST, MAX_IMAGE_BYTES, Pacer, ProviderBackpressure,
    WarmFlowerBoxDetector, sha_file,
)
from scripts.acquisition.mois_full93m_source_part_chunks import SCHEMA as CHUNK_SCHEMA
from scripts.acquisition.mois_full93m_signed_chunk_pixel_scan import run as scan_chunk

SCHEMA = "mois_full93m_long_lived_500px_photo_scanning_worker_v1"
EXPECTED = 93_209_351

def require(ok, msg):
    if not ok:
        raise ValueError("MOIS_PERSISTENT_IMAGE_WORKER: " + msg)

class PooledOpenDataDownloader:
    """Connection-pool reuse across 500px photos, with one slot's pace budget.

    All photo bytes remain transient in memory, never written to the result
    store. A provider 403/429/503 stops the entire worker rather than
    silently marking thousands of photos as biological absences.
    """
    def __init__(self, per_worker_rps: float, io_threads: int=12,
                 *, pool=None, sleep=time.sleep, retries=2):
        require(math.isfinite(per_worker_rps) and 0<per_worker_rps,
                "Source request budget must be positive and finite")
        require(1<=io_threads<=64 and 0<=retries<=4,
                "Invalid connection pool or retry limit")
        self.pacer = Pacer(per_worker_rps)
        self.sleep = sleep
        self.retries = retries
        self.io_threads = io_threads
        if pool is None:
            import urllib3
            self.pool = urllib3.PoolManager(
                num_pools=2, maxsize=io_threads, block=True,
                retries=False,
                timeout=urllib3.Timeout(connect=8, read=20),
                headers={"User-Agent":"MOIS-persistent-OpenData-licensed-flowers/1.0"},
            )
        else:
            self.pool = pool
        self.received=Counter()

    def __call__(self, url):
        require(isinstance(url,str) and url.startswith(HOST)
                and "/medium." in url, "Only signed iNaturalist Open Data medium images")
        last_error=None
        for attempt in range(self.retries+1):
            self.pacer.wait()
            response=None
            try:
                response=self.pool.request(
                    "GET",url,preload_content=False,redirect=False,retries=False,
                )
                status=int(response.status)
                self.received["http_"+str(status)]+=1
                if status in (403,429,503):
                    raise ProviderBackpressure(
                        "Open Data S3 "+str(status)+
                        ": stop this worker and diagnose aggregate request budget")
                if status in (404,410):
                    raise RuntimeError("original_photo_no_longer_available_http_"+str(status))
                if status in (500,502,504) and attempt<self.retries:
                    last_error=RuntimeError("transient_S3_http_"+str(status))
                elif status != 200:
                    raise RuntimeError("original_image_http_"+str(status))
                else:
                    header=response.headers.get("Content-Length")
                    if header and header.isdecimal() and int(header)>MAX_IMAGE_BYTES:
                        raise RuntimeError("medium_image_exceeds_limit")
                    data=response.read(MAX_IMAGE_BYTES+1)
                    if not data or len(data)>MAX_IMAGE_BYTES:
                        raise RuntimeError("empty_or_oversized_medium_original_image")
                    self.received["downloaded_photo_count"]+=1
                    self.received["downloaded_medium_bytes"]+=len(data)
                    return data
            except ProviderBackpressure:
                raise
            except Exception as exc:
                if isinstance(exc,RuntimeError):
                    # Permanent photo 404/410 and unexpected HTTP remain unknown.
                    raise
                last_error=exc
            finally:
                if response is not None:
                    try: response.release_conn()
                    except Exception: response.close()
            if attempt<self.retries:
                self.sleep(min(8, 2**attempt))
        raise RuntimeError("transient_medium_image_download_exhausted") from last_error

    def close(self):
        if hasattr(self.pool,"clear"):
            self.pool.clear()

def load_roots(root_paths):
    """Exact source metadata coverage, not actual decoded/flower counts."""
    seen=set()
    roots=[]
    for raw in root_paths:
        root=Path(raw)
        receipt_path=root/"chunk_receipt.json"
        receipt=json.loads(receipt_path.read_text(encoding="utf-8"))
        require(receipt.get("schema")==CHUNK_SCHEMA
                and receipt.get("source_fully_partitioned_once_into_nonoverlapping_photo_chains") is True
                and receipt.get("image_pixels_downloaded")==0,
                "Source chunk lineage must be signed and metadata-only")
        key=(receipt["original_tier"],receipt["original_part"])
        require(key not in seen,"Repeated original source tier/part across roots")
        seen.add(key)
        manifest=root/"chunk_manifest.csv"
        require(sha_file(manifest)==receipt["manifest_sha256"],
                "Original source chunk manifest SHA changed")
        with manifest.open(newline="",encoding="utf-8") as handle:
            rows=list(csv.DictReader(handle))
        require(len(rows)==receipt["chunks"] and
                sorted(int(x["chunk_index"]) for x in rows)==list(range(len(rows))),
                "Original chunk manifest missing or duplicated rows")
        require(sum(int(x["expected_photo_ids"]) for x in rows)
                ==receipt["actual_original_photo_IDs"],
                "Original photo ID source count drift across chunks")
        roots.append((root,receipt,rows))
    roots.sort(key=lambda entry:(entry[1]["original_tier"],entry[1]["original_part"]))
    return roots

def slot_for(tier,part,chunk,slots):
    """Stable SHA partition. All worker slots together cover every photo once."""
    require(1<=slots<=4096,"Worker slots outside operational supported range")
    digest=hashlib.sha256(f"mois93m-long-lived-v1|{tier}|{part}|{chunk}".encode()).digest()
    return int.from_bytes(digest[:8],"big")%slots

def plan(root_paths,*,slot=0,slots=1):
    require(isinstance(slot,int) and 0<=slot<slots<=4096,"Invalid worker slot partition")
    roots=load_roots(root_paths)
    all_records=sum(r["actual_original_photo_IDs"] for _,r,_ in roots)
    all_tasks=sum(len(rows) for _,_,rows in roots)
    assigned=[]
    for root,receipt,rows in roots:
        tier=receipt["original_tier"]; part=receipt["original_part"]
        for row in rows:
            chunk=int(row["chunk_index"])
            if slot_for(tier,part,chunk,slots)==slot:
                assigned.append((root,receipt,chunk,int(row["expected_photo_ids"])))
    return {
        "schema": SCHEMA,
        "source_photo_ids_in_provided_signed_roots":all_records,
        "signed_source_part_roots":len(roots),
        "all_photo_chunks_in_provided_roots":all_tasks,
        "assigned_chunks_to_this_slot":len(assigned),
        "assigned_original_photo_ids_to_this_slot":sum(x[3] for x in assigned),
        "worker_slot":slot,"total_worker_slots":slots,
        "metadata_covers_all_93209351_photo_ids":(
            len(roots)==14*128 and all_records==EXPECTED),
        "actual_photos_downloaded_by_planning":0,
        "source_photo_selection_probability":1,
        "assigned":assigned,
    }

def inspect_completed_chunk(sub, source_root, source_receipt, chunk_index,
                            expected, detector_identity):
    """Reuse a previously *fully* signed photo chunk, even at a new packet width.

    This removes costly rereads and avoids re-running original signed 500px
    photographs merely because the warm-worker batch size changed.
    """
    complete=sub/"complete_receipt.json"
    if not complete.exists():
        return None
    from scripts.acquisition.mois_full93m_medium_pixel_scan import SCHEMA as PIXEL_SCHEMA
    r=json.loads(complete.read_text(encoding="utf-8"))
    require(r.get("schema")==PIXEL_SCHEMA
            and r.get("source_role")=="full93m_chunk"
            and r.get("source_parent_sha256")==source_receipt["original_source_part_sha256"]
            and r.get("source_chunk_index")==chunk_index
            and r.get("source_tier")==source_receipt["original_tier"]
            and r.get("source_part")==source_receipt["original_part"]
            and r.get("detector_identity")==detector_identity
            and r.get("expected_photo_ids")==expected
            and r.get("actual_unique_source_photo_ids_attempted")==expected
            and r.get("all_original_input_photos_attempted") is True
            and r.get("entire_original_source_part_completed") is False
            and sum(r.get("screen_counts",{}).values())==expected,
            "Previously completed pixels differ from signed original photo/model")
    with (source_root/"chunk_manifest.csv").open(newline="",encoding="utf-8") as handle:
        rows=list(csv.DictReader(handle))
    part=next((row for row in rows if int(row["chunk_index"])==chunk_index),None)
    require(part is not None and
            r["source_gzip_sha256"]==part["chunk_gzip_sha256"],
            "Previously completed scan chunk source SHA mutated")
    seen=0
    packets=sorted(sub.glob("packet-*.csv.gz"))
    require(len(packets)==r["packets_sealed"] and packets,
            "Missing immutable photo packet for prior completed source chunk")
    for data in packets:
        proof=data.with_suffix("").with_suffix(".receipt.json")
        require(proof.exists(),"Previously completed source packet receipt missing")
        d=json.loads(proof.read_text(encoding="utf-8"))
        require(d.get("output_sha256")==sha_file(data)
                and d.get("source_gzip_sha256")==part["chunk_gzip_sha256"]
                and d.get("source_parent_sha256")==source_receipt["original_source_part_sha256"]
                and d.get("source_chunk_index")==chunk_index
                and d.get("detector_identity")==detector_identity
                and d.get("photo_ids",0)>0,
                "Prior completed original photo pixel/ROI packet cannot be reused")
        seen+=d["photo_ids"]
    require(seen==expected,"Prior photo-ID attempts missing from source packets")
    return {
        "actual_new_or_resumed_medium_photo_attempts":expected,
        "medium_decoded_photo_count":r["medium_image_pixels_decoded"],
        "screen_counts":r["screen_counts"],
        "source_chunk_all_photo_ids_scanned":True,
        "source_image_packet_reread_not_inferred":True,
    }

def execute(root_paths,output,*,weights,weight_sha,slot=0,slots=1,
            approved_aggregate_rps,download_workers=12,batch_size=512,
            max_chunks=None,device="cpu",detector=None,downloader=None,
            prefetch_next_packet=False,feature_mode="legacy_hsv"):
    """Process many source chunks in one model/network process.

    The aggregate limit is divided across exactly `slots` planned workers.
    This method never launches other workers or allocates cloud resources.
    Resume verifies every existing packet SHA, source/photo-ID and model pin.
    """
    require(math.isfinite(approved_aggregate_rps) and 0<approved_aggregate_rps,
            "An explicit aggregate Open Data provider request budget is required")
    require(1<=download_workers<=64 and 16<=batch_size<=1024,
            "Chunk processing bounds invalid")
    require(max_chunks is None or 1<=max_chunks<=1000000,
            "Invalid optional bounded maximum chunks")
    p=plan(root_paths,slot=slot,slots=slots)
    per_worker_rps=approved_aggregate_rps/slots
    require(per_worker_rps>0,"Aggregate request limit divided by zero")
    output=Path(output)
    output.mkdir(parents=True,exist_ok=True)
    is_real_model=detector is None
    if detector is None:
        detector=WarmFlowerBoxDetector(weights,weight_sha,0,"item",
                                      device=device,conf=.02)
    if downloader is None:
        downloader=PooledOpenDataDownloader(per_worker_rps,download_workers)
    selected=p["assigned"][:max_chunks] if max_chunks is not None else p["assigned"]
    completed=Counter(); observed=Counter()
    reused_completed=0
    start=time.monotonic()
    model_identity=(
        "yolo11:"+weight_sha+":item:0.02" if isinstance(detector,WarmFlowerBoxDetector)
        else "test_injected_detector:"+type(detector).__name__
    )
    try:
        with ThreadPoolExecutor(max_workers=download_workers) as executor:
            for k,(root,receipt,idx,expected) in enumerate(selected):
                sub=output/(f"tier-{receipt['original_tier']:02d}")/(
                    f"part-{receipt['original_part']:03d}")/(f"chunk-{idx:05d}")
                sub.parent.mkdir(parents=True,exist_ok=True)
                r=inspect_completed_chunk(
                    sub,root,receipt,idx,expected,model_identity)
                if r is not None:
                    reused_completed+=expected
                else:
                    r=scan_chunk(root,sub,idx,detector=detector,
                                 download=downloader,workers=download_workers,
                                 batch_size=batch_size,io_executor=executor,
                                 prefetch_next_packet=prefetch_next_packet,
                             feature_mode=feature_mode)
                require(r["source_chunk_all_photo_ids_scanned"] is True
                        and r["actual_new_or_resumed_medium_photo_attempts"]==expected,
                        "Missing actual source photo pixels in completed chunk")
                completed["original_photo_attempts"]+=expected
                completed["medium_decoded"]+=r["medium_decoded_photo_count"]
                completed["chunks_completed"]+=1
                observed.update(r["screen_counts"])
                if (k+1)%8==0 or k+1==len(selected):
                    result={
                        "schema":SCHEMA,
                        "source_original_photo_ids_in_signed_roots":
                            p["source_photo_ids_in_provided_signed_roots"],
                        "worker_slot":slot,"worker_slots":slots,
                        "approved_aggregate_source_requests_per_second":
                            approved_aggregate_rps,
                        "source_rate_per_active_worker_slot":per_worker_rps,
                        "assigned_original_photo_ids":p["assigned_original_photo_ids_to_this_slot"],
                        "chunks_selected_this_process":len(selected),
                        "actual_source_photo_attempts_in_verified_completed_chunks":
                            completed["original_photo_attempts"],
                        "medium_photo_pixels_decoded":completed["medium_decoded"],
                        "previously_finished_original_source_photo_pixels_reused":
                            reused_completed,
                        "new_source_photo_attempts_this_process":
                            completed["original_photo_attempts"]-reused_completed,
                        "chunks_source_sha_verified_complete":completed["chunks_completed"],
                        "technical_statuses":dict(observed),
                        "remaining_source_photo_ids_for_this_slot":
                            p["assigned_original_photo_ids_to_this_slot"]-completed["original_photo_attempts"],
                        "photo_pixels_acquired_are_automatically_verified_petals":False,
                        "full_93m_pixel_census_complete":False,
                        "detector_sha256":weight_sha,
                        "original_roi_classifications_changed":False,
                        "raw_source_images_written_to_output":False,
                        "next_packet_transfer_overlapped_with_yolo":prefetch_next_packet,
                        "colour_feature_mode":feature_mode,
                        "elapsed_process_seconds":round(time.monotonic()-start,3),
                    }
                    tmp=output/(f"slot-{slot:04d}-status.json.partial")
                    tmp.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
                    os.replace(tmp,output/(f"slot-{slot:04d}-status.json"))
                    print("MOIS_LONG_LIVED_SOURCE_PIXEL_PROGRESS",
                          json.dumps(result,sort_keys=True),flush=True)
        return result if selected else dict(p,scan_complete=False)
    finally:
        if hasattr(downloader,"close"):
            downloader.close()

def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("--chunks-roots-file",type=Path,required=True,
                   help="One already source-SHA-verified chunk-root directory per line")
    a.add_argument("--output",type=Path,required=True)
    a.add_argument("--weights",type=Path)
    a.add_argument("--weights-sha256")
    a.add_argument("--worker-slot",type=int,default=0)
    a.add_argument("--worker-slots",type=int,default=1)
    a.add_argument("--approved-aggregate-rps",type=float)
    a.add_argument("--download-workers",type=int,default=12)
    a.add_argument("--batch-size",type=int,default=512)
    a.add_argument("--max-chunks",type=int)
    a.add_argument("--device",default="cpu")
    a.add_argument("--colour-feature-mode",default="legacy_hsv",
                   choices=("legacy_hsv","floral_box_lab_v2"),
                   help="Use sRGB CIELAB floral-box major-colour descriptor in same photo pass")
    a.add_argument("--overlap-next-packet",action="store_true",
                   help="Prefetch next <=256 photo packet while detecting current; bounded RAM")
    a.add_argument("--execute",action="store_true",
                   help="Required to actually access public AWS S3 photos")
    x=a.parse_args()
    roots=[str((x.chunks_roots_file.parent / line.strip()).resolve())
           for line in x.chunks_roots_file.read_text().splitlines()
           if line.strip() and not line.lstrip().startswith("#")]
    require(bool(roots),"At least one signed source root required")
    p=plan(roots,slot=x.worker_slot,slots=x.worker_slots)
    if not x.execute:
        for k,value in p.items():
            if k=="assigned":continue
            print("MOIS_93M_PERSISTENT_WORKER_PLAN",k,value)
        return
    require(x.weights is not None and x.weights_sha256 is not None
            and x.approved_aggregate_rps is not None,
            "Full photo access requires signed model and explicit aggregate budget")
    execute(roots,x.output,weights=x.weights,weight_sha=x.weights_sha256,
            slot=x.worker_slot,slots=x.worker_slots,
            approved_aggregate_rps=x.approved_aggregate_rps,
            download_workers=x.download_workers,batch_size=x.batch_size,
            max_chunks=x.max_chunks,device=x.device,
            prefetch_next_packet=x.overlap_next_packet,
            feature_mode=x.colour_feature_mode)

if __name__=="__main__":
    main()
