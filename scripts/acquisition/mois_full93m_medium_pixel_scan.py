#!/usr/bin/env python3
"""Actually read medium-size pixels for every original MOIS 93m photo-ID.

One worker scans ONE signed disjoint tier/part gzip. Warm flower-box detection
runs on 500px images in batches; full-resolution EfficientSAM is NOT called.
Source observations are not prefiltered by taxa, native status or phenology.
A screen-negative is NOT botanically certified non-flowering. No screen result
is allowed to replace the original four-hue ROI labels. All download/decoder
and detection errors retain their photo ID and retryable UNKNOWN status.

Restart skips only SHA-verified immutable completed packets. This program does
not launch a 93m S3 crawl merely by import, file creation or unit testing.
"""
from __future__ import annotations
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from contextlib import nullcontext
import csv
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import tempfile
import threading
import time
import urllib.error
import urllib.request

from scripts.acquisition.mois_onepass_full93m_target_tiers import (
    COLUMNS, normalized_extension,
)
from scripts.acquisition.mois_flower_box_lab_colour import (
    PROTOCOL as BOX_LAB_PROTOCOL, DESCRIPTOR_COLUMNS as BOX_LAB_COLUMNS,
    flower_box_lab_descriptor,
)

SCHEMA = "mois_all_original_photo_medium_pixel_stage1_terminal_v1"
SOURCE_ARM = "full93m_observation_and_photo_complete_source"
COLUMNS_OUT = (
    "photo_id", "observation_id", "inat_taxon_id", "observed_on",
    "latitude", "longitude", "source_part", "source_tier", "photo_license",
    "screen_status", "technical_reason", "image_pixel_sha256",
    "downloaded_bytes", "image_width", "image_height",
    "mean_visible_brightness", "flower_box_score", "flower_box_xyxy",
    "box_proxy_white", "box_proxy_yellow_orange", "box_proxy_red_pink",
    "box_proxy_blue_purple", "box_proxy_other",
    "source_arm", "source_photo_inclusion_probability", "audit_stratum",
)
# v1 completed packets are immutable: Lab colours live in a SEPARATELY VERSIONED
# packet schema/output path. Never rewrite a legacy packet to add new columns.
FEATURE_MODES=("legacy_hsv","floral_box_lab_v2")
LAB_SCHEMA="mois_all_original_photo_medium_pixel_stage1_CIELAB_dominant_v2"
COLUMNS_LAB=COLUMNS_OUT+BOX_LAB_COLUMNS
ALLOWED = {"candidate", "negative", "uncertain", "error"}
HOST = "https://inaturalist-open-data.s3.amazonaws.com/photos/"
MAX_IMAGE_BYTES = 3_000_000

class ProviderBackpressure(RuntimeError):
    """Stop this photo shard on S3 rate-limit feedback; do not hammer provider."""

def require(ok, why):
    if not ok:
        raise ValueError("MOIS_FULL93M_PIXEL_SCAN: " + why)

def sha_bytes(data):
    return hashlib.sha256(data).hexdigest()

def sha_file(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(4 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def source_url(photo_id, extension):
    require(str(photo_id).isdecimal(), "Invalid photo ID")
    ext = normalized_extension(extension)
    return None if ext is None else HOST + str(photo_id) + "/medium." + ext

def audit_stratum(latitude, observation_date):
    try:
        lat = float(latitude)
        require(math.isfinite(lat) and -90 <= lat <= 90, "Bad latitude")
        zone = ("N" if lat >= 0 else "S") + str(min(2, int(abs(lat) // 30)))
    except (ValueError, TypeError):
        zone = "unknown_lat"
    d = str(observation_date)
    month = int(d[5:7]) if len(d) >= 7 and d[5:7].isdecimal() else 0
    season = (
        "DJF" if month in (12, 1, 2) else
        "MAM" if month in (3, 4, 5) else
        "JJA" if month in (6, 7, 8) else
        "SON" if month in (9, 10, 11) else "unknown_season"
    )
    return zone + "|" + season

class Pacer:
    """Per worker pace only; cluster must also enforce a shared S3 budget."""
    def __init__(self, requests_per_second):
        require(requests_per_second > 0, "Request rate must be positive")
        self.delay = 1 / requests_per_second
        self.lock = threading.Lock()
        self.next_at = 0.0

    def wait(self):
        with self.lock:
            now = time.monotonic()
            time_to_wait = max(0.0, self.next_at - now)
            self.next_at = max(self.next_at, now) + self.delay
        if time_to_wait:
            time.sleep(time_to_wait)

def download_medium(url, pacer, *, timeout=20, retries=2,
                    max_bytes=MAX_IMAGE_BYTES):
    """Photo comes from signed OPEN DATA S3, never the rate-limited iNat API."""
    last = None
    for attempt in range(retries + 1):
        pacer.wait()
        try:
            req = urllib.request.Request(url, headers={
                "User-Agent": "MOIS-public-licensed-open-data-research/1.0"
            })
            with urllib.request.urlopen(req, timeout=timeout) as stream:
                payload = stream.read(max_bytes + 1)
                require(len(payload) <= max_bytes,
                        "Media exceeds configured bounded image size")
                require(payload, "Empty original source image")
                return payload
        except (urllib.error.HTTPError, urllib.error.URLError,
                TimeoutError, OSError, ValueError) as exc:
            last = exc
            status = getattr(exc, "code", None)
            if status in (403, 429, 503):
                raise ProviderBackpressure(
                    f"HTTP {status} from Open Data S3: stop this shard "
                    "and diagnose provider access / global request rate "
                    "before more photographs are attempted") from exc
            # 404/410 object deletions stay explicit errors, not fake flower negatives.
            if status in (400, 401, 404, 410, 451):
                break
            if attempt < retries:
                time.sleep(min(15, 2 ** attempt))
    raise RuntimeError("image_access_failed: " + type(last).__name__)

def decode_medium(payload):
    from PIL import Image, ImageOps, ImageStat
    Image.MAX_IMAGE_PIXELS = 12_000_000
    with Image.open(io.BytesIO(payload)) as loaded:
        require(loaded.format in ("JPEG", "PNG"), "Unsupported decoded photo")
        im = ImageOps.exif_transpose(loaded)
        im.load()
        rgb = im.convert("RGB")
    # The S3 medium variant is max 500px; error rather than silently use huge originals.
    require(1 <= rgb.width <= 600 and 1 <= rgb.height <= 600,
            "Unexpected source medium resolution")
    small = rgb.resize((64, 64))
    brightness = sum(ImageStat.Stat(small).mean) / (3 * 255)
    return rgb, round(brightness, 6)

def proxy_from_box(image, xyxy):
    """VERY rough RGB/HSV distribution inside a detected box, NOT petal colour."""
    import numpy as np
    a, b, c, d = map(float, xyxy)
    l = max(0, min(image.width - 1, int(a)))
    t = max(0, min(image.height - 1, int(b)))
    r = max(l + 1, min(image.width, int(c)))
    bot = max(t + 1, min(image.height, int(d)))
    hsv = np.asarray(image.crop((l, t, r, bot)).resize((80, 80)).convert("HSV"))
    h = hsv[:, :, 0].astype(float) * 360 / 255
    s = hsv[:, :, 1]
    v = hsv[:, :, 2]
    white = (s < 60) & (v >= 175)
    pigment = (s >= 65) & (v >= 65)
    yellow_orange = pigment & (h >= 20) & (h < 75)
    red_pink = pigment & ((h < 20) | (h >= 305))
    blue_purple = pigment & (h >= 185) & (h < 305)
    fractions = [white.mean(), yellow_orange.mean(),
                 red_pink.mean(), blue_purple.mean()]
    fractions.append(max(0.0, 1.0 - sum(fractions)))
    return [format(float(p), ".6g") for p in fractions]

class WarmFlowerBoxDetector:
    """Reusable YOLO11 detector only; do NOT initialize EfficientSAM here."""
    def __init__(self, weights, expected_sha256, class_id, expected_class_name,
                 device="cpu", conf=0.02):
        require(sha_file(weights) == expected_sha256,
                "Flower detector weight SHA mismatch")
        require(0 < conf <= 0.2, "High-recall candidate screen requires low confidence")
        from ultralytics import YOLO
        self.model = YOLO(str(weights))
        names = self.model.names
        label = names[class_id] if isinstance(names, (list, tuple)) else names.get(class_id)
        require(label == expected_class_name, "YOLO flower class name mismatch")
        self.class_id, self.conf, self.device = class_id, conf, device
        self.model_sha = expected_sha256
        self.model_identity = expected_class_name

    def __call__(self, images):
        # Fixed-size inference microbatches prevent 512-image GPU OOM.
        output = []
        for offset in range(0, len(images), 32):
            subset = images[offset:offset + 32]
            predictions = self.model.predict(source=subset, imgsz=512,
                                             batch=len(subset),
                                             conf=self.conf, device=self.device,
                                             verbose=False, stream=False)
            require(len(predictions) == len(subset), "YOLO image counts changed")
            for prediction in predictions:
                best = None
                for box in prediction.boxes:
                    if int(box.cls.item()) != self.class_id:
                        continue
                    score = float(box.conf.item())
                    xyxy = list(map(float, box.xyxy[0].tolist()))
                    if best is None or score > best["score"]:
                        best = {"score": score, "xyxy": xyxy}
                output.append(best)
        return output

def packets(reader, n):
    from itertools import islice
    while True:
        part = list(islice(reader, n))
        if not part:
            return
        yield part

def fetch_decode_packet(rows, *, part, tier, download, workers,
                        source_role="full93m", io_executor=None):
    """Network and decode only. Pure source metadata enters; source images are ephemeral."""
    def one(row):
        base = {k: "" for k in COLUMNS_OUT}
        base.update(
            photo_id=row["photo_id"], observation_id=row["observation_uuid"],
            inat_taxon_id=row["inat_taxon_id"], observed_on=row["observed_on"],
            latitude=row["latitude"], longitude=row["longitude"],
            source_part=str(part), source_tier=str(tier),
            photo_license=row["photo_license"],
            source_arm=(SOURCE_ARM if source_role != "K1_calibration" else
                        "K1_signature_verified_calibration_only"),
            source_photo_inclusion_probability=("1" if source_role != "K1_calibration" else ""),
            audit_stratum=audit_stratum(row["latitude"], row["observed_on"]),
            screen_status="error",
        )
        url = source_url(row["photo_id"], row["image_extension"])
        if url is None:
            base["technical_reason"] = "unsupported_image_extension"
            return base, None
        try:
            data = download(url)
            image, brightness = decode_medium(data)
            base.update(image_pixel_sha256=sha_bytes(data),
                        downloaded_bytes=str(len(data)),
                        image_width=str(image.width), image_height=str(image.height),
                        mean_visible_brightness=str(brightness),
                        screen_status="uncertain", technical_reason="")
            return base, image
        except ProviderBackpressure:
            raise
        except Exception as exc:
            base["technical_reason"] = ("pixel_access_or_decode_error:" +
                                        type(exc).__name__ +
                                        ":" + str(exc)[:80].replace("\n", " "))
            return base, None
    began=time.monotonic()
    if io_executor is None:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            result=list(executor.map(one, rows))
    else:
        result=list(io_executor.map(one, rows))
    require(len(result)==len(rows),"Source photo rows lost during medium-pixel decode")
    return result, round(time.monotonic()-began,3)

def classify_fetched_packet(rows,pairs,detector,*,download_seconds,
                            performance=None,feature_mode="legacy_hsv"):
    """One immutable YOLO detector handles all fetched decoded photo candidates."""
    started=time.monotonic()
    out=[base for base,_ in pairs]
    valid=[(i,image) for i,(_,image) in enumerate(pairs) if image is not None]
    if valid and detector is not None:
        try:
            boxes=detector([image for _,image in valid])
            require(len(boxes)==len(valid),"Detector dropped input image")
            for (i,image),box in zip(valid,boxes):
                if box is None:
                    out[i]["screen_status"]="negative"
                    out[i]["technical_reason"]="model_screen_no_flower_candidate"
                else:
                    score,xyxy=float(box["score"]),box["xyxy"]
                    require(math.isfinite(score) and 0<=score<=1 and len(xyxy)==4,
                            "Invalid YOLO bounding box")
                    out[i]["screen_status"]="candidate"
                    out[i]["flower_box_score"]=format(score,".6g")
                    out[i]["flower_box_xyxy"]=",".join(format(float(v),".6g") for v in xyxy)
                    proxies=proxy_from_box(image,xyxy)
                    for key,value in zip(
                        ("box_proxy_white","box_proxy_yellow_orange",
                         "box_proxy_red_pink","box_proxy_blue_purple",
                         "box_proxy_other"),proxies):
                        out[i][key]=value
                    if feature_mode=="floral_box_lab_v2":
                        try:
                            out[i].update(flower_box_lab_descriptor(image,xyxy))
                        except Exception as lab_exc:
                            # Do NOT silently substitute a nine-palette or
                            # guessed petal value when image-Lab fails.
                            out[i]["screen_status"]="uncertain"
                            out[i]["technical_reason"]=(
                                "floral_box_Lab_proxy_error:"+
                                type(lab_exc).__name__)
        except Exception as exc:
            for i,_ in valid:
                out[i]["screen_status"]="error"
                out[i]["technical_reason"]="model_inference_error:"+type(exc).__name__
    elif valid:
        for i,_ in valid:
            out[i]["technical_reason"]="pixel_scanned_without_flower_detector"
    if performance is not None:
        performance.update(
            medium_image_access_and_decode_wall_seconds=round(download_seconds,3),
            flower_box_inference_and_proxy_wall_seconds=round(
                time.monotonic()-started,3))
    require(len(out)==len(rows) and all(row["screen_status"] in ALLOWED for row in out),
            "Stage1 photo row/status conservation failed")
    return out

def process_packet(rows, *, part, tier, detector, download, workers,
                   source_role="full93m", performance=None, io_executor=None,
                   preloaded_pairs=None, preloaded_transfer_seconds=None,
                   feature_mode="legacy_hsv"):
    """Backward-compatible source scan, optionally using prefetched 500px pixels."""
    if preloaded_pairs is None:
        require(preloaded_transfer_seconds is None,"Missing prefetched source photo pixels")
        pairs,download_seconds=fetch_decode_packet(
            rows,part=part,tier=tier,download=download,workers=workers,
            source_role=source_role,io_executor=io_executor)
    else:
        require(len(preloaded_pairs)==len(rows)
                and isinstance(preloaded_transfer_seconds,(float,int))
                and preloaded_transfer_seconds>=0,
                "Prefetched photo count or timing changed")
        pairs,download_seconds=preloaded_pairs,preloaded_transfer_seconds
    return classify_fetched_packet(
        rows,pairs,detector,download_seconds=download_seconds,
        performance=performance,feature_mode=feature_mode)

def write_packet_csv(path, rows, *, fieldnames=COLUMNS_OUT):
    tmp = path.with_suffix(path.suffix + ".partial")
    with tmp.open("wb") as raw:
        with gzip.GzipFile(filename="", fileobj=raw, mode="wb",
                           compresslevel=2, mtime=0) as zipped:
            with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as text:
                w = csv.DictWriter(text, fieldnames=fieldnames, lineterminator="\n")
                w.writeheader()
                w.writerows(rows)
    os.replace(tmp, path)

def scan_source(source, output, *, expected_source_sha256,
                expected_count, part, tier,
                detector=None, download=None, workers=8, batch_size=512,
                max_packets=None, source_role="full93m",
                parent_source_sha256=None, chunk_index=None, io_executor=None,
                prefetch_next_packet=False, feature_mode="legacy_hsv"):
    """Scan actual pixels; v2 adds direct floral BOX sRGB Lab without 4-hue gate."""
    source, output = Path(source), Path(output)
    require(sha_file(source) == expected_source_sha256,
            "Signed full93m original input photo gzip SHA changed")
    require(type(expected_count) is int and expected_count >= 0,
            "Signed expected source count is mandatory")
    require(0 <= part < 128 and 0 <= tier < 14,
            "Invalid original 93m disjoint source partition")
    require(1 <= batch_size <= 8192 and 1 <= workers <= 64,
            "Unbounded pixel packet or download concurrency")
    require(feature_mode in FEATURE_MODES,
            "An explicit supported immutable photo colour feature mode is required")
    record_schema=LAB_SCHEMA if feature_mode=="floral_box_lab_v2" else SCHEMA
    column_names=COLUMNS_LAB if feature_mode=="floral_box_lab_v2" else COLUMNS_OUT
    require(not prefetch_next_packet or
            (io_executor is not None and batch_size <= 256),
            "Prefetch must use one shared I/O executor and <=256 photos per packet")
    require(max_packets is None or max_packets > 0, "Invalid pilot limit")
    require(source_role in ("full93m", "full93m_chunk", "K1_calibration"),
            "A full source part, verified bounded source chunk, or K1 calibration is required")
    chunk_role = source_role == "full93m_chunk"
    require((chunk_role and type(chunk_index) is int and chunk_index >= 0
             and isinstance(parent_source_sha256,str)
             and len(parent_source_sha256)==64
             and all(ch in "0123456789abcdef" for ch in parent_source_sha256))
            or (not chunk_role and chunk_index is None
                and parent_source_sha256 is None),
            "A signed original parent SHA and chunk index must accompany every bounded chunk")
    require(source_role != "K1_calibration" or max_packets is not None,
            "K1 benchmark is never a full-source photo census")
    require(output.exists() is False or output.is_dir(), "Output is not a directory")
    output.mkdir(parents=True, exist_ok=True)
    if download is None:
        pace = Pacer(2.0)   # Limit is PER worker. Coordinate globally for S3 use.
        download = lambda url: download_medium(url, pace)
    stats = Counter()
    first = last = None
    completed_packets = decoded_count = payload_bytes = 0
    measured_packet_wall_seconds = 0.0
    download_decode_seconds = model_proxy_seconds = 0.0
    model_identity = (
        "pixel_only_pilot" if detector is None else
        "yolo11:" + detector.model_sha + ":" +
        detector.model_identity + ":" + str(detector.conf)
        if isinstance(detector, WarmFlowerBoxDetector) else
        "test_injected_detector:" + type(detector).__name__
    )
    if feature_mode=="floral_box_lab_v2":
        model_identity+="|"+BOX_LAB_PROTOCOL
    source_scan_started=time.monotonic()
    with (ThreadPoolExecutor(max_workers=1) if prefetch_next_packet
          else nullcontext(None)) as one_ahead, gzip.open(
              source, "rt", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        require(tuple(reader.fieldnames or ()) == COLUMNS,
                "Not an original photo-by-observation disjoint tier source")

        def pending_packets():
            items=iter(enumerate(packets(reader,batch_size)))
            if not prefetch_next_packet:
                for i,pkt in items:
                    yield i,pkt,None
                return
            def schedule(item):
                if item is None or (max_packets is not None
                                    and item[0]>=max_packets):
                    return None
                index,pkt=item
                dest=output/("packet-%07d.csv.gz"%index)
                proof=output/("packet-%07d.receipt.json"%index)
                if dest.exists() or proof.exists():
                    return None
                return one_ahead.submit(
                    fetch_decode_packet,pkt,part=part,tier=tier,
                    download=download,workers=workers,
                    source_role=source_role,io_executor=io_executor)
            current=next(items,None)
            pending=schedule(current)
            while current is not None:
                following=next(items,None)
                next_future=schedule(following)
                yield current[0],current[1],pending
                current,pending=following,next_future

        for idx, rows, fetch_future in pending_packets():
            if max_packets is not None and idx >= max_packets:
                break
            # Source is sorted by taxonomic source shard, not photo ID.
            ids = [r["photo_id"] for r in rows]
            require(len(set(ids)) == len(ids) and
                    all(pid.isdecimal() and r["tier"] == str(tier)
                        for pid, r in zip(ids, rows)),
                    "Invalid original photo-ID/tier within packet")
            dest = output / ("packet-%07d.csv.gz" % idx)
            receipt_path = output / ("packet-%07d.receipt.json" % idx)
            if dest.exists() or receipt_path.exists():
                require(dest.exists() and receipt_path.exists(),
                        "Incomplete previous packet; manual diagnosis required")
                proof = json.loads(receipt_path.read_text(encoding="utf-8"))
                require(proof.get("schema")==record_schema
                        and proof.get("source_gzip_sha256") == expected_source_sha256
                        and proof.get("output_sha256") == sha_file(dest)
                        and proof.get("photo_ids") == len(rows)
                        and proof.get("first_photo_id") == ids[0]
                        and proof.get("last_photo_id") == ids[-1]
                        and proof.get("source_part") == part
                        and proof.get("source_tier") == tier
                        and proof.get("detector_identity") == model_identity
                        and proof.get("source_role") == source_role
                        and (not chunk_role or
                             (proof.get("source_parent_sha256") == parent_source_sha256
                              and proof.get("source_chunk_index") == chunk_index)),
                        "Cannot resume a packet with different source identity")
                counts = proof["screen_counts"]
            else:
                started = time.monotonic()
                perf = {}
                fetched = fetch_future.result() if fetch_future is not None else None
                result = process_packet(
                    rows, part=part, tier=tier,
                    detector=detector, download=download, workers=workers,
                    source_role=source_role, performance=perf,
                    io_executor=io_executor,
                    preloaded_pairs=fetched[0] if fetched is not None else None,
                    preloaded_transfer_seconds=fetched[1] if fetched is not None else None,
                    feature_mode=feature_mode)
                write_packet_csv(dest, result,fieldnames=column_names)
                counts = dict(Counter(r["screen_status"] for r in result))
                proof = {
                    "schema": record_schema,
                    **({"CIELAB_box_source_protocol":BOX_LAB_PROTOCOL}
                       if feature_mode=="floral_box_lab_v2" else {}),
                    "source_gzip_sha256": expected_source_sha256,
                    "first_photo_id": ids[0], "last_photo_id": ids[-1],
                    "photo_ids": len(result), "source_part": part, "source_tier": tier,
                    "detector_identity": model_identity,
                    "source_role": source_role,
                    **({"source_parent_sha256": parent_source_sha256,
                        "source_chunk_index": chunk_index} if chunk_role else {}),
                    "screen_counts": counts,
                    "medium_image_pixels_attempted": len(result),
                    "medium_image_pixels_decoded": sum(bool(r["image_pixel_sha256"]) for r in result),
                    "medium_image_bytes_downloaded": sum(
                        int(r["downloaded_bytes"]) for r in result
                        if r["downloaded_bytes"].isdecimal()),
                    "flower_stage_ground_truth_verified": False,
                    "box_colour_proxy_not_verified_flower_petals": True,
                    "screen_negative_is_not_botanical_absence": True,
                    "original_full_resolution_ROI_modified": False,
                    "packet_wall_seconds": round(
                        time.monotonic() - started +
                        (fetched[1] if fetched is not None else 0), 3),
                    **perf,
                    "output_sha256": sha_file(dest),
                }
                temporary = receipt_path.with_suffix(".partial")
                temporary.write_text(json.dumps(proof, sort_keys=True, indent=2) + "\n")
                os.replace(temporary, receipt_path)
            stats.update(counts)
            decoded_count += int(proof["medium_image_pixels_decoded"])
            payload_bytes += int(proof["medium_image_bytes_downloaded"])
            measured_packet_wall_seconds += float(proof["packet_wall_seconds"])
            download_decode_seconds += float(
                proof.get("medium_image_access_and_decode_wall_seconds", 0))
            model_proxy_seconds += float(
                proof.get("flower_box_inference_and_proxy_wall_seconds", 0))
            completed_packets += 1
            first = ids[0] if first is None else first
            last = ids[-1]
    total = sum(stats.values())
    full = source_role in ("full93m", "full93m_chunk") and max_packets is None and total == expected_count
    if max_packets is None:
        require(total == expected_count,
                "Signed full source input photo count does not match true attempts")
    complete = {
        "schema": record_schema,
        **({"CIELAB_box_source_protocol":BOX_LAB_PROTOCOL,
            "floral_region_is_not_petals_semantically_separated":True,
            "principal_Lab_mode_is_not_calibrated_spectral_reflectance":True}
           if feature_mode=="floral_box_lab_v2" else {}),
        "source_gzip_sha256": expected_source_sha256,
        "source_part": part, "source_tier": tier,
        "source_role": source_role,
        "calibration_subset_not_global_census": source_role == "K1_calibration",
        **({"source_parent_sha256": parent_source_sha256,
            "source_chunk_index": chunk_index,
            "entire_original_source_part_completed": False,
            "all_source_photo_ids_scanned_only_in_this_bounded_chunk": True}
           if chunk_role else {}),
        "expected_photo_ids": expected_count,
        "actual_unique_source_photo_ids_attempted": total,
        "packets_sealed": completed_packets,
        "screen_counts": dict(stats),
        "detector_identity": model_identity,
        "medium_image_pixels_decoded": decoded_count,
        "medium_image_pixels_unavailable": total - decoded_count,
        "medium_image_payload_bytes_downloaded": payload_bytes,
        "aggregate_packet_work_seconds": round(measured_packet_wall_seconds, 3),
        "aggregate_medium_image_access_and_decode_wall_seconds":
            round(download_decode_seconds, 3),
        "aggregate_flower_box_inference_and_proxy_wall_seconds":
            round(model_proxy_seconds, 3),
        "images_per_worker_wall_second_if_sequential": (
            round(total / measured_packet_wall_seconds, 4)
            if measured_packet_wall_seconds > 0 else None),
        "all_original_input_photos_attempted": full,
        "screening_only_not_verified_flower_colour": True,
        "input_is_disjoint_upstream_by_photo_id": True,
    }
    if prefetch_next_packet:
        true_elapsed=time.monotonic()-source_scan_started
        complete["source_scan_prefetch_overlapped"]=True
        complete["actual_medium_pixel_source_wall_seconds"]=round(true_elapsed,3)
        complete["actual_photo_ids_per_second_after_prefetch"]=(
            round(total/true_elapsed,4) if true_elapsed>0 else None)
        complete["prefetched_photo_packet_size"]=batch_size
    if full:
        final_path = output / "complete_receipt.json"
        if final_path.exists():
            previous=json.loads(final_path.read_text(encoding="utf-8"))
            # Replays do not recompute original elapsed time, photo colours,
            # downloaded image bytes, or historical full ROI outcomes.
            if prefetch_next_packet:
                stable=[k for k in complete if k not in (
                    "source_scan_prefetch_overlapped",
                    "actual_medium_pixel_source_wall_seconds",
                    "actual_photo_ids_per_second_after_prefetch",
                    "prefetched_photo_packet_size")]
                require(all(previous.get(k)==complete[k] for k in stable),
                        "A completed photo source cannot be relabelled by resume")
                return previous
            require(previous == complete,
                    "Immutable complete source shard receipt changed")
        else:
            final_path.write_text(json.dumps(complete, sort_keys=True, indent=2) + "\n")
    return complete

def main():
    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument("--source", type=Path, required=True)
    a.add_argument("--source-sha256", required=True)
    a.add_argument("--expected-count", required=True, type=int)
    a.add_argument("--part", required=True, type=int)
    a.add_argument("--tier", required=True, type=int)
    a.add_argument("--output", required=True, type=Path)
    a.add_argument("--weights", type=Path)
    a.add_argument("--weights-sha256")
    a.add_argument("--flower-class-id", type=int, default=0)
    a.add_argument("--flower-class-name", default="item")
    a.add_argument("--device", default="cpu")
    a.add_argument("--model-conf", type=float, default=0.02)
    a.add_argument("--allow-pixel-only-pilot", action="store_true")
    a.add_argument("--workers", type=int, default=8)
    a.add_argument("--batch-size", type=int, default=512)
    a.add_argument("--max-packets", type=int)
    a.add_argument("--requests-per-second", type=float, default=2.0)
    x = a.parse_args()
    require((x.weights is not None and bool(x.weights_sha256)) or
            (x.weights is None and x.allow_pixel_only_pilot and
             x.max_packets is not None),
            "Full scan requires a SHA-pinned detector; pixel-only is pilot-only")
    model = (WarmFlowerBoxDetector(
        x.weights, x.weights_sha256, x.flower_class_id,
        x.flower_class_name, device=x.device, conf=x.model_conf)
        if x.weights else None)
    pacer = Pacer(x.requests_per_second)
    source_download = lambda url: download_medium(url, pacer)
    receipt = scan_source(
        x.source, x.output, expected_source_sha256=x.source_sha256,
        expected_count=x.expected_count, part=x.part, tier=x.tier,
        detector=model, download=source_download, workers=x.workers,
        batch_size=x.batch_size, max_packets=x.max_packets)
    print("MOIS_FULL93M_REAL_MEDIUM_PIXEL_SCAN",
          json.dumps(receipt, sort_keys=True))

if __name__ == "__main__":
    main()
