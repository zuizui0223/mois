#!/usr/bin/env python3
"""Bounded, actual medium-pixel throughput pilot on immutable MOIS K1 records.

Source outcome-stratified K1 selection is a CALIBRATION DATASET ONLY:
128 previously classifiable and 128 ROI/flip-gate unknown photos.
It is *not* an unbiased world-flora source or 93m-stage completion.
No provider calls until --execute. No flower model in this first transfer
benchmark; image-pixel decodes and wall times are genuinely measured.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
from pathlib import Path

from scripts.acquisition.mois_full93m_medium_pixel_scan import (
    scan_source, sha_file, download_medium, Pacer, WarmFlowerBoxDetector,
)
from scripts.acquisition.mois_onepass_full93m_target_tiers import COLUMNS

K1_SHA = "0d3de3374a06b1e7b9ad0281e2ab7844632ed139db871488524061191c8f34c1"
RECEIPT_SHA = "e86463bc61f9eebd8cbdbe01269f767c140bf243f9e374027008e9f095eb9044"
N_K1 = 76753
STATUS_CLASSIFIED = "classified_four_state_morph"
STATUS_GATE = "not_evaluable_roi_or_flip_gate"
SEED = "mois-93m-stage1-bounded-K1-real-pixels-throughput-v1"
N_PER = 128

def require(ok, why):
    if not ok:
        raise ValueError("MOIS_K1_MEDIUM_PILOT: " + why)

def rank(photo_id):
    return hashlib.sha256((SEED + "|" + str(photo_id)).encode()).digest()

def prepare(k1_parquet, k1_receipt, directory, *, n_per=N_PER,
            source_sha=K1_SHA, receipt_sha=RECEIPT_SHA,
            expected_count=N_K1):
    import pandas as pd
    directory = Path(directory)
    require(not directory.exists(), "Pilot source must not overwrite previous data")
    require(sha_file(k1_parquet) == source_sha and
            sha_file(k1_receipt) == receipt_sha,
            "Source K1 photo bytes do not match official SHA receipts")
    proof = json.loads(Path(k1_receipt).read_text(encoding="utf-8"))
    require(proof.get("schema") ==
            "mois_global93m_observation_baseline_real_RGB_photo_results_v1"
            and proof.get("attempted") == expected_count
            and proof.get("measured_photo_parquet_sha256") == source_sha,
            "K1 source terminal provenance changed")
    fields = ("photo_id", "observation_id", "inat_taxon_id", "observed_on",
              "latitude", "longitude", "photo_license", "photo_url_large",
              "measurement_status")
    records = pd.read_parquet(k1_parquet, columns=list(fields))
    require(len(records) == expected_count and records.photo_id.astype(str).is_unique,
            "Official K1 full photo-ID source count changed")
    selected = []
    for status in (STATUS_CLASSIFIED, STATUS_GATE):
        subset = records[records.measurement_status.eq(status)].copy()
        require(len(subset) >= n_per, "Not enough signed K1 original photo outcomes")
        subset["_rank"] = subset.photo_id.astype(str).map(rank)
        subset = subset.sort_values("_rank", kind="stable").head(n_per)
        selected.extend(subset.to_dict("records"))
    ids = [str(r["photo_id"]) for r in selected]
    obs = [str(r["observation_id"]) for r in selected]
    require(len(set(ids)) == len(ids) == n_per * 2 and len(set(obs)) == n_per * 2,
            "Benchmark must use distinct original photo AND observation IDs")
    selected.sort(key=lambda r: int(r["photo_id"]))
    directory.mkdir(parents=True)
    csv_path = directory / "K1_outcome_stratified_medium_256_original_photo_source.csv.gz"
    lookup = directory / "K1_signed_original_roi_status_calibration.csv"
    with (gzip.open(csv_path, "wt", encoding="utf-8", newline="") as csv_handle,
          lookup.open("w", encoding="utf-8", newline="") as lookup_handle):
        writer = csv.DictWriter(csv_handle, fieldnames=COLUMNS, lineterminator="\n")
        writer.writeheader()
        original = csv.DictWriter(lookup_handle,
                                  fieldnames=("photo_id", "K1_original_measurement_status"),
                                  lineterminator="\n")
        original.writeheader()
        for r in selected:
            pid = str(r["photo_id"])
            url = str(r["photo_url_large"])
            require(url.startswith(
                "https://inaturalist-open-data.s3.amazonaws.com/photos/" + pid +
                "/large."), "K1 source photo URL is not its original numeric photo ID")
            ext = url.rsplit(".", 1)[-1].lower()
            require(ext in ("jpg", "jpeg", "png"), "Unsupported K1 original source extension")
            row = {k: "" for k in COLUMNS}
            row.update(photo_id=pid,
                       observation_uuid=str(r["observation_id"]),
                       inat_taxon_id=str(r["inat_taxon_id"]),
                       observed_on=str(r["observed_on"]),
                       latitude=str(r["latitude"]), longitude=str(r["longitude"]),
                       photo_license=str(r["photo_license"]),
                       image_extension=ext, frozen_loader_supported="true",
                       tier="0")
            writer.writerow(row)
            original.writerow({
                "photo_id": pid,
                "K1_original_measurement_status": str(r["measurement_status"]),
            })
    report = {
        "schema": "mois_actual_K1_256_medium_pixel_calibration_source_v1",
        "is_calibration_not_93m_census": True,
        "K1_signed_source_sha256": source_sha,
        "K1_receipt_sha256": receipt_sha,
        "source_exact_K1_attempts": expected_count,
        "selected_original_photo_ids": len(selected),
        "signed_K1_four_hue_success_calibration": n_per,
        "signed_K1_ROI_flip_gate_unknown_calibration": n_per,
        "original_photo_source_sha256": sha_file(csv_path),
        "original_photo_K1_status_SHA256": sha_file(lookup),
        "new_photos_fetched_to_make_source": 0,
    }
    (directory / "calibration_source_receipt.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report

def benchmark(source_dir, result_dir, *, workers=8, rps=2.0,
              detector=None):
    require(workers <= 12 and 0 < rps <= 3.0,
            "Bounded provider access pilot only; no 93m bulk requests")
    source_dir, result_dir = Path(source_dir), Path(result_dir)
    meta = json.loads((source_dir / "calibration_source_receipt.json").read_text())
    require(meta["schema"] == "mois_actual_K1_256_medium_pixel_calibration_source_v1"
            and meta["is_calibration_not_93m_census"] is True
            and meta["selected_original_photo_ids"] == 256,
            "Unexpected original signed K1 benchmark source")
    source = source_dir / "K1_outcome_stratified_medium_256_original_photo_source.csv.gz"
    require(sha_file(source) == meta["original_photo_source_sha256"],
            "K1 benchmark source changed")
    pace = Pacer(rps)
    status = scan_source(
        source, result_dir, expected_source_sha256=meta["original_photo_source_sha256"],
        expected_count=256, part=0, tier=0, source_role="K1_calibration",
        detector=detector, download=lambda u: download_medium(u, pace),
        workers=workers, batch_size=256, max_packets=1)
    require(not status["all_original_input_photos_attempted"]
            and status["actual_unique_source_photo_ids_attempted"] == 256
            and status["calibration_subset_not_global_census"] is True
            and not (result_dir / "complete_receipt.json").exists(),
            "Calibration pilot accidentally claimed all-source completion")
    origin = {}
    with (source_dir / "K1_signed_original_roi_status_calibration.csv").open() as f:
        origin = {row["photo_id"]: row["K1_original_measurement_status"]
                  for row in csv.DictReader(f)}
    by_status = defaultdict(Counter)
    raw_bytes = []
    with gzip.open(result_dir / "packet-0000000.csv.gz", "rt") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        label = origin[row["photo_id"]]
        by_status[label][row["screen_status"]] += 1
        if row["downloaded_bytes"].isdecimal():
            raw_bytes.append(int(row["downloaded_bytes"]))
    require(len(rows) == 256 and len(origin) == 256,
            "Lost K1 source photo status during benchmarking")
    summary = {
        "schema": "mois_K1_256_actually_decoded_medium_pixel_benchmark_v1",
        "original_K1_photo_ROI_source_SHA256": meta["K1_signed_source_sha256"],
        "pilot_is_not_independent_photo_or_global_flora_prevalence": True,
        "pilot_images_attempted": 256,
        "medium_decoded": status["medium_image_pixels_decoded"],
        "medium_access_or_decode_unknown": status["medium_image_pixels_unavailable"],
        "medium_transfer_payload_bytes": status["medium_image_payload_bytes_downloaded"],
        "median_successful_medium_image_bytes": (
            sorted(raw_bytes)[len(raw_bytes) // 2] if raw_bytes else None),
        "stage1_wall_seconds": status["aggregate_packet_work_seconds"],
        "pixel_attempts_per_second": status["images_per_worker_wall_second_if_sequential"],
        "source_K1_roi_classification_by_medium_scan_status": {
            k: dict(v) for k, v in by_status.items()
        },
        "actual_flower_detection_model_used": detector is not None,
        "flower_box_model_original_weight_sha256": (
            detector.model_sha if detector is not None else None),
        "flower_box_candidate_count": int(status["screen_counts"].get("candidate", 0)),
        "flower_box_negative_screen_count": int(status["screen_counts"].get("negative", 0)),
        "flower_detection_recall_established": False,
        "no_raw_photo_bytes_published": True,
        "stage1_pixel_packet_SHA256": sha_file(result_dir / "packet-0000000.csv.gz"),
    }
    (result_dir / "K1_medium_pixel_benchmark_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary

def main():
    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument("--k1-parquet", type=Path, required=True)
    a.add_argument("--k1-receipt", type=Path, required=True)
    a.add_argument("--source-output", type=Path, required=True)
    a.add_argument("--pixel-output", type=Path, required=True)
    a.add_argument("--weights", type=Path)
    a.add_argument("--weights-sha256")
    a.add_argument("--flower-class-id", type=int, default=0)
    a.add_argument("--flower-class-name", default="item")
    a.add_argument("--device", default="cpu")
    a.add_argument("--detector-conf", type=float, default=0.02)
    a.add_argument("--execute", action="store_true",
                   help="Actually download no more than 256 medium images at <=3 req/s")
    a.add_argument("--workers", type=int, default=8)
    a.add_argument("--requests-per-second", type=float, default=2.0)
    x = a.parse_args()
    source = prepare(x.k1_parquet, x.k1_receipt, x.source_output)
    print("MOIS_K1_SIGNED_PIXEL_CALIBRATION_SOURCE", json.dumps(source, sort_keys=True))
    if x.execute:
        require((x.weights is None and x.weights_sha256 is None) or
                (x.weights is not None and x.weights_sha256 is not None),
                "Detector weights and SHA must be specified together")
        detector = (WarmFlowerBoxDetector(
            x.weights, x.weights_sha256,
            x.flower_class_id, x.flower_class_name,
            device=x.device, conf=x.detector_conf) if x.weights else None)
        result = benchmark(x.source_output, x.pixel_output,
                           workers=x.workers, rps=x.requests_per_second,
                           detector=detector)
        print("MOIS_K1_256_REAL_MEDIUM_PIXELS", json.dumps(result, sort_keys=True))

if __name__ == "__main__":
    main()
