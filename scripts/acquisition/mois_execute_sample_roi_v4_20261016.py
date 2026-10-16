#!/usr/bin/env python3
"""Execute ONE source-sealed Mois SAMPLE photo batch through original FCP ROI-v4.

The sample is <=40 original observations per taxon. This module NEVER treats it
as FCP's all-93M source, NEVER reads species, coordinates, environment, prior
colour outcomes or the private metadata join into the image classifier.

Production requires byte-verified FCP wrapper, measurement source and frozen
ROI-v4 YOLO/EfficientSAM model; exactly one bounded provider image request
sequence, with explicit on-host global pacing and prior access authorization.
An entire batch seals atomically or remains unfinished. Synthetic tests only
use injected fake image loaders and evaluators; no real pixels are retrieved.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import os
import re
import sys
import tempfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from scripts.acquisition.mois_stratified_photo_selection_20261016 import sha256_file
from scripts.acquisition.mois_build_sample_roi_handoff_20261016 import (
    WORKER, ACQUISITION, N_PARTITIONS,
)

# FCP PR #200, source-bound image-wrapper blob as read on 16 October 2026.
FCP_BACKEND_GIT_BLOB = "45de671444c0fac4fb4611375d19ed7e5a1608f3"
VALID_ID = re.compile(r"^FCPR-[0-9A-F]{24}$")
BATCH_CAP = 64
STATUS = "COMPLETE_MOIS_SAMPLE_ROI_WORK_UNIT_NOT_FULL_FCP_CENSUS"


def git_blob(path: Path) -> str:
    size = path.stat().st_size
    h = hashlib.sha1()
    h.update(b"blob " + str(size).encode() + b"\x00")
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(2 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def backend_from_fcp(fcp_root: Path):
    script = fcp_root / "scripts/acquisition/run_fcp_fresh_image_partition.py"
    if git_blob(script) != FCP_BACKEND_GIT_BLOB:
        raise ValueError("FCP provider/ROI adapter source blob drift; audit before pixel opening")
    sys.path.insert(0, str(fcp_root))
    spec = importlib.util.spec_from_file_location("_mois_pinned_fcp_roi_adapter", script)
    if spec is None or spec.loader is None:
        raise ValueError("Cannot load byte-frozen original FCP photo processing adapter")
    backend = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(backend)
    return backend


def read_packet(root: Path, partition: int, batch: int,
                batch_size: int) -> tuple[dict, list[dict], list[dict]]:
    if not 0 <= partition < N_PARTITIONS or batch < 0 or not 1 <= batch_size <= BATCH_CAP:
        raise ValueError("Source sample processing must use partition 0-127 and 1-64 photo units")
    receipt_file = root / "mois_roi_packet_receipt.json"
    receipt = json.loads(receipt_file.read_text(encoding="utf-8"))
    if (receipt.get("schema") != "mois_sample_scoped_original_roi_v4_blind_packets_v1"
        or receipt.get("status") != "SAMPLE_ONLY_BLIND_METADATA_SEALED_PIXELS_NOT_MEASURED"
        or receipt.get("worker_partitions") != N_PARTITIONS
        or receipt.get("only_mois_sample_not_all_fcp_source_photos") is not True
        or receipt.get("metadata_site_and_climate_hidden_from_worker") is not True
        or receipt.get("original_pixels_measured") is not False
        or receipt.get("historical_colour_outcomes_reused") is not False):
        raise ValueError("Source is not a sealed Mois-only blind sample")
    counts = receipt.get("partition_photo_counts")
    if (not isinstance(counts, list) or len(counts) != N_PARTITIONS
        or any(not isinstance(v, int) or v < 0 for v in counts)
        or sum(counts) != receipt.get("selected_photo_ids")):
        raise ValueError("Sample photo partition count integrity failure")
    selected = []
    for stem, headers in (("worker", WORKER), ("acquisition", ACQUISITION)):
        name = f"{stem}_partition_{partition:03}.csv"
        filepath = root / name
        if sha256_file(filepath) != receipt["file_sha256"].get(name):
            raise ValueError("A Mois sample worker/acquisition source packet changed")
        with filepath.open("r", newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames != headers:
                raise ValueError("Image-worker packet columns drifted from original location-blind method")
            entries = list(reader)
            if len(entries) != counts[partition]:
                raise ValueError("Sample partition photos missing or duplicated")
            selected.append(entries)
    worker, acquisition = selected
    if (len({w["measurement_id"] for w in worker}) != len(worker)
        or len({a["measurement_id"] for a in acquisition}) != len(acquisition)):
        raise ValueError("Duplicated Mois sample measurement identity")
    for w, a in zip(worker, acquisition):
        ident = w["measurement_id"]
        if (not VALID_ID.fullmatch(ident)
            or w["image_filename"] != ident + ".jpg"
            or not w["photo_license"]
            or any(w[key] != a[key] for key in WORKER)
            or not a["photo_url_large"]):
            raise ValueError("Sample worker/acquisition identity mismatch")
    left = batch * batch_size
    if left >= len(worker):
        raise ValueError("Requested sample photo batch is outside selected partition")
    return receipt, worker[left:left + batch_size], acquisition[left:left + batch_size]


def execute(packet_root: Path, sealed_root: Path, partition: int, batch: int,
            batch_size: int, *, backend, original_measure: Path,
            original_detector: Path, efficient_sam: Path,
            synthetic: bool = False) -> dict:
    receipt, workers, acq = read_packet(packet_root, partition, batch, batch_size)
    if not synthetic and not os.environ.get("FCP_PROVIDER_GLOBAL_RATE_FILE"):
        raise ValueError("Coordinated, same-host FCP_PROVIDER_GLOBAL_RATE_FILE required for provider access")
    fingerprint = backend.verify_frozen_method(original_measure, original_detector, efficient_sam)
    if synthetic and not fingerprint.get("synthetic_test_only"):
        raise ValueError("Test adapter needs explicit synthetic-only method attestation")
    if not synthetic and (fingerprint.get("location_blind") is not True
                          or fingerprint.get("no_old_measurement_csv_inputs") is not True):
        raise ValueError("Original ROI-v4 method did not certify location blindness")
    stem = f"terminal_p{partition:03}_b{batch:07}"
    target = sealed_root / stem
    inflight = sealed_root / (stem + ".INCOMPLETE")
    if target.exists() or inflight.exists():
        raise ValueError("Existing sealed or interrupted sample ROI batch cannot be overwritten")
    sealed_root.mkdir(parents=True, exist_ok=True)
    inflight.mkdir()
    with tempfile.TemporaryDirectory(prefix="mois_sample_roi_") as scratch:
        work = Path(scratch)
        images = work / "images"
        images.mkdir()
        good = []
        missing = []
        for row in acq:
            try:
                pixels = backend.image_loader(row["photo_url_large"])
                if not pixels:
                    raise ValueError("Empty image bytes from original adapter")
                (images / row["image_filename"]).write_bytes(pixels)
                good.append({key: row[key] for key in WORKER})
            except backend.PermanentOriginalPhotoError as exc:
                missing.append({
                    "measurement_id": row["measurement_id"],
                    "morph": "mixed_uncertain",
                    "measurement_status": "image_acquisition_failed",
                    "failure_reasons": type(exc).__name__ + ":" + str(exc)[:160],
                    **{"flower_fraction_" + c: 0.0 for c in backend.BIO},
                })
        columns = ["measurement_id", "morph", "measurement_status",
                   *["flower_fraction_" + c for c in backend.BIO]]
        measured = pd.DataFrame(columns=columns)
        if good:
            blind = work / "blind_worker.csv"
            pd.DataFrame(good, columns=WORKER).to_csv(blind, index=False)
            measured = backend.evaluate_with_frozen_source(
                blind, images, original_measure, original_detector,
                efficient_sam, work / "original_blind_roi_results.csv")
        combined = pd.concat([measured, pd.DataFrame(missing)], ignore_index=True, sort=False)
        expected = [v["measurement_id"] for v in workers]
        backend.validate_terminal(combined, expected)
        combined = combined.sort_values("measurement_id", kind="stable")
        outcome = inflight / "original_roi_sample_terminal.csv"
        combined.to_csv(outcome, index=False)
        num_classified = int((combined["measurement_status"] == "classified_four_state_morph").sum())
        seal = {
            "schema": "mois_sample_scoped_original_roi_v4_work_unit_v1",
            "status": STATUS,
            "batch_photo_count": len(expected),
            "photo_ids_classified": num_classified,
            "photo_ids_acquisition_failed": len(missing),
            "unclassifiable_other": len(expected) - num_classified - len(missing),
            "partition": partition,
            "batch": batch,
            "compute_batch_size_not_biological_inclusion_cap": batch_size,
            "full_sample_packets_sha256": sha256_file(packet_root / "mois_roi_packet_receipt.json"),
            "source_packet_worker_sha256": receipt["file_sha256"][f"worker_partition_{partition:03}.csv"],
            "source_packet_acquisition_sha256": receipt["file_sha256"][f"acquisition_partition_{partition:03}.csv"],
            "exact_selected_measurement_ids_sha256": hashlib.sha256(
                "".join(v + "\n" for v in expected).encode()).hexdigest(),
            "sealed_terminal_csv_sha256": sha256_file(outcome),
            "method": fingerprint,
            "original_fcp_roi_wrapper_git_blob": (
                "SYNTHETIC_TEST_ONLY" if synthetic else FCP_BACKEND_GIT_BLOB),
            "fcp_full_census_completed_by_this_batch": False,
            "mois_sample_complete": False,
            "original_pixels_persisted": False,
            "historical_colour_outcomes_reused": False,
            "test_only_synthetic": bool(synthetic),
        }
        (inflight / "work_unit_receipt.json").write_text(
            json.dumps(seal, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    # Rename ONE entire directory, so partial/terminal receipt pairs are never
    # misrepresented as an independent completed batch on interrupted writes.
    inflight.replace(target)
    return seal


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sample-locked-packets", type=Path, required=True)
    p.add_argument("--sealed-output-root", type=Path, required=True)
    p.add_argument("--fcp-repository", type=Path, required=True)
    p.add_argument("--original-measure-script", type=Path, required=True)
    p.add_argument("--original-detector", type=Path, required=True)
    p.add_argument("--efficient-sam", type=Path, required=True)
    p.add_argument("--partition", type=int, required=True)
    p.add_argument("--batch", type=int, required=True)
    p.add_argument("--batch-size", type=int, default=BATCH_CAP)
    p.add_argument("--provider-access-reviewed", action="store_true")
    a = p.parse_args()
    if not a.provider_access_reviewed:
        p.error("Explicit review of source-image provider terms, rights and request budget required")
    backend = backend_from_fcp(a.fcp_repository)
    result = execute(a.sample_locked_packets, a.sealed_output_root, a.partition,
                     a.batch, a.batch_size, backend=backend,
                     original_measure=a.original_measure_script,
                     original_detector=a.original_detector, efficient_sam=a.efficient_sam)
    print(json.dumps({k: result[k] for k in (
        "status", "partition", "batch", "batch_photo_count", "photo_ids_classified",
        "photo_ids_acquisition_failed", "unclassifiable_other"
    )}, sort_keys=True))


if __name__ == "__main__":
    main()
