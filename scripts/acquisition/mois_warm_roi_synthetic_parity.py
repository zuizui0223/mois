#!/usr/bin/env python3
"""Pair exact frozen CPU FCP main() cold vs warm on generated non-ecological pixels.

No iNaturalist downloads; no synthetic colour labels enter MOIS. Runs 2 x 4
blinded images, using exactly the original detector/ONNX/palette source in both
arms. Proves result parity and measures only model/evaluator process overhead.
A disjoint real-photo holdout is still required before production promotion.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time

import pandas as pd
from PIL import Image, ImageDraw

from scripts.acquisition.mois_warm_original_roi import WarmFrozenROI, git_blob, FROZEN_SOURCE_BLOB


def make_images(path: Path, n: int = 8):
    path.mkdir(parents=True, exist_ok=False)
    rows = []
    colours = [
        (245, 243, 235), (229, 192, 43), (215, 83, 107),
        (102, 96, 210), (250, 244, 230), (242, 179, 56),
        (180, 46, 112), (116, 91, 179),
    ]
    for i in range(n):
        im = Image.new("RGB", (320, 320), (30 + i, 68 + i * 3, 33 + i))
        d = ImageDraw.Draw(im)
        cx, cy = 150 + 10 * (i % 3), 148 + 8 * (i % 2)
        for j in range(9):
            angle = (2 * math.pi * j / 9) + i * 0.037
            px = cx + int(49 * math.cos(angle))
            py = cy + int(49 * math.sin(angle))
            d.ellipse((px - 36, py - 18, px + 36, py + 18), fill=colours[i])
        d.ellipse((cx - 26, cy - 26, cx + 26, cy + 26), fill=(135, 90, 45))
        filename = f"synthetic_{i:02d}.jpg"
        im.save(path / filename, format="JPEG", quality=95)
        rows.append({
            "measurement_id":f"SYNTH{str(i).zfill(3)}",
            "image_filename":filename,
            "photo_license":"synthetic_no_original_photo",
        })
    return rows


def benchmark(runtime: Path, sam: Path, output: Path):
    source = runtime / "scripts/analysis/measure_random_photo_first_blind_partition.py"
    detector = runtime / "data/atlas/qualification/roi_v4_training/jrc_yolo11n_last_v4.pt"
    if git_blob(source) != FROZEN_SOURCE_BLOB:
        raise ValueError("Pinned original ROI-v4 source mismatch")
    if not detector.is_file() or not sam.is_dir() or output.exists():
        raise ValueError("Frozen model not provisioned or benchmark output already exists")
    with tempfile.TemporaryDirectory(prefix="mois_real_fcp_warm_synthetic_") as tmp:
        root = Path(tmp)
        images = root / "images"
        rows = make_images(images)
        manifests = []
        for part in range(2):
            file = root / f"blind_{part}.csv"
            pd.DataFrame(rows[part*4:(part+1)*4],
                         columns=["measurement_id","image_filename","photo_license"]
            ).to_csv(file, index=False, lineterminator="\n")
            manifests.append(file)
        original_outputs, warm_outputs = [], []
        baseline_seconds = 0.0
        warm_seconds = 0.0
        image_sha = {
            f.name: hashlib.sha256(f.read_bytes()).hexdigest()
            for f in images.glob("*.jpg")
        }
        env = os.environ.copy()
        env["PYTHONPATH"] = os.pathsep.join(x for x in (
            str(runtime), str(Path(__file__).resolve().parents[2]),
            env.get("PYTHONPATH", "")
        ) if x)

        for part, worker in enumerate(manifests):
            out = root / f"cold_{part}.csv"
            start = time.perf_counter()
            subprocess.run([
                sys.executable, str(source),
                "--worker-manifest", str(worker),
                "--images-dir", str(images),
                "--detector-weight", str(detector),
                "--efficient-sam-dir", str(sam),
                "--output-csv", str(out),
            ], env=env, check=True)
            baseline_seconds += time.perf_counter() - start
            original_outputs.append(pd.read_csv(out, dtype={"measurement_id":str}))
        engine = WarmFrozenROI()
        for part, worker in enumerate(manifests):
            out = root / f"warm_{part}.csv"
            start = time.perf_counter()
            result = engine.evaluate(worker, images, source, detector, sam, out)
            warm_seconds += time.perf_counter() - start
            warm_outputs.append(result)
        mismatches = []
        for i, (cold, newer) in enumerate(zip(original_outputs, warm_outputs)):
            try:
                pd.testing.assert_frame_equal(cold, newer,
                    check_exact=True, check_dtype=True, check_like=False)
            except AssertionError as exc:
                mismatches.append({"packet":i, "reason":str(exc)[:500]})
        if engine.estimator_initializations != 1 or engine.calls != 2:
            raise RuntimeError("Warm evaluator initialized more than once or lost a packet")
        report = {
            "schema":"mois_same_frozen_roi_v4_synthetic_warm_cold_parity_v1",
            "date_utc":datetime.now(timezone.utc).isoformat(),
            "synthetic_non_ecological_images":8,
            "frozen_original_FCP_git_blob":FROZEN_SOURCE_BLOB,
            "two_blind_packets":True,
            "raw_iNaturalist_images_downloaded":0,
            "all_image_byte_SHA256":image_sha,
            "original_cold_model_initializations":2,
            "warm_model_initializations":engine.estimator_initializations,
            "cold_seconds":baseline_seconds,
            "warm_seconds":warm_seconds,
            "observed_synthetic_runtime_ratio":baseline_seconds/warm_seconds,
            "all_original_terminal_columns_exactly_equal":not mismatches,
            "mismatched_packet_columns":mismatches,
            "original_source_or_thresholds_modified":False,
            "this_is_not_a_real_photo_throughput_or_class_recall_estimate":True,
            "no_generated_pixels_or_hues_admitted_to_mois_data":True,
        }
        output.mkdir(parents=True, exist_ok=False)
        (output / "warm_frozen_synthetic_parity_receipt.json").write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n")
        print("MOIS_WARM_FROZEN_SYNTHETIC_PARITY",json.dumps(report, sort_keys=True),flush=True)
        if mismatches:
            raise AssertionError("Original cold/warm pixel outputs changed")
        return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("runtime", "sam", "output"):
        p.add_argument("--" + key, type=Path, required=True)
    a = p.parse_args()
    benchmark(a.runtime, a.sam, a.output)


if __name__ == "__main__":
    main()
