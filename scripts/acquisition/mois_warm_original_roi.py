#!/usr/bin/env python3
"""Opt-in warm frozen ROI-v4 evaluator for prospective original-photo workers.

Exactly the SHA-pinned FCP measurement main() runs on each blind 256-photo
packet, but the immutable YOLO11 + two EfficientSAM ONNX sessions can be
instantiated once across up to four photo packets in the SAME worker process.

This module neither changes deployed p0016/p0032/geo05 jobs nor downloads
photos as a side effect of import. A paired actual-photo parity/throughput
benchmark is required before using the opt-in work command in production.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
import time

FROZEN_SOURCE_BLOB = "b1d7091db44682559bbe4b48a9f3b8db648501a6"


def git_blob(path: Path) -> str:
    raw = Path(path).read_bytes()
    return hashlib.sha1(
        b"blob " + str(len(raw)).encode("ascii") + b"\0" + raw
    ).hexdigest()


class WarmFrozenROI:
    """Call unmodified source.main() repeatedly with one verified estimator."""

    def __init__(self, expected_source_blob: str = FROZEN_SOURCE_BLOB):
        # Override the original blob only in a controlled synthetic unit test.
        if len(expected_source_blob) != 40 or any(
            c not in "0123456789abcdef" for c in expected_source_blob
        ):
            raise ValueError("Invalid expected frozen source blob")
        self.expected_source_blob = expected_source_blob
        self.module = None
        self.source_path = None
        self.estimator = None
        self.estimator_key = None
        self.estimator_initializations = 0
        self.calls = 0
        self.total_evaluator_seconds = 0.0

    def _load(self, source: Path) -> None:
        source = Path(source).resolve(strict=True)
        if git_blob(source) != self.expected_source_blob:
            raise ValueError("Frozen FCP source Git blob changed")
        if self.module is not None:
            if source != self.source_path:
                raise ValueError("Source file changed inside warm ROI worker")
            return
        spec = importlib.util.spec_from_file_location("_mois_pinned_frozen_roi_main", source)
        if spec is None or spec.loader is None:
            raise RuntimeError("Cannot load the pinned original FCP evaluator")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        original_estimator = module.FrozenFlowerColourEstimator

        def estimator_factory(weight, sam, contract, *, torch_threads):
            # Explicit runtime freeze: never reuse sessions across new weights,
            # different ONNX paths, contracts or thread configuration.
            key = (
                str(Path(weight).resolve(strict=True)),
                str(Path(sam).resolve(strict=True)),
                json.dumps(contract, sort_keys=True, separators=(",", ":")),
                int(torch_threads),
            )
            if self.estimator is None:
                self.estimator = original_estimator(
                    weight, sam, contract, torch_threads=torch_threads
                )
                self.estimator_key = key
                self.estimator_initializations += 1
            elif key != self.estimator_key:
                raise ValueError("A warm ROI worker attempted to change frozen model settings")
            return self.estimator

        # Do NOT alter the saved FCP source, its main loop or ROI methodology.
        # Only replace the constructor *in this isolated imported module*.
        module.FrozenFlowerColourEstimator = estimator_factory
        self.module, self.source_path = module, source

    def evaluate(self, worker, images, source, detector, sam, out):
        """Drop-in run_batch(evaluator=) seam; input remains location-blind."""
        import pandas as pd

        self._load(source)
        # In addition to the run_batch() attestation on every photo work unit,
        # avoid continuing if the input code changed within this process.
        if git_blob(self.source_path) != self.expected_source_blob:
            raise ValueError("Frozen source edited while warm worker was running")
        args = SimpleNamespace(
            worker_manifest=Path(worker),
            images_dir=Path(images),
            detector_weight=Path(detector),
            efficient_sam_dir=Path(sam),
            output_csv=Path(out),
            measurement_contract=self.module.DEFAULT_MEASUREMENT,
            execution_contract=self.module.DEFAULT_EXECUTION,
            roi_contract=self.module.DEFAULT_ROI,
            locked_result=self.module.DEFAULT_LOCKED,
        )
        start = time.perf_counter()
        self.module.parse_args = lambda: args
        if self.module.main() != 0:
            raise RuntimeError("Unmodified frozen FCP main() failed")
        elapsed = time.perf_counter() - start
        self.calls += 1
        self.total_evaluator_seconds += elapsed
        return pd.read_csv(out, dtype={"measurement_id": str})

    def timing_receipt(self):
        return {
            "schema": "mois_prospective_frozen_roi_warm_evaluator_v1",
            "frozen_source_blob": self.expected_source_blob,
            "separate_original_main_calls": self.calls,
            "frozen_yolo_and_onnx_initializations": self.estimator_initializations,
            "total_frozen_evaluator_seconds": self.total_evaluator_seconds,
            "method_or_thresholds_modified": False,
            "old_photo_RGB_workers_modified": False,
            "warm_engine_production_parity_established": False,
        }


def prospective_work(source, runtime, sam, group, output):
    """Opt-in prospective photo group; **never** called by existing workflows."""
    from scripts.acquisition.mois_full_index_real_original_photo_measure import (
        plan, input_rows, measure_group,
    )
    from scripts.acquisition.mois_prefetch import PrefetchLoader
    from scripts.acquisition.run_fcp_fresh_image_partition import (
        image_loader, run_batch, PermanentOriginalPhotoError,
    )
    from scripts.analysis.mois_reassemble_original_new4096_measurements import METHOD

    info = plan(source)
    units = range(group * 4, min((group + 1) * 4, info["work_units"]))
    if not units:
        raise ValueError("Unknown prospective source group")
    urls = []
    for unit in units:
        receipt = Path(output) / f"shard_{unit:03d}" / "receipt.json"
        if not receipt.is_file():
            rows, _ = input_rows(source, unit)
            urls.extend(rows.photo_url_large.tolist())

    warm = WarmFrozenROI()
    def runner(*args, **kwargs):
        result = run_batch(*args, evaluator=warm.evaluate, **kwargs)
        if result["method"] != METHOD:
            raise ValueError("Frozen FCP photo-ROI method changed")
        return result

    start = time.perf_counter()
    with PrefetchLoader(
        urls, image_loader, Path(output) / "prefetch",
        permanent_errors=(PermanentOriginalPhotoError,),
    ) as prefetch:
        output_receipts = measure_group(
            source, group, runtime, sam, output,
            runner=runner, loader=prefetch,
        )
        stats = dict(prefetch.stats)
    result = warm.timing_receipt()
    result.update(
        group=group,
        photo_work_units=len(units),
        new_photos_downloaded=stats["download_calls"],
        photo_prefetch=stats,
        group_wall_seconds=time.perf_counter() - start,
        reused_completed_photo_units=len(units) - warm.calls,
        artificial_flower_labels_created=False,
    )
    Path(output).mkdir(parents=True, exist_ok=True)
    dest = Path(output) / f"warm_roi_group_{group:03d}_diagnostic.json"
    if dest.exists():
        raise ValueError("Refuse overwrite of previous warm diagnostic")
    dest.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print("MOIS_PROSPECTIVE_WARM_FROZEN_ROI", json.dumps(result, sort_keys=True), flush=True)
    return output_receipts


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("source", "runtime", "sam", "output"):
        p.add_argument("--" + key, type=Path, required=True)
    p.add_argument("--group", type=int, required=True)
    a = p.parse_args()
    prospective_work(a.source, a.runtime, a.sam, a.group, a.output)


if __name__ == "__main__":
    main()
