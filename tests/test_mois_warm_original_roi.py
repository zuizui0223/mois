"""Synthetic, no-network tests for isolated opt-in warm FCP ROI evaluator."""
from __future__ import annotations
import csv
import hashlib
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import pandas as pd
from scripts.acquisition import mois_warm_original_roi as warm


FAKE_SOURCE = '''
import csv
from pathlib import Path
DEFAULT_MEASUREMENT = Path("measurement.json")
DEFAULT_EXECUTION = Path("execution.json")
DEFAULT_ROI = Path("roi.json")
DEFAULT_LOCKED = Path("locked.json")
CREATED = 0
EXPECTED_VERSION = 1

class FrozenFlowerColourEstimator:
    def __init__(self, weight, sam, contract, *, torch_threads):
        global CREATED
        CREATED += 1
        self.scale = int(contract["version"])

    def measure(self, image):
        return image * self.scale

def parse_args():
    raise ValueError("Warm wrapper must provide location-blind args")

def main():
    args = parse_args()
    estimator = FrozenFlowerColourEstimator(
        args.detector_weight, args.efficient_sam_dir,
        {"version": EXPECTED_VERSION}, torch_threads=2,
    )
    with args.worker_manifest.open(newline="") as stream:
        original = list(csv.DictReader(stream))
    with args.output_csv.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["measurement_id","morph","measurement_status"])
        writer.writeheader()
        for row in original:
            n = estimator.measure(int(row["value"]))
            writer.writerow(dict(measurement_id=row["measurement_id"],
                                 morph="white" if n % 2 == 0 else "mixed_uncertain",
                                 measurement_status="classified_four_state_morph" if n % 2 == 0
                                                    else "not_evaluable_roi_or_flip_gate"))
    return 0
'''


class WarmOriginalROITests(unittest.TestCase):
    def setUp(self):
        self.tmp = TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.source = self.root / "original.py"
        self.source.write_text(FAKE_SOURCE)
        self.weight = self.root / "model.pt"
        self.weight.write_bytes(b"fake-test-only")
        self.sam = self.root / "sam"
        self.sam.mkdir()
        self.images = self.root / "images"
        self.images.mkdir()
        self.workers = [self.root / "blind0.csv", self.root / "blind1.csv"]
        for idx, worker in enumerate(self.workers):
            with worker.open("w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=["measurement_id","value"])
                w.writeheader()
                w.writerows([{"measurement_id":f"M{idx}{v}","value":str(v)} for v in (1,2,3)])
        self.engine = warm.WarmFrozenROI(expected_source_blob=warm.git_blob(self.source))

    def evaluate(self, idx):
        return self.engine.evaluate(
            self.workers[idx], self.images, self.source,
            self.weight, self.sam, self.root / f"output{idx}.csv",
        )

    def test_same_original_main_both_photo_work_units_single_yolo_onnx_init(self):
        frames = [self.evaluate(0), self.evaluate(1)]
        self.assertEqual(self.engine.estimator_initializations, 1)
        self.assertEqual(self.engine.module.CREATED, 1)
        self.assertEqual(self.engine.calls, 2)
        self.assertEqual(list(frames[0].morph), ["mixed_uncertain","white","mixed_uncertain"])
        self.assertEqual(list(frames[1].morph), ["mixed_uncertain","white","mixed_uncertain"])
        self.assertEqual(self.engine.timing_receipt()["frozen_yolo_and_onnx_initializations"], 1)

    def test_source_mutated_after_first_packet_fails_closed(self):
        self.evaluate(0)
        self.source.write_text(self.source.read_text() + "\n# modified\n")
        with self.assertRaisesRegex(ValueError, "source Git blob"):
            self.evaluate(1)
        self.assertEqual(self.engine.calls, 1)

    def test_different_model_weight_and_roi_contract_forbidden(self):
        self.evaluate(0)
        alternative = self.root / "other.pt"
        alternative.write_bytes(b"test")
        with self.assertRaisesRegex(ValueError, "change frozen model"):
            self.engine.evaluate(
                self.workers[1], self.images, self.source,
                alternative, self.sam, self.root / "changed.csv",
            )
        self.engine.module.EXPECTED_VERSION = 2
        with self.assertRaisesRegex(ValueError, "change frozen model"):
            self.evaluate(1)

    def test_wrong_pinned_source_refused_before_any_output(self):
        original = warm.WarmFrozenROI()
        with self.assertRaisesRegex(ValueError, "Frozen FCP source Git blob changed"):
            original.evaluate(
                self.workers[0], self.images, self.source,
                self.weight, self.sam, self.root / "none.csv",
            )
        self.assertFalse((self.root / "none.csv").exists())

    def test_repeat_original_main_results_have_matching_exact_rows(self):
        rows = [self.evaluate(0), self.evaluate(1)]
        # A fresh model creates the same result for each original blind packet.
        for i, worker in enumerate(self.workers):
            import importlib.util
            spec = importlib.util.spec_from_file_location(f"cold_fcp_{i}", self.source)
            cold = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(cold)
            from types import SimpleNamespace
            out = self.root / f"cold{i}.csv"
            cold.parse_args = lambda output=out, input=worker: SimpleNamespace(
                worker_manifest=input,images_dir=self.images,detector_weight=self.weight,
                efficient_sam_dir=self.sam,output_csv=output)
            self.assertEqual(cold.main(), 0)
            pd.testing.assert_frame_equal(
                rows[i], pd.read_csv(out, dtype={"measurement_id": str}),
                check_exact=True,
            )


if __name__ == "__main__":
    unittest.main()
