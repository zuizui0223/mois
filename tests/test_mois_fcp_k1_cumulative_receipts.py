"""Original finished K1 ROI receipt aggregation provenance and unknown safeguards."""
from __future__ import annotations

import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from scripts.acquisition import mois_fcp_k1_cumulative_receipts as k1


class TestCompletedK1CumulativeReceipts(unittest.TestCase):
    def setUp(self):
        self.tmp = TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.out = self.root / "published"

    def fixture(self, group, *, classified=400, fail=3):
        n = k1.group_expected(group)
        other = n - classified - fail
        self.assertGreaterEqual(other, 0)
        units = list(range(group * k1.GROUP_UNITS,
                           min((group + 1) * k1.GROUP_UNITS,
                               (k1.NEW_TARGET + k1.UNIT_SIZE - 1) // k1.UNIT_SIZE)))
        d = dict(
            schema=k1.SOURCE_SCHEMA,
            frozen_original_source_run=k1.SOURCE_RUN,
            full_93m_photo_target_catalogue=k1.SOURCE_PHOTOS,
            original_source_manifest_SHA256=k1.SOURCE_SHA,
            original_precolour_probability_source_SHA256=k1.PARENT_SHA,
            source_arm=k1.SOURCE_ARM,
            previously_measured_old_original_photo_IDs_reused_outside_this_group=k1.REUSED_OLD_ORIGINALS,
            original_new_selected_photo_targets_all_groups=k1.NEW_TARGET,
            original_worker_group=group,
            original_completed_source_units=units,
            original_photo_terminal_source_receipts_SHA256=[
                dict(unit=u, receipt_sha256="a" * 64,
                     photos_sha256="b" * 64, terminal_sha256="c" * 64)
                for u in units],
            actual_new_original_photo_attempts=n,
            actual_automatic_four_hue_ROI_classifications=classified,
            original_image_acquisition_failed=fail,
            original_acquired_but_four_hue_unclassifiable=other,
            original_distinct_observations=n,
            distinct_conditional_joint_photo_probabilities=11,
            original_hue_counts={
                "white": 100, "yellow_orange": 100,
                "red_pink": 100, "blue_purple": classified - 300,
                "mixed_uncertain": fail + other,
            },
            frozen_original_FCP_roi_v4_source_method_attested=True,
            original_photo_identity_and_image_terminal_SHA_audited=True,
            source_inclusion_probability_retained_at_each_photo=True,
            original_warm_FCP_model_inits=1,
            original_full_roi_wall_seconds=2000.5,
            original_full_roi_evaluator_seconds=1900.0,
            original_photo_prefetch={"consumed": n, "download_calls": n,
                                    "permanent_failures": fail},
            image_downloads_OR_ROI_reinferences_in_this_sealer=0,
            old_historical_original_photo_parquets_rewritten=False,
            technical_unknown_ROIs_not_turned_into_other_hues=True,
            true_botanical_plant_flower_or_native_status_validated=False,
            probability_weighted_global_flora_white_rate_estimated=False,
            artifact_SHA256={
                "original_photo_ROIs.parquet": "d" * 64,
                "original_photo_source_coverage.json": "e" * 64,
                "original_photo_hue_counts.csv": "f" * 64,
                "warm_roi_group_000_diagnostic.json": "1" * 64,
            })
        path = self.root / f"g{group}.json"
        path.write_text(json.dumps(d), encoding="utf-8")
        return path

    @staticmethod
    def alter(path, key, value):
        obj = json.loads(path.read_text())
        obj[key] = value
        path.write_text(json.dumps(obj))

    def test_real_worker_progress_is_order_invariant_and_unknown_preserving(self):
        a = self.fixture(0)
        b = self.fixture(2, classified=450, fail=2)
        result = k1.aggregate([b, a], self.out)
        self.assertEqual(result["completed_original_worker_group_IDs"], [0, 2])
        self.assertEqual(result["completed_real_new_original_photo_attempts"], 2048)
        self.assertEqual(result["completed_automatic_four_hue_classifications"], 850)
        self.assertEqual(result["completed_original_acquisition_failures"], 5)
        self.assertEqual(result["completed_acquired_but_automatic_four_hue_unclassifiable"], 1193)
        self.assertEqual(result["completed_mixed_unknown_photo_count"], 1198)
        self.assertEqual(result["completed_hue_photo_counts_CONDITIONAL_ON_AUTOMATIC_CLASSIFICATION"]["white"], 200)
        self.assertTrue(result["this_is_NOT_a_crossarm_photo_ID_deduplicated_global_total"])
        self.assertTrue(result["global_wild_flora_or_blooming_petal_status_not_independently_verified"])
        self.assertEqual(len(result["not_yet_sealed_worker_groups"]), 73)
        path = self.out / "k1_finished_original_group_counts.csv"
        with path.open(newline="") as stream:
            self.assertEqual(len(list(csv.DictReader(stream))), 2)
        self.assertTrue((self.out / "k1_cumulative_original_ROI_receipt.json").is_file())

    def test_final_group_is_partial_not_a_missing_47_photos(self):
        last = self.fixture(74, classified=330, fail=1)
        row = k1.audit_receipt(last)
        self.assertEqual(row["new_original_photo_attempts"], 977)
        self.assertEqual(row["original_worker_group"], 74)

    def test_duplicate_group_and_immutable_output_rejected(self):
        a = self.fixture(0)
        with self.assertRaisesRegex(ValueError, "Duplicate signed K1"):
            k1.aggregate([a, a], self.out)
        self.assertFalse(self.out.exists())
        k1.aggregate([a], self.out)
        with self.assertRaisesRegex(ValueError, "Refuse to rewrite immutable"):
            k1.aggregate([a], self.out)

    def test_wrong_source_SHA_and_missing_unknown_must_fail_closed(self):
        a = self.fixture(0)
        self.alter(a, "original_source_manifest_SHA256", "0" * 64)
        with self.assertRaisesRegex(ValueError, "source SHA"):
            k1.aggregate([a], self.out)
        self.assertFalse(self.out.exists())
        a = self.fixture(0)
        self.alter(a, "original_acquired_but_four_hue_unclassifiable", 0)
        with self.assertRaisesRegex(ValueError, "nonresponse"):
            k1.aggregate([a], self.out)
        self.assertFalse(self.out.exists())

    def test_pending_group_cannot_be_pretended_completed(self):
        a = self.fixture(0)
        self.alter(a, "actual_new_original_photo_attempts", 800)
        with self.assertRaisesRegex(ValueError, "nonresponse counts"):
            k1.aggregate([a], self.out)
        self.assertFalse(self.out.exists())

    def test_missing_terminal_sums_and_false_botanical_confirmation_fail(self):
        a = self.fixture(0)
        self.alter(a, "original_photo_terminal_source_receipts_SHA256", [])
        with self.assertRaisesRegex(ValueError, "terminal/source"):
            k1.aggregate([a], self.out)
        self.assertFalse(self.out.exists())
        a = self.fixture(0)
        self.alter(a, "true_botanical_plant_flower_or_native_status_validated", True)
        with self.assertRaisesRegex(ValueError, "provenance"):
            k1.aggregate([a], self.out)
        self.assertFalse(self.out.exists())

    def test_receipt_aggregate_digest_is_content_based(self):
        a = self.fixture(0)
        b = self.fixture(1, classified=450)
        first = k1.aggregate([a, b], self.out)
        with TemporaryDirectory() as place:
            second = k1.aggregate([b, a], Path(place) / "another")
        self.assertEqual(first["source_original_group_receipt_aggregate_SHA256"],
                         second["source_original_group_receipt_aggregate_SHA256"])


if __name__ == "__main__":
    unittest.main()
