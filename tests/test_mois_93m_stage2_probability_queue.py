"""Stage-2 93m photo queue only samples screen negatives with known probability."""
import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from scripts.acquisition import mois_93m_stage2_probability_queue as q


class Stage2ProbabilityQueueTests(unittest.TestCase):
    def setUp(self):
        self.tmp = TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.input = self.root / "signed_stage1.csv"
        self.policy = self.root / "policy.json"
        self.out = self.root / "decision"
        self.write_rows([
            ("1", "O1", "global_probability", "0.0032", "tropical", "negative"),
            ("2", "O1", "global_probability", "0.0032", "tropical", "negative"),
            ("3", "O2", "geo05", "0.0008", "southern", "candidate"),
            ("4", "O3", "geo05", "0.0008", "southern", "uncertain"),
            ("5", "O4", "geo05", "0.0008", "southern", "error"),
        ])
        self.policy.write_text(json.dumps({
            "schema": "mois_93m_prespecified_negative_audit_policy_v1",
            "negative_audit_rates": {"tropical": 0.2, "southern": 0.5}
        }))

    def write_rows(self, rows):
        with self.input.open("w", newline="") as out:
            w = csv.writer(out, lineterminator="\n")
            w.writerow(q.INPUT_COLUMNS)
            w.writerows(rows)

    def signed(self, *, mode="shadow", review=None):
        return q.plan(self.input, self.policy, self.out, mode=mode,
                      review_receipt=review,
                      expected_input_sha256=q.sha256(self.input))

    def decisions(self):
        with (self.out / "all_photo_screen_decisions.csv").open(newline="") as stream:
            return list(csv.DictReader(stream))

    def test_shadow_never_discards_without_botanical_validation(self):
        receipt = self.signed()
        self.assertEqual(receipt["counts"]["all_screened"], 5)
        self.assertEqual(receipt["counts"]["queued_full_roi"], 5)
        self.assertEqual(receipt["counts"]["negative_not_sampled"], 0)
        self.assertTrue(all(x["full_roi_selected"] == "true" for x in self.decisions()))
        self.assertEqual(self.decisions()[0]["joint_design_probability"], "0.0032000000000000002")
        self.assertIsNone(receipt["independent_review_sha256"])

    def test_stage_mode_requires_independent_photo_specific_review(self):
        with self.assertRaisesRegex(ValueError, "independent botanical validation"):
            self.signed(mode="staged")
        self.assertFalse(self.out.exists())

    def test_authenticated_independent_review_unlocks_randomized_negatives(self):
        review = self.root / "review.json"
        review.write_text(json.dumps({
            "schema": "mois_independent_botanical_screen_review_v1",
            "screen_input_sha256": q.sha256(self.input),
            "audit_policy_sha256": q.sha256(self.policy),
            "review_basis": "independent_photo_specific_open_flower_review",
            "independent_review_complete": True,
            "review_by_stratum": {
                "tropical": {"reviewed_negative_count": 100, "white_flower_recall_lower95": 0.96},
                "southern": {"reviewed_negative_count": 100, "white_flower_recall_lower95": 0.99}
            }
        }))
        receipt = self.signed(mode="staged", review=review)
        rows = self.decisions()
        self.assertEqual(rows[0]["full_roi_selected"], rows[1]["full_roi_selected"])
        self.assertEqual(rows[0]["stage2_selection_probability"], "0.20000000000000001")
        self.assertEqual(rows[0]["joint_design_probability"], format(0.0032 * 0.2, ".17g"))
        self.assertTrue(all(x["full_roi_selected"] == "true" for x in rows[2:]))
        self.assertEqual(receipt["counts"]["all_screened"], 5)
        self.assertEqual(receipt["counts"]["queued_full_roi"] +
                         receipt["counts"]["negative_not_sampled"], 5)
        self.assertEqual(receipt["output_sha256"]["all_photo_screen_decisions.csv"],
                         q.sha256(self.out / "all_photo_screen_decisions.csv"))

    def test_validation_with_low_hue_recall_cannot_enable_skip(self):
        r = self.root / "review.json"
        r.write_text(json.dumps({
            "schema": "mois_independent_botanical_screen_review_v1",
            "screen_input_sha256": q.sha256(self.input),
            "audit_policy_sha256": q.sha256(self.policy),
            "review_basis": "independent_photo_specific_open_flower_review",
            "independent_review_complete": True,
            "review_by_stratum": {
                "tropical": {"reviewed_negative_count": 100, "white_flower_recall_lower95": 0.90},
                "southern": {"reviewed_negative_count": 100, "white_flower_recall_lower95": 0.99}
            }
        }))
        with self.assertRaisesRegex(ValueError, "recall lower bound"):
            self.signed(mode="staged", review=r)

    def test_stale_source_sha_denied_before_assignment(self):
        with self.assertRaisesRegex(ValueError, "checksum differs"):
            q.plan(self.input, self.policy, self.out,
                   expected_input_sha256="0" * 64)
        self.assertFalse(self.out.exists())

    def test_photo_duplication_or_unsorted_id_fails_closed(self):
        self.write_rows([
            ("2", "O1", "base", "0.0032", "tropical", "candidate"),
            ("1", "O2", "base", "0.0032", "tropical", "candidate")
        ])
        with self.assertRaisesRegex(ValueError, "sorted and unique"):
            self.signed()
        self.assertFalse(self.out.exists())

    def test_unknown_state_and_invalid_source_probability_denied(self):
        self.write_rows([("1", "O1", "base", "0.0032", "tropical", "certified_nonflower")])
        with self.assertRaisesRegex(ValueError, "Unknown screen state"):
            self.signed()
        self.write_rows([("1", "O1", "base", "NaN", "tropical", "candidate")])
        with self.assertRaisesRegex(ValueError, "Probability"):
            self.signed()

    def test_existing_outputs_are_never_overwritten(self):
        self.out.mkdir()
        with self.assertRaisesRegex(ValueError, "Refuse to overwrite"):
            self.signed()


if __name__ == "__main__":
    unittest.main()
