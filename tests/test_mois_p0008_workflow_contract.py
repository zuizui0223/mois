"""Static integration gates for actual p0004 -> p0008 original-image acquisition.

Only verifies wired source filenames, exact p0004 schema, previous SHA,
nested source probability and the right upstream workflow trigger.
"""
from pathlib import Path
import unittest


WORKFLOW = (
    Path(__file__).resolve().parents[1] / ".github" / "workflows" /
    "mois-full93m-third-nested-p0008-RGB.yml"
)


class TestNestedP0008InputReceipt(unittest.TestCase):
    def setUp(self):
        self.workflow = WORKFLOW.read_text(encoding="utf-8")

    def test_only_real_p0004_pixel_union_is_accepted(self):
        s = self.workflow
        self.assertIn(
            'mois-genuine630602-original-pixels-first-two-93m-waves-v1', s)
        self.assertIn(
            '--pattern mois_p0004_actual_photo_append_receipt.json', s)
        self.assertIn(
            "root/'mois_p0004_actual_photo_append_receipt.json'", s)
        self.assertIn(
            "doc['schema']=='mois_FCP_and_full93m_second_p0004_real_original_photo_append_v1'",
            s)
        self.assertIn("doc['output_sha256'][file.name]", s)
        self.assertNotIn('mois_p0008_actual_photo_append_receipt.json', s)
        self.assertNotIn('sha256_original_pixel_outputs', s)

    def test_p0008_uses_previous_completed_measurement_without_reinference(self):
        s = self.workflow
        self.assertIn(
            '"MOIS source-photo append genuine second p0004 18119 to 612483 original ROI results"', s)
        self.assertIn("assert len(old)==630602", s)
        self.assertIn("assert doc['old_four_hue_classified']==295801", s)
        self.assertIn("mois-complete-93209351-original-photo-target-tiers-v1",s)
        self.assertIn("tier-02-group-*.tar",s)
        self.assertIn("mois_queue_from_full93m_target_tiers",s)
        self.assertIn("'--tier','2'",s)
        self.assertIn("'--compute-bucket-count','1'",s)
        self.assertIn("d['source_tier_total_photo_candidates']==36583",s)
        self.assertIn("d['source_tier_compute_bucket_total_photos']==36583",s)
        self.assertNotIn("gh release download mois-full-93209351-indexed-photo-observation-source-v1",s)
        self.assertIn("assert d['selection_probability']==0.0008", s)
        self.assertIn(
            "mois-full93m-observation-third-original-RGB-p0008-v1", s)
        self.assertIn("if: needs.freeze-source-queue.outputs.ready == 'true'", s)


if __name__ == "__main__":
    unittest.main()
