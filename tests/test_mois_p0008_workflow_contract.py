"""Source identity and nonblocking p0008 acquisition workflow contracts."""
from pathlib import Path
import unittest

WORKFLOW = Path(__file__).resolve().parents[1] / '.github/workflows/mois-full93m-third-nested-p0008-RGB.yml'

class TestNestedP0008InputReceipt(unittest.TestCase):
    def setUp(self):
        self.workflow = WORKFLOW.read_text(encoding='utf-8')

    def test_only_completed_prior_actual_pixels_are_called_measured(self):
        s = self.workflow
        self.assertIn('mois-genuine-FCP447500-plus93m-first18462-photo-ID-union-v1', s)
        self.assertIn("doc['all_genuinely_attempted_unique_original_photo_IDs']==612483", s)
        self.assertIn("doc['all_original_photo_four_hue_classified']==295801", s)
        self.assertIn("doc['completed_output_sha256'][file.name]", s)
        self.assertIn('assert len(old)==612483', s)
        self.assertNotIn('mois-genuine630602-original-pixels-first-two-93m-waves-v1', s)

    def test_parallel_work_requires_photo_and_observation_disjointness(self):
        s = self.workflow
        self.assertIn('mois_verify_parallel_photo_queues import verify_disjoint', s)
        self.assertIn("verify_disjoint('current','active-'+name)", s)
        self.assertIn("[('p0004',18119),('shadow',569)]", s)
        self.assertIn('mois-93m-p0004-precolour-observation --dir active-p0004', s)
        self.assertIn('mois-tier03-shadow-b0-source --dir active-shadow', s)
        self.assertIn('cp source/parallel_*_receipt.json genuine_photos/', s)
        self.assertIn("github.event_name == 'push' && github.ref == 'refs/heads/main'", s)
        self.assertIn('cancel-in-progress: false', s)

    def test_reads_only_new_full_tier_and_keeps_the_frozen_measurement_model(self):
        s = self.workflow
        self.assertIn('tier-02-group-*.tar', s)
        self.assertIn("'--tier','2','--compute-bucket','0','--compute-bucket-count','1'", s)
        self.assertIn("d['source_tier_total_photo_candidates']==36583", s)
        self.assertIn("d['source_tier_compute_bucket_total_photos']==36583", s)
        self.assertIn("d['source_tier_probability_lower_bound']==0.0004", s)
        self.assertIn("d['selection_probability']==0.0008", s)
        self.assertNotIn('gh release download mois-full-93209351-indexed-photo-observation-source-v1', s)
        self.assertIn('74bf69a52cd68bfc7dbe6158f4a681235caa6ca0', s)
        self.assertIn('max-parallel: 4', s)
        self.assertIn("if: needs.freeze-source-queue.outputs.ready == 'true'", s)
        self.assertNotIn('mois_lowres_screen_shadow work', s)

if __name__ == '__main__': unittest.main()
