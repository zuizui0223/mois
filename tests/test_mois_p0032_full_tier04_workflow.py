"""Contract-only tests for the next disjoint original-photo RGB wave.

All source identities and colour counters are checked again on real pixels by
the frozen source queue, existing batch runner, and terminal collector.
"""
import re
import subprocess
import unittest
from pathlib import Path

import yaml

PATH = (Path(__file__).resolve().parents[1] / ".github/workflows" /
        "mois-full93m-fifth-tier04-p0032-RGB.yml")
SOURCE_TIER = "4"
COUNT = 148983
PRIOR_COUNT = 653699
PRIOR_CLASSIFIED = 311781
ACTIVE = (("p0004", 39353343403, 18119),
          ("p0008", 39371199127, 36375),
          ("p0016", 39374811084, 73872),
          ("shadow", 39369517433, 569))


class TestFullOriginalTier04Acquisition(unittest.TestCase):
    def setUp(self):
        self.assertTrue(PATH.is_file(), "Missing tier04 original-photo acquisition workflow")
        self.text = PATH.read_text(encoding="utf-8")
        self.flow = yaml.safe_load(self.text)

    def test_new_target_tier_only_not_entire_93m_reread(self):
        s = self.text
        for expected in (
            "mois-complete-93209351-original-photo-target-tiers-v1",
            "tier-04-group-*.tar",
            "'--tier','4'",
            "'--compute-bucket-count','1'",
            "full['new_target_counts_by_tier']['4']==148983",
            "q['source_tier_total_photo_candidates']==148983",
            "q['source_tier_compute_bucket_total_photos']==148983",
            "q['source_tier_probability_lower_bound']==0.0016",
            "q['selection_probability']==0.0032",
            "mois_queue_from_full93m_target_tiers",
        ):
            self.assertIn(expected, s)
        self.assertNotIn(
            "gh release download mois-full-93209351-indexed-photo-observation-source-v1",
            s)

    def test_previous_measured_pixel_ids_are_authenticated_and_not_reacquired(self):
        s = self.text
        for expected in (
            "mois-real-pixel-checkpoint-39400235870-v1",
            "r['combined_unique_photo_attempts']==653699",
            "r['combined_four_hue_classified']==311781",
            "sha(file)==r['output_sha256'][file.name]",
            "'--previous','prior.csv'",
            "q['previously_measured_photo_ids_skipped']",
            "q['new_photo_attempts_planned']",
            "q['native_photo_status_not_inferred_or_hard_filtered']",
        ):
            self.assertIn(expected, s)

    def test_all_four_preexisting_queues_have_zero_photo_and_observation_overlap(self):
        s = self.text
        for name, run, count in ACTIVE:
            self.assertIn(str(run), s)
            self.assertIn(f"('{name}',{count})", s)
        for expected in (
            "verify_disjoint('current','reserved-'+name)",
            "proof['overlapping_photo_ids']==0",
            "proof['overlapping_observation_ids']==0",
            "source/parallel_*_receipt.json",
        ):
            self.assertIn(expected, s)

    def test_frozen_full_colour_roi_retained_and_no_lowres_or_native_filter(self):
        measure = self.flow["jobs"]["measure"]
        self.assertLessEqual(measure["strategy"]["max-parallel"], 2)
        self.assertEqual(measure["strategy"]["fail-fast"], False)
        s = "\n".join(z.get("run", "") for z in measure["steps"])
        self.assertIn("74bf69a52cd68bfc7dbe6158f4a681235caa6ca0", s)
        self.assertIn("mois_full_index_real_original_photo_measure work", s)
        self.assertIn("FCP_PROVIDER_GLOBAL_RATE_FILE", s)
        self.assertNotIn("mois_lowres_screen_shadow", s)
        self.assertNotIn("native_only", s)
        self.assertIn("inputs.resume_run", self.text)
        self.assertIs(self.flow["concurrency"]["cancel-in-progress"], False)

    def test_fail_closed_collection_and_original_photo_provenance(self):
        collection = self.flow["jobs"]["collect"]
        self.assertEqual(set(collection["needs"]), {"source", "measure"})
        s = "\n".join(z.get("run", "") for z in collection["steps"])
        self.assertIn("mois_full_index_real_original_photo_measure collect", s)
        self.assertIn("source/parallel_*_receipt.json", s)
        self.assertIn("mois-full93m-observation-fifth-original-RGB-p0032-v1", s)
        self.assertNotIn("mois_lowres_screen_shadow", s)

    def test_all_embedded_python_and_shell_statements_compile(self):
        inline_blocks = 0
        for job in self.flow["jobs"].values():
            for step in job["steps"]:
                if "run" not in step:
                    continue
                script = re.sub(r"\$\{\{.*?\}\}", "0", step["run"])
                result = subprocess.run(["bash", "-n"], input=script,
                                        text=True, capture_output=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                for code in re.findall(r"python - <<'PY'\n(.*?)\nPY(?:\n|$)",
                                       script, re.S):
                    compile(code, "frozen-acquisition-workflow", "exec")
                    inline_blocks += 1
        self.assertGreaterEqual(inline_blocks, 1)


if __name__ == "__main__":
    unittest.main()
