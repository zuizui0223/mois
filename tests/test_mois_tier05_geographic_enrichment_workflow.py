"""Static contracts for separate tropical/southern photo sourcing, not RGB data."""
import re
import subprocess
import unittest
from pathlib import Path
import yaml

PATH=Path(__file__).resolve().parents[1]/".github/workflows/mois-tier05-geographic-enrichment-RGB.yml"


class TestGeographicEnrichmentOriginalRGB(unittest.TestCase):
    def setUp(self):
        self.assertTrue(PATH.is_file(), "Separate geographic original RGB workflow missing")
        self.text=PATH.read_text(encoding="utf-8")
        self.workflow=yaml.safe_load(self.text)

    def test_only_next_signed_source_tier_not_full_93m_is_used(self):
        for expected in (
            "tier-05-group-*.tar",
            "full['new_target_counts_by_tier']['5']==297967",
            "'--tier','5'",
            "'--geographic-enrichment'",
            "source_tier_total_photo_candidates']==297967",
            "source_tier_probability_lower_bound']==0.0032",
            "q['selection_probability'] is None",
            "geographic_selection_excluded_photo_candidates",
            "mois_queue_from_full93m_target_tiers"
        ):
            self.assertIn(expected,self.text)
        self.assertNotIn("gh release download mois-full-93209351-indexed-photo-observation-source-v1",
                         self.text)

    def test_authenticate_already_measured_and_verify_zero_full_active_queues(self):
        for item in (
            "mois-real-pixel-checkpoint-39416137341-v1",
            "r['combined_unique_photo_attempts']==670083",
            "r['combined_four_hue_classified']==318062",
            "sha(file)==r['output_sha256'][file.name]",
            "('p0004',18119)","('p0008',36375)",
            "('p0016',73872)","('p0032',148071)","('shadow',569)",
            "39353343403","39371199127","39374811084",
            "39401484641","39369517433",
            "verify_disjoint('current','reserved-'+name)",
            "overlapping_photo_ids']==0",
            "overlapping_observation_ids']==0",
        ):
            self.assertIn(item,self.text)

    def test_never_call_hard_native_or_photo_screen_filter(self):
        s=self.text
        self.assertIn("native_photo_status_not_inferred_or_hard_filtered",s)
        self.assertIn("geographic_enrichment_strictly_separate_from_complete_probability_arm",s)
        measure=self.workflow["jobs"]["measure"]
        self.assertEqual(measure["strategy"]["max-parallel"],1)
        self.assertIs(measure["strategy"]["fail-fast"],False)
        shell="\n".join(step.get("run","") for step in measure["steps"])
        self.assertIn("mois_full_index_real_original_photo_measure work",shell)
        self.assertIn("74bf69a52cd68bfc7dbe6158f4a681235caa6ca0",shell)
        self.assertNotIn("mois_lowres_screen_shadow",shell)
        self.assertIn("FCP_PROVIDER_GLOBAL_RATE_FILE",shell)

    def test_only_signed_full_photo_results_can_be_published(self):
        s=self.text
        collect=self.workflow["jobs"]["collect"]
        self.assertEqual(set(collect["needs"]),{"source","measure"})
        self.assertIn("mois_full_index_real_original_photo_measure collect",s)
        self.assertIn("source/parallel_*_receipt.json",s)
        self.assertIn("mois-full93m-observation-sixth-geo-enrichment-original-RGB-tier05-v1",s)
        self.assertIn("resume_run",s)
        self.assertFalse(self.workflow["concurrency"]["cancel-in-progress"])

    def test_embedded_python_and_shell_compile(self):
        n=0
        for job in self.workflow["jobs"].values():
            for step in job["steps"]:
                if "run" not in step:continue
                shell=re.sub(r"\$\{\{.*?\}\}","0",step["run"])
                checked=subprocess.run(["bash","-n"],input=shell,
                    capture_output=True,text=True)
                self.assertEqual(checked.returncode,0,checked.stderr)
                for code in re.findall(r"python - <<'PY'\n(.*?)\nPY(?:\n|$)",shell,re.S):
                    compile(code,"signed-geographic-workflow","exec");n+=1
        self.assertGreaterEqual(n,1)

if __name__=="__main__":unittest.main()
