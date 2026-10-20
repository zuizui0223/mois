"""Execution wiring checks; no synthetic value is a measured photo outcome."""
import re
import subprocess
import unittest
from pathlib import Path
import yaml

PATH = Path(__file__).resolve().parents[1]/'.github/workflows/mois-full93m-fourth-tier03-p0016-RGB.yml'

class TestNextOriginalPhotoWave(unittest.TestCase):
    def setUp(self):
        self.assertTrue(PATH.is_file(), 'Missing full tier03 original-photo acquisition workflow')
        self.text=PATH.read_text()
        self.flow=yaml.safe_load(self.text)

    def test_only_new_tier_and_completed_checkpoint_are_read(self):
        s=self.text
        for expected in ('tier-03-group-*.tar', "'--tier','3'", "'--compute-bucket-count','1'",
                         "full['new_target_counts_by_tier']['3']==74883",
                         'mois-real-pixel-checkpoint-39373763787-v1',
                         "r['combined_unique_photo_attempts']==625340",
                         "r['combined_four_hue_classified']==300828",
                         "sha(file)==r['output_sha256'][file.name]"):
            self.assertIn(expected,s)
        self.assertNotIn('gh release download mois-full-93209351-indexed-photo-observation-source-v1',s)

    def test_source_excludes_all_three_reservations(self):
        s=self.text
        for expected in ('39353343403','39371199127','39369517433',
                         "('p0004',18119)", "('p0008',36375)", "('shadow',569)",
                         "proof=verify_disjoint('current','reserved-'+name)",
                         "q['source_tier_probability_lower_bound']==0.0008",
                         "q['selection_probability']==0.0016"):
            self.assertIn(expected,s)

    def test_full_model_unchanged_and_bounded_workers(self):
        j=self.flow['jobs']['measure']
        self.assertEqual(j['strategy']['max-parallel'],2)
        s='\n'.join(step.get('run','') for step in j['steps'])
        self.assertIn('74bf69a52cd68bfc7dbe6158f4a681235caa6ca0',s)
        self.assertIn('mois_full_index_real_original_photo_measure work',s)
        self.assertIn('FCP_PROVIDER_GLOBAL_RATE_FILE',s)
        self.assertNotIn('mois_lowres_screen_shadow',s)
        self.assertIn('inputs.resume_run',s)
        self.assertFalse(self.flow['concurrency']['cancel-in-progress'])

    def test_only_complete_reconciled_results_can_publish(self):
        j=self.flow['jobs']['collect']
        self.assertEqual(set(j['needs']),{'source','measure'})
        s='\n'.join(step.get('run','') for step in j['steps'])
        self.assertIn('mois_full_index_real_original_photo_measure collect',s)
        self.assertIn('source/parallel_*_receipt.json',s)
        self.assertIn('mois-full93m-observation-fourth-original-RGB-p0016-v1',s)

    def test_embedded_python_and_shell_compile_without_execution(self):
        count=0
        for job in self.flow['jobs'].values():
            for step in job['steps']:
                if 'run' not in step:continue
                shell=re.sub(r'\$\{\{.*?\}\}', '0', step['run'])
                checked=subprocess.run(['bash','-n'],input=shell,text=True,capture_output=True)
                self.assertEqual(checked.returncode,0,checked.stderr)
                for code in re.findall(r"python - <<'PY'\n(.*?)\nPY(?:\n|$)",shell,re.S):
                    compile(code,'workflow-inline','exec');count+=1
        self.assertGreaterEqual(count,1)

if __name__=='__main__':unittest.main()
