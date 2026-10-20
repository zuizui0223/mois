"""Check recovery wiring; synthetic fixtures are not botanical measurements."""
import ast
import csv
import importlib.util
import json
from pathlib import Path
import re
import tempfile
import textwrap
import unittest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT/'scripts/acquisition/mois_completed_pixel_checkpoint.py'
WORKFLOW = ROOT/'.github/workflows/mois-completed-pixel-checkpoint.yml'
P0016 = (39374811084, 73872, '16baf617396dd7b388debc65fd25eb8f9041f49ca5e6f78d0b3e2898f6dfeed5')
P0032 = (39401484641, 148071, '3628794fbb70f9e1515dabea6db731862265c0ad6c251732cfd9ab4df4badce8')
GEO05 = (39417990675, 22278, '3b80efe07287c026a34eec7636038e61fa5aa4e16d49aeb139a533bb46c8feda')
RUNS = [39353343403, 39371199127, 39369517433, 39374811084, 39401484641, 39417990675]
NAMES = [
    'MOIS second nested 93m observation RGB p0004 after first source',
    'MOIS third nested 93m observation RGB p0008 after second full integration',
    'MOIS new tier03 pixels and lowres screen shadow without image exclusion',
    'MOIS full remaining tier03 original RGB p0016 without lowres exclusion',
    'MOIS fifth disjoint tier04 original RGB p0032 without lowres exclusion',
    'MOIS sixth geographic enrichment tier05 original RGB (tropics and south)',
]

class TestFourWaveRecovery(unittest.TestCase):
    def setUp(self):
        spec=importlib.util.spec_from_file_location('checkpoint_wiring', SCRIPT)
        self.module=importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.module)
        self.workflow=WORKFLOW.read_text(encoding='utf-8')

    def test_new_wave_has_exact_source_manifest_and_real_planned_count(self):
        self.assertEqual(self.module.WAVES.get('p0016'), P0016)
        self.assertEqual(self.module.WAVES.get('p0032'), P0032)
        self.assertEqual(self.module.WAVES.get('geo05'), GEO05)
        self.assertEqual(set(self.module.WAVES), {'p0004','p0008','shadow','p0016','p0032','geo05'})
        self.assertEqual(self.module.WAVES['p0004'][0:2],(39353343403,18119))
        self.assertEqual(self.module.WAVES['p0008'][0:2],(39371199127,36375))
        self.assertEqual(self.module.WAVES['shadow'][0:2],(39369517433,569))

    def test_downloader_and_collector_cover_the_same_four_waves(self):
        blocks=re.findall(r"python - <<'PY'\n(.*?)\n\s+PY(?:\n|$)",self.workflow,re.S)
        specs=None
        for block in blocks:
            tree=ast.parse(textwrap.dedent(block))
            for node in ast.walk(tree):
                if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='specs' for t in node.targets):
                    specs=ast.literal_eval(node.value)
        self.assertIsNotNone(specs)
        self.assertEqual({s[0] for s in specs}, {'p0004','p0008','shadow','p0016','p0032','geo05'})
        self.assertEqual({s[0]:s[1] for s in specs}, {k:v[0] for k,v in self.module.WAVES.items()})
        fourth=next(s for s in specs if s[0]=='p0016')
        self.assertEqual(fourth[2:4], ('mois-93m-p0016-precolour-observation','mois-93m-p0016-RGB-group-'))
        fifth=next(s for s in specs if s[0]=='p0032')
        self.assertEqual(fifth[2:4], ('mois-93m-p0032-precolour-observation','mois-93m-p0032-RGB-group-'))
        sixth=next(s for s in specs if s[0]=='geo05')
        self.assertEqual(sixth[2:4], ('mois-93m-geo05-precolour-observation','mois-93m-geo05-RGB-group-'))
        self.assertIn("'--paginate'",self.workflow)
        self.assertIn('source_artifact_snapshots', SCRIPT.read_text(encoding='utf-8'))

    def test_only_pinned_main_runs_can_trigger_automatic_recovery(self):
        self.assertIn('  workflow_run:', self.workflow)
        for name in NAMES:
            self.assertIn('      - "'+name+'"',self.workflow)
        self.assertIn('    types: [completed]',self.workflow)
        gate=self.workflow.split('  recover:',1)[1].split('    runs-on:',1)[0]
        self.assertIn("github.event_name != 'pull_request'",gate)
        self.assertIn("github.event_name != 'workflow_run'",gate)
        self.assertIn("github.event.workflow_run.head_repository.full_name == 'zuizui0223/mois'",gate)
        self.assertIn("github.event.workflow_run.head_branch == 'main'",gate)
        literal='['+','.join(map(str,RUNS))+']'
        self.assertIn("contains(fromJSON('"+literal+"'), github.event.workflow_run.id)",gate)

    def test_source_backed_milestone_push_can_refresh_checkpoint_without_inference(self):
        s=self.workflow
        self.assertIn("  push:\n    branches: [main]\n    paths:",s)
        self.assertIn("      - docs/MOIS_CHECKPOINT_REQUEST_20261020.md",s)
        self.assertIn("  workflow_dispatch:",s)
        self.assertIn("      - name: Validate actual terminals and append only measured IDs to baseline projection",s)
        self.assertNotIn("mois_full_index_real_original_photo_measure work",s)

    def test_synthetic_waiting_wave_is_not_measured_or_nonwhite(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'source';workers=root/'workers'
            (source/'shards').mkdir(parents=True);workers.mkdir()
            file=source/'shards/shard-000.csv'
            with file.open('w',newline='') as handle:
                writer=csv.DictWriter(handle,fieldnames=['photo_id','measurement_id','shard'])
                writer.writeheader();writer.writerow(dict(photo_id='901',measurement_id='fixture',shard='0'))
            manifest=dict(schema='mois_full_source_observation_sampling_cohort_v1',
                colour_measured=False,full_frame_photo_ids=93209351,worker_ids=[0],work_units=1,
                new_photo_attempts_planned=1,photo_shard_sha256={'0':self.module.sha(file)})
            (source/'manifest.json').write_text(json.dumps(manifest))
            rows,result=self.module.collect_wave(source,workers,{'synthetic_fixture_only':True})
            self.assertEqual(rows,[])
            self.assertEqual(result['attempted'],0)
            self.assertEqual(result['downloaded'],0)
            self.assertEqual(result['classifiable'],0)
            self.assertEqual(result['pending_photo_attempts'],1)
            self.assertFalse(result['wave_complete'])

    def test_recovery_never_invokes_image_inference_or_changes_acquisition(self):
        self.assertNotIn('mois_full_index_real_original_photo_measure work',self.workflow)
        self.assertNotIn('bootstrap_fcp_original_roi_v4',self.workflow)
        self.assertIn('--prerelease',self.workflow)
        self.assertIn('cancel-in-progress: false',self.workflow)
        self.assertIn('partial_completion_is_not_representative_probability_sample',SCRIPT.read_text())

if __name__=='__main__':unittest.main()
