"""Synthetic queue-identity tests; no image access or measured colour claims."""
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

MODULE = Path(__file__).resolve().parents[1] / 'scripts/acquisition/mois_verify_parallel_photo_queues.py'

class TestParallelPhotoQueues(unittest.TestCase):
    def setUp(self):
        self.assertTrue(MODULE.is_file(), 'Parallel source-queue verification is not implemented')
        spec = importlib.util.spec_from_file_location('parallel_queues', MODULE)
        self.m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.m)
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def source(self, name, rows):
        root = self.root / name
        (root / 'shards').mkdir(parents=True)
        file = root / 'shards/shard-000.csv'
        with file.open('w', newline='') as out:
            w = csv.writer(out)
            w.writerow(['photo_id', 'observation_id', 'measurement_id', 'shard'])
            w.writerows((p, o, 'MOIS-'+str(p), 0) for p, o in rows)
        receipt = {'schema': 'mois_full_source_observation_sampling_cohort_v1',
                   'full_frame_photo_ids': 93209351, 'colour_measured': False,
                   'new_photo_attempts_planned': len(rows), 'worker_ids': [0],
                   'work_units': 1,
                   'photo_shard_sha256': {'0': hashlib.sha256(file.read_bytes()).hexdigest()}}
        (root / 'manifest.json').write_text(json.dumps(receipt))
        return root

    def test_disjoint_in_progress_queue_needs_no_measured_release(self):
        new = self.source('new', [(4, 'obs-c'), (5, 'obs-c')])
        active = self.source('active', [(1, 'obs-a'), (2, 'obs-b')])
        r = self.m.verify_disjoint(new, active)
        self.assertEqual(r['new_photo_targets'], 2)
        self.assertEqual(r['active_photo_targets'], 2)
        self.assertEqual(r['overlapping_photo_ids'], 0)
        self.assertEqual(r['overlapping_observation_ids'], 0)
        self.assertFalse(r['active_targets_counted_as_measured_photos'])
        self.assertFalse(r['image_pixels_read'])

    def test_same_photo_id_rejected_even_with_other_observation(self):
        new = self.source('new', [(1, 'different')])
        active = self.source('active', [(1, 'original')])
        with self.assertRaisesRegex(ValueError, 'photo ID overlap'):
            self.m.verify_disjoint(new, active)

    def test_same_observation_other_photo_rejected(self):
        new = self.source('new', [(2, 'same')])
        active = self.source('active', [(1, 'same')])
        with self.assertRaisesRegex(ValueError, 'observation ID overlap'):
            self.m.verify_disjoint(new, active)

    def test_modified_shard_rejected_before_using_ids(self):
        new = self.source('new', [(2, 'new')])
        active = self.source('active', [(1, 'old')])
        with (active / 'shards/shard-000.csv').open('a') as out:
            out.write('3,extra,MOIS-3,0\n')
        with self.assertRaisesRegex(ValueError, 'SHA'):
            self.m.verify_disjoint(new, active)

    def test_count_mismatch_rejected(self):
        new = self.source('new', [(2, 'new')])
        active = self.source('active', [(1, 'old')])
        p = active / 'manifest.json'
        r = json.loads(p.read_text()); r['new_photo_attempts_planned'] = 9
        p.write_text(json.dumps(r))
        with self.assertRaisesRegex(ValueError, 'count'):
            self.m.verify_disjoint(new, active)

    def test_duplicate_photo_within_one_queue_rejected(self):
        new = self.source('new', [(2, 'new'), (2, 'other')])
        active = self.source('active', [(1, 'old')])
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            self.m.verify_disjoint(new, active)

    def test_empty_observation_id_rejected(self):
        new = self.source('new', [(2, '')])
        active = self.source('active', [(1, 'old')])
        with self.assertRaisesRegex(ValueError, 'Missing'):
            self.m.verify_disjoint(new, active)

    def test_measured_manifest_cannot_masquerade_as_source_queue(self):
        new = self.source('new', [(2, 'new')])
        active = self.source('active', [(1, 'old')])
        p = active / 'manifest.json'
        r = json.loads(p.read_text()); r['colour_measured'] = True
        p.write_text(json.dumps(r))
        with self.assertRaisesRegex(ValueError, 'pre-colour'):
            self.m.verify_disjoint(new, active)

if __name__ == '__main__': unittest.main()
