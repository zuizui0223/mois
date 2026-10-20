"""Synthetic contracts; these tests are not measured botanical evidence."""
from pathlib import Path
import importlib.util
import unittest
import pandas as pd

PATH = Path(__file__).resolve().parents[1] / 'scripts/acquisition/mois_lowres_screen_shadow.py'

class TestShadowScreen(unittest.TestCase):
    def setUp(self):
        if not PATH.is_file():
            self.fail('Missing shadow-screen implementation')
        spec = importlib.util.spec_from_file_location('shadow_test_module', PATH)
        self.m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.m)

    def test_negative_screen_does_not_skip_or_change_baseline(self):
        calls = []
        result = {'flower_mask': object(), 'automated_colour_state_status': 'automated_colour_state_admitted'}
        class Full:
            def measure(self, image):
                calls.append('full')
                return result
        def screen(image):
            calls.append('screen')
            return 0.0, 0
        records = []
        wrapper = self.m.ShadowEstimator(Full(), screen, records)
        self.assertIs(wrapper.measure(object()), result)
        self.assertEqual(calls, ['full', 'screen'])
        self.assertEqual(records[0]['screen_status'], 'ok')
        self.assertEqual(records[0]['screen_max_confidence'], 0.0)
        self.assertFalse(records[0]['production_filter_enabled'])

    def test_screen_failure_cannot_destroy_original_measured_output(self):
        expected = {'flower_mask': 'immutable'}
        class Full:
            def measure(self, image): return expected
        def broken(image): raise RuntimeError('temporary model failure')
        records = []
        self.assertIs(self.m.ShadowEstimator(Full(), broken, records).measure(object()), expected)
        self.assertEqual(records[0]['screen_status'], 'error')
        self.assertIn('RuntimeError', records[0]['screen_error'])

    def test_baseline_error_propagates_and_gets_a_screen_receipt(self):
        class Full:
            def measure(self, image): raise ValueError('baseline failed')
        records = []
        with self.assertRaisesRegex(ValueError, 'baseline failed'):
            self.m.ShadowEstimator(Full(), lambda image: (0.8, 1), records).measure(object())
        self.assertEqual(len(records), 1)
        self.assertTrue(records[0]['baseline_runtime_error'])

    def test_invalid_confidence_is_screen_error_not_a_negative(self):
        class Full:
            def measure(self, image): return {'ok': True}
        for score in (float('nan'), -0.1, 1.1):
            records = []
            self.m.ShadowEstimator(Full(), lambda image: (score, 1), records).measure(object())
            self.assertEqual(records[0]['screen_status'], 'error')

    def fixture(self):
        return pd.DataFrame([
            dict(measurement_id='a', morph='white', screen_status='ok', screen_max_confidence=0.0,
                 baseline_roi_seconds=8., screen_seconds=1.),
            dict(measurement_id='b', morph='red_pink', screen_status='ok', screen_max_confidence=0.5,
                 baseline_roi_seconds=10., screen_seconds=1.),
            dict(measurement_id='c', morph='mixed_uncertain', screen_status='ok', screen_max_confidence=0.0,
                 baseline_roi_seconds=2., screen_seconds=1.),
            dict(measurement_id='d', morph='blue_purple', screen_status='error', screen_max_confidence=None,
                 baseline_roi_seconds=10., screen_seconds=1.),
        ])

    def test_report_separates_white_and_chromatic_misses_unknowns_not_negative_flowers(self):
        r = self.m.summarize(self.fixture())
        t = r['thresholds']['0.05']
        self.assertEqual(r['photos_with_baseline_and_shadow'], 4)
        self.assertEqual(t['hypothetical_skipped_images'], 2)
        self.assertEqual(t['baseline_classifiable_lost'], 1)
        self.assertEqual(t['baseline_unknown_in_skipped'], 1)
        self.assertEqual(t['baseline_classifiable_retention'], 2 / 3)
        self.assertEqual(t['by_hue']['white']['lost'], 1)
        self.assertEqual(t['by_hue']['red_pink']['lost'], 0)
        self.assertIsNone(t['by_hue']['yellow_orange']['retention'])
        self.assertFalse(r['production_filter_enabled'])
        self.assertFalse(r['independent_botanical_flowering_validation_completed'])
        self.assertEqual(t['projected_net_roi_seconds_saved_NOT_realized'], 6.0)

    def test_all_shadow_failures_are_kept_by_hypothetical_filter(self):
        f = self.fixture()
        f['screen_status'] = 'error'
        f['screen_max_confidence'] = float('nan')
        r = self.m.summarize(f)
        self.assertEqual(r['thresholds']['0.05']['hypothetical_skipped_images'], 0)
        self.assertEqual(r['thresholds']['0.05']['baseline_classifiable_lost'], 0)

    def test_duplicate_or_missing_measurement_ids_rejected(self):
        f = self.fixture()
        f.loc[1, 'measurement_id'] = 'a'
        with self.assertRaisesRegex(ValueError, 'identity'):
            self.m.summarize(f)
        f.loc[1, 'measurement_id'] = ''
        with self.assertRaisesRegex(ValueError, 'identity'):
            self.m.summarize(f)

    def test_claimed_bad_or_nonfinite_timings_rejected(self):
        for bad in (float('nan'), -1., float('inf')):
            f = self.fixture()
            f.loc[0, 'screen_seconds'] = bad
            with self.assertRaisesRegex(ValueError, 'timing'):
                self.m.summarize(f)

    def test_no_biological_colour_can_be_filled_from_screen(self):
        f = self.fixture()
        f.loc[0, 'morph'] = 'unknown_renamed'
        with self.assertRaisesRegex(ValueError, 'hue'):
            self.m.summarize(f)

    def test_empty_acquired_shadow_is_valid_without_inventing_colour(self):
        f = self.fixture().iloc[:0]
        r = self.m.summarize(f)
        self.assertEqual(r['photos_with_baseline_and_shadow'], 0)
        self.assertIsNone(r['thresholds']['0.05']['baseline_classifiable_retention'])

    def test_evaluator_calls_original_main_and_returns_original_frame(self):
        import hashlib, tempfile
        from PIL import Image
        fake_source = """from pathlib import Path
import hashlib
import pandas as pd
from PIL import Image
DEFAULT_MEASUREMENT=DEFAULT_EXECUTION=DEFAULT_ROI=DEFAULT_LOCKED=Path('.')
class FrozenFlowerColourEstimator:
    def __init__(self,*args,**kwargs): pass
    def measure(self,image): return {'baseline':'unchanged'}
def parse_args(): raise AssertionError('must inject only paths')
def main():
    a=parse_args()
    estimator=FrozenFlowerColourEstimator(a.detector_weight,a.efficient_sam_dir,{},torch_threads=1)
    rows=[]
    for row in pd.read_csv(a.worker_manifest,dtype=str).itertuples(index=False):
        image=a.images_dir/row.image_filename
        with Image.open(image) as im: estimator.measure(im)
        rows.append({'measurement_id':row.measurement_id,'morph':'white',
                     'image_sha256':hashlib.sha256(image.read_bytes()).hexdigest()})
    pd.DataFrame(rows).to_csv(a.output_csv,index=False)
    return 0
"""
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            original=root/'original.py'
            original.write_text(fake_source)
            raw=original.read_bytes()
            self.m.SOURCE_BLOB=hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
            self.m.make_screen=lambda weight,contract: lambda image:(0.,0)
            Image.new('RGB',(8,8),(240,240,240)).save(root/'sample.png')
            pd.DataFrame([dict(measurement_id='abc',image_filename='sample.png',photo_license='cc-by')]).to_csv(root/'worker.csv',index=False)
            frame=self.m.evaluate_shadow(root/'worker.csv',root,original,root/'weight',root,
                                         root/'out.csv',sidecar=root/'sidecar')
            self.assertEqual(frame.morph.tolist(),['white'])
            sidecar=pd.read_csv(root/'sidecar/screen_shadow.csv')
            self.assertEqual(sidecar.morph.tolist(),['white'])
            self.assertEqual(sidecar.screen_max_confidence.tolist(),[0.])
            self.assertEqual(sidecar.image_sha256.tolist(),frame.image_sha256.tolist())
            self.assertEqual(original.read_text(),fake_source)

if __name__ == '__main__':
    unittest.main()
