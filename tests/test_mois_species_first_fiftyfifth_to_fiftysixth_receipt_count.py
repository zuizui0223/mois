"""Guard the number of actually acquired wave55 original photos against stale wave50 constants."""
import inspect
import unittest
from scripts.analysis.mois_species_first_fiftyfifth_to_fiftysixth import read_fiftyfifth


class FiftyfifthSignedReceiptCountTest(unittest.TestCase):
    def test_one_original_matches_manifest_worker_and_audit_contract(self):
        body = inspect.getsource(read_fiftyfifth)
        self.assertIn('manifest.get("new_photo_attempts_planned") == 1', body)
        self.assertIn('signed.get("actual_independent_original_photo_attempts") == 1', body)
        self.assertIn('len(raw) == 1', body)
        self.assertIn('len(audit) == 1', body)
        self.assertIn('signed["actual_genuine_strict9_numeric_success"]', body)
        self.assertIn('signed["actual_technical_unknowns"]', body)
        self.assertIn('y.loc[~flags, PIX].eq("").all().all()', body)
        self.assertIn('np.allclose(measured.sum(axis=1), 1, atol=0.001)', body)


if __name__ == "__main__":
    unittest.main()
