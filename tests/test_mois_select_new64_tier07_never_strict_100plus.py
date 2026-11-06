"""MOIS tiny tier07 supplementary sampling rules are outcome-aware and bounded."""
import unittest
from scripts.acquisition.mois_select_new64_tier07_never_strict_100plus import (
    band4,effort_class,rank,STRATA,TARGET_PER_CELL,
    TIER,SCHEMA
)

class RealOriginalPhotoSourceSamplingRules(unittest.TestCase):
    def test_separate_both_tropical_hemispheres_and_extratropical(self):
        self.assertEqual(band4(67),"north_extratropical")
        self.assertEqual(band4(24),"north_extratropical")
        self.assertEqual(band4(0),"north_tropical")
        self.assertEqual(band4(-22),"south_tropical")
        self.assertEqual(band4(-24),"south_extratropical")
        self.assertEqual(band4(-68),"south_extratropical")
        self.assertEqual(band4("missing"),"unknown")
        self.assertEqual(len(STRATA),4)
        self.assertEqual(TIER,7)

    def test_prior_one_to_four_attempts_only(self):
        self.assertEqual(effort_class(1),"one_attempt")
        self.assertEqual(effort_class(2),"two_to_four_attempts")
        self.assertEqual(effort_class(4),"two_to_four_attempts")
        self.assertIsNone(effort_class(0))
        self.assertIsNone(effort_class(5))

    def test_salt_reproducible_and_independent_original_source_ids(self):
        self.assertEqual(rank("taxon","north_tropical","one_attempt","71"),
                         rank("taxon","north_tropical","one_attempt","71"))
        self.assertNotEqual(rank("taxon","north_tropical","one_attempt","71"),
                            rank("taxon","north_tropical","one_attempt","72"))
        self.assertNotEqual(rank("taxon","north_tropical","one_attempt","71"),
                            rank("photo","north_tropical","one_attempt","71"))
        self.assertEqual(TARGET_PER_CELL,8)
        self.assertIn("minimal_64",SCHEMA)

if __name__=="__main__":
    unittest.main()
