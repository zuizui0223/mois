"""Source taxon/actual photo outcomes, rare species and unknown safeguards."""
from __future__ import annotations
import unittest

from scripts.analysis import mois_full42111_taxon_pixel_coverage as m


class TestActualWorldOriginalImageTaxonomicBreadth(unittest.TestCase):
    def test_full_original_93m_unbiased_source_denominator_is_not_a_pixel_count(self):
        self.assertEqual(m.CENSUS_SHA,"6adab400c72932ab1079cdb66632018e90b14f20af2aaf53d6a8521dd8eb8441")
        self.assertEqual(len(m.CENSUS_HEADER),8)
        self.assertEqual(sum(x[0] for x in m.OLD_COUNTS.values()),929138)
        self.assertEqual(sum(x[1] for x in m.OLD_COUNTS.values()),416915)
        self.assertEqual(929138+76753+3141,1009032)
        self.assertEqual(416915+27689+1201,445805)
        self.assertEqual(m.K1_PARQUET_SHA,
            "0d3de3374a06b1e7b9ad0281e2ab7844632ed139db871488524061191c8f34c1")
        self.assertEqual(m.FCP_PARQUET_SHA,
            "ee79c60671dda65a3d2f4f9319af952dbe81ab805ce269fe32b3bb3448b0d6b8")

    def test_rare_taxon_below100_must_not_be_false_FCP_eligible(self):
        self.assertFalse(m.source_taxon_eligible(0,2,"false"))
        self.assertFalse(m.source_taxon_eligible(1,2,"false"))
        self.assertFalse(m.source_taxon_eligible(99,300,"false"))
        self.assertTrue(m.source_taxon_eligible(100,101,"true"))
        self.assertTrue(m.source_taxon_eligible(999999,2000000,"true"))
        with self.assertRaisesRegex(ValueError,"eligibility mutated"):
            m.source_taxon_eligible(99,99,"true")
        with self.assertRaisesRegex(ValueError,"eligibility mutated"):
            m.source_taxon_eligible(100,100,"false")

    def test_no_photo_source_count_can_be_smaller_than_original_observation_count(self):
        for a,b in ((1,0),(100,99),(-1,0)):
            with self.subTest(a=a,b=b),self.assertRaisesRegex(ValueError,"eligibility mutated"):
                m.source_taxon_eligible(a,b,"false")

    def test_technical_unclassifiable_kept_in_photo_attempt_count(self):
        x=m.taxon_coverage_record("100","plant",15,20,"false",5,2)
        self.assertEqual(x["genuine_original_photo_attempts"],5)
        self.assertEqual(x["automatic_four_hue_ROIs"],2)
        self.assertEqual(x["acquired_or_evaluated_but_hue_unknown"],3)
        self.assertEqual(x["source_taxon_under100_observations"],"true")
        self.assertEqual(x["has_any_four_hue_ROI"],"true")

    def test_rare_taxon_without_genuine_photo_results_is_prioritized(self):
        x=m.taxon_coverage_record("101","plant",1,1,"false",0,0)
        self.assertEqual(x["has_any_photo_attempt"],"false")
        self.assertEqual(x["has_any_four_hue_ROI"],"false")
        self.assertEqual(x["acquired_or_evaluated_but_hue_unknown"],0)
        self.assertEqual(x["has_independent_original_observation_to_resample"],"true")
        noobs=m.taxon_coverage_record("102","plant",0,1,"false",0,0)
        self.assertEqual(noobs["has_independent_original_observation_to_resample"],"false")

    def test_no_unmeasured_image_can_be_assigned_automatic_hue(self):
        with self.assertRaisesRegex(ValueError,"denominators invalid"):
            m.taxon_coverage_record("101","plant",1,1,"false",0,1)
        with self.assertRaisesRegex(ValueError,"denominators invalid"):
            m.taxon_coverage_record("101","plant",1,1,"false",1,2)
        with self.assertRaisesRegex(ValueError,"denominators invalid"):
            m.taxon_coverage_record("101","plant",1,1,"false",1,-1)
        with self.assertRaisesRegex(ValueError,"denominators invalid"):
            m.taxon_coverage_record("not-a-numeric-id","plant",1,1,"false",1,0)

    def test_source_taxon_DIL_names_are_not_claimed_as_native_blooming_petals(self):
        row=m.taxon_coverage_record("8","plant",102,200,"true",2,0)
        self.assertEqual(row["source_taxon_under100_observations"],"false")
        self.assertEqual(row["acquired_or_evaluated_but_hue_unknown"],2)
        self.assertEqual(row["has_any_four_hue_ROI"],"false")


if __name__=="__main__":
    unittest.main()
