"""Synthetic fixtures only: not ecological observations."""
import unittest
import pandas as pd

from scripts.analysis.mois_p0008_real_photo_coverage import summarize, latitude_band

HUES={"white","yellow_orange","red_pink","blue_purple"}

class TestFullP0008RealPhotoCoverage(unittest.TestCase):
    def setUp(self):
        self.frame=pd.DataFrame([
            dict(photo_id="1",observation_id="A",inat_taxon_id="31",observed_on="2024-06-01",
                latitude="42",longitude="140",morph="white",measurement_status="classified_four_state_morph",
                classified_four_hue="true",source_arm="global93m_onepass_target_tier_02",selection_probability="0.0008"),
            dict(photo_id="2",observation_id="A",inat_taxon_id="31",observed_on="2024-06-01",
                latitude="42",longitude="140",morph="mixed_uncertain",measurement_status="unclassifiable",
                classified_four_hue="false",source_arm="global93m_onepass_target_tier_02",selection_probability="0.0008"),
            dict(photo_id="3",observation_id="B",inat_taxon_id="32",observed_on="2023-11-05",
                latitude="-15",longitude="70",morph="blue_purple",measurement_status="classified_four_state_morph",
                classified_four_hue="true",source_arm="global93m_onepass_target_tier_02",selection_probability="0.0008"),
            dict(photo_id="4",observation_id="C",inat_taxon_id="33",observed_on="2025-02-15",
                latitude="-40",longitude="150",morph="mixed_uncertain",measurement_status="image_acquisition_failed",
                classified_four_hue="false",source_arm="global93m_onepass_target_tier_02",selection_probability="0.0008"),
        ])

    def test_actual_terminal_outcomes_and_observation_dependence(self):
        r=summarize(self.frame)
        self.assertEqual(r["overall"]["photo_attempts"],4)
        self.assertEqual(r["overall"]["classified_four_hue"],2)
        self.assertEqual(r["overall"]["acquisition_failed"],1)
        self.assertEqual(r["overall"]["acquired_unclassifiable"],1)
        self.assertEqual(r["distinct_observations"],3)
        self.assertEqual(r["multi_photo_observations"],1)
        self.assertEqual(r["taxon_ids_with_classified_photos"],2)
        self.assertEqual(sum(x["photo_attempts"] for x in r["latitude_bands"]),4)
        self.assertEqual(sum(x["photo_attempts"] for x in r["months"]),4)
        self.assertEqual(sum(x["photo_attempts"] for x in r["latitude_by_month"]),4)
        self.assertFalse(r["globally_representative_white_fraction_claimed"])

    def test_tropical_latitude_bounds(self):
        for v in ("-23.44","23.44","0"):
            self.assertEqual(latitude_band(v),"tropical")
        self.assertEqual(latitude_band("-23.441"),"south_extratropical")
        self.assertEqual(latitude_band("23.441"),"north_extratropical")

    def test_duplicate_original_image_id_rejected(self):
        a=self.frame.copy();a.loc[3,"photo_id"]="1"
        with self.assertRaisesRegex(ValueError,"photo ID"):
            summarize(a)

    def test_invalid_hue_or_false_classification_rejected(self):
        a=self.frame.copy();a.loc[1,"classified_four_hue"]="true"
        with self.assertRaisesRegex(ValueError,"classification"):
            summarize(a)
        a=self.frame.copy();a.loc[3,"morph"]="white"
        with self.assertRaisesRegex(ValueError,"classification"):
            summarize(a)

    def test_geographic_enrichment_must_not_enter_original_probability_arm(self):
        a=self.frame.copy();a.loc[2,"source_arm"]="global93m_geographic_enrichment_tier_05"
        with self.assertRaisesRegex(ValueError,"source arm"):
            summarize(a)

    def test_original_observation_selection_probability_must_be_fixed(self):
        a=self.frame.copy();a.loc[2,"selection_probability"]="0.0032"
        with self.assertRaisesRegex(ValueError,"sampling probability"):
            summarize(a)

    def test_unknown_is_not_a_rose_or_a_white_flower(self):
        r=summarize(self.frame)
        self.assertEqual(r["overall"]["white_classified"],1)
        self.assertAlmostEqual(r["overall"]["white_fraction_conditional_classified"],0.5)

if __name__=="__main__":
    unittest.main()
