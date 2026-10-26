"""Real original ROI9 dominant-colour-as-primary-proxy, preserving all unknowns."""
import unittest
import numpy as np
import pandas as pd

from scripts.analysis.mois_main_petal_colour_proxy_from_roi9 import (
    assign_main_colour,COLUMNS
)

class PhotoMajorFlowerColourTests(unittest.TestCase):
    def fixture(self):
        d=pd.DataFrame({
            "photo_id":["100","101","102","103"],
            "source_arm_original":["A","A","A","A"],
            "original_four_hue":["red_pink","white","mixed_uncertain","mixed_uncertain"],
            "continuous9_strict":[True,True,False,False],
            "continuous9_exploratory_ambiguous":[False,False,True,False],
        })
        for col in COLUMNS:d[col]=np.nan
        d.loc[0,list(COLUMNS)]=[.51,0,0,0,.30,.19,0,0,0]
        d.loc[1,list(COLUMNS)]=[.70,0,0,0,0,.30,0,0,0]
        d.loc[2,list(COLUMNS)]=[.20,0,0,.60,.20,0,0,0,0]
        return d

    def test_dominant_color_does_not_depend_on_petals_or_4_hue_gates(self):
        r=assign_main_colour(self.fixture())
        self.assertEqual(len(r),4)
        self.assertEqual(r.main_petal_colour_proxy.tolist(),
                         ["white","white","red",""])
        self.assertEqual(r.original_four_hue.iloc[0],"red_pink")
        self.assertEqual(r.main_colour_dominance_tag.tolist(),
                         ["mixed_or_low_dominance_below_0.70",
                          "dominant_fraction_ge_0.70",
                          "mixed_or_low_dominance_below_0.70",
                          "not_measured"])
        self.assertAlmostEqual(r.main_colour_fraction.iloc[0],.51)
        self.assertEqual(r.second_colour_component.iloc[0],"pink")
        self.assertAlmostEqual(r.second_colour_fraction.iloc[0],.30)
        self.assertEqual(r.main_colour_provenance.iloc[2],
                         "ambiguous_ROI9_sensitivity_only")
        self.assertTrue(r.ambiguous_main_colour_sensitivity_only.iloc[2])
        self.assertFalse(r.strict_main_colour_primary.iloc[2])
        self.assertFalse(r.main_colour_is_a_botanically_verified_petal.any())
        self.assertEqual(r.main_petal_colour_proxy.iloc[3],"")
        self.assertTrue(np.isnan(r.main_colour_fraction.iloc[3]))

    def test_equal_major_fraction_fixed_order_and_no_minimum_purity_exclusion(self):
        d=self.fixture().iloc[:1].copy()
        d.loc[0,list(COLUMNS)]=[.50,0,0,0,.50,0,0,0,0]
        r=assign_main_colour(d)
        self.assertEqual(r.main_petal_colour_proxy.iloc[0],"white")
        self.assertEqual(r.second_colour_component.iloc[0],"pink")
        self.assertAlmostEqual(r.dominant_minus_second_fraction.iloc[0],0)
        self.assertEqual(r.main_colour_dominance_tag.iloc[0],
                         "mixed_or_low_dominance_below_0.70")

    def test_strict_and_ambiguous_cannot_overlap(self):
        d=self.fixture()
        d.loc[0,"continuous9_exploratory_ambiguous"]=True
        with self.assertRaisesRegex(ValueError,"overlap"):
            assign_main_colour(d)

    def test_technical_unknown_cannot_have_imputed_major_colour(self):
        d=self.fixture()
        d.loc[3,"continuous9_white"]=.9
        with self.assertRaisesRegex(ValueError,"unaccepted"):
            assign_main_colour(d)

if __name__=="__main__":
    unittest.main()
