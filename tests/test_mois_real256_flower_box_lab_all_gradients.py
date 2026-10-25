"""Source independent circular hue and no-feature/pilot polar QA."""
import unittest
import numpy as np
import pandas as pd
from scripts.analysis.mois_real256_flower_box_lab_all_gradients import (
    circular_stats,med,show
)

class PilotBoxLabTests(unittest.TestCase):
    def test_hue_near_zero_crosses_360_without_arithmetic_error(self):
        h=circular_stats([359,1],weights=[30,30])
        self.assertTrue(h['hue_circular_mean_degrees']<3
                        or h['hue_circular_mean_degrees']>357)
        self.assertGreater(h['hue_resultant_length'],.99)
    def test_opposite_hue_vectors_have_undefined_mean(self):
        h=circular_stats([90,270])
        self.assertIsNone(h['hue_circular_mean_degrees'])
        self.assertLess(h['hue_resultant_length'],.00001)
    def test_zero_chroma_hue_remains_missing(self):
        h=circular_stats([np.nan,np.nan],weights=[0,0])
        self.assertIsNone(h['hue_circular_mean_degrees'])
    def test_empty_photo_pilot_reporting_retains_n_zero(self):
        frame=pd.DataFrame({"screen_status":[],
         "floral_box_main_L_star":[],"floral_box_main_C_star":[]})
        a=show("S_polar_66.56_90","latitude_zone",frame)
        self.assertEqual(a["signed_original_K1_photo_ID_attempts_in_pilot"],0)
        self.assertEqual(a["actual_Lab_original_pixel_assays"],0)
        self.assertIsNone(a["median_main_L_star"])
        self.assertFalse(a["N_10_for_minimum_descriptive_support"])

if __name__=="__main__":
    unittest.main()
