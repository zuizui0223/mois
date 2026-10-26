"""Synthetic source-sealed 2048 ROI9 output summaries preserve technical missingness."""
import unittest
import pandas as pd
import numpy as np
from scripts.analysis.mois_collect_balanced_real2048_original_roi9 import summarize,SCHEMA

class MeasurementSeparateBiasCohortsTests(unittest.TestCase):
    def test_reference_vs_targeted_seasonal_roi_failure_remains_unknown(self):
        d=pd.DataFrame({
            "sampling_design_cohort":["reference","reference","targeted","targeted"],
            "latitude_region":["N_temperate","N_temperate","S_tropical","S_tropical"],
            "recorded_local_season":["summer","summer","winter","winter"],
            "photo_id":["101","102","103","104"],
            "observation_id":["o1","o2","o3","o4"],
            "inat_taxon_id":["1","2","3","4"],
            "continuous9_strict":[True,False,True,False],
            "continuous9_white":[.85,np.nan,.10,np.nan],
            "original_four_hue":["white","mixed_uncertain","red_pink","mixed_uncertain"],
            "image_acquisition_failed":[False,False,False,True]
        })
        groups=summarize(d)
        self.assertEqual(len(groups),2)
        north=next(v for v in groups if v["design_cohort"]=="reference")
        south=next(v for v in groups if v["design_cohort"]=="targeted")
        self.assertEqual(north["genuinely_attempted_original_photos"],2)
        self.assertEqual(north["new_authentic_strict9_ROIs"],1)
        self.assertAlmostEqual(north["strict_roi9_mean_white_fraction"],.85)
        self.assertEqual(north["new_nonclassifiable_image_Rois"],1)
        self.assertEqual(south["new_image_acquisition_failures"],1)
        self.assertEqual(south["new_authentic_strict9_ROIs"],1)
        self.assertEqual(south["distinct_photo_observations"],2)
        self.assertTrue(south["all_photo_attempts_kept_even_when_ROI_fails"])
        self.assertIn("actual2048",SCHEMA)

if __name__=="__main__":
    unittest.main()
