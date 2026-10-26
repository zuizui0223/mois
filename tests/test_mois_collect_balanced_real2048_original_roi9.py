"""Synthetic source-sealed 2048 ROI9 output summaries preserve technical missingness."""
import unittest
import pandas as pd
import numpy as np
from scripts.analysis.mois_collect_balanced_real2048_original_roi9 import (
    summarize, SCHEMA, normalize_source_design_provenance,
)

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

    def test_actual_signed_precolour_design_cohort_header_survives_photo_id_join(self):
        """Real tier07 source CSV calls the field design_cohort, not sampling_design_cohort."""
        source = pd.DataFrame({
            "photo_id": ["102", "101"],
            "original_observation_id": ["b", "a"],
            "inat_taxon_id": ["22", "11"],
            "design_cohort": ["targeted", "reference"],
        })
        validated = normalize_source_design_provenance(source)
        self.assertNotIn("sampling_design_cohort", source.columns)
        self.assertNotIn("design_cohort", validated.columns)
        self.assertEqual(validated["sampling_design_cohort"].tolist(),
                         ["targeted", "reference"])
        actual_pixel_rows = pd.DataFrame({
            "photo_id": ["101", "102"],
            "observation_id": ["a", "b"],
            "inat_taxon_id": ["11", "22"],
        })
        joined = actual_pixel_rows.merge(validated, on="photo_id",
                                         how="left", validate="one_to_one")
        self.assertEqual(joined["sampling_design_cohort"].tolist(),
                         ["reference", "targeted"])
        self.assertEqual(joined["taxon_id_provenance"].tolist(), ["11", "22"])
        self.assertEqual(joined["original_observation_id"].tolist(), ["a", "b"])

    def test_missing_or_ambiguous_signed_source_design_is_not_silently_guessed(self):
        input_rows = pd.DataFrame({
            "photo_id": ["100"],
            "inat_taxon_id": ["12"],
            "design_cohort": ["targeted"],
        })
        with self.assertRaisesRegex(ValueError, "header changed"):
            normalize_source_design_provenance(input_rows.drop(columns=["design_cohort"]))
        with self.assertRaisesRegex(ValueError, "unrecognized or missing"):
            normalize_source_design_provenance(
                input_rows.assign(design_cohort="unknown")
            )
        with self.assertRaisesRegex(ValueError, "header changed"):
            normalize_source_design_provenance(
                input_rows.assign(sampling_design_cohort="reference")
            )

if __name__=="__main__":
    unittest.main()
