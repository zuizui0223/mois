"""Source-authentic nine solar×rain ROI failure sensitivity; no missing hues imputed."""
import unittest
import pandas as pd

from scripts.analysis.mois_geb_real313535_mnar_white_tipping import (
    calculate,validate,extract,mean_white_under_unknown
)

# Actual original photo/valid original four-hue/white and other visible
# colour counts copied from the signed nine-cell source release.
# The production workflow still requires original CSV SHA256.
SRC=[
    ("low","low",32271,17446,8260,5490,1348,2348),
    ("low","middle",22387,11691,5356,3757,1091,1487),
    ("low","high",50084,23551,9964,8346,2796,2445),
    ("middle","low",41783,21268,10570,6290,1851,2557),
    ("middle","middle",36466,18229,8574,5623,2065,1967),
    ("middle","high",25813,11695,4996,3822,1589,1288),
    ("high","low",30348,14822,7600,4083,1450,1689),
    ("high","middle",45556,20651,10561,5779,2398,1913),
    ("high","high",28496,11754,5487,3756,1601,910),
]
def frame():
    return pd.DataFrame([
      {
       "rainfall_tertile":rain,"solar_tertile":solar,
       "original_processed_photo_attempts":total,
       "real_original_classifiable_ROI_photos":known,
       "technically_unknown_hue_photo_attempts":total-known,
       "white_photo_ROIs":white,
       "yellow_orange_photo_ROIs":yellow,
       "red_pink_photo_ROIs":red,
       "blue_purple_photo_ROIs":blue,
      } for rain,solar,total,known,white,yellow,red,blue in SRC
    ])

class TrueRealSourceMNARTests(unittest.TestCase):
    def test_signed_source_counts_and_nonresponse_by_solar(self):
        rows=validate(frame())
        self.assertEqual(int(rows.original_processed_photo_attempts.sum()),313204)
        self.assertEqual(int(rows.real_original_classifiable_ROI_photos.sum()),151107)
        self.assertEqual(int(rows.white_photo_ROIs.sum()),71368)
        self.assertEqual(int(rows.technically_unknown_hue_photo_attempts.sum()),162097)

    def test_empirical_missing_label_tipping_points(self):
        strat,scenarios,weighted,summary=calculate(frame())
        actual=dict(zip(strat.solar_tertile,
                        strat.photographic_white_wet_minus_dry_completecase_pp))
        self.assertAlmostEqual(actual["low"],3.9290350347684244,places=7)
        self.assertAlmostEqual(actual["middle"],5.3273620498219545,places=7)
        self.assertAlmostEqual(actual["high"],4.373798359555902,places=7)
        drops=dict(zip(strat.solar_tertile,
                       strat.wet_missing_minus_wet_observed_probability_needed_to_tie_pp))
        self.assertAlmostEqual(drops["low"],-7.6799146744,places=5)
        self.assertAlmostEqual(drops["middle"],-9.7447623185,places=5)
        self.assertAlmostEqual(drops["high"],-7.4444963597,places=5)
        self.assertAlmostEqual(summary[
            "completecase_attempt_weighted_solar_stratified_wet_minus_dry_white_pp"],
            4.550412063034302,places=6)
        self.assertAlmostEqual(summary[
            "one_sided_assumed_white_fraction_depletion_of_wet_unmeasured_to_erase_contrast_pp"],
            8.25060612687695,places=6)
        for r in strat.itertuples():
            self.assertLess(r.unidentified_wet_minus_dry_low_all_wet_unknown_0_dry_unknown_1_pp,0)
            self.assertGreater(r.unidentified_wet_minus_dry_high_all_wet_unknown_1_dry_unknown_0_pp,0)
        subset=weighted[(weighted.counterfactual_MNAR_scenario==
            "symmetric_wet_depletion_dry_enrichment")&
            (weighted.assumed_unmeasured_delta_percentage_points==5)]
        self.assertEqual(len(subset),1)
        self.assertAlmostEqual(
            subset.iloc[0].photo_attempt_equal_solar_support_standardized_wet_minus_dry_pp,
            -0.6662467668018599,places=6)
        self.assertFalse(subset.iloc[0].direction_remains_positive_under_assumption)
        self.assertTrue(summary[
            "MNAR_scenario_results_do_not_correct_genus_adjusted_309439_worldclim_model"])
        self.assertTrue(summary["NO_unknown_original_photo_assigned_white_or_chromatic"])

    def test_impossible_source_rows_fail_without_creating_white_photos(self):
        sample=frame();sample.loc[0,"real_original_classifiable_ROI_photos"]+=1
        with self.assertRaisesRegex(ValueError,"equal measured"):
            validate(sample,verify_global=False)
        sample=frame();sample.loc[0,"white_photo_ROIs"]=999999
        with self.assertRaisesRegex(ValueError,"do not exhaust"):
            validate(sample,verify_global=False)
        cell=extract(frame().iloc[0])
        with self.assertRaisesRegex(ValueError,"must be in"):
            mean_white_under_unknown(cell,-.01)

if __name__=="__main__":
    unittest.main()
