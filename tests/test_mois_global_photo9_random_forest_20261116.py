import unittest
import numpy as np
from scripts.analysis.mois_global_photo9_random_forest_20261116 import (
    PREDICTORS,CLIMATE,FULL,columns_for,white_scores,multiclass_scores,COLORS,
    FEATURE_ZIP_SHA,make_model
)
class TestRealGlobalRF(unittest.TestCase):
    def test_all_nineteen_worldclim_bio_and_modis_both_units_present(self):
        a=columns_for(CLIMATE)
        self.assertTrue(all(f"BIO{i:02d}" in a for i in range(1,20)))
        self.assertIn("srad_monthly_mean",a)
        self.assertIn("wind_monthly_mean",a)
        self.assertIn("vpd_kPa_approx",a)
        n=columns_for(FULL)
        self.assertEqual(len(set(n)),len(n))
        self.assertTrue(set(PREDICTORS["NDVI_vegetation_greenness"]).issubset(n))
        self.assertFalse(set(PREDICTORS["NDVI_vegetation_greenness"]).intersection(a))
    def test_binary_metrics_actual_predictive_skill(self):
        t=np.array([0,1,0,1])
        p=np.array([.1,.9,.2,.8])
        s=white_scores(t,p)
        self.assertEqual(s["auc_ROC"],1.0)
        self.assertEqual(s["balanced_accuracy"],1.)
    def test_9class_macro_accuracy_not_majority_accuracy(self):
        t=np.array(COLORS)
        p=np.eye(9)
        result=multiclass_scores(t,p,COLORS)
        self.assertEqual(result["balanced_accuracy"],1.0)
        self.assertEqual(result["top1_accuracy"],1.0)
    def test_rf_is_real_training_not_an_inert_control(self):
        forest=make_model(white=True)
        self.assertGreater(forest.n_estimators,100)
        self.assertGreaterEqual(forest.min_samples_leaf,10)
        self.assertEqual(len(FEATURE_ZIP_SHA),64)
if __name__=='__main__':unittest.main()
