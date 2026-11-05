"""Identifiability and same-photo input contracts: no synthetic science claimed."""
import unittest
import pandas as pd
from scripts.analysis import mois_latent_error_identifiability as h

class TestHWSourceContract(unittest.TestCase):
    def test_two_tests_two_populations_standard_just_identified(self):
        a=h.identifiability(2,2,"shared")
        self.assertEqual((a["observable_degrees"],a["free_parameters"]),(6,6))
        self.assertEqual(a["count_status"],"JUST_IDENTIFIED_COUNT_ONLY")
        self.assertFalse(a["conditional_independence_verified_by_counts"])
    def test_two_tests_two_populations_environment_specific_underidentified(self):
        a=h.identifiability(2,2,"population_specific")
        self.assertEqual((a["observable_degrees"],a["free_parameters"]),(6,10))
        self.assertEqual(a["count_status"],"UNDERIDENTIFIED")
    def test_three_tests_two_populations_specific_requires_conditional_independence(self):
        a=h.identifiability(3,2,"population_specific")
        self.assertEqual((a["observable_degrees"],a["free_parameters"]),(14,14))
        self.assertEqual(a["count_status"],"JUST_IDENTIFIED_COUNT_ONLY")
    def test_two_tests_three_populations_differential_still_underidentified(self):
        a=h.identifiability(2,3,"population_specific")
        self.assertEqual((a["observable_degrees"],a["free_parameters"]),(9,15))
        self.assertEqual(a["count_status"],"UNDERIDENTIFIED")
    def test_bad_design_rejected(self):
        with self.assertRaises(ValueError):h.identifiability(1,2,"shared")
        with self.assertRaises(ValueError):h.identifiability(2,1,"shared")
        with self.assertRaises(ValueError):h.identifiability(2,2,"bogus")
    def sample(self):
        return pd.DataFrame([
            {"photo_id":str(i),"observation_id":"o"+str(i),"inat_taxon_id":"t"+str(i%5),
             "source_arm":"reference","environment_stratum":"dry" if i<4 else "wet",
             "test1_label":("white" if i%2 else "nonwhite"),
             "test2_label":("white" if i%3 else "nonwhite"),
             "test1_method":"original_FCP_ROI9_white","test2_method":"independent_VLM_white",
             "test1_model_sha256":"a"*64,"test2_model_sha256":"b"*64,
             "source_photo_sha256":"c"*64}
            for i in range(8)])
    def test_two_separate_binary_white_tests_generate_cell_counts(self):
        z=h.paired_count_audit(self.sample())
        self.assertEqual(z["same_photo_paired_complete"],8)
        self.assertEqual(sum(v["n_complete"] for v in z["strata"]),8)
        self.assertEqual(len(z["strata"]),2)
        self.assertFalse(z["conditional_independence_verified"])
        self.assertFalse(z["latent_white_prevalence_estimated"])
    def test_unknown_is_not_negative(self):
        d=self.sample()
        d.loc[0,"test2_label"]="unknown"
        z=h.paired_count_audit(d)
        self.assertEqual(z["same_photo_paired_complete"],7)
        self.assertEqual(z["unknown_abstentions"],1)
    def test_flower_detector_or_box_lab_is_not_white_test(self):
        d=self.sample()
        d["test2_method"]="YOLO_500px_BOX_Lab"
        with self.assertRaisesRegex(ValueError,"same botanical white"):
            h.paired_count_audit(d)
    def test_different_photos_fail_same_photo_identity(self):
        d=self.sample()
        d.loc[0,"photo_id"]="1"
        with self.assertRaisesRegex(ValueError,"duplicated"):
            h.paired_count_audit(d)
    def test_same_observation_two_photos_must_be_accounted_for(self):
        d=self.sample()
        d.loc[0,"observation_id"]="o1"
        with self.assertRaisesRegex(ValueError,"observation"):
            h.paired_count_audit(d)
    def test_no_provenance_for_second_classifier_fails_closed(self):
        d=self.sample()
        d.loc[0,"test2_model_sha256"]=""
        with self.assertRaisesRegex(ValueError,"provenance"):
            h.paired_count_audit(d)
    def test_no_forced_prevalence_difference(self):
        z=h.paired_count_audit(self.sample())
        self.assertEqual(z["strata_prevalence_distinct_proven"],False)
        self.assertEqual(z["observed_white_label_rates_can_prove_prevalence"],False)
    def test_test_dependence_and_open_flower_not_solved_by_many_photos(self):
        a=h.identifiability(2,2,"shared")
        self.assertFalse(a["conditional_independence_verified_by_counts"])
        self.assertFalse(a["flower_petal_truth_verified_by_counts"])

if __name__=="__main__":unittest.main()
