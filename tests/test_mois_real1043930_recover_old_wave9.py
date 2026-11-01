"""Old original ROI9 recovery refuses 4-colour-only guesses and duplicate photo IDs."""
import unittest
import numpy as np
import pandas as pd
from scripts.analysis.mois_real1043930_recover_old_wave9 import (
    candidate_patch,FRACTIONS
)

class GenuineExistingPhotoFractionTests(unittest.TestCase):
    def atlas(self):
        d=pd.DataFrame({
           "photo_id":["100","101","102","103"],
           "source_arm_original":["historical_main"]*4,
           "original_four_hue":["white","blue_purple","red_pink","mixed_uncertain"],
           "continuous9_strict":[False,False,True,False],
        })
        return d.set_index("photo_id")

    def source(self):
        d=pd.DataFrame({
           "photo_id":["100","101","102","103"],
           "morph":["white","white","red_pink","mixed_uncertain"],
           "measurement_status":["classified_four_state_morph",
                 "classified_four_state_morph",
                 "classified_four_state_morph","not_evaluable_roi_or_flip_gate"]
        })
        for f in FRACTIONS:d[f]=0.
        d.loc[0,FRACTIONS[0]]=.65
        d.loc[0,FRACTIONS[1]]=.35
        d.loc[1,FRACTIONS[0]]=1.
        d.loc[2,FRACTIONS[3]]=1.
        return d

    def test_only_matching_original_id_and_frozen_hue_recovers(self):
        a,selfaudit=candidate_patch(self.atlas(),self.source(),"wave")
        self.assertEqual(a.photo_id.tolist(),["100"])
        self.assertEqual(selfaudit["source_label_agreements_with_latest"],2)
        self.assertEqual(selfaudit["strict_new_9vectors_photoID_hue_and_source_verified"],1)
        self.assertAlmostEqual(float(a.iloc[0][FRACTIONS[0]]),.65)
        self.assertAlmostEqual(float(a.iloc[0][FRACTIONS[1]]),.35)

    def test_no_invention_from_legacy_four_class_when_raw_9_is_missing(self):
        d=self.source()
        d.loc[0,list(FRACTIONS)]=np.nan
        a,_=candidate_patch(self.atlas(),d,"wave")
        self.assertTrue(a.empty)

    def test_original_photo_wrong_acquisition_arm_is_fatal(self):
        d=self.atlas().copy()
        d.loc["100","source_arm_original"]="K1"
        with self.assertRaisesRegex(ValueError,"wrong source arm"):
            candidate_patch(d,self.source(),"wave")

    def test_original_classifier_status_must_match_original_categorical_hue(self):
        d=self.source()
        d.loc[0,"measurement_status"]="not_evaluable_roi_or_flip_gate"
        with self.assertRaisesRegex(ValueError,"terminal status/hue"):
            candidate_patch(self.atlas(),d,"wave")

if __name__=="__main__":
    unittest.main()
