"""Original signed 2048 regional photo acquisitions must improve source-taxon coverage without imputing unknowns."""
import unittest
import pandas as pd
from scripts.analysis.mois_audit2048_taxon_gap_gain_after400710 import apply_gain,priority

class MeasuredRealPhotoTaxonGapTests(unittest.TestCase):
    def frame(self):
        return pd.DataFrame([
            ("1",0,0,0,priority(0,0)),
            ("2",1,0,0,priority(1,0)),
            ("3",100,1,0,priority(100,0)),
            ("4",100,4,4,priority(100,4)),
            ("5",100,9,9,priority(100,9)),
            ("6",1000,10,10,priority(1000,10)),
        ],columns=["inat_taxon_id","n_unique_source_photo_ids",
                   "n_image_photo_attempts","n_strict_roi9","priority"])
    def newphotos(self):
        return pd.DataFrame({
            "photo_id":["101","102","103","104","105","106"],
            "observation_id":["obs1","obs2","obs3","obs4","obs5","obs6"],
            "inat_taxon_id":["2","3","4","5","6","888"],
            "continuous9_strict":[True,True,True,True,False,True],
            "sampling_design_cohort":["targeted","targeted","reference",
                                      "targeted","reference","targeted"],
            "latitude_region":["N_tropical"]*6,
            "recorded_local_season":["winter"]*6,
            "original_four_hue":["white","yellow_orange","red_pink","blue_purple",
                                 "mixed_uncertain","white"],
            "image_acquisition_failed":[False]*6,
        })

    def test_source_frame_priority_transitions_and_outside_frame_are_separate(self):
        before=self.frame()
        sample=self.newphotos()
        after,outside,_=apply_gain(before,sample)
        indexed=after.set_index("inat_taxon_id")
        self.assertEqual(indexed.loc["2","priority_after"],"P3_thin_colour_1to4")
        self.assertEqual(indexed.loc["3","priority_after"],"P3_thin_colour_1to4")
        self.assertEqual(indexed.loc["4","priority_after"],"P4_thin_colour_5to9")
        self.assertEqual(indexed.loc["5","priority_after"],"P5_existing_10plus_strict9")
        self.assertEqual(indexed.loc["6","priority_after"],"P5_existing_10plus_strict9")
        self.assertEqual(indexed.loc["1","priority_after"],"P0_no_original_source_photo")
        self.assertEqual(int(after.newly_colour_supported.sum()),2)
        self.assertEqual(int(after.newly_reach_5.sum()),1)
        self.assertEqual(int(after.newly_reach_10.sum()),1)
        self.assertEqual(int(outside.new_strict9.sum()),1)
        self.assertEqual(int(outside.new_photo_attempts.sum()),1)
        self.assertEqual(int(before.n_strict_roi9.sum()),23)
        self.assertEqual(int(after.n_strict_roi9.sum()),23)
        self.assertEqual(int(after.n_strict9_after.sum()),27)

    def test_reject_pseudoreplication_and_unmeasured_colours(self):
        s=self.newphotos()
        duplicated=s.copy()
        duplicated.loc[1,"observation_id"]=duplicated.loc[0,"observation_id"]
        with self.assertRaisesRegex(ValueError,"duplicated"):
            apply_gain(self.frame(),duplicated)
        fake=s.copy()
        fake.loc[4,"continuous9_strict"]=True
        with self.assertRaisesRegex(ValueError,"must never gain"):
            apply_gain(self.frame(),fake)
        invalid=s.copy()
        invalid.loc[0,"image_acquisition_failed"]=True
        with self.assertRaisesRegex(ValueError,"Unacquired"):
            apply_gain(self.frame(),invalid)

    def test_bad_frame_and_missing_probability_role_fail_closed(self):
        s=self.newphotos()
        s.loc[0,"sampling_design_cohort"]="unknown"
        with self.assertRaisesRegex(ValueError,"cohort"):
            apply_gain(self.frame(),s)
        s=self.newphotos()
        invalid=self.frame()
        invalid.loc[1,"n_unique_source_photo_ids"]=0
        with self.assertRaisesRegex(ValueError,"zero upstream"):
            apply_gain(invalid,s)

if __name__=="__main__":
    unittest.main()
