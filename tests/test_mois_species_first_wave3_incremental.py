"""Synthetic tests: wave3 adds only genuinely new source taxa, no false fixation."""
import unittest
import pandas as pd
from scripts.analysis.mois_species_first_wave3_incremental import add_to_catalogue


class ThirdWaveIncrementalTests(unittest.TestCase):
    def test_colour_success_new_taxon_and_unknown_stays_missing(self):
        palette=("white","yellow","orange","red","pink","magenta","purple","blue","bronze")
        rows=[]
        for x in range(42111):
            r={"inat_taxon_id":str(x+1),"hue_status":"no_strict9",
               "dominant_hue9":"","observed_hue9_set":"",
               "after_strict9":0,"after_photo_attempts":0,
               "independent_roi9_observations":0,"distinct_palette_votes":0,
               "ambiguous_observation_rois":0,"next_measure_priority":0}
            r.update({"observed_mean_"+c:None for c in palette})
            rows.append(r)
        src=pd.DataFrame(rows)
        a=pd.DataFrame([{"inat_taxon_id":"1","photo_id":"100",
                         "observation_id":"one","strict9":True,
                         "dominant_label9":"purple"},
                        {"inat_taxon_id":"2","photo_id":"200",
                         "observation_id":"two","strict9":True,
                         "dominant_label9":""},
                        {"inat_taxon_id":"3","photo_id":"300",
                         "observation_id":"three","strict9":False,
                         "dominant_label9":""}])
        for c in palette:
            a["strict9_"+c]=[1.0 if c=="purple" else 0.0,
                             1.0/9 if c in palette else 0, float("nan")]
        revised=add_to_catalogue(src,a,expected_photos=3).set_index("inat_taxon_id")
        self.assertEqual(revised.loc["1","hue_status"],"provisional_single")
        self.assertEqual(revised.loc["1","dominant_hue9"],"purple")
        self.assertEqual(revised.loc["2","hue_status"],"mixed_image_unresolved")
        self.assertEqual(revised.loc["3","hue_status"],"no_strict9")
        self.assertEqual(revised.loc["3","after_strict9"],0)
        self.assertEqual(int(revised.after_strict9.sum()),2)
        self.assertEqual(int(revised.after_photo_attempts.sum()),3)
        self.assertTrue(bool(revised.loc["1","routine_hue_assigned"]))
        self.assertFalse(bool(revised.loc["3","routine_hue_assigned"]))

    def test_does_not_relabel_previously_colour_supported(self):
        a=pd.DataFrame([{"inat_taxon_id":str(i),
                         "hue_status":"no_strict9" if i!=1 else "provisional_single",
                         "dominant_hue9":"" if i!=1 else "white",
                         "after_photo_attempts":1,"after_strict9":int(i==1)}
                         for i in range(1,42112)])
        b=pd.DataFrame({"inat_taxon_id":["1"],"photo_id":["100"],
                        "observation_id":["u"],"strict9":[True],
                        "dominant_label9":["white"]})
        with self.assertRaisesRegex(ValueError,"already had"):
            add_to_catalogue(a,b,expected_photos=1)


    def test_csv_empty_hue_nan_remains_an_eligible_zero_colour_taxon(self):
        src=pd.DataFrame([
            {"inat_taxon_id":str(i+1),"hue_status":"no_strict9",
             "dominant_hue9":float("nan"),"observed_hue9_set":float("nan"),
             "after_photo_attempts":0,"after_strict9":0,
             "independent_roi9_observations":0,"distinct_palette_votes":0,
             "ambiguous_observation_rois":0,"next_measure_priority":0}
             for i in range(42111)])
        third=pd.DataFrame([{"inat_taxon_id":"1","photo_id":"100",
            "observation_id":"independent-1","strict9":False,"dominant_label9":""}])
        updated=add_to_catalogue(src,third,expected_photos=1)
        zeros=updated.loc[updated.after_strict9.eq(0)]
        self.assertEqual(len(zeros),42111)
        self.assertTrue(zeros.dominant_hue9.eq("").all())
        self.assertTrue(zeros.observed_hue9_set.eq("").all())


if __name__=="__main__":
    unittest.main()
