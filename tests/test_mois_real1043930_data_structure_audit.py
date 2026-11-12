"""Synthetic MOIS source-by-taxon concentration, observation cluster and ROI accounting."""
import unittest
import numpy as np
import pandas as pd
from scripts.analysis import mois_real1043930_data_structure_audit as m

class PhotoDataStructureTests(unittest.TestCase):
    def test_photo_population_concentration_not_equal_photographic_count(self):
        self.assertAlmostEqual(m.top_share([5,5,5,20],1),100*20/35)
        self.assertAlmostEqual(m.top_share([5,5,5,20],2),100*25/35)
        self.assertAlmostEqual(m.concentration([10,10,10,10]),4)
        self.assertAlmostEqual(m.concentration([90,10]),1/(.9**2+.1**2))
        self.assertIsNone(m.concentration([0,0]))
        self.assertEqual(m.finite_pct(3,4),75.)
        self.assertIsNone(m.finite_pct(1,0))

    def test_source_group_distinguishes_unclassifiable_and_ambiguous(self):
        f=pd.DataFrame({
          "collection_source":["FCP","FCP","OTHER","OTHER"],
          "source_arm_original":["historical_main"]*4,
          "inat_taxon_id":["13","14","13","14"],
          "genus":["Campanula","Cirsium","Campanula","Cirsium"],
          "observation_id":["obs1","obs1","obs2",""],
          "original_four_hue":["white","mixed_uncertain","red_pink","mixed_uncertain"],
          "continuous9_strict":[True,False,True,False],
          "continuous9_exploratory_ambiguous":[False,True,False,False]})
        z=m.source_rows(f,"collection_source")
        self.assertEqual(len(z),2)
        a=next(x for x in z if x["source_group"]=="FCP")
        self.assertEqual(a["photo_attempts"],2)
        self.assertEqual(a["strict_nine_component_success"],1)
        self.assertEqual(a["ambiguous_nine_component_sensitivity"],1)
        self.assertEqual(a["strict9_colour_missing"],1)
        self.assertEqual(a["distinct_original_observations"],1)
        self.assertEqual(a["median_photos_per_distinct_observation"],2.)
        b=next(x for x in z if x["source_group"]=="OTHER")
        self.assertEqual(b["source_observation_ID_missing"],1)
        self.assertEqual(b["strict9_colour_missing"],1)

    def test_no_unverified_white_when_source_photo_ROI_unknown(self):
        d=pd.DataFrame({"strict":[True,False],"white_fraction":[.7,np.nan]})
        self.assertTrue(d.loc[~d.strict,"white_fraction"].isna().all())

if __name__=="__main__":
    unittest.main()
