"""No guessed flowers in 128-shard bulk data-only receipt; source arm remains explicit."""
import unittest
import numpy as np
import pandas as pd
from scripts.analysis.mois_bulk128_collect_original_roi9 import classify_nine,summary,PALETTE
from scripts.acquisition.mois_bulk128_tier10_source import GAP_ARM,REFERENCE_ARM,NO_COLOUR_PRIORITY,P5

class Bulk128PhotographicOutcomeTests(unittest.TestCase):
    def test_protect_technical_unknown_and_two_estimator_designs(self):
        d=pd.DataFrame({
            "photo_id":["1","2","3","4"],
            "inat_taxon_id":["101","202","303","404"],
            "observation_id":["o1","o2","o3","o4"],
            "source_arm":[GAP_ARM,GAP_ARM,REFERENCE_ARM,REFERENCE_ARM],
            "acquisition_design_cohort":["targeted","targeted","reference","reference"],
            "original_source_priority":[NO_COLOUR_PRIORITY[1],
                                        NO_COLOUR_PRIORITY[1],P5,P5],
            "latitude_region":["S_tropical","N_tropical","N_temperate","N_temperate"],
            "recorded_local_season":["summer","winter","spring","summer"],
            "morph":["white","mixed_uncertain","red_pink","mixed_uncertain"],
            "measurement_status":["classified_four_state_morph",
                                  "not_evaluable_ambiguous_palette_composition",
                                  "classified_four_state_morph",
                                  "image_acquisition_failed"],
        })
        for hue in PALETTE:d["flower_fraction_"+hue]=0.0
        d.loc[0,"flower_fraction_white"]=.95
        d.loc[0,"flower_fraction_yellow"]=.05
        d.loc[1,"flower_fraction_bronze"]=1.0
        d.loc[2,"flower_fraction_pink"]=1.0
        for hue in PALETTE:d.loc[3,"flower_fraction_"+hue]=np.nan
        output=classify_nine(d)
        self.assertEqual(output.strict9.tolist(),[True,False,True,False])
        self.assertEqual(int(output.original_nine_component_available_but_not_strict.sum()),1)
        self.assertTrue(output.loc[~output.strict9,"strict9_white"].isna().all())
        report=summary(output,["acquisition_design_cohort","original_source_priority"])
        targeted=[v for v in report if v["acquisition_design_cohort"]=="targeted"]
        reference=[v for v in report if v["acquisition_design_cohort"]=="reference"]
        self.assertEqual(len(targeted),1)
        self.assertEqual(targeted[0]["attempted_original_photos"],2)
        self.assertEqual(targeted[0]["strict9_newly_classified_zero_colour_taxa"],1)
        self.assertEqual(targeted[0]["technical_or_ambiguous_unknown"],1)
        self.assertEqual(reference[0]["photo_acquisition_failed"],1)
        self.assertEqual(reference[0]["strict9_newly_classified_zero_colour_taxa"],0)

    def test_technical_unclassifiable_image_must_not_gain_a_pigment_trait(self):
        d=pd.DataFrame({
            "morph":["mixed_uncertain"],
            "measurement_status":["image_acquisition_failed"],
        })
        for hue in PALETTE:d["flower_fraction_"+hue]=.1
        out=classify_nine(d)
        self.assertFalse(out.loc[0,"strict9"])
        self.assertTrue(pd.isna(out.loc[0,"strict9_white"]))

if __name__=="__main__":
    unittest.main()
