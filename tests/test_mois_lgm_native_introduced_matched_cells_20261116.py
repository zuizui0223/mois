import unittest
import pandas as pd
from scripts.analysis.mois_lgm_native_introduced_matched_cells_20261116 import (
    same_cell_cohort
)

class TestMatchedWCVPNative(unittest.TestCase):
    def test_same_grid_support_and_no_unmatched_join(self):
        local=pd.DataFrame({
            "tdwg3_assignment_mode":["strict","strict","moderate"],
            "grid_x":[1,2,3],"grid_y":[1,1,1],
            "native_species_n":[70,30,90],
            "introduced_species_n":[45,12,40],
            "native_white_n":[35,20,40],
            "introduced_white_n":[18,5,10],
            "GBIF_occurrence_count":[100,50,60],
            "photo_label_species_n":[125,42,130],
            "white_photo_label_species_n":[53,25,50]
        })
        climate=pd.DataFrame({
            "grid_x":[1,2,3],"grid_y":[1,1,1],
            "photo_label_species_n":[125,42,130],
            "white_photo_label_species_n":[53,25,50],
            "block_500km":["0|0","0|0","0|0"],
            "block_1000km":["0|0","0|0","0|0"],
            "present_temperature_C":[20,21,22],
            "present_precip_mm":[1000,1200,1300],
            "threeGCM_delta_temperature_C_median":[5,6,7],
            "threeGCM_delta_log_precip_median":[.1,.2,.3],
            "ESA_AGB_woody_Mg_ha_2020":[50,40,60],
            "lon":[10,20,30],"lat":[30,40,50]})
        a=same_cell_cohort(local,climate,"strict",30)
        self.assertEqual(len(a),1)
        self.assertEqual(a.grid_x.iloc[0],1)
        self.assertAlmostEqual(a.native_minus_intro_white.iloc[0],.5-.4)
        self.assertEqual(a.GBIF_occurrence_count.iloc[0],100)

if __name__=="__main__":
    unittest.main()
