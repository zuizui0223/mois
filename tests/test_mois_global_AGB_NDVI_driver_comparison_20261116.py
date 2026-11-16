import unittest
import numpy as np
import pandas as pd
from scripts.analysis.mois_global_AGB_NDVI_driver_comparison_20261116 import (
    ENV_GROUPS,GROUPS,BIOMASS,NDVI,MOISTURE,
    columns,genus_folds,white_scores,SEED
)
class TestBiomassCandidateDrivers(unittest.TestCase):
    def test_biomass_is_independent_from_vegetation_greenness(self):
        self.assertIn("log_ESACCI_AGB_woody_Mg_ha_2020",GROUPS[BIOMASS])
        self.assertNotIn("MODIS_NDVI_2001_2021_mean",GROUPS[BIOMASS])
        self.assertEqual(set(GROUPS[BIOMASS])&set(GROUPS[NDVI]),set())
        self.assertIn("vpd_kPa_approx",GROUPS[MOISTURE])
    def test_group_ablation_preserves_other_correlated_groups(self):
        full=set(columns(ENV_GROUPS))
        no_biomass=set(columns([g for g in ENV_GROUPS if g!=BIOMASS]))
        self.assertEqual(full-no_biomass,set(GROUPS[BIOMASS]))
        self.assertTrue(set(GROUPS[MOISTURE]).issubset(no_biomass))
    def test_3fold_genus_disjoint_9photo_colour(self):
        colors=("white","yellow","orange","red","pink","magenta","purple","blue","bronze")
        rows=[]
        for gi in range(360):
            for j in range(3):
                color=colors[(gi*5+j)%len(colors)]
                rows.append({"genus":f"g{gi}","source_hue9":color,
                             "white_photo":int(color=="white")})
        d=pd.DataFrame(rows)
        splits=genus_folds(d)
        self.assertEqual(len(splits),3)
        self.assertEqual(sum(len(te) for _,te in splits),len(d))
        for tr,te in splits:
            self.assertTrue(set(d.genus.iloc[tr]).isdisjoint(set(d.genus.iloc[te])))
    def test_AUC_not_probably_biomass_causality(self):
        y=np.array([0,1,0,1])
        p=np.array([0.1,0.9,0.3,0.8])
        self.assertEqual(white_scores(y,p)["AUC"],1.0)
if __name__=="__main__":unittest.main()
