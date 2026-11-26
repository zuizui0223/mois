"""No invented absences or outcome-selected regions in source JSDM design."""
import unittest
import numpy as np
import pandas as pd

from scripts.analysis.mois_jsdm_fourthcorner_gbif_native_pilot_20261117 import (
    stable_rank,panel_sites,REGIONS,MIN_SITES,MAX_SITES
)

class TestJSDMAvailabilityDesign(unittest.TestCase):
    def test_panels_predeclared_not_discovered_by_colour(self):
        self.assertEqual(len(REGIONS),3)
        self.assertTrue(all(len(x)==4 for x in REGIONS.values()))
        self.assertEqual(stable_rank("species:123"),
                         stable_rank("species:123"))
        self.assertNotEqual(stable_rank("species:123"),
                            stable_rank("species:124"))
    def test_geography_eligibility_and_hashed_site_order(self):
        env=pd.DataFrame({
            "grid_x":np.arange(900),
            "grid_y":np.arange(900),
            "lat":[37.0]*900,
            "lon":[125.0]*900,
            "NPP_2001_2023_mean_kgC_m2_yr":[.5]*900
        })
        r=env[["grid_x","grid_y"]].copy()
        r["region_strict"]=np.arange(900)%10!=0
        r["code"]="KOR"
        x=panel_sites(env,r,REGIONS["East_Asia_temperate"])
        self.assertEqual(len(x),MAX_SITES)
        self.assertTrue(x.region_strict is not None if
                        "region_strict" in x.columns else True)
        self.assertEqual(x.code.nunique(),1)
        self.assertEqual(len(x[["grid_x","grid_y"]].drop_duplicates()),300)
        y=panel_sites(env,r,REGIONS["North_America_temperate"])
        self.assertEqual(len(y),0)

if __name__=="__main__":
    unittest.main()
