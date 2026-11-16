"""Synthetic proof: genus-free floral hue regression has no taxonomic fields."""
import unittest
import numpy as np
import pandas as pd
from scripts.analysis.mois_genus_free_photo_white_NPP_LGM_20261117 import (
    design,estimate,sample_frame,NPP_NAME,NPP_COVER,LGM_TEMP,LGM_RAIN
)

class TestGenusFreeFlowerColour(unittest.TestCase):
    def source(self):
        rng=np.random.default_rng(20261117)
        n=600
        NPP=rng.uniform(.1,1.5,n)
        z=pd.DataFrame({
            "key":[str(x) for x in range(n)],
            "genus":np.repeat([f"T{i}" for i in range(60)],10),
            "source_hue9":np.where(NPP+rng.normal(0,.15,n)>.75,"white","purple"),
            "mean9_hue":np.where(NPP+rng.normal(0,.2,n)>.70,"white","yellow"),
            "hue_status":np.where(rng.random(n)<.5,
                               "sampling_consensus_single","provisional_single"),
            NPP_NAME:NPP,NPP_COVER:rng.uniform(.8,1,n),
            LGM_TEMP:rng.normal(4,2,n),
            LGM_RAIN:rng.normal(0,.1,n),
            "lgm_climate_coverage":rng.uniform(.9,1,n),
            "BIO01":rng.normal(16,5,n),
            "BIO12":rng.uniform(200,2400,n),
            "MODIS_NDVI_2001_2021_mean":rng.uniform(.1,.85,n),
            "ESA_AGB_woody_Mg_ha_2020":rng.uniform(0,200,n),
            "mean_abs_latitude_over_occupied_cells":rng.uniform(10,70,n),
            "records_passed_spatial_QA":rng.integers(10,8000,n),
            "occupied_grid_cells":rng.integers(1,200,n),
        })
        return z
    def test_no_genus_in_data_matrix_and_species_once(self):
        d=sample_frame(self.source())
        x,sd=design(d)
        self.assertTrue("NPP" in x and "LGM_warming" in x)
        self.assertFalse(any("genus" in str(k).lower() for k in x.columns))
        self.assertEqual(len(d),len(x))
        self.assertEqual(d.key.nunique(),len(d))
    def test_white_photo_NPP_direction_without_genus(self):
        d=self.source()
        m=estimate(d)
        self.assertEqual(m["n_species"],600)
        self.assertGreater(m["effects"]["NPP"]["beta_white_probability_per_species_SD"],0)
        self.assertEqual(m["model"],
             "genus_free_one_accepted_species_one_vote_OLS_HC3")
        r=estimate(d,label="mean9")
        self.assertEqual(r["status"],"fitted")
        g=estimate(d,adjustment="unadjusted")
        self.assertGreater(g["effects"]["NPP"]["beta_white_probability_per_species_SD"],0)
        logit=estimate(d,logistic=True)
        self.assertGreater(logit["effects"]["NPP"]["odds_ratio_per_species_SD"],1.)
    def test_low_quality_not_assumed_true_petal(self):
        d=self.source()
        d.loc[d.index[:3],"lgm_climate_coverage"]=.1
        d.loc[d.index[3:10],NPP_COVER]=.1
        filtered=sample_frame(d)
        self.assertEqual(len(filtered),590)
        self.assertEqual(estimate(filtered)["status"],"fitted")

if __name__=="__main__":
    unittest.main()
