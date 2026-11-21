"""Independent unit-level guards: species photo-colour, NPP, paleoclimate."""
import unittest
import numpy as np
import pandas as pd

from scripts.analysis.mois_flower_colour_LGM_NPP_GBIF_20261117 import (
    species_model_advanced,cell_white_model,cell_cohort,NPP_NAME,
    LGM_TEMP,LGM_RAIN
)

class TestFlowerColourCurrentProductivity(unittest.TestCase):
    def synthetic_species(self):
        rng=np.random.default_rng(17112026)
        n=960
        genera=np.repeat([f"G_{i}" for i in range(80)],12)
        npp=rng.normal(.6,.2,n)
        colour_white=(npp+.2*rng.normal(size=n)>.57)
        return pd.DataFrame({
          "genus":genera,
          "source_hue9":np.where(colour_white,"white","purple"),
          "hue_status":np.where(rng.random(n)>.5,"sampled_multicolour",
                                 "sampling_consensus_single"),
          "BIO01":rng.normal(12,4,n),
          "BIO12":rng.uniform(200,2500,n),
          "MODIS_NDVI_2001_2021_mean":rng.uniform(.1,.9,n),
          "ESA_AGB_woody_Mg_ha_2020":rng.uniform(0,150,n),
          "mean_abs_latitude_over_occupied_cells":rng.uniform(4,68,n),
          "records_passed_spatial_QA":rng.integers(10,10000,n),
          LGM_TEMP:rng.normal(5,1.3,n),
          LGM_RAIN:rng.normal(.1,.12,n),
          NPP_NAME:npp,
          "NPP_2020_GBIF_kgC_m2_yr":npp*.9,
        })
    def test_within_genus_source_colour_and_NPP(self):
        d=self.synthetic_species()
        x=species_model_advanced(d,"white")
        self.assertEqual(x["status"],"fitted")
        self.assertEqual(x["n_genera"],80)
        for w in ("species_equal","genus_equal"):
            model=x["models"][w]
            self.assertIn("NPP",model["with_NPP"]["effects"])
            self.assertNotIn("NPP",model["without_NPP"]["effects"])
            self.assertGreater(model["with_NPP"]["effects"]["NPP"][
                "beta_probability_points_per_sd"],0)
            self.assertGreaterEqual(model["delta_R2_add_NPP"]+1e-9,0)
    def test_synthetic_cell_block_and_photo_confounder(self):
        rng=np.random.default_rng(20261117)
        n=1000
        b=np.repeat([f"B{i}" for i in range(100)],10)
        npp=rng.uniform(.05,1.6,n)
        x=pd.DataFrame({
          "block_500km":b,"block_1000km":np.repeat([
              f"C{i}" for i in range(50)],20),
          "threeGCM_delta_temperature_C_median":rng.normal(5,2,n),
          "threeGCM_delta_log_precip_median":rng.normal(0,.2,n),
          "present_temperature_C":rng.normal(13,4,n),
          "present_precip_mm":rng.uniform(100,2500,n),
          "MODIS_NDVI_2001_2021_mean":rng.uniform(.1,.85,n),
          "ESA_AGB_woody_Mg_ha_2020":rng.uniform(0,180,n),
          "NPP_2001_2023_mean_kgC_m2_yr":npp,
          "GBIF_occurrence_count":rng.integers(100,25000,n),
          "lat":rng.uniform(-70,70,n),
          "colour_species_n":rng.integers(45,250,n),
          "white_photo_species_fraction":np.clip(.45+.16*npp
                   +rng.normal(0,.04,n),.02,.99)
        })
        m=cell_white_model(x)
        n0=cell_white_model(x,npp=False)
        self.assertEqual(m["n_cells"],1000)
        self.assertEqual(m["n_spatial_blocks_500km"],100)
        self.assertIn("NPP",m["effects"])
        self.assertGreater(m["effects"]["NPP"]["beta_fraction_per_SD"],0)
        self.assertNotIn("NPP",n0["effects"])
        self.assertTrue(m["no_local_flower_colour_measured"])
        w=cell_white_model(x,weight=True,block="1000")
        self.assertEqual(w["n_cells"],1000)

if __name__=="__main__":
    unittest.main()
