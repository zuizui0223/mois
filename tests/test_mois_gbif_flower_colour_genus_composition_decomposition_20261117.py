"""Synthetic validation: colour genus-composition partition exactly reconstructs."""
import unittest
import numpy as np
import pandas as pd
from scripts.analysis.mois_gbif_flower_colour_genus_composition_decomposition_20261117 import (
    local_genus_expected, matched_spatial_slopes
)

class TestFlowerColourGenusDecomp(unittest.TestCase):
    def test_observed_white_expected_genus_plus_residual(self):
        source=pd.DataFrame({
            "specieskey":["1","2","3","4","5"],
            "grid_x":[1]*5,"grid_y":[2]*5,
            "genus":["A","A","B","B","C"],
            "status":["native_at_TDWG3"]*3+["introduced_at_TDWG3"]*2,
            "white":[True,False,True,False,True]
        })
        ref=pd.DataFrame({
            "genus":["A","B","C"],
            "genus_global_white_rate":[.5,.25,1.],
            "source_catalogue_species_n":[2,4,1]
        })
        d=local_genus_expected(source,ref)
        self.assertEqual(len(d),2)
        native=d[d.status.eq("native_at_TDWG3")].iloc[0]
        self.assertEqual(native.n_status_source_species,3)
        self.assertAlmostEqual(native.observed_white,2/3)
        self.assertAlmostEqual(native.expected_white_from_genus_composition,(.5+.5+.25)/3)
        self.assertAlmostEqual(native.within_genus_colour_residual,2/3-1.25/3)

    def test_spatial_slope_exact_reconstruction_without_bias(self):
        rng=np.random.default_rng(20261117)
        n=800
        group=np.repeat([f"g{i}" for i in range(80)],10)
        npp=rng.uniform(.05,1.5,n)
        pgen=rng.uniform(.2,.8,n)
        expected=np.clip(.3+.03*npp+.1*pgen,0,1)
        extra=np.clip(-.1+.04*npp+rng.normal(0,.02,n),-.3,.3)
        d=pd.DataFrame({
            "block_500km":group,
            "block_1000km":np.repeat([f"G{i}" for i in range(40)],20),
            "colour_species_n":rng.integers(30,500,n),
            "GBIF_occurrence_count":rng.integers(100,60000,n),
            "NPP_2001_2023_mean_kgC_m2_yr":npp,
            "present_temperature_C":rng.normal(12,5,n),
            "present_precip_mm":rng.uniform(100,2500,n),
            "MODIS_NDVI_2001_2021_mean":rng.uniform(.12,.86,n),
            "ESA_AGB_woody_Mg_ha_2020":rng.uniform(0,160,n),
            "threeGCM_delta_temperature_C_median":rng.normal(5,2,n),
            "threeGCM_delta_log_precip_median":rng.normal(0,.2,n),
            "lat":rng.uniform(-68,70,n),
            "observed_white":expected+extra,
            "expected_white_from_genus_composition":expected,
            "within_genus_colour_residual":extra
        })
        result=matched_spatial_slopes(d)
        self.assertLess(result["additivity_max_abs"],1e-12)
        for key in ("NPP","LGM_warming"):
            v=result["models"]["raw_observed_white"]["effects"][key]["beta_fraction_per_SD"]
            x=result["models"]["expected_from_genus_pool"]["effects"][key]["beta_fraction_per_SD"]
            y=result["models"]["residual_within_genus"]["effects"][key]["beta_fraction_per_SD"]
            self.assertAlmostEqual(v,x+y,places=12)

if __name__=="__main__":
    unittest.main()
