"""Independent implementation check: same-cell family null vs AGB/NDVI."""
import unittest
import numpy as np
import pandas as pd
from scripts.analysis.mois_gbif_family_excess_biomass_ndvi_20261116 import (
    normalized_model_matrix,fit
)

class TestBiomassAndGenusNull(unittest.TestCase):
    def sample(self):
        rng=np.random.default_rng(20261116)
        n=600
        b=np.repeat([f"B{i}" for i in range(60)],10)
        a=np.maximum(0,rng.normal(35,15,n))
        v=rng.normal(.52,.13,n)
        v=np.clip(v,0,1)
        y=.015*np.log1p(a)+rng.normal(0,.012,n)
        f=pd.DataFrame({
            "block_500km":b,"block_1000km":b,
            "threeGCM_delta_temperature_C_median":rng.normal(5,2,n),
            "threeGCM_delta_log_precip_median":rng.normal(0,.05,n),
            "present_temperature_C":rng.normal(12,5,n),
            "present_precip_mm":np.exp(rng.normal(6,1,n)),
            "ESA_AGB_woody_Mg_ha_2020":a,
            "MODIS_NDVI_2001_2021_mean":v,
            "GBIF_occurrence_count":rng.integers(100,100000,n),
            "native_species_n":rng.integers(40,1000,n),
            "introduced_species_n":rng.integers(40,350,n),
            "family_union":rng.integers(20,150,n),
            "genus_union":rng.integers(60,900,n),
            "lat":rng.uniform(-80,80,n),
            "residual_after_family_conditioning":y,
            "observed_jaccard":np.clip(.08+y,0,1)
        })
        return f
    def test_positive_woody_AGB_signal_is_fittable(self):
        d=self.sample()
        x,sd=normalized_model_matrix(d)
        self.assertEqual(len(x),len(d))
        self.assertTrue("log_woody_AGB" in x)
        r=fit(d)
        self.assertGreater(r["n_500km_clusters"],25)
        self.assertGreater(r["effect"]["log_woody_AGB"]["beta_per_global_SD"],0)
        self.assertFalse(r["woody_AGB_is_NPP_or_flower_biomass"])
    def test_ablation_and_species_support_weight(self):
        d=self.sample()
        no_agb=fit(d,include_agb=False)
        no_ndvi=fit(d,include_ndvi=False)
        weighted=fit(d,weighted=True)
        self.assertNotIn("log_woody_AGB",no_agb["effect"])
        self.assertNotIn("NDVI_mean",no_ndvi["effect"])
        self.assertEqual(weighted["weights"],"min_native_intro_species_count")

if __name__=="__main__":
    unittest.main()
