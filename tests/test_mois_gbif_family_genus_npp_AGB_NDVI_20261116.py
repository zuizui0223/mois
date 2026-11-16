"""GBIF + 23-year NPP regression synthesis tests on synthetic spatial blocks."""
import unittest
import numpy as np
import pandas as pd
from scripts.analysis.mois_gbif_family_genus_npp_AGB_NDVI_20261116 import (
    matrix_with_NPP,fit
)

class TestNPPIncrement(unittest.TestCase):
    def fixture(self,n=600):
        rng=np.random.default_rng(20261116)
        npp=rng.uniform(.02,1.4,n)
        agb=rng.uniform(0,150,n)
        ndvi=np.clip(.5*npp+rng.normal(.1,.1,n),0,1)
        d=pd.DataFrame({
          "block_500km":np.repeat([f"B{i:03d}" for i in range(n//10)],10),
          "block_1000km":np.repeat([f"U{i:03d}" for i in range(n//20)],20),
          "threeGCM_delta_temperature_C_median":rng.normal(5,2,n),
          "threeGCM_delta_log_precip_median":rng.normal(.02,.1,n),
          "present_temperature_C":rng.normal(14,4,n),
          "present_precip_mm":np.exp(rng.normal(6,.7,n)),
          "ESA_AGB_woody_Mg_ha_2020":agb,
          "MODIS_NDVI_2001_2021_mean":ndvi,
          "GBIF_occurrence_count":rng.integers(30,20000,n),
          "native_species_n":rng.integers(30,300,n),
          "introduced_species_n":rng.integers(30,180,n),
          "family_union":rng.integers(10,80,n),
          "genus_union":rng.integers(50,250,n),
          "lat":rng.uniform(-65,65,n),
          "NPP_2001_2023_mean_kgC_m2_yr":npp,
          "NPP_2020_kgC_m2_yr":npp+rng.normal(0,.05,n)
        })
        residual=.025*npp+rng.normal(0,.008,n)
        d["residual_after_family_conditioning"]=residual
        d["observed_jaccard"]=np.clip(.12+residual,0,1)
        return d
    def test_npp_support_independent_from_agb_ndvi(self):
        d=self.fixture()
        x,sds=matrix_with_NPP(d)
        self.assertEqual(x.shape[0],600)
        self.assertIn("NPP_annual_carbon_flux",x.columns)
        full=fit(d)
        no_npp=fit(d,npp=False)
        self.assertGreater(full["effects"]["NPP_annual_carbon_flux"]["beta_per_global_SD"],0)
        self.assertNotIn("NPP_annual_carbon_flux",no_npp["effects"])
        self.assertGreaterEqual(full["R2_in_sample"]+1e-10,no_npp["R2_in_sample"])
        self.assertTrue(full["NPP_is_annual_carbon_flux_not_biomass"] if
                        "NPP_is_annual_carbon_flux_not_biomass" in full
                        else full["NPP_is_annual_carbon_flux_not_biomass"])
    def test_sensitivity_all_three_vegetation_groups(self):
        d=self.fixture()
        rows=[
          fit(d,agb=False),
          fit(d,ndvi=False),
          fit(d,agb=False,ndvi=False),
          fit(d,block="1000"),
          fit(d,support_weighted=True),
          fit(d,year2020_only=True)
        ]
        self.assertTrue(all(x["n_cells"]==600 for x in rows))
        self.assertNotIn("log_woody_AGB",rows[0]["effects"])
        self.assertNotIn("NDVI_mean",rows[1]["effects"])
        self.assertIn("NPP_annual_carbon_flux",rows[2]["effects"])
        self.assertEqual(rows[4]["weights"],"min_WCVP_status_group_species_count")
        self.assertEqual(rows[5]["npp_years"],"2020")
if __name__=="__main__":
    unittest.main()
