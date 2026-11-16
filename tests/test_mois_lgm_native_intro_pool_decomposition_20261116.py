"""Independent algebra/guard tests for signed GBIF WCVP species-pool partition."""
import unittest
import numpy as np
import pandas as pd
from scripts.analysis.mois_lgm_native_intro_pool_decomposition_20261116 import (
    decompose_cells, fit_shared_design, COMPONENTS
)

class TestPoolPartition(unittest.TestCase):
    def test_arithmetic_is_exact_and_unshared_genera_are_accounted(self):
        d=pd.DataFrame({
            "native_species_n":[100,80],
            "introduced_species_n":[80,120],
            "native_white_n":[60,45],
            "introduced_white_n":[40,90],
            "shared_native_species_n":[70,30],
            "shared_intro_species_n":[45,55],
            "shared_native_white_n":[42,20],
            "shared_intro_white_n":[18,30],
            "native_minus_intro_white":[.6-.5,45/80-90/120],
        })
        z=decompose_cells(d)
        self.assertTrue(np.allclose(z.reconstructed_white_difference,
                                    z.native_minus_intro_white,atol=1e-12))
        self.assertAlmostEqual(z.iloc[0].native_species_share_in_common_genera,.7)
        self.assertAlmostEqual(z.iloc[0].introduced_species_share_in_common_genera,.5625)
        self.assertAlmostEqual(z.iloc[0].contrast_inside_common_genus_pool,.12625)
        self.assertAlmostEqual(z.iloc[0].representation_of_common_genera,.06875)
        self.assertAlmostEqual(z.iloc[0].exclusive_genus_pool_white_contribution,-.095)

    def test_impossible_native_white_counts_rejected(self):
        d=pd.DataFrame({
            "native_species_n":[20],"introduced_species_n":[15],
            "native_white_n":[3],"introduced_white_n":[8],
            "shared_native_species_n":[10],"shared_intro_species_n":[8],
            "shared_native_white_n":[4],"shared_intro_white_n":[4],
            "native_minus_intro_white":[3/20-8/15]})
        with self.assertRaises(ValueError):
            decompose_cells(d)

    def test_model_coefficients_sum_exactly(self):
        rng=np.random.default_rng(20261116)
        n=900
        blocks=np.repeat([f"B{i}" for i in range(90)],10)
        warming=rng.normal(5,2,n)
        rain=rng.normal(0,1,n)
        temp=rng.normal(13,4,n)
        agb=np.abs(rng.normal(30,10,n))
        d=pd.DataFrame({
          "block_500km":blocks,"block_1000km":blocks,
          "threeGCM_delta_temperature_C_median":warming,
          "threeGCM_delta_log_precip_median":rain,
          "present_temperature_C":temp,
          "present_precip_mm":np.exp(rng.normal(6,1,n)),
          "ESA_AGB_woody_Mg_ha_2020":agb,
          "lat":rng.normal(40,8,n),
          "lon":rng.normal(10,12,n),
          "native_species_n":rng.integers(90,300,n),
          "introduced_species_n":rng.integers(70,180,n),
          "GBIF_occurrence_count":rng.integers(500,100000,n),
          "n_shared_genera":rng.integers(5,65,n),
          "genus_equal_native_minus_intro":rng.normal(0,.06,n)
        })
        d["shared_native_species_n"]=np.minimum(d.native_species_n,
                       np.maximum(10,(d.native_species_n*.65).astype(int)))
        d["shared_intro_species_n"]=np.minimum(d.introduced_species_n,
                       np.maximum(10,(d.introduced_species_n*.7).astype(int)))
        d["shared_native_white_n"]=np.floor(d.shared_native_species_n*.6).astype(int)
        d["shared_intro_white_n"]=np.floor(d.shared_intro_species_n*.45).astype(int)
        # Strictly no more whites in unmatched genera than unmatched species.
        d["native_white_n"]=d.shared_native_white_n+np.floor(
            (d.native_species_n-d.shared_native_species_n)*.5).astype(int)
        d["introduced_white_n"]=d.shared_intro_white_n+np.floor(
            (d.introduced_species_n-d.shared_intro_species_n)*.4).astype(int)
        d["native_minus_intro_white"]=(d.native_white_n/d.native_species_n
                    -d.introduced_white_n/d.introduced_species_n)
        z=decompose_cells(d)
        fit=fit_shared_design(z, "500")
        m=fit["component_models"]
        self.assertLess(abs(m[COMPONENTS[0]]["beta_fraction_per_SD"]-sum(
            m[name]["beta_fraction_per_SD"] for name in COMPONENTS[1:4])),1e-10)
        self.assertEqual(m["native_minus_intro_white"]["n_cells"],900)

if __name__=="__main__":
    unittest.main()
