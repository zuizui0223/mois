"""Exact family-constrained genus overlap under WCVP class-label exchange."""
import unittest
import numpy as np
import pandas as pd
from scripts.analysis.mois_lgm_gbif_family_conditional_genus_null_20261116 import (
    genus_null_both_species_counts
)

class TestFamilyConditionedGBIF(unittest.TestCase):
    def test_family_status_segregation_is_explained_by_family_null(self):
        # Two families have only one status each. A naive global shuffle
        # expects shared genera; a within-family shuffle cannot produce any.
        d=pd.DataFrame({
            "grid_x":[1]*4,"grid_y":[2]*4,
            "specieskey":["a","b","c","d"],
            "WCVP_family":["F1","F1","F2","F2"],
            "genus":["A","A","B","B"],
            "status":["native_at_TDWG3","native_at_TDWG3",
                      "introduced_at_TDWG3","introduced_at_TDWG3"]
        })
        z=genus_null_both_species_counts(d).iloc[0]
        self.assertEqual(z.family_union,2)
        self.assertEqual(z.native_families,1)
        self.assertEqual(z.introduced_families,1)
        self.assertEqual(z.shared_families,0)
        self.assertAlmostEqual(z.observed_jaccard,0)
        self.assertAlmostEqual(z.family_conditioned_expected_jaccard,0)
        self.assertAlmostEqual(z.unstratified_expected_jaccard,2/3)
        self.assertAlmostEqual(z.residual_after_family_conditioning,0)

    def test_hypergeometric_family_and_global_expected_differ(self):
        # F1 has native A+B and introduced A+C; F2 has native D,D and
        # introduced E. Species status is evaluated independently per cell.
        d=pd.DataFrame({
            "grid_x":[4]*7,"grid_y":[8]*7,
            "specieskey":list("abcdefg"),
            "WCVP_family":["F1"]*4+["F2"]*3,
            "genus":["A","A","B","C","D","D","E"],
            "status":["native_at_TDWG3","introduced_at_TDWG3",
                      "native_at_TDWG3","introduced_at_TDWG3",
                      "native_at_TDWG3","native_at_TDWG3",
                      "introduced_at_TDWG3"]
        })
        z=genus_null_both_species_counts(d).iloc[0]
        self.assertEqual(z.native_total,4)
        self.assertEqual(z.genus_union,5)
        self.assertEqual(z.observed_shared_genera,1)
        self.assertEqual(z.shared_families,2)
        self.assertAlmostEqual(z.observed_jaccard,1/5)
        self.assertAlmostEqual(z.family_conditioned_expected_jaccard,4/15)
        self.assertAlmostEqual(z.unstratified_expected_jaccard,8/35)
        self.assertAlmostEqual(z.residual_after_family_conditioning,1/5-4/15)
        self.assertTrue(np.isfinite(z.residual_after_unstratified_conditioning))

    def test_repeat_species_cell_is_forbidden(self):
        d=pd.DataFrame({
            "grid_x":[1,1],"grid_y":[2,2],
            "specieskey":["x","x"],
            "WCVP_family":["F","F"],
            "genus":["A","A"],
            "status":["native_at_TDWG3"]*2
        })
        with self.assertRaises(ValueError):
            genus_null_both_species_counts(d)

if __name__=="__main__":
    unittest.main()
