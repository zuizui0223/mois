import unittest
import numpy as np
import pandas as pd
from scripts.analysis.mois_lgm_gbif_genus_overlap_conditional_null_20261116 import (
    p_all_members_same_status, expected_genus_overlap
)

class TestGBIFConditionalGenusNull(unittest.TestCase):
    def test_exact_hypergeom_small_groups(self):
        # From four species, split two native, two introduced;
        # the chance a particular two-species genus is all native is 1/6.
        x=p_all_members_same_status(np.array([4]),np.array([2]),np.array([2]))
        self.assertAlmostEqual(x[0],1/6,places=12)
        zero=p_all_members_same_status(np.array([4]),np.array([1]),np.array([2]))
        self.assertEqual(zero[0],0)
        whole=p_all_members_same_status(np.array([4]),np.array([4]),np.array([2]))
        self.assertAlmostEqual(whole[0],1.)

    def test_jaccard_expected_both_species_in_one_genus(self):
        q=pd.DataFrame({
            "specieskey":list("abcd"),
            "grid_x":[2]*4,"grid_y":[3]*4,
            "genus":["A","A","B","C"],
            "status":["native_at_TDWG3","introduced_at_TDWG3",
                      "native_at_TDWG3","introduced_at_TDWG3"]
        })
        a=expected_genus_overlap(q).iloc[0]
        self.assertEqual(a.species_total,4)
        self.assertEqual(a.native_total,2)
        self.assertEqual(a.genus_union,3)
        self.assertEqual(a.observed_shared_genera,1)
        self.assertAlmostEqual(a.observed_jaccard,1/3)
        # Genus A has m=2 -> P(both)=1-1/6-1/6=2/3.
        # B,C have m=1 -> zero expected shared.
        self.assertAlmostEqual(a.expected_jaccard_if_status_random,2/9)
        self.assertAlmostEqual(a.excess_jaccard_vs_status_random,1/9)

    def test_disjoint_single_species_genera_null_zero(self):
        d=pd.DataFrame({
            "specieskey":list("abcd"),
            "grid_x":[1]*4,"grid_y":[5]*4,
            "genus":list("ABCD"),
            "status":["native_at_TDWG3"]*2+["introduced_at_TDWG3"]*2
        })
        a=expected_genus_overlap(d).iloc[0]
        self.assertAlmostEqual(a.expected_jaccard_if_status_random,0)
        self.assertAlmostEqual(a.observed_jaccard,0)
        self.assertAlmostEqual(a.excess_jaccard_vs_status_random,0)

if __name__=="__main__":
    unittest.main()
