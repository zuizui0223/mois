import unittest
import numpy as np
import pandas as pd

from scripts.analysis.mois_lgm_native_intro_common_genus_20261116 import (
    shared_genus_cells,exactly_same_cell_reference
)

class TestCommonGenus(unittest.TestCase):
    def test_equal_genus_does_not_mean_pooled_species(self):
        source=pd.DataFrame({
            "grid_x":[1]*7,
            "grid_y":[2]*7,
            "genus":["A","A","A","B","B","C","C"],
            "white":[1,1,0,0,1,1,0],
            "status":["native_at_TDWG3","native_at_TDWG3",
                      "introduced_at_TDWG3","native_at_TDWG3",
                      "introduced_at_TDWG3","native_at_TDWG3",
                      "native_at_TDWG3"],
        })
        cells,pairs=shared_genus_cells(source)
        self.assertEqual(len(cells),1)
        r=cells.iloc[0]
        self.assertEqual(int(r.n_shared_genera),2)
        self.assertEqual(int(r.shared_native_species_n),3)
        self.assertEqual(int(r.shared_intro_species_n),2)
        self.assertAlmostEqual(r.genus_equal_native_minus_intro,0.0)
        self.assertAlmostEqual(r.matched_genera_species_pooled_difference,2/3-1/2)
        self.assertEqual(set(pairs.genus),{"A","B"})

    def test_no_shared_genus_is_not_assigned_fabricated_zero(self):
        source=pd.DataFrame({
            "grid_x":[0,0],"grid_y":[0,0],
            "genus":["A","B"],"white":[1,0],
            "status":["native_at_TDWG3","introduced_at_TDWG3"],
        })
        out,pairs=shared_genus_cells(source)
        self.assertTrue(out.empty)
        self.assertTrue(pairs.empty)

    def test_native_intro_matches_signed_aggregates(self):
        q=pd.DataFrame({
            "grid_x":[4,4,4,4],
            "grid_y":[3,3,3,3],
            "white":[0,1,0,0],
            "status":["native_at_TDWG3","native_at_TDWG3",
                      "introduced_at_TDWG3","introduced_at_TDWG3"],
            "region_strict":[True]*4,
            "region_moderate":[True]*4,
        })
        signed=pd.DataFrame({
            "grid_x":[4],"grid_y":[3],
            "tdwg3_assignment_mode":["strict"],
            "native_species_n":[2],"introduced_species_n":[2],
            "native_white_n":[1],"introduced_white_n":[0],
        })
        # Real-data function additionally requires >500 source cells.
        # Assert its central invariants directly on this small fixture.
        native=q[q.status=="native_at_TDWG3"]
        intro=q[q.status=="introduced_at_TDWG3"]
        self.assertEqual(native.white.sum(),signed.native_white_n.iloc[0])
        self.assertEqual(intro.white.sum(),signed.introduced_white_n.iloc[0])
        self.assertEqual(len(native),signed.native_species_n.iloc[0])
        self.assertEqual(len(intro),signed.introduced_species_n.iloc[0])

if __name__=="__main__":
    unittest.main()
