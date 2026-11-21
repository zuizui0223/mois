"""Native-introduced genus Jaccard is defined without photographic hue."""
import unittest
import pandas as pd
from scripts.analysis.mois_lgm_gbif_native_intro_genus_turnover_20261116 import (
    genus_overlap
)

class TestGBIFGenusTurnover(unittest.TestCase):
    def test_species_counts_and_genus_union(self):
        d=pd.DataFrame({
            "grid_x":[1]*6+[2]*3,
            "grid_y":[3]*6+[4]*3,
            "specieskey":list("abcdefghi"),
            "genus":["A","A","B","A","C","D","E","E","F"],
            "status":["native_at_TDWG3","native_at_TDWG3",
                      "native_at_TDWG3","introduced_at_TDWG3",
                      "introduced_at_TDWG3","introduced_at_TDWG3",
                      "native_at_TDWG3","introduced_at_TDWG3",
                      "introduced_at_TDWG3"]
        })
        o=genus_overlap(d)
        self.assertEqual(len(o),2)
        a=o[o.grid_x.eq(1)].iloc[0]
        self.assertEqual(int(a.native_genus_n),2)
        self.assertEqual(int(a.introduced_genus_n),3)
        self.assertEqual(int(a.shared_genus_n),1)
        self.assertAlmostEqual(a.jaccard_native_introduced_genus,1/4)
        self.assertAlmostEqual(a.genus_turnover_one_minus_jaccard,.75)
        self.assertAlmostEqual(a.native_exclusive_genus_fraction,.5)
        self.assertAlmostEqual(a.introduced_exclusive_genus_fraction,2/3)
        b=o[o.grid_x.eq(2)].iloc[0]
        self.assertAlmostEqual(b.jaccard_native_introduced_genus,.5)
    def test_repeated_species_cell_id_is_not_an_independent_taxon(self):
        d=pd.DataFrame({
            "grid_x":[1,1],"grid_y":[2,2],
            "specieskey":["a","a"],
            "genus":["A","A"],"status":["native_at_TDWG3"]*2
        })
        with self.assertRaises(ValueError):
            genus_overlap(d)

if __name__=="__main__":
    unittest.main()
