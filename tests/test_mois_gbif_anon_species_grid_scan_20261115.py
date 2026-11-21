import unittest
import numpy as np
from scripts.analysis.mois_gbif_anon_species_grid_scan_20261115 import (
    canonical_pair,flags,allowed_coords,BAD_ISSUES
)
class TestAnonGBIFSpatial(unittest.TestCase):
    def test_hyphenated_accepted_species_preserved(self):
        self.assertEqual(canonical_pair("Aphelandra lingua-bovis (L.)"),"Aphelandra lingua-bovis")
        self.assertEqual(canonical_pair("Trigonella foenum-graecum"),"Trigonella foenum-graecum")
        self.assertEqual(canonical_pair("Solanum polyacanthos Lam."),"Solanum polyacanthos")
        self.assertIsNone(canonical_pair("Solanum"))
    def test_issue_nested_array_dict_flattened(self):
        self.assertEqual(flags([{"array":"ZERO_COORDINATE"},
                                {"array":{"array_element":"COUNTRY_COORDINATE_MISMATCH"}}]),
                         ["ZERO_COORDINATE","COUNTRY_COORDINATE_MISMATCH"])
        self.assertEqual(flags(None),[])
    def test_0_0_and_geospatial_issues_rejected(self):
        self.assertFalse(allowed_coords(0,0,float("nan"),None))
        self.assertFalse(allowed_coords(130,35,float("nan"),["COUNTRY_COORDINATE_MISMATCH"]))
        self.assertFalse(allowed_coords(130,35,90000,None))
        self.assertFalse(allowed_coords(181,35,0,None))
        self.assertTrue(allowed_coords(130,35,300,"GEODETIC_DATUM_ASSUMED_WGS84"))
    def test_exact_geospatial_issue_flags_pinned(self):
        self.assertEqual(len(BAD_ISSUES),4)
        self.assertIn("COUNTRY_COORDINATE_MISMATCH",BAD_ISSUES)
if __name__=="__main__":unittest.main()
