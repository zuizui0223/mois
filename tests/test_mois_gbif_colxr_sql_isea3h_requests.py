"""No-request tests for COL XR global ISEA3H species-key SQL cube packaging."""
import unittest
import pandas as pd
from scripts.analysis.mois_gbif_colxr_sql_isea3h_requests import (
    sql_for_species_keys,make_requests
)


class GbifSQLCubeTests(unittest.TestCase):
    def test_sql_uses_colxr_specieskey_global_equal_area_grid(self):
        q=sql_for_species_keys(["MKYB","34RMK","MKYB2"],resolution=8)
        self.assertIn("GBIF_ISEA3HCode(8, decimalLatitude, decimalLongitude, 0)",q)
        self.assertIn("specieskey IN",q)
        self.assertIn("hasgeospatialissues = false",q)
        self.assertIn("GROUP BY specieskey",q)
        self.assertNotIn("gbif_classification",q)
        with self.assertRaisesRegex(ValueError,"valid alpha-numeric"):
            sql_for_species_keys(["'; DELETE FROM occurrence; --"])

    def test_requests_do_not_lose_species_or_synonym_conflict(self):
        rows=pd.DataFrame({
            "status":["colxr_exact_accepted_species"]*3+
                     ["synonym_or_doubtful_review"],
            "key":["ABC","DEF","GHI","OLD"]})
        chunks=make_requests(rows,chunk_size=2)
        self.assertEqual([len(x["source_species_keys"]) for x in chunks],[2,1])
        s=str([x["GBIF_download_request"] for x in chunks])
        self.assertNotIn("'OLD'",s)
        self.assertIn("ABC",s)


if __name__=="__main__":
    unittest.main()
