"""Offline guards for 2025 true-year global original-photo island colour × ALAN."""
import unittest
import pandas as pd
from scripts.analysis import mois_2025_sameyear_island_ALAN_white_20261107 as m
from scripts.analysis import mois_2019_2023_sameyear_island_ALAN_white_20261107 as earlier

class Year2025Tests(unittest.TestCase):
    def test_source_FCP_original_geo_hash_and_original_photos_pinned(self):
        self.assertEqual(len(m.GEO_SHA),64)
        self.assertEqual(len(m.SITES_SHA),64)
        self.assertEqual(m.PHOTO_SHA,earlier.PHOTO_SHA)
        self.assertEqual(m.WCVP_SHA,earlier.WCVP_SHA)
    def test_photo_derived_island_slopes_require_supported_island_cells(self):
        d=pd.DataFrame([dict(region5=str(i),land_class="island" if i<3 else "mainland",
             land_id=str(i) if i<3 else "bigcontinent",white_pp=i)
             for i in range(100)])
        result=m.fit_if_supported(d,"white_pp~1","region5",interaction=True)
        self.assertFalse(result["estimable"])
        self.assertEqual(result["independent_islands"],3)
    def test_source_year_raster_contract_is_explicit(self):
        import inspect
        s=inspect.getsource(m.run)
        self.assertIn("population_raster_year",s)
        self.assertIn("VIIRS_year",s)
        self.assertIn("source_photo_year",s)
        self.assertIn("population_density_per_km2",s)
if __name__=="__main__":
    unittest.main()
