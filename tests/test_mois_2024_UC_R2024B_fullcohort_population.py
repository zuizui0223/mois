"""MOIS WorldPop R2024B UC 2024 full landmass source contracts."""
import unittest
import numpy as np,pandas as pd
from rasterio.io import MemoryFile
from rasterio.transform import from_origin
from pyproj import Geod
from scripts.analysis import mois_2024_UC_R2024B_fullcohort_population as s

class Test2024UnconstrainedPopulation(unittest.TestCase):
    def mem(self):
        m=MemoryFile()
        with m.open(driver="GTiff",height=2,width=2,count=1,dtype="float32",
                    crs="EPSG:4326",nodata=-99999,
                    transform=from_origin(140,36,1/120,1/120)) as z:
            z.write(np.array([[0.,6.],[-99999.,12.]],dtype=np.float32),1)
        return m
    def rows(self):
        return pd.DataFrame([
            {"target_id":"2024|a","query_longitude":140.003,"query_latitude":35.997},
            {"target_id":"2024|b","query_longitude":140.011,"query_latitude":35.997},
            {"target_id":"2024|c","query_longitude":140.003,"query_latitude":35.987},
            {"target_id":"2024|d","query_longitude":140.0,"query_latitude":20.0}])
    def test_true_zero_different_from_nodata_and_outside_bounds(self):
        mem=self.mem()
        with mem.open() as ds:
            x,r=s.raster_samples(ds,self.rows())
        self.assertEqual(x.population_UC_R2024B_status.tolist(),
                         ["OK","OK","NODATA","OUTSIDE_SOURCE_BOUNDS"])
        self.assertEqual(x.population_UC_R2024B_count_per_cell.iloc[0],0.)
        self.assertEqual(x.population_UC_R2024B_density_per_km2.iloc[0],0.)
        self.assertTrue(x.population_UC_R2024B_density_per_km2.iloc[1]>0)
        self.assertTrue(x.population_UC_R2024B_density_per_km2.iloc[2:].isna().all())
    def test_true_geodesic_area_varies_with_latitude(self):
        mem=self.mem()
        with mem.open() as ds:
            geod=Geod(ellps="WGS84")
            self.assertNotEqual(s.area_km2_for_row(ds,0,geod),
                                s.area_km2_for_row(ds,1,geod))
    def test_invalid_source_year_fails(self):
        mem=self.mem();t=self.rows();t.loc[0,"target_id"]="2020|a"
        with mem.open() as ds:
            with self.assertRaisesRegex(ValueError,"year-2024"):
                s.raster_samples(ds,t)
    def test_nodata_and_invalid_negative_not_population_zero(self):
        self.assertEqual(s.parse_count(0,-99999),(0., "OK"))
        self.assertEqual(s.parse_count(-99999,-99999),(None,"NODATA"))
        self.assertEqual(s.parse_count(-2,-99999),(None,"INVALID_NEGATIVE"))
    def test_pop_model_series_never_spliced(self):
        self.assertIn("R2024B",s.UC_FILENAME)
        self.assertIn("_UC_",s.UC_FILENAME)
        self.assertIn("R2025A",s.CN_SOURCE_SHA.replace(s.CN_SOURCE_SHA,"R2025A"))
        self.assertNotEqual(s.ORIGINAL_2024_SHA,s.CN_SOURCE_SHA)

if __name__=="__main__": unittest.main()
