"""Source guards for 2024 actual source-photo spatial association and island cluster bootstrap."""
import unittest
import numpy as np
import pandas as pd
from scripts.analysis import mois_2024_sameyear_global_island_photo_white_20261107 as m

class SourceGuardTests(unittest.TestCase):
    def test_no_filter_on_high_lat_or_small_island_in_cell_grid(self):
        lon=pd.Series([-179.99,179.99,139.7,-122.4])
        lat=pd.Series([-89.99,89.99,35.7,68.2])
        key=m.cell_id(lon,lat)
        self.assertEqual(len(key),4)
        self.assertTrue(key.notna().all())
    def test_bootstrap_refuses_false_low_number_of_islands(self):
        fake=pd.DataFrame([
           dict(land_id=str(i),log_light=j,log_pop=j+.01,
              bio5=float(i),log_rain=.3,elevation=i*10,log_attempts=j+.5,
              white_pp=20+i+j)
           for i in range(3) for j in (.1,.4)])
        with self.assertRaisesRegex(ValueError,"island number"):
            m.island_bootstrap(fake)
    def test_all_data_source_hashes_present(self):
        self.assertEqual(set(m.SOURCE_SHA),
          {'canonical','wcvp','photo_landmass','sameyear_2024_environment'})
        self.assertTrue(all(len(s)==64 for s in m.SOURCE_SHA.values()))
    def test_model_rejects_insufficient_independent_geo_blocks(self):
        fake=pd.DataFrame([dict(land_id='1',log_light=.1,white_pp=30)]*4)
        with self.assertRaisesRegex(ValueError,"Too few"):
            m.cluster_fit(fake,"white_pp~log_light","land_id")
if __name__=="__main__":unittest.main()
