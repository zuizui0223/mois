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
    def test_island_light_interaction_is_a_separate_estimand(self):
        # Source-free synthetic design: verify the model actually estimates
        # a separate island-minus-mainland nightlight slope, not just an
        # island mean offset. No ecological effect is inferred from these data.
        rows=[]
        for region in range(8):
            for is_island in (0,1):
                for j in range(4):
                    light=0.2 + 0.17*j + 0.015*region
                    rows.append(dict(region5=str(region),is_island=is_island,
                        log_light=light,
                        white_pp=41 + 0.3*region + 1.5*is_island
                                 - 1.0*light - 3.2*is_island*light
                                 + 0.02*((region*7+j*11+is_island*13)%5)))
        d=pd.DataFrame(rows)
        out=m.cluster_fit(d,
            "white_pp ~ log_light * is_island + C(region5)","region5")
        self.assertIn("log_light:is_island",out)
        self.assertLess(out["log_light:is_island"]["beta"],-2.0)
        self.assertIn("is_island",out)
        self.assertTrue(np.isfinite(out["log_light:is_island"]["SE"]))

    def test_model_rejects_insufficient_independent_geo_blocks(self):
        fake=pd.DataFrame([dict(land_id='1',log_light=.1,white_pp=30)]*4)
        with self.assertRaisesRegex(ValueError,"Too few"):
            m.cluster_fit(fake,"white_pp~log_light","land_id")
if __name__=="__main__":unittest.main()
