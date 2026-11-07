"""Same genuine photo cell + GSHHG landmass, time FE and measured ROI contracts."""
import unittest
import numpy as np
import pandas as pd
import statsmodels.formula.api as sf
from scripts.analysis import mois_island_ALAN_repeated_sites_photo_year_20261108 as m

class OriginalRepeatedPhotoSites(unittest.TestCase):
    def test_site_has_original_cell_and_original_landmass_but_not_year(self):
        d=pd.DataFrame([
            dict(target_id="2019|6_7|10",year=2019,land_id="10",
                 land_class="island",region5="A",is_island=1),
            dict(target_id="2021|6_7|10",year=2021,land_id="10",
                 land_class="island",region5="A",is_island=1),
            dict(target_id="2023|6_7|20",year=2023,land_id="20",
                 land_class="mainland",region5="A",is_island=0),
            dict(target_id="2025|7_7|10",year=2025,land_id="10",
                 land_class="island",region5="A",is_island=1),
        ])
        a=m.authentic_site_identity(d)
        self.assertEqual(int(a.repeat_site_id.nunique()),3)
        full,repeat=m.repeat_site_frame(d)
        self.assertEqual(len(repeat),2)
        self.assertEqual(int(repeat.year.nunique()),2)
        self.assertSetEqual(set(repeat.repeat_site_id),{"6_7|10"})
    def test_same_source_photo_site_landmass_conflict_rejected(self):
        d=pd.DataFrame([
            dict(target_id="2019|6_7|10",year=2019,land_id="99",
                 land_class="island",region5="A",is_island=1)])
        with self.assertRaisesRegex(ValueError,"landmass ID"):
            m.authentic_site_identity(d)
    def test_year_and_site_demean_equivalent_to_dense_ols(self):
        rng=np.random.default_rng(123456)
        rows=[]
        for i in range(36):
            island=int(i%2==0)
            site=f"{i+150}_{i+240}|{11 if island else 25}"
            region=f"g{i%15}"
            site_baseline=float(rng.normal(0,1))
            years=m.YEARS if i%3==0 else (2019,2021,2024,2025)
            for year in years:
                light=.7+float(rng.normal(0,.25))+(year-2019)*.04
                population=1.6+float(rng.normal(0,.35))
                effort=1.5+float(rng.normal(0,.25))
                value=(43.0+site_baseline+2.*light
                    +.2*(year-2019)+.8*island*light-.4*population+
                    .5*effort+float(rng.normal(0,.02)))
                rows.append(dict(repeat_site_id=site,
                    year=year,log_light=light,is_island=island,
                    log_pop=population,log_attempts=effort,
                    white_pp=value,region5=region))
        d=pd.DataFrame(rows)
        ordered,x,y,names=m.site_demean_with_year_FE(d)
        original=sf.ols(
           "white_pp ~ log_light + log_light:is_island + log_pop + log_attempts"
           " + C(repeat_site_id) + C(year)",data=d).fit()
        value=np.linalg.lstsq(x,y,rcond=None)[0]
        self.assertAlmostEqual(float(value[0]),float(original.params["log_light"]),places=6)
        self.assertAlmostEqual(float(value[1]),float(original.params["log_light:is_island"]),places=6)
        self.assertAlmostEqual(float(value[2]),float(original.params["log_pop"]),places=6)
        self.assertAlmostEqual(float(value[3]),float(original.params["log_attempts"]),places=6)
    def test_identical_original_photo_years_do_not_create_change_information(self):
        source=pd.DataFrame([
          dict(target_id="2019|1_2|1",year=2019,land_id="1",land_class="island",
               region5="x",is_island=1,log_light=.3,white_pp=51.,
               log_pop=0.,log_attempts=1.2),
          dict(target_id="2019|2_3|2",year=2019,land_id="2",land_class="mainland",
               region5="x",is_island=0,log_light=.5,white_pp=42.,
               log_pop=0.,log_attempts=1.1),
        ])
        result=m.estimate_within_site(source)
        self.assertFalse(result["estimable"])
        self.assertEqual(result["reason"],"source_repeated_site_or_temporal_exposure_support_gate")
        self.assertEqual(result["source_support"]["same_geographic_cell_and_GSHHG_landmass_repeated_2plus_years"]["photo_cell_years"],0)
    def test_observed_panel_uses_only_true_years_and_no_false_species_evolution(self):
        self.assertEqual(m.YEARS,(2019,2021,2023,2024,2025))
        self.assertEqual(m.SOURCE_N,{"strict_5p3r3s":7982,
                                     "measured_3p2r2s":20543})
        import inspect
        source=inspect.getsource(m.run)
        self.assertIn('"no_MNAR_or_old_tipping_sensitivity":True',source)
        self.assertIn('"no_missing_site_year_source_value_imputed":True',source)
        self.assertIn('"no_new_original_photographs_or_ROI_inference":True',source)

if __name__=="__main__":
    unittest.main()
