"""Offline no-network guards for 5-year original photographic island × ALAN panel."""
import unittest
import pandas as pd
import numpy as np
from scripts.analysis import mois_island_ALAN_multiyear_sameyear_source_panel_20261107 as m

class MultiyearPanelContracts(unittest.TestCase):
    def test_true_years_only_and_no_2022_or_2020_source_splicing(self):
        self.assertEqual(m.PHOTO_YEARS,(2019,2021,2023,2024,2025))
        self.assertNotIn(2022,m.PHOTO_YEARS)
        self.assertNotIn(2020,m.PHOTO_YEARS)
        self.assertEqual(sum(m.YEAR_COUNTS.values()),7982)
    def test_same_spatial_photo_cell_in_other_year_is_not_another_island(self):
        d=pd.DataFrame([
            {"land_id":"12","cell_id_without_year":"1_1","year":2019,"log_light":1.2},
            {"land_id":"12","cell_id_without_year":"1_1","year":2021,"log_light":1.7},
            {"land_id":"12","cell_id_without_year":"2_1","year":2023,"log_light":1.3}])
        s=d.groupby("land_id").agg(
            sites=("cell_id_without_year","nunique"),
            years=("year","nunique"))
        self.assertEqual(int(s.sites.loc["12"]),2)
        self.assertEqual(int(s.years.loc["12"]),3)
    def test_source_target_year_not_inferred_from_radiance(self):
        frame=pd.DataFrame({
            "target_id":["2021|1_2|12"],
            "land_id":["12"],"land_class":["island"],
            "is_island":[1],"region5":["region"],"white_pp":[45.0],
            "log_light":[1.1],"log_pop":[0.2],"bio5":[17.0],
            "log_rain":[2.0],"elevation":[20.0],"log_attempts":[2.0],
            "photo_attempts":[6],"accepted_species_n":[3],
            "strict9_success":[4]})
        with self.assertRaisesRegex(ValueError,"original.*2023|2023.*changed"):
            m.year_safeguards(2023,frame)
    def test_common_geographic_support_never_uses_island_only_blocks(self):
        a=pd.DataFrame([
            {"region5":"both","is_island":1,"year":2019},
            {"region5":"both","is_island":0,"year":2019},
            {"region5":"cross_year_only","is_island":1,"year":2019},
            {"region5":"cross_year_only","is_island":0,"year":2021},
            {"region5":"only_island","is_island":1,"year":2019},
            {"region5":"only_island","is_island":1,"year":2025},
            {"region5":"only_mainland","is_island":0,"year":2024}])
        self.assertEqual(m.common_geographic_support_blocks(a),{(2019,"both")})
    def test_generalisability_checks_leave_year_out_and_count_island_concentration(self):
        import inspect
        source=inspect.getsource(m.test_generalizability)
        self.assertIn("for omitted in PHOTO_YEARS",source)
        self.assertIn("largest_island_share_of_island_cell_years",source)
        self.assertIn("common_geographic_support_blocks",source)

    def test_panel_model_has_year_AND_regional_geographic_effects(self):
        import inspect
        s=inspect.getsource(m.run)
        self.assertIn("C(region5) + C(year)",s)
        self.assertIn("C(land_id) + C(year)",s)
        self.assertIn('cluster_fit(allcells,reg,"region5")',s)
        self.assertIn('cluster_fit(within,within_form,"land_id")',s)
if __name__=="__main__":
    unittest.main()
