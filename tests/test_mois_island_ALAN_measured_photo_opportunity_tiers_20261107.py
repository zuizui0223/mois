"""Tests never infer colour from unclassified ROI and verify quality tier nesting."""
import unittest, inspect
import pandas as pd
from scripts.analysis import mois_island_ALAN_measured_photo_opportunity_tiers_20261107 as m

class OpportunityTierTests(unittest.TestCase):
    def test_source_quality_tiers_nested_and_not_imputed(self):
        self.assertEqual(m.TIERS["strict_5p3r3s"],(5,3,3))
        self.assertEqual(m.TIERS["measured_3p2r2s"],(3,2,2))
        self.assertEqual(m.TIERS["measured_1p1r1s_support_only"],(1,1,1))
        f=pd.DataFrame({
          "target_id":["2021|1_1|2","2021|2_2|3","2021|3_3|4"],
          "photo_attempts":[5,3,1],"strict9_success":[3,2,1],
          "accepted_species_n":[3,2,1],
          "population_status":["OK","OK","OK"],
          "population_density_per_km2":[0.0,2.0,8.0],
          "land_class":["island","mainland","island"],
          "land_id":["2","3","4"],"viirs_valid_months":[4,3,4],
          "viirs_median":[1.3,2.2,0.4],"bio12":[120,120,120],
          "bio5":[17,18,15],"elevation":[15,20,30],
          "white_equal":[0.2,0.5,0.7],
          "query_longitude":[130,130,130],
          "query_latitude":[35,35,35]})
        a=m.verified_colour_frame(2021,f,(5,3,3))
        b=m.verified_colour_frame(2021,f,(3,2,2))
        c=m.verified_colour_frame(2021,f,(1,1,1))
        self.assertEqual((len(a),len(b),len(c)),(1,2,3))
        self.assertEqual(c.loc[c.target_id.eq("2021|1_1|2"),"log_pop"].iloc[0],0)
        self.assertEqual(c.white_pp.tolist(),[20.,50.,70.])
    def test_unknown_no_numerical_colour_is_dropped(self):
        d=pd.DataFrame({
            "target_id":["2023|2_2|4"],"photo_attempts":[2],
            "strict9_success":[1],"accepted_species_n":[1],
            "population_status":["OK"],"population_density_per_km2":[2.0],
            "land_class":["island"],"land_id":["4"],
            "viirs_valid_months":[4],"viirs_median":[0.2],
            "bio12":[120],"bio5":[15],"elevation":[5],
            "white_equal":[None],"query_longitude":[120],
            "query_latitude":[30]})
        self.assertEqual(len(m.verified_colour_frame(2023,d,(1,1,1))),0)
    def test_sameyear_and_geoblock_both_necessary(self):
        df=pd.DataFrame([
           dict(year=2019,region5="A",is_island=1),
           dict(year=2021,region5="A",is_island=0),
           dict(year=2023,region5="B",is_island=1),
           dict(year=2023,region5="B",is_island=0)])
        subset,groups=m.exact_year_block_support(df)
        self.assertEqual(groups,{(2023,"B")})
        self.assertEqual(len(subset),2)
    def test_additional_measured_pairs_separated_from_existing_strict_groups(self):
        # Original strict cohort has one island+mainland pair in 2019/A;
        # newly-measured rows create a new 2021/B paired year×region.
        strict=pd.DataFrame([
          dict(target_id="2019|1_1|1",year=2019,region5="A",is_island=1),
          dict(target_id="2019|1_2|2",year=2019,region5="A",is_island=0),
          dict(target_id="2021|1_3|1",year=2021,region5="B",is_island=1)])
        measured=pd.concat([strict,pd.DataFrame([
          dict(target_id="2021|1_4|2",year=2021,region5="B",is_island=0),
          dict(target_id="2019|1_5|1",year=2019,region5="A",is_island=1)])],
          ignore_index=True)
        result=m.provenance_of_expanded_co_support(strict,measured)
        self.assertEqual(result["old_strict_paired_year_region_groups"],1)
        self.assertEqual(result["new_measured_paired_year_region_groups"],2)
        self.assertEqual(result["newly_paired_year_region_groups"],1)
        self.assertEqual(result["new_original_strict9_measured_added_cell_years"],2)
        self.assertEqual(result["newly_added_cells_in_previously_paired_groups"],1)
        self.assertFalse(result["measured_2species_only_new_paired_blocks"]["estimable"])
    def test_auxiliary_tier_cohort_not_primary_strict_claim(self):
        source=inspect.getsource(m.run)
        self.assertIn("strict_primary_preserved",source)
        self.assertIn("measured_1_species_is_support_counts_only",source)
        self.assertIn("no_2022_provider_400_month_replacement",source)
        fit=inspect.getsource(m.fit_tier)
        self.assertIn('C(year_region5)',fit)
        self.assertIn('cluster_fit(df,f+" + C(region5) + C(year)","region5")',fit)
if __name__=="__main__":
    unittest.main()
