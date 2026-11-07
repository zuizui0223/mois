"""MOIS observed photo-opportunity × accepted species factorial invariants."""
import unittest
import pandas as pd
from scripts.analysis import mois_island_ALAN_photo_opportunity_species_factorial_20261107 as m

class PhotoOpportunityFactorialTests(unittest.TestCase):
    def test_true_2x2_design_thresholds(self):
        self.assertEqual(m.DESIGN,{
          "strict_5photos_3species":(5,3,3),
          "photos_only_3photos_3species":(3,3,3),
          "species_only_5photos_2species":(5,2,2),
          "both_3photos_2species":(3,2,2)})
        for photo,suc,species in m.DESIGN.values():
            self.assertGreaterEqual(suc,species)
    def test_original_accepted_species_cannot_outnumber_strict9(self):
        d=pd.DataFrame([
          {"year":2019,"target_id":"2019|1|1","is_island":1,
           "land_id":"1","region5":"r","accepted_species_n":3,
           "strict9_success":2}])
        with self.assertRaisesRegex(ValueError,"outnumber"):
            m.diagnostic_counts(d)
    def test_observed_photo_and_species_strata_bins(self):
        d=pd.DataFrame([
          dict(year=2019,region5="A",photo_attempts=3,accepted_species_n=2,
               is_island=1,target_id="2019|1|1"),
          dict(year=2019,region5="A",photo_attempts=4,accepted_species_n=2,
               is_island=0,target_id="2019|2|2"),
          dict(year=2021,region5="A",photo_attempts=5,accepted_species_n=3,
               is_island=1,target_id="2021|3|3"),
          dict(year=2021,region5="A",photo_attempts=19,accepted_species_n=3,
               is_island=0,target_id="2021|4|4")])
        shared,info=m.exact_observed_opportunity_overlap(d)
        self.assertEqual(info["sameyear_region_photo_species_paired_strata"],1)
        self.assertEqual(info["balanced_pairs_possible_from_real_strata"],1)
        self.assertEqual(len(shared),2)
        self.assertSetEqual(set(shared.year),{2019})
        self.assertEqual(shared.attempt_bin.astype(str).tolist(),["3-4","3-4"])
    def test_unsupported_matched_island_photo_opportunities_are_explicit(self):
        d=pd.DataFrame([
          dict(year=2019,region5="A",photo_attempts=3,accepted_species_n=2,
               is_island=1,target_id="2019|1|1"),
          dict(year=2019,region5="A",photo_attempts=3,accepted_species_n=2,
               is_island=0,target_id="2019|2|2")])
        estimate=m.matched_opportunity_estimate(d)
        self.assertFalse(estimate["estimable"])
        self.assertEqual(estimate["balanced_pairs_possible_from_real_strata"],1)
    def test_only_both_relaxed_observations_are_disjoint_original_photos(self):
        source=pd.DataFrame([
           dict(target_id="2019|s|1",year=2019,is_island=1,
                white_pp=20.,viirs_median=1.,photo_attempts=5,accepted_species_n=3),
           dict(target_id="2019|p|1",year=2019,is_island=1,
                white_pp=30.,viirs_median=2.,photo_attempts=3,accepted_species_n=3),
           dict(target_id="2019|c|2",year=2019,is_island=0,
                white_pp=40.,viirs_median=3.,photo_attempts=5,accepted_species_n=2),
           dict(target_id="2019|j|2",year=2019,is_island=0,
                white_pp=50.,viirs_median=4.,photo_attempts=3,accepted_species_n=2)])
        byid=lambda ids: source.loc[source.target_id.isin(ids)].copy()
        a=m.disjoint_original_measured_source_classes({
            "strict_5photos_3species":byid(["2019|s|1"]),
            "photos_only_3photos_3species":byid(["2019|s|1","2019|p|1"]),
            "species_only_5photos_2species":byid(["2019|s|1","2019|c|2"]),
            "both_3photos_2species":source})
        self.assertEqual({k:v["cell_years"] for k,v in a.items()},{
            "strict_original":1,
            "photo_requirement_only_relaxed":1,
            "species_requirement_only_relaxed":1,
            "joint_photo_and_species_relaxation_only":1})
        self.assertEqual(a["joint_photo_and_species_relaxation_only"]["mainland"]["white_ROI_mean_pp"],50.)
    def test_not_mnar_and_not_photo_reinference(self):
        import inspect
        s=inspect.getsource(m.run)
        self.assertIn('"no_MNAR_tipping_or_IPW_scenarios":True',s)
        self.assertIn('"no_original_photo_RGB_downloads_or_reinference":True',s)
        self.assertIn('"strict_primary_reference_unchanged":True',s)

if __name__=="__main__":
    unittest.main()
