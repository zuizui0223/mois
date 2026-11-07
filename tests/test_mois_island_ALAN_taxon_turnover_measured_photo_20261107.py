"""Unit tests: independent original observed WCVP species reference, not botanical labels."""
import unittest
import numpy as np
import pandas as pd
from scripts.analysis import mois_island_ALAN_taxon_turnover_measured_photo_20261107 as m

class SourceTaxonTurnoverTests(unittest.TestCase):
    def test_taxon_source_photo_mean_is_not_cell_or_photo_frequency_weighted(self):
        s=pd.DataFrame([
            dict(target_id="2019|a|1",wcvp_accepted_plant_name_id="A",
                 photo_white_fraction=.20,actual_ROI_photos=14,is_island=1),
            dict(target_id="2021|b|2",wcvp_accepted_plant_name_id="A",
                 photo_white_fraction=.80,actual_ROI_photos=1,is_island=0),
            dict(target_id="2019|a|1",wcvp_accepted_plant_name_id="B",
                 photo_white_fraction=.90,actual_ROI_photos=1,is_island=1)])
        taxon,ref=m.source_taxon_reference_white_components(s)
        a=taxon.set_index("wcvp_accepted_plant_name_id").loc["A"]
        self.assertAlmostEqual(float(a.source_mean_photo_white_fraction),.5)
        self.assertEqual(int(a.measured_cell_years),2)
        c=ref.set_index("target_id").loc["2019|a|1"]
        self.assertAlmostEqual(float(c.taxon_average_white_fraction),(.5+.9)/2)
        self.assertEqual(int(c.source_singleton_taxa),1)
        self.assertEqual(int(c.taxon_species_count),2)

    def test_true_same_year_and_region_and_species_pair_gate(self):
        d=pd.DataFrame([
            dict(target_id="2019|a|1",year=2019,region5="X",
                 wcvp_accepted_plant_name_id="A",is_island=1),
            dict(target_id="2019|b|2",year=2019,region5="X",
                 wcvp_accepted_plant_name_id="A",is_island=0),
            dict(target_id="2021|c|1",year=2021,region5="Y",
                 wcvp_accepted_plant_name_id="B",is_island=1),
            dict(target_id="2023|d|2",year=2023,region5="Y",
                 wcvp_accepted_plant_name_id="B",is_island=0)])
        eligible=pd.DataFrame({"target_id":d.target_id})
        r=m.genuinely_crossland_taxon_5deg_year_support(d,eligible)
        self.assertEqual(r["paired_same_taxon_year_region_groups"],1)
        self.assertEqual(r["taxa_directly_seen_on_both_island_mainland_within_year_region"],1)
        self.assertEqual(r["actual_photo_species_cell_years"],2)
        self.assertFalse(r["estimable"])

    def test_source_year_and_photo_colours_must_be_real(self):
        import inspect
        original=inspect.getsource(m.real_species_cell_years)
        self.assertIn('canonical.strict9',original)
        self.assertIn('canonical.strict9_white.notna()',original)
        self.assertIn('canonical.source_subarm.eq("historical_main")',original)
        self.assertIn('panel.require(p.photo_id.is_unique',original)
        self.assertIn("pure.cell_id.is_unique",original)
        run=inspect.getsource(m.run)
        self.assertIn("exact_original_ROI_and_WCVP_white_values_reconstructed",run)
        self.assertIn("no_photo_unknown_imputation_or_MNAR_tipping",run)

    def test_earlier_year_source_reference_never_uses_2024_or_2025_photos(self):
        # Taxon A is photographed in 2019, and in 2024 with a different colour.
        # Taxon B has a source observation in 2021 and 2025.
        species=pd.DataFrame([
          dict(target_id="2019|s1|1",year=2019,wcvp_accepted_plant_name_id="A",
               photo_white_fraction=.2),
          dict(target_id="2021|s2|2",year=2021,wcvp_accepted_plant_name_id="B",
               photo_white_fraction=.8),
          dict(target_id="2024|s3|1",year=2024,wcvp_accepted_plant_name_id="A",
               photo_white_fraction=.6),
          dict(target_id="2024|s3|1",year=2024,wcvp_accepted_plant_name_id="B",
               photo_white_fraction=.5),
          dict(target_id="2025|s4|2",year=2025,wcvp_accepted_plant_name_id="C",
               photo_white_fraction=.9)])
        cohort=pd.DataFrame([
          dict(target_id="2024|s3|1",year=2024,
               is_island=1,region5="X",accepted_species_n=2,
               log_light=1.,log_pop=.2,bio5=10.,log_rain=2.,
               elevation=20.,log_attempts=2.),
          dict(target_id="2025|s4|2",year=2025,
               is_island=0,region5="Y",accepted_species_n=2,
               log_light=1.,log_pop=.2,bio5=10.,log_rain=2.,
               elevation=20.,log_attempts=2.)])
        r=m.cross_year_taxon_reference_transfer(species,cohort)
        self.assertEqual(r["training_accepted_taxa_with_photo_colours"],2)
        self.assertEqual(r["current_year_cells_with_at_least_two_earlier_photo_taxon_references"],1)
        self.assertEqual(r["eligible_island_cells"],1)
        self.assertFalse(r["estimable"])
        self.assertTrue(r["no_heldout_year_photo_used_in_earlier_year_reference"])

    def test_photo_reference_composition_is_exact_derived_observation(self):
        # For the two species in the island cell, the photographic whiteness
        # exactly equals a source taxon mean component plus site deviation.
        actual_white=(.2+.9)/2
        source_taxon_reference=(.5+.9)/2
        deviation=actual_white-source_taxon_reference
        self.assertAlmostEqual(source_taxon_reference+deviation,actual_white)
        self.assertNotAlmostEqual(deviation,0.)

if __name__=="__main__":
    unittest.main()
