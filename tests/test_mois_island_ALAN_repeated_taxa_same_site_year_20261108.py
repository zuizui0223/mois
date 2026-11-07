"""Offline original-source accepted-taxon overlap and repeat-site-year gates."""
import unittest
import pandas as pd
from scripts.analysis import mois_island_ALAN_repeated_taxa_same_site_year_20261108 as m

class GenuinePhotoTaxonSiteYearTests(unittest.TestCase):
    def test_observed_original_taxon_overlap_is_not_interpreted_as_absence(self):
        source=pd.DataFrame([
          dict(target_id="2019|12_13|9",year=2019,land_id="9",
               land_class="island",is_island=1,region5="X",
               white_pp=50.,log_light=.4,log_pop=.2,log_attempts=1.4),
          dict(target_id="2021|12_13|9",year=2021,land_id="9",
               land_class="island",is_island=1,region5="X",
               white_pp=55.,log_light=.5,log_pop=.2,log_attempts=1.5),
          dict(target_id="2019|14_15|10",year=2019,land_id="10",
               land_class="mainland",is_island=0,region5="X",
               white_pp=41.,log_light=.3,log_pop=.2,log_attempts=1.4),
          dict(target_id="2024|14_15|10",year=2024,land_id="10",
               land_class="mainland",is_island=0,region5="X",
               white_pp=44.,log_light=.6,log_pop=.2,log_attempts=1.8)
        ])
        t="wcvp_accepted_plant_name_id"
        rows=[]
        for target,taxa_list in [
            ("2019|12_13|9",["A","B"]),
            ("2021|12_13|9",["B","C"]),
            ("2019|14_15|10",["D"]),
            ("2024|14_15|10",["E"])]:
            for species in taxa_list:
                rows.append(dict(target_id=target,
                     wcvp_accepted_plant_name_id=species,
                     photo_white_fraction=.3,actual_ROI_photos=1,
                     year=int(target[:4]),region5="X",
                     is_island=int(target.endswith("|9")),
                     land_id=target.split("|")[-1],
                     log_light=.5,log_pop=.2,log_attempts=1.4))
        sp=pd.DataFrame(rows)
        cells,all_taxa=m.true_site_species_years(sp,source)
        _,repeat=m.same_site.repeat_site_frame(cells)
        overlaps=m.observed_taxon_overlap(repeat,all_taxa)
        self.assertEqual(len(overlaps),2)
        self.assertAlmostEqual(overlaps.loc[overlaps.land_class.eq("island"),"jaccard"].iloc[0],1/3)
        self.assertTrue(bool(overlaps.loc[overlaps.land_class.eq("mainland"),"no_same_taxon_both_photographic_years"].iloc[0]))
        r=m.species_site_repeats(all_taxa)
        self.assertEqual(set(r[t]),{"B"})
        self.assertEqual(len(r),2)
        self.assertEqual(int(r.species_repeat_site_id.nunique()),1)
        self.assertTrue(r.white_pp.eq(30).all())
    def test_nonrepeated_source_taxa_must_never_influence_year_change_fit(self):
        d=pd.DataFrame([
            dict(species_repeat_site_id="A|s",year=2019,region5="X",is_island=1,
                 land_id="9",log_light=.5,log_pop=.1,log_attempts=2.,
                 photo_white_fraction=.5,white_pp=50.),
            dict(species_repeat_site_id="A|s",year=2024,region5="X",is_island=1,
                 land_id="9",log_light=.7,log_pop=.2,log_attempts=2.,
                 photo_white_fraction=.4,white_pp=40.)])
        out=m.fit_within_real_species_site(d)
        self.assertFalse(out["estimable"])
        self.assertIn("support",out["reason"])
    def test_reference_is_observed_same_species_not_new_plant_colour(self):
        import inspect
        s=inspect.getsource(m.run)
        self.assertIn('"observation_taxon_Jaccard_is_not_field_plant_turnover":True',s)
        self.assertIn('"no_technical_unknown_colour_label_imputation":True',s)
        self.assertIn('"no_MNAR_tipping_or_IPW_scenarios":True',s)
        self.assertIn('"no_new_photographs_or_ROI_inference":True',s)
        self.assertEqual(m.EXPECTED,{"strict_5p3r3s":7982,
                                     "measured_3p2r2s":20543})

if __name__=="__main__":
    unittest.main()
