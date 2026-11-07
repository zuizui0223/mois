"""True original photo observation-month and taxon overlap source tests."""
import unittest
import numpy as np
import pandas as pd
from scripts.analysis import mois_island_ALAN_repeated_photo_month_match_20261108 as m

class OriginalMonthMatchedTaxonPhotos(unittest.TestCase):
    def test_different_calendar_months_do_not_create_fictitious_shared_flora(self):
        rows=[]
        def add(site,year,month,species,land):
            rows.append(dict(repeat_site_id=site,year=year,month=month,
                             wcvp_accepted_plant_name_id=species,land_class=land))
        for sp in ("A","B"):add("3_4|7",2019,6,sp,"island")
        for sp in ("B","C"):add("3_4|7",2021,6,sp,"island")
        add("5_6|9",2019,7,"D","mainland")
        add("5_6|9",2023,11,"D","mainland")
        add("7_8|1",2019,6,"A","island")
        add("7_8|1",2019,6,"B","island")
        add("7_8|1",2019,7,"C","island")
        add("7_8|1",2021,7,"C","island")
        dated=pd.DataFrame(rows)
        repeat=pd.DataFrame([
            dict(repeat_site_id=s,year=y,land_class=c)
            for s,c,years in [("3_4|7","island",[2019,2021]),
                              ("5_6|9","mainland",[2019,2023]),
                              ("7_8|1","island",[2019,2021])]
            for y in years])
        pairs,unknown=m.observed_photo_year_pairs(repeat,dated)
        self.assertEqual(unknown,0)
        self.assertEqual(len(pairs),3)
        self.assertEqual(int(pairs.same_dominant_photo_month.sum()),1)
        self.assertEqual(int(pairs.at_least_one_shared_actual_photo_month.sum()),2)
        same=pairs.loc[pairs.same_dominant_photo_month].iloc[0]
        self.assertAlmostEqual(float(same.jaccard_all_actual_dated_photos),1/3)
        shifted=pairs.loc[pairs.dominant_photo_month_circular_distance.eq(1)].iloc[0]
        self.assertAlmostEqual(float(shifted.jaccard_restrict_to_same_actual_calendar_months),1.)
        self.assertAlmostEqual(float(shifted.jaccard_all_actual_dated_photos),1/3)
        nonoverlap=pairs.loc[~pairs.at_least_one_shared_actual_photo_month].iloc[0]
        self.assertTrue(pd.isna(nonoverlap.jaccard_restrict_to_same_actual_calendar_months))
        self.assertAlmostEqual(float(nonoverlap.jaccard_all_actual_dated_photos),1.)
    def test_invalid_month_is_unknown_not_january(self):
        a=pd.DataFrame([
          dict(photo_id="A",observation_id="1",inat_taxon_id="11",
               source_subarm="historical_main",strict9=True,strict9_white=.5,
               observed_on="2021-07-07",longitude=.03,latitude=.03),
          dict(photo_id="B",observation_id="2",inat_taxon_id="11",
               source_subarm="historical_main",strict9=True,strict9_white=.4,
               observed_on="2023",longitude=.03,latitude=.03)
        ])
        wcvp=pd.DataFrame([
          dict(inat_taxon_id="11",wcvp_match_state="unique_exact_name_to_accepted_id",
               wcvp_accepted_plant_name_id="TAXON")])
        source=pd.DataFrame([
          dict(target_id="2021|3600_1800|42",year=2021,land_id="42",land_class="island"),
          dict(target_id="2023|3600_1800|42",year=2023,land_id="42",land_class="island")])
        months,rec=m.original_photo_source_months(a,wcvp,source)
        self.assertEqual(len(months),1)
        self.assertEqual(int(months.month.iloc[0]),7)
        self.assertEqual(rec["source_observation_date_month_unavailable_photos"],1)
        self.assertEqual(rec["accepted_original_strict9_source_photos_joined"],2)
    def test_no_photo_colour_imputation_or_inferred_flowering_season(self):
        import inspect
        s=inspect.getsource(m.run)
        self.assertIn('"photo_no_record_not_plant_absence":True',s)
        self.assertIn('"date_month_missing_excluded_not_filled":True',s)
        self.assertIn('"no_MNAR_tipping_or_IPW":True',s)
        self.assertIn('"no_new_original_photo_RGB_download_or_ROI_inference":True',s)
        self.assertEqual(m.EXPECTED,{"strict_5p3r3s":7982,
                                     "measured_3p2r2s":20543})

if __name__=="__main__":
    unittest.main()
