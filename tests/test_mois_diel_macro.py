import unittest
import pandas as pd
from scripts.analysis import mois_diel_macro as m

class SourceContracts(unittest.TestCase):
    def test_ambiguous_alias_is_not_assigned(self):
        names=pd.DataFrame({'source_species_exact_binomial':['Old alba','Old alba','Good alba'], 'wcvp_accepted_plant_name_id':['1','2','3'], 'wcvp_match_state':['unique_exact_name_to_accepted_id']*3})
        accepted=pd.DataFrame({'wcvp_accepted_name':['New alba','Other alba','Good alba'],'wcvp_accepted_plant_name_id':['1','2','3']})
        lookup=m.name_lookup(names,accepted)
        self.assertNotIn('Old alba',lookup)
        self.assertEqual(lookup['Good alba'],'3')
    def test_conflicting_synonyms_are_not_resolved_by_first_match(self):
        self.assertIsNone(m.resolve_species(['Old alba','New alba'],{'Old alba':'1','New alba':'2'}))
        self.assertEqual(m.resolve_species(['Old alba','New alba'],{'Old alba':'1','New alba':'1'}),'1')
    def test_night_share_keeps_both_zero_unknown(self):
        x=m.night_share(pd.Series([0,2,3]),pd.Series([0,2,1]))
        self.assertTrue(pd.isna(x.iloc[0]));self.assertEqual(x.iloc[1],.5);self.assertEqual(x.iloc[2],.75)
    def test_same_year_join_rejects_mismatch(self):
        x=pd.DataFrame({'site_id':['a'],'year':[2020]})
        e=pd.DataFrame({'site_id':['a'],'year':[2019],'radiance':[1]})
        joined=m.join_exposure(x,e)
        self.assertTrue(joined.radiance.isna().all())
    def test_duplicate_exposure_fails(self):
        x=pd.DataFrame({'site_id':['a'],'year':[2020]})
        e=pd.DataFrame({'site_id':['a','a'],'year':[2020,2020],'radiance':[1,2]})
        with self.assertRaises(ValueError):m.join_exposure(x,e)
    def test_species_equal_weight_is_not_photo_weight(self):
        x=pd.DataFrame({'species':['a','a','b']})
        w=m.species_weights(x)
        self.assertEqual(w.tolist(),[.5,.5,1.])

if __name__=='__main__':unittest.main()
