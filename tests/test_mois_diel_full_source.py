import unittest
import pandas as pd
from scripts.acquisition import mois_diel_full_source as m

class FullSourceTests(unittest.TestCase):
    def test_extract_has_no_per_species_cap_or_roi_gate(self):
        d=pd.DataFrame({'observation_uuid':['a','b','c'],'inat_taxon_id':['1','1','2'],'observed_on':['2020-01-01']*3,'latitude':[20]*3,'longitude':[10]*3,'positional_accuracy':[10]*3,'observer_id':['u','u','v']})
        selected,bg=m.select_chunk(d,{'1'})
        self.assertEqual(len(selected),2);self.assertEqual(bg.n_observations.sum(),3)
    def test_wrong_date_and_bad_coordinate_do_not_enter(self):
        d=pd.DataFrame({'observation_uuid':['a','b','c'],'inat_taxon_id':['1']*3,'observed_on':['2020-01-01','2020-99-99','2019-01-01'],'latitude':[20,20,100],'longitude':[10]*3,'positional_accuracy':[10]*3,'observer_id':['u']*3})
        s,b=m.select_chunk(d,{'1'});self.assertEqual(s.observation_uuid.tolist(),['a'])
    def test_missing_accuracy_is_preserved_not_precise(self):
        d=pd.DataFrame({'observation_uuid':['a'],'inat_taxon_id':['1'],'observed_on':['2020-01-01'],'latitude':[20],'longitude':[10],'positional_accuracy':[None],'observer_id':['u']})
        s,b=m.select_chunk(d,{'1'});self.assertFalse(s.coordinate_precision_5km.iloc[0])

if __name__=='__main__':unittest.main()
