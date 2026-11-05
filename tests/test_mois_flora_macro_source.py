import unittest
import pandas as pd
from scripts.acquisition import mois_flora_macro_source as m

class BroadSourceContracts(unittest.TestCase):
    def test_multistate_not_forced_white_or_nonwhite(self):
        d=pd.DataFrame({'accepted_species':['A alba','B alba','C alba'],'normalized_value':['white','white|red_pink','red_pink'],'trait_name':['flower_primary_color']*3,'quality':['high']*3,'resolution_status':['resolved']*3,'classification':['single']*3})
        self.assertEqual(m.direct_colours(d).accepted_species.tolist(),['A alba','C alba'])
    def test_photo_repetition_does_not_inflate_species_richness(self):
        d=pd.DataFrame({'cell025':['a']*4,'year':[2020]*4,'wcvp_id':['1','1','2','2'],'white':[1,1,0,0],'genus':['A','A','B','B']})
        c,g=m.aggregate_occurrence_species(d)
        self.assertEqual(c.n_species.iloc[0],2);self.assertEqual(c.n_white_species.iloc[0],1)
    def test_conflicting_species_colour_fails(self):
        d=pd.DataFrame({'cell025':['a','a'],'year':[2020,2020],'wcvp_id':['1','1'],'white':[1,0],'genus':['A','A']})
        with self.assertRaises(ValueError):m.aggregate_occurrence_species(d)

if __name__=='__main__':unittest.main()
