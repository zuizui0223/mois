"""Synthetic data contracts only: never botanical field evidence."""
import unittest
import pandas as pd
import numpy as np
from scripts.analysis import mois_viirs_independent_white_2020 as s

class TestSourceDirectWhite(unittest.TestCase):
    def test_direct_species_single_state_only(self):
        d=pd.DataFrame({'accepted_species':['a','b','c','d','e','f'],
            'trait_name':['flower_primary_color']*6,
            'normalized_value':['white','red_pink','white|blue_purple','white','white','other_described'],
            'quality':['high','medium','high','low','high','high'],
            'resolution_status':['resolved']*6,
            'classification':['species_direct']*5+['cultivar_contamination']})
        x=s.exact_direct_species_labels(d)
        self.assertEqual(x.accepted_species.tolist(),['a','b','e'])
        self.assertEqual(x.white_direct.tolist(),[1.,0.,1.])
    def test_one_accepted_id_multiple_nominal_synonyms(self):
        name=pd.DataFrame({'inat_taxon_id':['1','2','3'],
            'wcvp_match_state':['unique_exact_name_to_accepted_id']*2+['multiple_accepted_ids_ambiguous'],
            'wcvp_accepted_plant_name_id':['A','A','B']})
        acc=pd.DataFrame({'wcvp_accepted_plant_name_id':['A','B'],
                          'wcvp_accepted_name':['Genusa alpha','Genusb beta']})
        direct=pd.DataFrame({'accepted_species':['Genusa alpha','Genusb beta'],'white_direct':[1.,0.]})
        mapping=s.source_taxonomy(name,acc,direct)
        self.assertEqual(len(mapping),2)
        self.assertEqual(set(mapping.wcvp_accepted_plant_name_id),{'A'})
    def test_paired_support_requires_three_documented_species_both_sides(self):
        pairs=pd.DataFrame({'bright_cell_id':['x','z'],'dark_cell_id':['y','w'],
            'match_block':['block','block'],'white_diff_bright_minus_dark_pp':[4.,8.]})
        cells=pd.DataFrame({'cell_id':['x','y','z','w'],
            'n_documented_species':[3,5,2,6],'direct_white_share':[.7,.3,.9,.2]})
        out=s.paired_direct(pairs,cells)
        self.assertEqual(len(out),1)
        self.assertAlmostEqual(float(out.direct_white_difference_pp.iloc[0]),40)
    def test_no_reuse_fails(self):
        pairs=pd.DataFrame({'bright_cell_id':['x','x'],'dark_cell_id':['y','z'],
                            'match_block':['a','a']})
        with self.assertRaisesRegex(ValueError,'reused'):
            s.paired_direct(pairs,pd.DataFrame({'cell_id':['x','y','z'],
                  'n_documented_species':[5]*3,'direct_white_share':[.2]*3}))
    def test_block_bootstrap_reproducibility_and_low_support(self):
        d=pd.DataFrame({'match_block':['a']*10+['b']*10,'effect':[2.]*20})
        self.assertEqual(s.block_bootstrap(d,'effect',iters=100),[2.,2.])
        self.assertIsNone(s.block_bootstrap(d.iloc[:19],'effect'))
    def test_exact_six_signed_sources(self):
        self.assertEqual(set(s.PINS),{'canonical','wcvp_lookup','wcvp_names',
            'island_direct','candidate_cells','matched_pairs'})
        self.assertTrue(all(len(x)==64 for x in s.PINS.values()))
if __name__=='__main__':unittest.main()
