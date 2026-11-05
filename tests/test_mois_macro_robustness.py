import unittest
import pandas as pd
from scripts.analysis.mois_macro_robustness import remove_genus

class GenusRemoval(unittest.TestCase):
    def test_counts_and_support_are_recomputed_after_removal(self):
        d=pd.DataFrame({'cell025':['a','b'],'year':[2020,2020],'n_species':[12,10],'n_white_species':[5,3],'n_genera':[6,5],'radiance':[1,1],'population_density':[1,1],'precipitation_mm':[1,1],'n_precise_observations':[10,10]})
        g=pd.DataFrame({'cell025':['a','b'],'year':[2020,2020],'genus':['G','G'],'n_species':[2,1],'n_white_species':[1,1]})
        v=remove_genus(d,g,'G')
        self.assertEqual(v.cell025.tolist(),['a'])
        self.assertEqual(v.iloc[0].n_genera,5)
        self.assertEqual(v.iloc[0].white_species_fraction,.4)

if __name__=='__main__':unittest.main()
