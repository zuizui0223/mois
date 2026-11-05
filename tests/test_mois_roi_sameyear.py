import unittest
import pandas as pd
import numpy as np
from scripts.analysis.mois_roi_sameyear_supplement import aggregate

class OriginalEstimand(unittest.TestCase):
    def test_repeated_photos_do_not_change_taxon_balanced_mean(self):
        d=pd.DataFrame({'cell025':['1_2']*11,'year':[2020]*11,'observation_id':range(11),'inat_taxon_id':['A']*9+['B','unknown'],'strict9':[True]*10+[False],'strict9_white':[0.]*9+[1.,np.nan]})
        r=aggregate(d).iloc[0]
        self.assertEqual(r.roi_white_fraction,.5)
        self.assertEqual(r.roi_white_photo_weighted,.1)
        self.assertEqual(r.n_roi_taxa,2)
        self.assertEqual(r.n_valid_roi,10)

if __name__=='__main__':unittest.main()
