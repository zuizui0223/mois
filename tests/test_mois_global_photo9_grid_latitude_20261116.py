import unittest
import numpy as np
import pandas as pd
from scripts.analysis.mois_global_photo9_grid_latitude_20261116 import (
    lat_band_floor,within_group,photo_mean_top,BANDS
)
class TestGlobalPhotoHuesGrid(unittest.TestCase):
    def test_latitude_band_half_open_including_negative_lat(self):
        self.assertEqual(lat_band_floor(35.0),30)
        self.assertEqual(lat_band_floor(-32.1),-40)
        self.assertEqual(lat_band_floor(-0.1),-10)
        self.assertEqual(lat_band_floor(0.0),0)
        self.assertEqual(lat_band_floor(89.0),80)
        self.assertEqual(len(BANDS),18)
    def test_genus_effect_demeaning_not_biased_by_between_genus_colour(self):
        d=pd.DataFrame({'genus':['A','A','B','B'],'x':[1.,0.,1.,1.]})
        r=within_group(d[['x']],d.genus)
        self.assertTrue(np.allclose(r.x,[0.5,-0.5,0,0]))
    def test_modal_vs_mean_uses_real_numeric_nine(self):
        self.assertEqual(photo_mean_top([0.2,0.8,0,0,0,0,0,0,0]),'yellow')
        with self.assertRaises(ValueError):photo_mean_top([1]+[0]*7)
    def test_missing_latitude_stays_nan_not_tropical(self):
        self.assertTrue(np.isnan(lat_band_floor(np.nan)))
if __name__=='__main__':unittest.main()
