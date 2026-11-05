import unittest
import numpy as np
import pandas as pd
from scripts.acquisition.mois_macro_climate import annual

class AnnualClimate(unittest.TestCase):
    def test_leap_year_precipitation_integrates_seconds(self):
        d=pd.date_range('2020-01-01',periods=12,freq='MS')
        self.assertAlmostEqual(annual(np.ones((12,1))/86400,d,'precipitation')[0],366.)
    def test_temperature_and_missingness(self):
        d=pd.date_range('2020-01-01',periods=12,freq='MS');v=np.ones((12,2))*283.15;v[2,1]=np.nan
        r=annual(v,d,'temperature');self.assertAlmostEqual(r[0],10.);self.assertTrue(np.isnan(r[1]))
    def test_cannot_replace_missing_month_with_duplicate(self):
        d=list(pd.date_range('2020-01-01',periods=12,freq='MS'));d[-1]=d[0]
        with self.assertRaises(ValueError):annual(np.ones((12,1)),d,'temperature')

if __name__=='__main__':unittest.main()
