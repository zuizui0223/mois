import unittest
import numpy as np
import pandas as pd
from scripts.analysis.mois_global_worldclim_hue_20261116 import (
    extract_tif_samples, demean, climate_coverage_status
)
class TestMOISGlobalWorldClim(unittest.TestCase):
    def test_valid_2d_lonlat_raster_indices_and_nodata(self):
        values=np.array([[1.,2.,-9999.],[4.,5.,6.]])
        # West=-180, north=90, res=120 degrees x and 90 degrees y
        lons=np.array([-179,-59,61,179,np.nan])
        lats=np.array([89,89,89,89,0])
        v=extract_tif_samples(values,lons,lats,-180,90,120,-90,nodata=-9999.)
        self.assertEqual(v[:2].tolist(),[1.,2.])
        self.assertTrue(np.isnan(v[2]))
        self.assertTrue(np.isnan(v[3]))
        self.assertTrue(np.isnan(v[4]))
    def test_nodata_does_not_turn_into_zero_rainfall(self):
        values=np.array([[-9999.,0.],[1.,2.]])
        x=extract_tif_samples(values,np.array([-100,100]),np.array([45,45]),
                              -180,90,180,-90,nodata=-9999.)
        self.assertTrue(np.isnan(x[0]))
        self.assertEqual(x[1],0.)
    def test_within_genus_demean_removes_between_genus_white_effect(self):
        d=pd.DataFrame({'x':[1.,0.,1.,1.],'genus':['A','A','B','B']})
        z=demean(d[['x']],d.genus)
        self.assertTrue(np.allclose(z.x,[0.5,-0.5,0.,0.]))
    def test_climate_coverage_missing_not_zero(self):
        self.assertEqual(climate_coverage_status(0,10),'not_sampled_coastal_or_missing')
        self.assertEqual(climate_coverage_status(5,10),'partial')
        self.assertEqual(climate_coverage_status(10,10),'complete')
        self.assertEqual(climate_coverage_status(0,0),'no_eligible_GBIF_grid')
if __name__=='__main__':unittest.main()
