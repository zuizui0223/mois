"""Synthetic QA of official 0.5-degree yearly NPP extraction to GBIF grid."""
import unittest
import numpy as np

from scripts.analysis.mois_gbif_mod17_npp_annual_100km_20261116 import (
    nearest_05degree_indices,nc_sample
)

class TestMOD17NPP(unittest.TestCase):
    def test_regular_geographic_coordinate_in_ascending_and_descending_order(self):
        lon=np.arange(-179.75,180,.5)
        lat=np.arange(89.75,-90,-.5)
        self.assertEqual(len(lon),720)
        self.assertEqual(len(lat),360)
        x=nearest_05degree_indices(lon,np.array([-179.75,0.25,179.75,180.3]))
        y=nearest_05degree_indices(lat,np.array([89.75,-0.25,-89.75,100.]))
        self.assertEqual(x.tolist(),[0,360,719,-1])
        self.assertEqual(y.tolist(),[0,180,359,-1])
    def test_QC_zero_valid_500m_is_not_zero_NPP(self):
        lon=np.arange(-179.75,180,.5)
        lat=np.arange(89.75,-90,-.5)
        npp=np.ones((360,720),dtype=float)*.6
        count=np.ones((360,720),dtype=float)*1000
        ix=360
        iy=180
        npp[iy,ix]=0.   # true ecological zero valid
        count[iy,ix]=10
        count[iy,ix+1]=0 # missing / technically invalid rather than 0 NPP
        coords=(np.array([[.25,.75,180.5], [.25,.75,180.5]]),
                np.array([[-.25,-.25,-.25],[-.25,-.25,-.25]]))
        v,n,q=nc_sample(npp,count,lat,lon,coords)
        self.assertEqual(n.tolist(),[2,0,0])
        self.assertAlmostEqual(v[0],0.)
        self.assertTrue(np.isnan(v[1]) and np.isnan(v[2]))
        self.assertEqual(q[0],10)
    def test_source_qualifies_negative_NPP_and_disallows_masked(self):
        lon=np.arange(-179.75,180,.5)
        lat=np.arange(89.75,-90,-.5)
        npp=np.full((360,720),np.nan)
        q=np.zeros((360,720))
        npp[180,360]=-0.04;q[180,360]=10.
        coords=(np.array([[0.25]]),np.array([[-.25]]))
        result,n,quality=nc_sample(npp,q,lat,lon,coords)
        self.assertAlmostEqual(result[0],-.04)
        self.assertEqual(n[0],1)
if __name__=="__main__":
    unittest.main()
