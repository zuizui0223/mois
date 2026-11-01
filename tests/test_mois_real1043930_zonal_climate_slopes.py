"""Synthetic within-genus regional climate/photo slope and unsupported polar controls."""
import unittest
import numpy as np
import pandas as pd

from scripts.analysis import mois_real1043930_zonal_climate_slopes as m

class PolarHemispherePhotoSlopesTests(unittest.TestCase):
    def test_geographic_sectors_polar_and_southern_are_not_conflated(self):
        lat=np.asarray([75,60,40,10,-10,-40,-60,-75],float)
        self.assertEqual(m.sector(lat,"north_arctic").tolist(),
                         [True,False,False,False,False,False,False,False])
        self.assertEqual(m.sector(lat,"south_antarctic").tolist(),
                         [False,False,False,False,False,False,False,True])
        self.assertEqual(m.sector(lat,"south_tropical").tolist(),
                         [False,False,False,False,True,False,False,False])
        self.assertEqual(m.sector(lat,"north_extratropics").tolist(),
                         [True,True,True,False,False,False,False,False])

    def test_named_taxon_local_white_slope_removes_genus_cell_mean(self):
        data=[]
        for genus in ("Campanula","Cirsium"):
            for j in range(2):
                data.append(dict(latitude=40.,longitude=140.,
                     photo_id=str(len(data)+1),taxon=str(len(data)+1),
                     genus=genus,source_cohort="TEST",
                     lat5=26,lon5=64,
                     clim_vapour=10.+j,clim_rain=3.+j,
                     clim_elevation=100+j*100,
                     white_fraction=float(j),dominant_white=float(j),
                     entropy9=.2+j*.3,bronze_fraction=.9-j*.6))
        df=pd.DataFrame(data)
        previous=(m.MIN_CELLS,m.MIN_GENUS)
        try:
            m.MIN_CELLS=2;m.MIN_GENUS=1
            z=m.compute_region(df,"north_extratropics",{
                "atmospheric_vapour":1,
                "log_annual_precipitation":1,
                "elevation":100},draws=40)
        finally:
            m.MIN_CELLS,m.MIN_GENUS=previous
        self.assertEqual(len(z),len(m.OUTCOMES)*len(m.FEATURES))
        x=next(v for v in z if v["outcome"]=="white_fraction"
               and v["predictor"]=="atmospheric_vapour")
        self.assertAlmostEqual(x["slope_pp_per_global_source_SD"],100)
        self.assertEqual(x["eligible_taxon_geocells"],4)
        zero=next(v for v in z if v["outcome"]=="bronze_fraction_NOT_brown_petals"
                  and v["predictor"]=="elevation")
        self.assertLess(zero["slope_pp_per_global_source_SD"],0)

    def test_empty_antarctic_biotic_photo_source_has_no_false_estimate(self):
        df=pd.DataFrame({"latitude":np.asarray([],dtype=float),
           "genus":pd.Series([],dtype=str),"taxon":pd.Series([],dtype=str)})
        out=m.compute_region(df,"south_antarctic",{v:1. for v in m.FEATURES})
        self.assertEqual(len(out),12)
        self.assertTrue(all(v["slope_pp_per_global_source_SD"] is None for v in out))
        self.assertTrue(all(v["reason"]=="no_actual_original_photo_ROI_in_geographic_zone"
                            for v in out))

if __name__=="__main__":
    unittest.main()
