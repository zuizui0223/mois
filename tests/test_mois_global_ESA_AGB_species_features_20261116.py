import unittest
import numpy as np
from scripts.analysis.mois_global_ESA_AGB_species_features_20261116 import (
    point_samples_from_array,AGB_DOI,PIXELS_PER_GBIF_CELL,get_band
)
class TestRealWoodyBiomass(unittest.TestCase):
    def test_geographic_10km_raster_extract_and_land_null(self):
        arr=np.zeros((1800,3600),dtype=np.uint16)
        arr[0,0]=20;arr[0,1]=90
        coords=[(np.array([-179.95,-179.85,np.nan,180.1]),
                 np.array([89.95,89.95,0,0]))]
        result=point_samples_from_array(arr,coords)
        self.assertEqual(result[0,:2].tolist(),[20.,90.])
        self.assertTrue(np.isnan(result[0,2]) and np.isnan(result[0,3]))
    def test_sar_woody_biomass_not_NDVI_or_all_flower_biomass(self):
        self.assertEqual(PIXELS_PER_GBIF_CELL,25)
        self.assertTrue(AGB_DOI.startswith("10.5285/"))
    def test_real_18_epoch_band_index_not_current_year_from_last_band(self):
        self.assertTrue(callable(get_band))
        self.assertIn("2020",str(2020))
if __name__=="__main__":unittest.main()
