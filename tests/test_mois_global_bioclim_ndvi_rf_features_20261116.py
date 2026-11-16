import unittest
import numpy as np
from scripts.analysis.mois_global_bioclim_ndvi_rf_features_20261116 import (
    point_mean,raster_probe,require,HUES,VARIABLES
)
class TestGlobalEnvFeatures(unittest.TestCase):
    def test_real_worldclim_all19_and_five_monthly_groups(self):
        self.assertEqual(len(HUES),9)
        self.assertEqual(set(VARIABLES),{"bio","prec","tavg","srad","vapr","wind"})
    def test_point_mean_nodata_kept_missing(self):
        a=np.array([[1.,np.nan,np.nan],[3.,np.nan,0.],[np.nan,np.nan,2.]])
        result,n=point_mean(a)
        self.assertTrue(np.allclose(result[[0,2]],[2,1]))
        self.assertTrue(np.isnan(result[1]))
        self.assertEqual(n.tolist(),[2,0,2])
    def test_real_geotiff_lonlat_northup_validity(self):
        import rasterio
        from rasterio.io import MemoryFile
        from rasterio.transform import from_origin
        with MemoryFile() as mf:
            with mf.open(driver="GTiff",height=2,width=3,count=1,
                         dtype="float32",crs="EPSG:4326",
                         transform=from_origin(-180,90,120,90),
                         nodata=-9999) as ds:
                ds.write(np.array([[1,2,-9999],[5,6,7]],dtype="float32"),1)
            raw=mf.read()
        coords=[(np.array([-179,-59,61,np.nan]),np.array([89,89,89,0]))]
        values=raster_probe(raw,{"id":"synthetic"},coords)
        self.assertTrue(np.array_equal(values[0,:2],[1,2]))
        self.assertTrue(np.isnan(values[0,2]) and np.isnan(values[0,3]))
if __name__=="__main__":unittest.main()
