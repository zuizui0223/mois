"""Small synthetic functional tests for the real-input MOIS LGM pipeline."""
import tempfile
import unittest
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis.mois_lgm_climate_history_species_20261116 import (
    layer_number, paired_change, extract_at_points, extract_two_bioclim,
    genus_model,
)

class TestLGMClimate(unittest.TestCase):
    def test_worldclim_bi_names(self):
        for name in ["bio1.bil", "bio12.bil", "cclgmbi1.tif",
                     "cclgmbi12.tif", "mrlgmbi_12.tif", "melgmbi1_10m.tif"]:
            self.assertIn(layer_number(name), {1,12},name)
        self.assertEqual(layer_number("cclgmbi2.tif"),2)
        self.assertIsNone(layer_number("pr_10m.tif"))

    def test_extract_layer_and_bil_header(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"example.zip"
            with zipfile.ZipFile(p,"w") as z:
                for k in (1,12):
                    z.writestr(f"bio{k}.bil", b"\x00\x00")
                    z.writestr(f"bio{k}.hdr", "BYTEORDER I\n")
                z.writestr("bio5.bil", b"\x00\x00")
            found=extract_two_bioclim(p,Path(d)/"contents","test")
            self.assertEqual(set(found),{1,12})
            self.assertTrue(found[1][0].with_suffix(".hdr").exists())
            self.assertTrue(found[12][0].exists())

    def test_paired_differences_do_not_treat_missing_as_zero(self):
        ct=np.array([[10,14,1],[20,18,5]],float)
        cp=np.array([[100,200,0],[120,300,0]],float)
        pt=np.array([[5,np.nan,0],[12,8,0]],float)
        pp=np.array([[50,0,0],[60,120,np.nan]],float)
        n,t,p=paired_change((ct,cp),(pt,pp))
        # 0 mm of precipitation remains a real observation, not missing.
        self.assertEqual(n.tolist(),[2,1,1])
        self.assertAlmostEqual(t[0],6.5)
        self.assertAlmostEqual(t[1],10.)
        self.assertAlmostEqual(t[2],1.)
        self.assertAlmostEqual(p[2],0.)
        # A separate pair of wholly missing columns should become NaN.
        bad=np.array([[np.nan],[np.nan]])
        miss, dt, dp=paired_change((bad,bad),(bad,bad))
        self.assertEqual(miss.tolist(),[0])
        self.assertTrue(np.isnan(dt[0]) and np.isnan(dp[0]))

    def test_fast_raster_index_and_outside_missing(self):
        # 4x4 synthetic raster, x -2..2, y 2..-2, north-up
        arr=np.arange(16,dtype=float).reshape(4,4)
        coords=(np.array([[0.,-1.,10.]]),np.array([[1.,-1.,1.]]))
        got=extract_at_points(arr,coords,(-2.,2.,1.,-1.))
        self.assertEqual(got.shape,(1,3))
        self.assertEqual(got[0,0],6)
        self.assertEqual(got[0,1],13)
        self.assertTrue(np.isnan(got[0,2]))

    def test_genus_cluster_effect_and_units(self):
        rng=np.random.default_rng(20261116)
        n=640
        g=np.repeat([f"G{i}" for i in range(80)],8)
        climate=rng.normal(size=n)
        history=rng.normal(size=n)
        d=pd.DataFrame({
            "genus":g,
            "source_hue9":np.where(history+0.6*rng.normal(size=n)>0,"white","yellow"),
            "hue_status":"sampling_consensus_single",
            "BIO01":15+4*climate,
            "BIO12":np.exp(6+0.3*climate),
            "MODIS_NDVI_2001_2021_mean":rng.uniform(.2,.85,n),
            "ESA_AGB_woody_Mg_ha_2020":rng.uniform(1,100,n),
            "mean_abs_latitude_over_occupied_cells":rng.uniform(10,70,n),
            "records_passed_spatial_QA":rng.integers(10,2000,n),
            "threeGCM_delta_temperature_C_median":history+5
        })
        a=genus_model(d,["threeGCM_delta_temperature_C_median"],equal_genus=True)
        self.assertEqual(a["n_genera"],80)
        self.assertGreater(a["historic_effects"]["threeGCM_delta_temperature_C_median"]["beta_per_species_SD"],0.)
        self.assertEqual(a["species_equal_vs_genus_equal"],"genus_equal")

if __name__=="__main__":
    unittest.main()
