"""WorldClim raster location, nodata and coarse topography source QA.

All small rasters below are artificial code fixtures; never biological results.
"""
import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin
from scripts.analysis.mois_sample_248312_worldclim_static_and_monthnormals import (
    sample_tiff, raster_terrain, SOURCE_PREVIOUS_SHA, BIO_IDS,
)

def fake(path,grid,nodata=-9999.0):
    grid=np.array(grid,dtype=np.float32)
    with rasterio.open(path,"w",driver="GTiff",height=len(grid),width=len(grid[0]),
                       count=1,crs="EPSG:4326",transform=from_origin(0,40,1/24,1/24),
                       dtype="float32",nodata=nodata) as dst:
        dst.write(grid,1)

def test_original_worldclim_coordinates_include_land_and_missing(tmp_path):
    p=tmp_path/"test.tif"
    fake(p,[[10,20,30],[40,-9999,60],[70,80,90]])
    x=sample_tiff(p,np.array([.02,.06,0.06,400,np.nan]),
                  np.array([39.98,39.98,39.94,39.98,39.9]))
    assert x[0]==10
    assert x[1]==20
    assert np.isnan(x[2])  # Earth coordinates into explicit nodata
    assert np.isnan(x[3]) and np.isnan(x[4])

def test_flat_terrain_has_valid_height_zero_slope_but_unknown_aspect(tmp_path):
    p=tmp_path/"elev.tif"
    fake(p,[[100]*5 for _ in range(5)])
    z=raster_terrain(p,np.array([2/24+.001]),np.array([40-2/24-.001]))
    assert z["wc2_elev_2p5m"][0]==100
    assert z["terrain_slope_2p5m_degree"][0]==pytest.approx(0)
    assert z["terrain_roughness_2p5m_m"][0]==pytest.approx(0)
    assert np.isnan(z["terrain_eastness_2p5m"][0])
    assert np.isnan(z["terrain_northness_2p5m"][0])

def test_terrain_is_missing_if_three_by_three_original_dem_contains_nodata(tmp_path):
    p=tmp_path/"elev.tif"
    a=[[float(i*10+j) for j in range(5)] for i in range(5)]
    a[1][1]=-9999
    fake(p,a)
    z=raster_terrain(p,np.array([2/24+.001]),np.array([40-2/24-.001]))
    assert np.isfinite(z["wc2_elev_2p5m"][0])
    assert np.isnan(z["terrain_slope_2p5m_degree"][0])
    assert np.isnan(z["terrain_TPI_2p5m_m"][0])

def test_original_fcp_frozen_worldclim_inputs_not_redefined():
    assert set(SOURCE_PREVIOUS_SHA)=={"bio","srad"}
    assert len(BIO_IDS)==19
    assert all(len(x)==64 for x in SOURCE_PREVIOUS_SHA.values())
