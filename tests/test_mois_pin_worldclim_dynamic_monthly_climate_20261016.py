"""Synthetic-only 48-raster original WorldClim T/S/P/vapour source admission.

Never downloads provider media, botanical photo pixels or published original
source rasters. Confirm strictly original scales, 12 months, SHA and grids.
"""
from __future__ import annotations
import importlib.util
import json
import zipfile
from pathlib import Path
import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location(
    "mois_wc_dynamic",ROOT/"scripts/acquisition/mois_pin_worldclim_dynamic_monthly_climate_20261016.py")
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def fixture(tmp, *, wrong_grid=False, bad_extra=None):
    originals={}
    stamps={}
    for k in m.VARIABLES:
        zpath=tmp/(m.HEADER+k+".zip")
        with zipfile.ZipFile(zpath,"w",compression=zipfile.ZIP_DEFLATED) as arch:
            arch.writestr("readme.txt","synthetic, not an original WorldClim raster")
            for month in range(1,13):
                default={"tavg":17+month*.2,"srad":12000+month*100,
                         "prec":25+month,"vapr":1.2+month*.01}[k]
                p=tmp/f"{k}_{month}.tif"
                transform=from_origin(-180,90,60,45)
                if wrong_grid and month==12 and k=="tavg":
                    transform=from_origin(-179,90,60,45)
                with rasterio.open(p,"w",driver="GTiff",count=1,width=6,height=4,
                                   crs="EPSG:4326",transform=transform,
                                   dtype="float32",nodata=-9999) as ds:
                    ar=np.full((4,6),default,dtype="float32")
                    ar[3,0]=-9999
                    ds.write(ar,1)
                arch.write(p,f"wc2.1_10m_{k}_{month}.tif")
            if bad_extra and k=="srad":
                arch.writestr(bad_extra,"wrong")
        originals[k]=zpath
        stamps[k]=m.sha256_file(zpath)
    return originals,stamps


def test_four_complete_original_climate_axes_share_grid_and_are_sha_pinned(tmp_path):
    archives,pins=fixture(tmp_path)
    result=m.validate_originals(archives,pins,tmp_path/"verified",synthetic=True)
    assert result["status"]==m.STATUS
    assert result["n_georeferenced_original_raster_layers"]==48
    assert result["climate_grid_shape"]==[6,4]
    assert result["iNaturalist_observation_date_joined"] is False
    assert result["biological_photo_pixels_measured"] is False
    assert result["global_white_fraction_estimated"] is False
    manifest=tmp_path/"verified/dynamic_monthly_manifest.json"
    d=json.loads(manifest.read_text())
    assert d["units"]=={
        "tavg":"degC","srad":"kJ m-2 day-1",
        "prec":"mm/month","vapr":"kPa"
    }
    assert d["n_original_monthly_rasters"]==48
    assert result["manifest_sha256"]==m.sha256_file(manifest)
    assert all(len(d["months"][k])==12 for k in m.VARIABLES)
    with pytest.raises(ValueError,match="Cannot overwrite"):
        m.validate_originals(archives,pins,tmp_path/"verified",synthetic=True)


def test_wrong_grid_hard_fails_before_promoting_climate_source(tmp_path):
    archives,pins=fixture(tmp_path,wrong_grid=True)
    with pytest.raises(ValueError,match="not identical"):
        m.validate_originals(archives,pins,tmp_path/"out",synthetic=True)
    assert not (tmp_path/"out").exists()


def test_malicious_archive_entries_or_mismatched_original_sha_hard_fail(tmp_path):
    archives,pins=fixture(tmp_path,bad_extra="../../wrong.jpg")
    with pytest.raises(ValueError,match="unsafe or non-monthly"):
        m.validate_originals(archives,pins,tmp_path/"out",synthetic=True)
    assert not (tmp_path/"out").exists()
    archives,pins=fixture(tmp_path/"other") if False else (archives,pins)
    bad=dict(pins,tavg="a"*64)
    with pytest.raises(ValueError,match="archive SHA256 mismatch"):
        m.validate_originals(archives,bad,tmp_path/"otherout",synthetic=True)
