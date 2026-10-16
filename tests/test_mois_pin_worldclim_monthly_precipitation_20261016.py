"""Synthetic 12-GeoTIFF-only WorldClim original archive admission tests."""
from __future__ import annotations
import importlib.util
import json
import zipfile
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]

def load(name,path):
    spec=importlib.util.spec_from_file_location(name,ROOT/path)
    obj=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj

clim=load("mois_real_monthly_rainfall", "scripts/analysis/mois_monthly_rainfall_flower_colour_20261016.py")
m=load("mois_worldclim_zip_admission","scripts/acquisition/mois_pin_worldclim_monthly_precipitation_20261016.py")
raster_test=load("mois_monthly_raster_fixture","tests/test_mois_monthly_rainfall_flower_colour_20261016.py")


def zip_fixture(tmp_path,extra=None):
    _,_,rasters=raster_test.raster_fixture(tmp_path)
    archive=tmp_path/m.ARCHIVE_BASENAME
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED) as z:
        for p in rasters:
            z.write(p,p.name)
        if extra:
            z.writestr(extra,'unexpected')
    return archive


def test_twelve_geotiff_source_sealed_and_can_reuse_monthly_season_classifier(tmp_path):
    archive=zip_fixture(tmp_path)
    out=tmp_path/'sealed_precip'
    result=m.seal(archive,m.sha256_file(archive),out,_synthetic_small=True)
    assert result['status']==m.STATUS
    assert result['n_byte_verified_monthly_rasters']==12
    assert result['synthetic_test_only'] is True
    assert result['worldwide_flowering_season_or_colour_result'] is False
    assert set(result['source_original_geotiff_sha256'])=={f'{i:02}' for i in range(1,13)}
    manifest=out/'monthly_precipitation_manifest.json'
    module=clim.MonthlyNormals(manifest,clim.sha256_file(manifest),test_small=True)
    vals,_=module.extract('40','35')
    assert clim.rainfall_label(vals,10)['label']=='WET'
    with pytest.raises(ValueError,match='Refuse overwrite'):
        m.seal(archive,m.sha256_file(archive),out,_synthetic_small=True)


def test_reject_zip_member_path_escape(tmp_path):
    archive=zip_fixture(tmp_path,extra='../../unsafe.txt')
    with pytest.raises(ValueError,match='Unexpected'):
        m.seal(archive,m.sha256_file(archive),tmp_path/'output',_synthetic_small=True)
    assert not (tmp_path/'output').exists()


def test_wrong_zip_digest_and_incomplete_or_extra_months_are_not_sealed(tmp_path):
    archive=zip_fixture(tmp_path)
    with pytest.raises(ValueError,match='SHA256'):
        m.seal(archive,'f'*64,tmp_path/'out',_synthetic_small=True)
    with zipfile.ZipFile(archive,'a') as z:
        z.writestr('wc2.1_10m_prec_13.tif','extra')
    with pytest.raises(ValueError,match='Duplicate/out-of-range'):
        m.seal(archive,m.sha256_file(archive),tmp_path/'out',_synthetic_small=True)


def test_synthetic_tiny_rasters_cannot_impersonate_production_global_raster(tmp_path):
    archive=zip_fixture(tmp_path)
    with pytest.raises(ValueError,match='10m precip'):
        m.seal(archive,m.sha256_file(archive),tmp_path/'prod')
