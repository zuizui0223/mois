"""Test incremental genuine-photo WorldClim annual append and physical solar units."""
import numpy as np
import pandas as pd
import pytest

from scripts.analysis import mois_original313535_full19_bio_vapr_terrain_lake as ref
from scripts.analysis import mois_global313535_water_light_heat_habitat_partition as parent
from scripts.analysis.mois_real350399_incremental_all19_climate_20261019 import (
    append_immutable, new_sites_climate,
)


def sample_photos():
    return pd.DataFrame({
        "photo_id": ["10", "11"], "latitude": [35.0, -22.0],
        "longitude": [135.0, 130.0]
    })


def mock_worldclim(root):
    for name in ("bio", "srad", "vapr", "elev"):
        (root/name).mkdir(parents=True, exist_ok=True)
    for i in range(1,20):
        (root/"bio"/f"wc2.1_10m_bio_{i}.tif").touch()
    for i in range(1,13):
        (root/"vapr"/f"wc2.1_10m_vapr_{i:02d}.tif").touch()
        (root/"srad"/f"wc2.1_10m_srad_{i:02d}.tif").touch()
    (root/"elev"/"wc2.1_2.5m_elev.tif").touch()


def fake_sample(path, lon, lat):
    assert len(lon)==len(lat)==2
    if path.parent.name=="srad":
        if path.stem.endswith("_12"):
            return np.array([86.4, np.nan],dtype=float)
        return np.array([86.4,172.8],dtype=float)
    if path.parent.name=="vapr":
        if path.stem.endswith("_07"):
            return np.array([2.0,np.nan],dtype=float)
        return np.array([2.0,3.0],dtype=float)
    return np.array([1.0,2.0],dtype=float)


def fake_terrain(path, lon, lat):
    return {key:np.array([100.,200.]) for key in ref.TERRAIN}


def test_new_source_only_sampling_original_solars_and_missing_months(tmp_path):
    mock_worldclim(tmp_path)
    actual = new_sites_climate(sample_photos(),tmp_path,
                               sampler=fake_sample,
                               terrain_sampler=fake_terrain)
    assert actual.photo_id.tolist()==["10","11"]
    assert len(actual)==2
    assert np.isclose(actual[parent.SOLAR].iloc[0],1.0)
    assert np.isnan(actual[parent.SOLAR].iloc[1])
    assert np.isclose(actual[ref.VAPR].iloc[0],2.0)
    assert np.isnan(actual[ref.VAPR].iloc[1])
    assert len([name for name in actual.columns if name.startswith("wc2_bio")])==19
    assert set(ref.TERRAIN).issubset(actual.columns)


def test_immutable_append_old_original_raster_values_unchanged(tmp_path):
    mock_worldclim(tmp_path)
    fresh=new_sites_climate(sample_photos(),tmp_path,
                            sampler=fake_sample,
                            terrain_sampler=fake_terrain)
    core=fresh.iloc[:1].copy()
    core.loc[:,"photo_id"]=["9"]
    new=fresh.iloc[1:].copy()
    full=append_immutable(core,
            pd.DataFrame({"photo_id":["9"]}),
            pd.DataFrame({"photo_id":["11"]}),
            new,expected_base=1,expected_new=1)
    assert len(full)==2
    assert full.photo_id.tolist()==["9","11"]
    assert full.iloc[:1].reset_index(drop=True).equals(
        core.reset_index(drop=True))


def test_reject_existing_photo_or_changed_abiotic_schema(tmp_path):
    mock_worldclim(tmp_path)
    fresh=new_sites_climate(sample_photos(),tmp_path,
                            sampler=fake_sample,
                            terrain_sampler=fake_terrain)
    same=fresh.iloc[[0]].copy()
    with pytest.raises(ValueError,match="New photo site already"):
        append_immutable(same, pd.DataFrame({"photo_id":["10"]}),
                         pd.DataFrame({"photo_id":["10"]}), same,
                         expected_base=1,expected_new=1)
    dropped=fresh.iloc[[1]].drop(columns=[ref.VAPR])
    with pytest.raises(ValueError,match="ABIOTIC variable"):
        append_immutable(same,pd.DataFrame({"photo_id":["10"]}),
                         pd.DataFrame({"photo_id":["11"]}),dropped,
                         expected_base=1,expected_new=1)


def test_fails_on_missing_month_source_not_infers_weather(tmp_path):
    mock_worldclim(tmp_path)
    (tmp_path/"srad"/"wc2.1_10m_srad_03.tif").unlink()
    with pytest.raises(ValueError,match="monthly VAPR/SRAD"):
        new_sites_climate(sample_photos(),tmp_path,
                          sampler=fake_sample,
                          terrain_sampler=fake_terrain)
