"""Tiny source-integrity tests: real input contract, not mock ecological inference."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from scripts.analysis import mois_full350399_real_site_annual_climate as target
from scripts.analysis import mois_original313535_full19_bio_vapr_terrain_lake as original
from scripts.analysis import mois_global313535_water_light_heat_habitat_partition as climate


def row(pid, morph, arm):
    return dict(photo_id=str(pid), inat_taxon_id=str(100 + pid),
                species="Campanula punctata", morph=morph,
                observation_id=f"obs-{pid}", observed_on="2025-06-15",
                latitude=32., longitude=131.,
                classified=morph in ("white", "red_pink", "blue_purple", "yellow_orange"),
                white=morph=="white", arm=arm, technical_failure=None)


def fixed(tmp_path, monkeypatch):
    for name, value in (("N_CORE", 3), ("N_NEW", 4), ("N_TOTAL", 7),
                        ("N_CLASSIFIED", 5), ("N_WHITE", 3),
                        ("OLD_CLASSIFIED", 2), ("OLD_WHITE", 1)):
        monkeypatch.setattr(target, name, value)
    core=pd.DataFrame([row(1,"white","core"),row(2,"blue_purple","core"),
                       row(3,"mixed_uncertain","core")])
    extras=[
        pd.DataFrame([row(4,"white","broad")]),
        pd.DataFrame([row(5,"yellow_orange","wave2")]),
        pd.DataFrame([row(6,"mixed_uncertain","cached8192")]),
        pd.DataFrame([row(7,"white","tail16384")]),
    ]
    monkeypatch.setattr(target, "assemble",
                        lambda *args,**kw: (
                            pd.concat((core,*extras),ignore_index=True), [core,*extras]))
    input_files=[]
    for i in range(5):
        f=tmp_path / f"input_{i}.parquet";f.write_bytes(b"original"+bytes([i]))
        input_files.append(f)
    original_env=pd.DataFrame({"photo_id":core.photo_id})
    for k in (*original.BIO, original.VAPR, climate.SOLAR, *original.TERRAIN):
        original_env[k]=np.array([20.,21.,np.nan],dtype=np.float32)
    path=tmp_path/"oldenv.parquet"
    original_env.to_parquet(path,index=False)
    raster=tmp_path/"worldclim"
    source_digests={}
    for name, file in target.RASTER_FILES.items():
        folder=raster/name;folder.mkdir(parents=True,exist_ok=True)
        f=folder/file;f.write_bytes(name.encode())
        source_digests[name]=target.sha256(f)
    monkeypatch.setattr(original,"PIN_BIO",source_digests["bio"])
    monkeypatch.setattr(original,"PIN_SRAD",source_digests["srad"])
    receipt={
        "schema":"mois_real317631_fifth4096_incremental_RGB_annual_bio19_terrain_v1",
        "all_original_photo_attempts":3,
        "all_classified_four_hue_photo_ROIs":2,
        "all_white_photo_ROIs":1,
        "photo_ID_aligned_Zstd_files_sha256":{path.name:target.sha256(path)},
        "source_original_raster_sha256":source_digests,
    }
    doc=tmp_path/"receipt.json";doc.write_text(json.dumps(receipt))
    def synthetic_sample(photos, _):
        d=pd.DataFrame({"photo_id":photos.photo_id.to_numpy()})
        for k in (*original.BIO,original.VAPR,climate.SOLAR,*original.TERRAIN):
            d[k]=np.full(len(photos),11,dtype=np.float32)
        return d
    return input_files,path,doc,raster,synthetic_sample


def test_real_photo_lake_append_retains_climate_and_roi_unknown(tmp_path,monkeypatch):
    p,env,receipt,raster,sampler=fixed(tmp_path,monkeypatch)
    result=target.assemble_lake(p,env,receipt,raster,tmp_path/"out",sampler=sampler)
    assert result["actually_processed_original_photo_pixels"]==7
    assert result["genuinely_new_photo_sites_sampled_now"]==4
    photo=pd.read_parquet(tmp_path/"out"/target.PHOTO)
    cov=pd.read_parquet(tmp_path/"out"/target.ENVIRONMENT)
    assert photo.photo_id.astype(str).tolist()==cov.photo_id.astype(str).tolist()
    assert len(photo)==7 and photo.classified.sum()==5 and photo.white.sum()==3
    assert photo.loc[photo.photo_id.eq("3"),"classified"].eq(False).all()
    assert cov.iloc[:3].reset_index(drop=True).equals(pd.read_parquet(env))
    assert cov.loc[cov.photo_id.eq("7"),climate.RAIN].iloc[0]==11
    assert result["soil_unmeasured_for_new_taxa_not_imputed"]
    assert result["output_sha256"][target.ENVIRONMENT] == target.sha256(tmp_path/"out"/target.ENVIRONMENT)
    with pytest.raises(ValueError,match="Cannot overwrite"):
        target.assemble_lake(p,env,receipt,raster,tmp_path/"out",sampler=sampler)


def test_fails_if_old_abiotic_source_has_drifted(tmp_path,monkeypatch):
    p,env,receipt,raster,sampler=fixed(tmp_path,monkeypatch)
    e=pd.read_parquet(env);e.loc[0,climate.RAIN]=-1
    e.to_parquet(env,index=False)
    with pytest.raises(ValueError,match="environment Parquet modified"):
        target.assemble_lake(p,env,receipt,raster,tmp_path/"out",sampler=sampler)


def test_new_raster_sampling_keeps_units_and_missing(tmp_path,monkeypatch):
    monkeypatch.setattr(target,"N_NEW",4)
    raster=tmp_path/"wc"
    for i in range(1,20):
        f=raster/"bio"/f"wc2.1_10m_bio_{i}.tif";f.parent.mkdir(parents=True,exist_ok=True);f.touch()
    for m in range(1,13):
        for folder in ("vapr","srad"):
            f=raster/folder/f"wc2.1_10m_{folder}_{m:02d}.tif"
            f.parent.mkdir(parents=True,exist_ok=True);f.touch()
    f=raster/"elev"/"wc2.1_2.5m_elev.tif";f.parent.mkdir(parents=True,exist_ok=True);f.touch()
    def mock_tif(p,lon,lat):
        if "/srad/" in str(p):return np.full(len(lon),8640,dtype=np.float32)
        if "/vapr/" in str(p):return np.full(len(lon),22,dtype=np.float32)
        return np.arange(len(lon),dtype=np.float32)
    def mock_dem(p,lon,lat):
        return {key:np.full(len(lon),120,dtype=np.float32)
                for key in original.TERRAIN}
    monkeypatch.setattr(target,"sample_tiff",mock_tif)
    monkeypatch.setattr(target,"raster_terrain",mock_dem)
    d=pd.DataFrame([row(i,"white","tail16384") for i in range(1,5)])
    result=target.raster_site_variables(d,raster)
    assert len(result)==4
    assert result[climate.SOLAR].tolist()==[100.]*4
    assert result[original.VAPR].tolist()==[22.]*4
    assert result[original.TERRAIN[0]].eq(120).all()
    assert len(result.columns)==1+19+1+1+7
