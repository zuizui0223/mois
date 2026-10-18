"""Additional original 4096 photo site climate is sampled only after signed pixels."""
import numpy as np
import pandas as pd
import pytest

from scripts.analysis import mois_postwave_real317631_incremental_feature_lake as m


def test_true_extra_original_pixel_counts_and_no_missing_hue_to_colour():
    assert m.NEW==4096 and m.TOTAL==317631
    assert m.PRIOR_SOURCE.startswith("mois_full313535")
    assert m.INPUT_SOURCE.startswith("mois_full317631")
    assert m.base.N_HUES==151239


def test_raster_extract_new_photos_only(monkeypatch,tmp_path):
    monkeypatch.setattr(m,"NEW",3)
    x=pd.DataFrame({"photo_id":["10","11","12"],
                    "longitude":["138.2","133.4","152.0"],
                    "latitude":["34.3","42.7","12.3"]})
    for i in range(1,20):
        folder=tmp_path/"bio";folder.mkdir(exist_ok=True)
        (folder/f"wc2.1_10m_bio_{i}.tif").touch()
    for i in range(1,13):
        for var in ("vapr","srad"):
            folder=tmp_path/var;folder.mkdir(exist_ok=True)
            (folder/f"wc2.1_10m_{var}_{i:02d}.tif").touch()
    def fake_sample(p,lon,lat):
        if "vapr" in p.name:return np.array([10.,20.,np.nan])
        if "srad" in p.name:return np.array([8640.,17280.,np.nan])
        n=int(p.stem.split("_")[-1])
        return np.array([float(n),float(n*2),np.nan])
    monkeypatch.setattr(m,"sample_tiff",fake_sample)
    monkeypatch.setattr(m,"raster_terrain",
        lambda path,lon,lat:{k:np.array([1.,2.,np.nan])
                              for k in m.expanded.TERRAIN})
    out=m.new_raster_variables(x,tmp_path)
    assert len(out)==3 and out.photo_id.tolist()==["10","11","12"]
    assert np.allclose(out[m.expanded.BIO[0]].iloc[:2],[1,2])
    assert np.allclose(out[m.expanded.VAPR].iloc[:2],[10,20])
    assert np.allclose(out[m.base.SOLAR].iloc[:2],[100,200])
    assert np.isnan(out[m.base.SOLAR].iloc[2])
    assert np.isnan(out[m.expanded.VAPR].iloc[2])


def test_previous_photo_rows_immutable_append_only(monkeypatch):
    monkeypatch.setattr(m,"TOTAL",4)
    monkeypatch.setattr(m,"NEW",2)
    monkeypatch.setattr(m.base,"N_PHOTOS",2)
    old=pd.DataFrame({"photo_id":["1","2"],"morph":["white","red_pink"]})
    new=pd.DataFrame({"photo_id":["3","4"],"morph":["mixed_uncertain","white"]})
    combo=m.append(old,new,"foo")
    assert combo.photo_id.tolist()==["1","2","3","4"]
    assert combo.iloc[:2].reset_index(drop=True).equals(old)
    with pytest.raises(ValueError,match="Cannot append"):
        m.append(old,pd.DataFrame({"photo_id":["2","4"],
                                  "morph":["white","white"]}),"foo")
