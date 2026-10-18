"""WorldClim full nineteen BIO / VAPR / DEM original photograph source gates."""
import numpy as np
import pandas as pd
import pytest

from scripts.analysis import mois_original313535_full19_bio_vapr_terrain_lake as m

def test_exactly_19_climate_bio_plus_vapour_and_annual_solar():
    assert len(m.BIO)==19 and m.BIO[0]=="wc2_bio01_10m"
    assert m.BIO[-1]=="wc2_bio19_10m"
    assert m.VAPR=="wc2_annual_vapr_10m"
    assert len(m.CLIMATE)==21 and len(m.TERRAIN)==7
    assert m.PIN_BIO=="00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf"

def test_photo_source_and_actual_vapour_average_no_pseudoweather(monkeypatch,tmp_path):
    monkeypatch.setattr(m.previous,"N_PHOTOS",3)
    p=pd.DataFrame({"photo_id":["1","2","3"],
                    "latitude":[34.,-41.,15.],"longitude":[135.,149.,98.],
                    m.previous.SOLAR:[120.,200.,240.]})
    for i in range(1,20):
        t=tmp_path/"bio"/f"wc2.1_10m_bio_{i}.tif"
        t.parent.mkdir(parents=True,exist_ok=True)
        t.touch()
    for i in range(1,13):
        t=tmp_path/"vapr"/f"wc2.1_10m_vapr_{i:02d}.tif"
        t.parent.mkdir(parents=True,exist_ok=True)
        t.touch()
    def sample(path,lon,lat):
        name=path.name
        if "vapr" in name:
            return np.array([2.,4.,np.nan])
        index=int(name.split("_")[-1].split(".")[0])
        return np.array([index, index+1,index+2.],dtype=float)
    monkeypatch.setattr(m,"sample_tiff",sample)
    monkeypatch.setattr(m,"raster_terrain",
        lambda path,lon,lat:{k:np.array([100.,250.,50.]) for k in m.TERRAIN})
    out=m.climatology(p,tmp_path)
    assert out.photo_id.tolist()==["1","2","3"]
    assert np.allclose(out[m.BIO[0]],[1,2,3])
    assert np.allclose(out[m.BIO[-1]],[19,20,21])
    assert np.allclose(out[m.VAPR].iloc[:2],[2.,4.])
    assert np.isnan(out[m.VAPR].iloc[2])
    assert out[m.previous.SOLAR].tolist()==[120.,200.,240.]

def test_previously_recorded_original_bio_cannot_be_silently_replaced():
    p=pd.DataFrame({k:[v+1.,v+2.,v+3.] for k,v in [
        (m.previous.HEAT,20.),
        (m.previous.RAIN,100.),
        (m.previous.DRY,9.),
        (m.previous.SEASON,2.)]})
    same=p.copy()
    monkeypatch=None
    old=m.need
    def permissive_coverage(condition,why):
        if "unavailable" in why:return None
        return old(condition,why)
    # Exercise equality without needing full 313535-photo synthetic data.
    m.need=permissive_coverage
    try:
        r=m.verify_existing_photo_environment(p,same)
        assert set(r)==set([m.previous.HEAT,m.previous.RAIN,m.previous.DRY,m.previous.SEASON])
        bad=same.copy()
        bad[m.previous.RAIN]+=10
        with pytest.raises(ValueError,match="disagree"):
            m.verify_existing_photo_environment(p,bad)
    finally:
        m.need=old
