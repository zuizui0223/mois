"""Four visible hues are compositional, not anthocyanin-only states."""
import numpy as np
import pandas as pd
from scripts.analysis import mois_original313535_four_visible_hues_climate_mixture as m


def photos(n=450):
    rng=np.random.default_rng(53)
    rows=[]
    for i in range(n):
        rain=float(100+(i%3)*600+i*.02)
        solar=float(70+((i//3)%3)*100)
        hue=m.HUES[(i*13+i//9)%4]
        if i%17==0:hue="mixed_uncertain"
        rows.append({
            "photo_id":str(i+1),"morph":hue,"classified":hue in m.HUES,
            "inat_taxon_id":str(i//3),
            "genus":f"G{i//90}","species":f"G{i//90} sp{i//3}",
            m.RAIN:rain,m.SOLAR:solar,
            "latitude":float(i%9*8-36),
            "abs_lat":float(abs(i%9*8-36)),
            "longitude_sin":float(np.sin(i/300)),
            "longitude_cos":float(np.cos(i/300)),
            "source_new":0.,
        })
    return pd.DataFrame(rows)


def test_hue_grid_explicit_technical_unknown_denominator():
    p=photos()
    b,cut=m.classify_climate(p)
    cells=m.four_hue_source_cells(b)
    assert len(cells)==9
    assert cells.original_processed_photo_attempts.sum()==len(p)
    assert cells.real_original_classifiable_ROI_photos.sum()==int(p.classified.sum())
    assert cells.technically_unknown_hue_photo_attempts.sum()==int((~p.classified).sum())
    for row in cells.itertuples():
        assert sum(getattr(row,h+"_photo_ROIs") for h in m.HUES)==row.real_original_classifiable_ROI_photos
        assert abs(sum(getattr(row,h+"_among_4class_photos") for h in m.HUES)-1)<1e-12


def test_wet_dry_hue_effect_sums_zero_not_anthocyanin_proxy():
    p=photos()
    t,_=m.classify_climate(p)
    summary=m.four_hue_source_cells(t)
    delta=m.wet_dry_changes(summary)
    assert len(delta)==12
    assert delta.visible_hue.nunique()==4
    for solar,group in delta.groupby("solar_tertile"):
        assert abs(group.wet_minus_dry_photo_percentage_points.sum())<1e-10
    assert "mixed_uncertain" not in m.HUES


def test_simplex_normalizes_even_out_of_range_model_outputs():
    p=m.simplex(np.array([[-.2,1.5,.2,-.4],[.1,.1,.5,.8]]))
    assert np.all(p>=0)
    assert np.allclose(p.sum(axis=1),1)
    assert np.isfinite(m.brier(np.ones((2,4))*.25,p,np.array([1.,2.])))


def test_train_predict_photo_hue_brier_holdout(monkeypatch):
    rng=np.random.default_rng(310)
    orig=photos(540)
    # Override genus and source taxon to produce enough distinct folds.
    orig["genus"]=[f"G{i//30}" for i in range(len(orig))]
    orig["inat_taxon_id"]=[str(i//3) for i in range(len(orig))]
    for block in m.CLIMATE:
        for covar in m.parent.BLOCKS[block]:
            if covar==m.RAIN or covar==m.SOLAR:continue
            orig[covar]=rng.normal(size=len(orig))+np.arange(len(orig))*.003
    # 180 eligible taxa /18 genera; test the taxon table and multi-target Ridge.
    monkeypatch.setattr(m,"FOLDS",3)
    # Named source requirement is intended for the real ~150k cases.
    # Instead directly exercise predict_fold with a small verified four-hue
    # fraction panel whose taxa and features are constructed independently.
    original=m.taxon_table
    try:
        # temporary neutral length gate through rows replicated over enough
        # nominal taxa with distinct IDs; avoid modifying production bounds.
        bulk=pd.concat([
            orig.assign(inat_taxon_id=orig.inat_taxon_id.astype(str)+"_k"+str(k),
                        genus=orig.genus.astype(str)+"_k"+str(k))
            for k in range(24)],ignore_index=True)
        t,labels,meta=original(bulk,tuple(
            v for b in m.CLIMATE for v in m.parent.BLOCKS[b]))
    finally:
        pass
    assert len(t)>1000
    assert np.allclose(t[labels].sum(axis=1),1)
    train=np.arange(0,len(t),2);test=np.arange(1,len(t),2)
    result=m.predict_fold(t,train,test,t[labels].to_numpy(float),m.CLIMATE)
    assert set(result)=={"full","geo_only","drop_temperature",
                         "drop_moisture","drop_solar"}
    assert all(np.isfinite(v["weighted_brier"]) for v in result.values())
