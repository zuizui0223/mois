"""PCA is not selected by outcome; every group comparison uses matched source photos."""
import numpy as np
import pandas as pd
from scripts.analysis import mois_global313535_all_environment_pca_ablation as m

def fake_panel(n_genera=205,taxa_per_genus=6,seed=20261018):
    rng=np.random.default_rng(seed)
    rows=[];number=0
    for g in range(n_genera):
        genus=f"G{g:04d}"
        la=-70+(g%20)*7.0
        lo=-155+(g%17)*17.0
        for t in range(taxa_per_genus):
            number+=1
            x=rng.normal(size=len(m.ABIOTIC))
            row={k:float(x[i]+0.03*(g%5)) for i,k in enumerate(m.ABIOTIC)}
            # Temperature/precipitation correlation and solar interactions
            row[m.env.BIO[0]]=float(row[m.env.BIO[11]]*.75+rng.normal(0,.4))
            white=float(rng.uniform()<.45+.12*np.tanh(row[m.env.BIO[11]]))
            row.update({"photo_id":str(number),
                "inat_taxon_id":str(number),
                "species":genus+" species",
                "genus":genus,
                "white":bool(white),
                "classified":True,
                "latitude":la+rng.normal(0,.3),
                "abs_lat":abs(la+rng.normal(0,.3)),
                "longitude_sin":np.sin(np.deg2rad(lo)),
                "longitude_cos":np.cos(np.deg2rad(lo)),
                "source_new":float(g%2),
            })
            rows.append(row)
    return pd.DataFrame(rows)

def test_all_19_bio_and_static_predictors_count_and_identity():
    assert len(m.env.BIO)==19
    assert len(m.ABIOTIC)==33
    assert len(set(m.ABIOTIC))==33
    assert set(m.BLOCKS)=={"temperature","moisture","solar","terrain","soil"}
    assert m.VARIANCE==.9

def test_photographic_species_panel_gates_soil_before_species_aggregation():
    p=fake_panel()
    expected=len(p)
    p.loc[p.index[:8],m.prev.SOIL[0]]=np.nan
    s,r=m.species_panel(p)
    assert r["source_classified_original_photo_rows_exact_complete_case"]==expected-8
    assert len(s)>1000
    assert r["source_nominal_genera_after_gate"]>=100
    assert r["complete_case_photo_subset_identical_for_all_predictor_ablation"]

def test_train_only_scaling_and_pca_no_outcome_variable():
    rng=np.random.default_rng(14)
    tr=pd.DataFrame(rng.normal(size=(100,6)),columns=list("abcdef"))
    te=pd.DataFrame([[100]*6],columns=list("abcdef"))
    a,b,n=m.fit_block_transform(tr,te,tr.columns)
    assert a.shape==(100,n)
    assert b.shape==(1,n)
    assert 1<=n<=6
    assert abs(float(a.mean()))<1e-6
    assert not np.isclose(float(np.max(np.abs(b))),0)

def test_genuine_group_comparison_without_leaking_heldout_test_group():
    s,receipt=m.species_panel(fake_panel())
    test=m.cv_abiotic_importance(s)
    assert len(test)==2*(2+2*len(m.BLOCKS))
    for holdout,rows in test.groupby("heldout_design"):
        assert len(rows)==12
        assert rows.no_training_or_PCA_parameter_leakage_into_test.all()
        assert set(rows.model).issuperset(
            {"geo_only","ALL_five_environment_groups",
             "drop_soil","drop_moisture","only_terrain"})
        assert np.isfinite(rows.whole_holdout_group_cv_R2).all()

def test_pca_PC1_environment_variance_not_assumed_white_importance():
    s,_=m.species_panel(fake_panel())
    axes,loads,corr=m.pca_descriptives(s)
    p1=axes[(axes.pca_block=="all_33_abiotic_full_community_descriptive") &
            (axes.axis=="PC1")]
    assert 0<float(p1.explained_abiotic_variance_fraction.iloc[0])<1
    assert len(loads[loads.axis=="PC1"])>=33
    assert len(corr.columns)>0
