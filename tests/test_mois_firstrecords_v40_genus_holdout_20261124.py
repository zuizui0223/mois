import importlib.util
from pathlib import Path
import os
import numpy as np
import pandas as pd

FILE=Path(__file__).resolve().parents[1]/'scripts'/'analysis'/'mois_firstrecords_v40_genus_holdout_predictive_audit_20261124.py'
spec=importlib.util.spec_from_file_location('mois_v40_cv',FILE)
audit=importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def test_pinned_source_and_record_grain():
    source=os.environ.get('MOIS_V40_SOURCE_ZIP', '/mnt/data/MOIS_FirstRecords_v40_50year_new_regions_999_20261122.zip')
    data,meta=audit.load_source(source)
    assert len(data)==2813
    assert data.key.nunique()==2813
    assert data.genus.nunique()==1294
    assert data.any_new_regions_50yr.sum()==793
    assert set(data.first_record_continent)=={'Africa','Asia','Europe','Oceania','North America','MULTIPLE_CONTINENTS'}
    assert data.first_year.between(1800,1973).all()


def test_genus_held_out_never_shared_and_each_species_tested_once():
    genera=np.array(['a']*12+['b']*4+['c']*9+['d']*5+['e']*5+['f']*2+['g']*2)
    test=np.concatenate([idx for _,idx in audit.make_folds(genera,20261124,5)])
    assert len(test)==len(genera)
    assert np.array_equal(np.sort(test),np.arange(len(genera)))
    for tr,te in audit.make_folds(genera,20261124,5):
        assert not set(genera[tr])&set(genera[te])


def test_logloss_orientation_and_categorical_ood_safety():
    y=np.array([0,1],dtype=int)
    good=audit.prob_loss(y,np.array([.1,.9]))
    bad=audit.prob_loss(y,np.array([.9,.1]))
    assert (good<bad).all()
    assert (audit.prob_loss(y,np.array([0,1]))>=0).all()
    d=pd.DataFrame({'genus':['A']*3+['B']*3+['C']*3+['D']*3+['E']*3,
                    'first_cohort25':list('000111222333444'),
                    'first_record_continent':['Asia']*6+['Europe']*9,
                    'cohort_continent':['0|Asia']*3+['1|Asia']*3+['2|Europe']*3+['3|Europe']*3+['4|Europe']*3,
                    'white_label':['white','white','nonwhite']*5,
                    'source_hue9':['white','yellow','white']*5,
                    'any_new_regions_50yr':[0,1,0,1,0,1,0,1,0,1,0,1,0,1,0]})
    predictions=audit.fit_oof(d,20261124)
    for pred in predictions.values():
        assert len(pred)==15
        assert np.isfinite(pred).all()
        assert ((pred>0)&(pred<1)).all()
