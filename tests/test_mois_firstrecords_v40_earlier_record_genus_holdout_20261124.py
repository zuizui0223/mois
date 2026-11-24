"""Contracts for earliest-any-source historical-record provenance and leakage controls."""
import numpy as np
import pandas as pd
from scripts.analysis.mois_firstrecords_v40_earlier_record_genus_holdout_20261124 import (
    MODELS,NUM,fit_models,input_rows,evaluate
)
from scripts.analysis.mois_firstrecords_v40_genus_holdout_predictive_audit_20261124 import make_folds


def test_historical_predictor_references_only_prebaseline_or_atbaseline_source_fields():
    restricted={'first_cohort25','first_record_continent','cohort_continent',
                'prior25_events_log1p','prior25_locations_log1p','earlier_record_flag',
                'white_label','source_hue9'}
    for cols in MODELS.values():
        assert set(cols)<=restricted
        assert 'any_new_regions_50yr' not in cols
        assert 'new_regions_50yr' not in cols
        assert 'first_any_source_year' not in cols
    assert NUM==['prior25_events_log1p','prior25_locations_log1p']


def test_folds_keep_all_species_in_one_genus_out_of_training():
    genera=np.array(['A']*4+['B']*6+['C']*4+['D']*2+['E']*2+['F']*2)
    for tr,te in make_folds(genera,20261124,n_folds=5):
        assert not set(genera[tr])&set(genera[te])
        assert len(tr)+len(te)==len(genera)


def test_small_dataset_refuses_false_crossvalidated_metric():
    d=pd.DataFrame({'genus':['A','B'],'any_new_regions_50yr':[0,1]})
    result,preds=evaluate(d)
    assert result['status'].startswith('insufficient')
    assert preds is None
