#!/usr/bin/env python3
"""MOIS: test whether photo9 has predictive value once pre-index source recording intensity is known.

The earlier FirstRecords v4.0 first observation of an alien species and its
later appearance in other first-record regions are NOT biological introductions,
colonization rates, floral colour evolution, or adaptive radiation.

The historical recording-opportunity index is only the number of recorded alien
plant first records in a previous 25-year window in the initial continent. It
is a reporting-intensity PROXY, not actual observer effort.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import platform
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd
import sklearn
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss, log_loss, roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from scripts.analysis.mois_firstrecords_v40_genus_holdout_predictive_audit_20261124 import (
    load_source, make_folds, prob_loss, bootstrap_diff_by_genus,
    SOURCE_ZIP_SHA256, SEEDS, guard,
)

FILENAME_MD5 = {
    'FirstRecords_dataset_public_v4.0.csv':'cbba8f2c2c518519674fd4f0540809d4',
    'FirstRecords_taxonomy_table_public_v4.0.csv':'d3e21741364b33fabb3cef3b249543bd',
    'FirstRecords_location_table_public_v4.0.csv':'cb69dd1ad613cbfe24db316b89f5884a',
}
FIRST_YEAR=1500
LAST_YEAR=2023
WINDOW=25
REPEATS=5
MODELS = {
    'history_only': ('first_cohort25','first_record_continent','cohort_continent'),
    'history_prior25': ('first_cohort25','first_record_continent','cohort_continent','prior25_events_log1p','prior25_locations_log1p'),
    'history_prior25_white': ('first_cohort25','first_record_continent','cohort_continent','prior25_events_log1p','prior25_locations_log1p','white_label'),
    'history_prior25_all9': ('first_cohort25','first_record_continent','cohort_continent','prior25_events_log1p','prior25_locations_log1p','source_hue9'),
}
NUM_COLS = ['prior25_events_log1p','prior25_locations_log1p']


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_v4_sources(root):
    root=Path(root)
    for name,md5 in FILENAME_MD5.items():
        path=root/name
        guard(path.exists(),'missing canonical source '+name)
        guard(hashlib.md5(path.read_bytes()).hexdigest()==md5,'v4.0 source checksum drift '+name)
    raw=pd.read_csv(root/'FirstRecords_dataset_public_v4.0.csv',sep=';',encoding='utf-8-sig',
                    encoding_errors='replace',dtype=str,keep_default_na=False)
    tax=pd.read_csv(root/'FirstRecords_taxonomy_table_public_v4.0.csv',encoding='utf-8-sig',
                    encoding_errors='replace',dtype=str,keep_default_na=False)
    loc=pd.read_csv(root/'FirstRecords_location_table_public_v4.0.csv',encoding='utf-8-sig',
                    encoding_errors='replace',dtype=str,keep_default_na=False)
    guard((len(raw),len(tax),len(loc))==(90023,56971,557),'source table shape changed')
    return raw,tax,loc


def canonical_v4_plant_records(raw,tax,loc):
    # Preserve source-native grain: one canonical accepted taxon per source location.
    accepted=tax.loc[(tax.taxaGroup=='Vascular plants')&
                     (tax.GBIFstatus=='ACCEPTED')&
                     (tax.GBIFmatchType=='EXACT')&
                     (tax.GBIFtaxonRank=='SPECIES')&
                     (tax.species.ne('')),
                     ['taxonID','species','genus']].copy()
    guard(not accepted.isna().any().any(),'invalid species identifiers')
    conflicts=accepted.groupby('taxonID').species.nunique()
    accepted=accepted.loc[~accepted.taxonID.isin(conflicts[conflicts>1].index)]
    accepted=accepted.drop_duplicates('taxonID')
    location=loc[['locationID','location','continent']].drop_duplicates()
    conflict=location.groupby('locationID')[['location','continent']].nunique().gt(1).any(axis=1)
    location=location.loc[~location.locationID.isin(conflict[conflict].index)].drop_duplicates('locationID')
    guard(location.locationID.is_unique,'ambiguous locationID')
    indexed=raw.merge(accepted,on='taxonID',how='inner',validate='many_to_one',suffixes=('','_accepted'))
    indexed=indexed.merge(location[['locationID','continent']],on='locationID',how='inner',validate='many_to_one')
    indexed=indexed.loc[~(indexed.species.str.contains('\ufffd',regex=False)|
                          indexed.continent.str.contains('\ufffd',regex=False))].copy()
    indexed['year']=pd.to_numeric(indexed.firstRecordEvent,errors='coerce')
    indexed=indexed.loc[indexed.year.between(FIRST_YEAR,LAST_YEAR)&
                        indexed.year.eq(np.floor(indexed.year))].copy()
    indexed['year']=indexed.year.astype(int)
    all_record=indexed.groupby(['species','locationID'],as_index=False).year.min()
    first_unfiltered=all_record.groupby('species').year.min().rename('first_any_source_year')
    exact_orig=indexed.verbatimFirstRecordEvent.str.fullmatch(r'(1[5-9]\d{2}|20[0-2]\d)')
    match_orig=indexed.verbatimFirstRecordEvent.eq(indexed.year.astype(str))
    strict=indexed.loc[
        exact_orig & match_orig &
        indexed.confidenceFirstRecordEvent.eq('high confidence') &
        indexed.occurrenceStatus.eq('present') &
        indexed.establishmentMeans.str.casefold().eq('introduced') &
        indexed.degreeOfEstablishment.str.casefold().isin(['established','invasive']) &
        indexed.accessRights.str.startswith('Free')
    ].copy()
    strict=strict.sort_values(['species','locationID','year']).drop_duplicates(['species','locationID'])
    guard(not strict.duplicated(['species','locationID']).any(),'strict raw duplicate taxon-location found')
    return strict,first_unfiltered,{
        'source_accepted_exact_species':int(accepted.species.nunique()),
        'source_all_accepted_taxon_location_valid_year':len(all_record),
        'source_strict_recorded_alien_taxon_location':len(strict),
        'source_strict_alien_distinct_taxa':int(strict.species.nunique()),
        'source_strict_alien_distinct_locations':int(strict.locationID.nunique()),
        'source_strict_continents':{str(k):int(v) for k,v in strict.continent.value_counts().items()},
        'earlier_source_definition':'earliest 1500–2023 source FirstRecordEvent from any status/confidence for accepted species; may be imputed or nonestablished',
        'effort_definition':'accepted vascular high confidence, original exact 4digit first year, introduced established/invasive, distinct species × source locationID events in first continent at years t-25...t-1',
    }


def source_opportunity_for_species(dataset,strict):
    # Precompute full-year calendar per continent once to avoid repeated filtering.
    all_years=np.arange(FIRST_YEAR,LAST_YEAR+2,dtype=int)
    index={}
    for continent,block in strict.groupby('continent'):
        yr=block.year.to_numpy(dtype=int)
        loc=block.locationID.to_numpy(dtype=str)
        # Number of historical species×location records in 25 years PRIOR to index.
        bincount=np.bincount(yr-FIRST_YEAR,minlength=LAST_YEAR-FIRST_YEAR+1)
        cumsum=np.r_[0,bincount.cumsum()]
        low=np.maximum(all_years-WINDOW,FIRST_YEAR)
        past_events=cumsum[all_years-FIRST_YEAR]-cumsum[low-FIRST_YEAR]
        past_regions=np.zeros(len(all_years),dtype=int)
        for t in range(FIRST_YEAR,LAST_YEAR+2):
            past_regions[t-FIRST_YEAR]=np.unique(loc[(yr>=max(t-WINDOW,FIRST_YEAR))&(yr<t)]).size
        index[continent]=(past_events,past_regions)
    values=np.zeros((len(dataset),2),dtype=int)
    present=np.zeros(len(dataset),dtype=bool)
    for i,(continent,t) in enumerate(zip(dataset.first_record_continent,dataset.first_year)):
        if continent in index and FIRST_YEAR<=int(t)<=LAST_YEAR:
            a,b=index[continent];j=int(t)-FIRST_YEAR
            values[i]=a[j],b[j]
            present[i]=True
    return values,present


def prepare_sample(data,strict,first_any):
    frame=data.copy().reset_index(drop=True)
    opportunity,present=source_opportunity_for_species(frame,strict)
    frame['past25_source_records']=opportunity[:,0]
    frame['past25_reporting_locations']=opportunity[:,1]
    frame['effort_available']=present
    frame['prior25_events_log1p']=np.log1p(frame.past25_source_records)
    frame['prior25_locations_log1p']=np.log1p(frame.past25_reporting_locations)
    frame=frame.join(first_any,on='species')
    frame['source_has_earlier_any_year']=frame.first_any_source_year.lt(frame.first_year)
    frame['source_missing_any_year']=frame.first_any_source_year.isna()
    guard(not frame.duplicated('key').any(),'sample no longer one species per row')
    guard(frame.source_missing_any_year.sum()==0,'source first-record check missing for focal species')
    guard((frame.first_any_source_year<=frame.first_year).all(),'first-any date cannot exceed strict first date')
    # Nonunique first year continents cannot be assigned a single historical
    # exposure. Keep as a different completeness class, not zero source effort.
    return frame


def probability_loss(y,p):
    pp=np.clip(p,1e-10,1-1e-10)
    return -(y*np.log(pp)+(1-y)*np.log1p(-pp))


def fit_oof(folds_data,seed):
    y=folds_data.any_new_regions_50yr.to_numpy(int)
    out={k:np.empty(len(folds_data)) for k in MODELS}
    for tr,te in make_folds(folds_data.genus.to_numpy(),seed):
        for name,cols in MODELS.items():
            nums=[v for v in cols if v in NUM_COLS]
            cats=[v for v in cols if v not in NUM_COLS]
            transformers=[('cat',OneHotEncoder(handle_unknown='ignore'),cats)]
            if nums:
                transformers.append(('historic_proxy',StandardScaler(),nums))
            pre=ColumnTransformer(transformers)
            model=make_pipeline(pre,LogisticRegression(max_iter=1500,C=1.0,solver='lbfgs'))
            model.fit(folds_data.iloc[tr][list(cols)],y[tr])
            out[name][te]=model.predict_proba(folds_data.iloc[te][list(cols)])[:,1]
    for key,values in out.items():
        guard(np.isfinite(values).all() and ((values>=0)&(values<=1)).all(),
              'invalid out-of-genus score '+key)
    return out


def score_sample(sample,repeats=5):
    sample=sample.reset_index(drop=True)
    if len(sample)<300 or sample.genus.nunique()<100:
        return {'status':'insufficient sample for out-of-genus 5fold',
                'n_species':len(sample),'n_genera':int(sample.genus.nunique())},None
    y=sample.any_new_regions_50yr.to_numpy(int)
    preds={key:[] for key in MODELS}
    repeated=[]
    for seed in SEEDS[:repeats]:
        predictions=fit_oof(sample,seed)
        for key,p in predictions.items():
            preds[key].append(p)
            repeated.append({'seed':seed,'model':key,
                             'logloss':float(log_loss(y,p)),
                             'brier':float(brier_score_loss(y,p)),
                             'AUC':float(roc_auc_score(y,p))})
    means={k:np.mean(np.stack(v),axis=0) for k,v in preds.items()}
    loss={k:np.mean(np.stack([probability_loss(y,p) for p in v]),axis=0)
          for k,v in preds.items()}
    metrics={k:{'logloss':float(log_loss(y,p)),
                'brier':float(brier_score_loss(y,p)),
                'AUC':float(roc_auc_score(y,p))} for k,p in means.items()}
    compares={}
    base='history_prior25'
    for model in ['history_prior25_white','history_prior25_all9']:
        d=loss[base]-loss[model]
        compares[model]={
            'history_effort_baseline_minus_model_positive_means_hue_improves':True,
            'delta_logloss':float(d.mean()),
            'delta_logloss_genus_cluster_bootstrap':bootstrap_diff_by_genus(
                sample,d,seed=20261124,draws=3000),
            'delta_AUC':float(metrics[model]['AUC']-metrics[base]['AUC']),
        }
    increment=loss['history_only']-loss[base]
    compares['added_historical_opportunity']={
        'delta_logloss':float(increment.mean()),
        'positive_means_historical_opportunity_improves':True,
        'delta_logloss_genus_cluster_bootstrap':bootstrap_diff_by_genus(
            sample,increment,seed=20261125,draws=3000),
        'delta_AUC':float(metrics[base]['AUC']-metrics['history_only']['AUC'])}
    frame=sample[['key','species','genus','first_year','first_record_continent',
                  'source_hue9','any_new_regions_50yr','first_any_source_year',
                  'source_has_earlier_any_year','past25_source_records',
                  'past25_reporting_locations']].copy()
    for k,v in means.items():frame[k+'_OOF_mean_probability']=v
    for k,v in loss.items():frame[k+'_meanOOF_logloss']=v
    return {'status':'completed','n_species':len(sample),'n_genera':int(sample.genus.nunique()),
            'n_observed_firstrecords_with_added_region_in_50y':int(y.sum()),
            'event_fraction':float(y.mean()),
            'historical_effort_prior25_proxy_variation':{
                'n_distinct_events':int(sample.past25_source_records.nunique()),
                'n_distinct_reporting_locations':int(sample.past25_reporting_locations.nunique()),
                'median_prior25_events':float(sample.past25_source_records.median()),
                'median_prior25_reporting_locations':float(sample.past25_reporting_locations.median()),
                'fraction_zero_historical_records':float(sample.past25_source_records.eq(0).mean())},
            'CV_repeats':repeats,
            'CV_group':'genus 5-fold, no genus in both training and validation',
            'meanOOF_metrics':metrics,
            'paired_CV_deltas':compares,
            'repeat_metrics':repeated},frame


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--source',required=True)
    parser.add_argument('--v4',required=True)
    parser.add_argument('--out',required=True)
    args=parser.parse_args()
    out=Path(args.out)
    guard(not out.exists(),'do not overwrite frozen historical opportunity source audit')
    data,metadata=load_source(args.source)
    raw,tax,loc=load_v4_sources(args.v4)
    strict,any_first,source_audit=canonical_v4_plant_records(raw,tax,loc)
    frame=prepare_sample(data,strict,any_first)
    mainframe=frame.loc[frame.effort_available].copy()
    earlier=mainframe.loc[~mainframe.source_has_earlier_any_year].copy()
    guard(len(mainframe)>=2500,'Unexpected loss of primary species due to continent opportunity index')
    guard((mainframe.first_year.between(1800,1973)).all(), 'followup years leaked')
    all_result,main_predictions=score_sample(mainframe,repeats=5)
    guard(all_result['status']=='completed','Main predictive test did not execute')
    print('MOIS_HISTORICAL_EFFORT_CV_RESULT full_same_continent',
          json.dumps({'n_species':len(mainframe),'scores':all_result['meanOOF_metrics'],
                      'delta':all_result['paired_CV_deltas']},ensure_ascii=False),flush=True)
    sensitivity,_=score_sample(earlier,repeats=3)
    print('MOIS_HISTORICAL_EFFORT_CV_RESULT no_earlier_any_first_record',
          json.dumps({'n_species':len(earlier),'scores':sensitivity.get('meanOOF_metrics'),
                      'delta':sensitivity.get('paired_CV_deltas')},ensure_ascii=False),flush=True)
    out.mkdir(parents=True)
    main_predictions.to_csv(out/'species_level_prior_reporting_opportunity_and_OOF.csv.gz',
                            index=False,compression={'method':'gzip','mtime':0})
    pd.DataFrame(all_result['repeat_metrics']).to_csv(out/'five_repeat_genus_OOF_history_effort_scores.csv',index=False)
    by_continent=mainframe.groupby('first_record_continent',dropna=False).agg(
        n_species=('key','size'),observed_additional_record_frac=('any_new_regions_50yr','mean'),
        mean_first_year=('first_year','mean'),
        mean_prior25_events=('past25_source_records','mean'),
        median_prior25_events=('past25_source_records','median'),
        mean_prior25_reporting_locations=('past25_reporting_locations','mean'),
        n_earlier_any_year=('source_has_earlier_any_year','sum')).reset_index()
    by_continent.to_csv(out/'continent_recording_opportunity_profile.csv',index=False)
    mainframe['cohort25']=((mainframe.first_year-1800)//25).astype(int)
    by_cohort=mainframe.groupby(['cohort25','first_record_continent'],dropna=False).agg(
        n_species=('key','size'),observed_added_record_frac=('any_new_regions_50yr','mean'),
        median_prior25_event_count=('past25_source_records','median'),
        median_prior25_location_count=('past25_reporting_locations','median')).reset_index()
    by_cohort.to_csv(out/'continent_year_cohort_recording_proxy.csv',index=False)
    summary={
        'source_firstrecords_v40_original_MD5':FILENAME_MD5,
        'source_50y_species_archive_SHA256':SOURCE_ZIP_SHA256,
        'source_audit':source_audit,
        'n_previous_50y_species':len(data),
        'n_primary_with_one_first_continent':len(mainframe),
        'n_multi_first_continent_no_effort_mapping':int((~frame.effort_available).sum()),
        'n_earlier_any_source_year':int(mainframe.source_has_earlier_any_year.sum()),
        'n_earlier_any_source_year_among_whole_previous_sample':int(frame.source_has_earlier_any_year.sum()),
        'main_CV':all_result,
        'restricted_to_no_earlier_ANY_first_record':sensitivity,
        'earlier_any_firstrecord_is_nonstrict_and_can_be_uncertain_or_imputed':True,
        'source_event_counts_are_not_observer_effort':True,
        'strict_firstrecord_is_not_true_introduction':True,
        'photo9_is_not_independently_verified_perianth_hue':True,
        'predictive_nonimprovement_is_not_biological_equivalence':True,
        'software':{'python':platform.python_version(),'numpy':np.__version__,
                    'pandas':pd.__version__,'scikit_learn':sklearn.__version__}
    }
    (out/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    (out/'receipt.json').write_text(json.dumps({
        'script_SHA256':digest(__file__),
        'source_firstrecords_v40_original_MD5':FILENAME_MD5,
        'source_50y_species_archive_SHA256':SOURCE_ZIP_SHA256,
        'outputs_sha256':{p.name:digest(p) for p in sorted(out.iterdir())}
    },indent=2)+'\n')
    print('MOIS_HISTORICAL_RECORD_OPPORTUNITY_COMPLETE',json.dumps({
        'source_audit':source_audit,
        'species_prior_record_qc':{
            'n':len(mainframe),'earlier_firstrecords':int(mainframe.source_has_earlier_any_year.sum())},
        'CV_scores':all_result['meanOOF_metrics'],
        'CV_deltas':all_result['paired_CV_deltas'],
        'earlier_excluded_scores':sensitivity.get('meanOOF_metrics')},ensure_ascii=False),flush=True)

if __name__=='__main__': main()