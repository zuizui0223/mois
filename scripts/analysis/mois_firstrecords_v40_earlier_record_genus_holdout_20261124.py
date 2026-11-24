#!/usr/bin/env python3
"""Second-stage signed-data diagnostic: earliest-source record completeness vs photo9.

The historical earliest-any-source event can be low-confidence, imputed or
casual; its flag measures heterogeneous archival eligibility NOT actual
introduction date or confirmed colonization. All included predictors are from
year < baseline or baseline metadata (not future first-record outcomes).
"""
from __future__ import annotations
import argparse,io,hashlib,json,platform
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd
import sklearn
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OneHotEncoder,StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import log_loss,brier_score_loss,roc_auc_score
from scripts.analysis.mois_firstrecords_v40_genus_holdout_predictive_audit_20261124 import (
    make_folds,prob_loss,bootstrap_diff_by_genus,SEEDS,guard
)
SOURCE_ZIP_SHA='ee212ebd94d374e31d5511f22f0042f20f45618b47961faf8bce8fd45124e26e'
EVENT_SOURCE='species_level_prior_reporting_opportunity_and_OOF.csv.gz'
NUM=['prior25_events_log1p','prior25_locations_log1p']
COMMON=['first_cohort25','first_record_continent','cohort_continent']+NUM
MODELS={
  'historical_cohort_continent':COMMON[:3],
  'history_prior25':COMMON,
  'history_prior25_earlier_record':COMMON+['earlier_record_flag'],
  'history_prior25_earlier_record_white':COMMON+['earlier_record_flag','white_label'],
  'history_prior25_earlier_record_all9':COMMON+['earlier_record_flag','source_hue9']
}

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def input_rows(source):
    guard(sha(source)==SOURCE_ZIP_SHA,'signed original recording opportunity ZIP changed')
    with ZipFile(source) as z:
        receipt=json.load(z.open('receipt.json'))
        for name,h in receipt['outputs_sha256'].items():
            guard(hashlib.sha256(z.read(name)).hexdigest()==h,'signed source member changed '+name)
        d=pd.read_csv(io.BytesIO(z.read(EVENT_SOURCE)),compression='gzip',dtype={'key':str})
        s=json.load(z.open('summary.json'))
    guard(len(d)==2799 and d.key.is_unique and int(d.any_new_regions_50yr.sum())==786,
          'fixed sample 2799 species / 786 recorded new-region events changed')
    guard(int(d.source_has_earlier_any_year.sum())==1486,'archive eligibility flag changed')
    d['first_cohort25']=((d.first_year-1800)//25).astype(str)
    d['cohort_continent']=d.first_cohort25+'|'+d.first_record_continent
    d['white_label']=np.where(d.source_hue9=='white','white','nonwhite')
    d['earlier_record_flag']=np.where(d.source_has_earlier_any_year,'earlier-any-source','no-earlier-source')
    for col,src in [('prior25_events_log1p','past25_source_records'),
                    ('prior25_locations_log1p','past25_reporting_locations')]:
        d[col]=np.log1p(d[src])
    guard(d[NUM].notna().all().all(),'pre-event source proxies missing')
    return d,s

def fit_models(dataset,seed):
    y=dataset.any_new_regions_50yr.to_numpy(int)
    out={k:np.empty(len(dataset)) for k in MODELS}
    for tr,te in make_folds(dataset.genus.to_numpy(),seed):
        for label,cols in MODELS.items():
            cat=[c for c in cols if c not in NUM]
            num=[c for c in cols if c in NUM]
            transforms=[('categories',OneHotEncoder(handle_unknown='ignore'),cat)]
            if num:transforms.append(('preceding25y',StandardScaler(),num))
            clf=make_pipeline(ColumnTransformer(transforms),
                              LogisticRegression(max_iter=1500,C=1.,solver='lbfgs'))
            clf.fit(dataset.iloc[tr][cols],y[tr])
            out[label][te]=clf.predict_proba(dataset.iloc[te][cols])[:,1]
    for label,p in out.items():
        guard(np.isfinite(p).all() and ((p>0)&(p<1)).all(),label+' OOF invalid')
    return out

def evaluate(dataset,n_repeats=5):
    data=dataset.reset_index(drop=True)
    y=data.any_new_regions_50yr.to_numpy(int)
    if len(data)<300 or data.genus.nunique()<100:
        return {'status':'insufficient grouped-CV evaluation','n_species':len(data)},None
    preds={k:[] for k in MODELS}
    rows=[]
    for seed in SEEDS[:n_repeats]:
        oof=fit_models(data,seed)
        for model,p in oof.items():
            preds[model].append(p)
            rows.append({'seed':seed,'model':model,'n_species':len(data),
                'logloss':float(log_loss(y,p)),'brier':float(brier_score_loss(y,p)),
                'AUC':float(roc_auc_score(y,p))})
    means={k:np.mean(np.stack(v),axis=0) for k,v in preds.items()}
    losses={k:np.mean(np.stack([prob_loss(y,p) for p in v]),axis=0) for k,v in preds.items()}
    metrics={k:{'OOF_mean_probability_logloss':float(log_loss(y,p)),
                'OOF_mean_probability_brier':float(brier_score_loss(y,p)),
                'OOF_mean_probability_AUC':float(roc_auc_score(y,p))} for k,p in means.items()}
    paired={}
    comparisons={
        'historical_prior25_vs_no_effort':('historical_cohort_continent','history_prior25'),
        'earlier_record_vs_history_prior25':('history_prior25','history_prior25_earlier_record'),
        'white_vs_earlier_record_baseline':('history_prior25_earlier_record','history_prior25_earlier_record_white'),
        'all9_vs_earlier_record_baseline':('history_prior25_earlier_record','history_prior25_earlier_record_all9')}
    for key,(b,a) in comparisons.items():
        d=losses[b]-losses[a]
        paired[key]={'baseline':b,'added':a,
                     'positive_delta_means_additional_feature_improves':True,
                     'delta_logloss':float(d.mean()),
                     'delta_AUC_on_meanOOF':metrics[a]['OOF_mean_probability_AUC']-metrics[b]['OOF_mean_probability_AUC'],
                     'genus_cluster_bootstrap_conditional_on_fitted_OOF':bootstrap_diff_by_genus(
                          data,d,draws=3000)}
    by_record=data.groupby('earlier_record_flag').agg(
        n_species=('key','size'),n_additional_record_species=('any_new_regions_50yr','sum'),
        added_region_record_fraction=('any_new_regions_50yr','mean'),
        first_year_median=('first_year','median')).reset_index()
    pred=data[['key','genus','species','first_year','first_record_continent','source_hue9',
               'earlier_record_flag','first_any_source_year','any_new_regions_50yr',
               'past25_source_records','past25_reporting_locations']].copy()
    for key,v in means.items():pred[key+'_meanOOF_probability']=v
    for key,v in losses.items():pred[key+'_meanOOF_logloss']=v
    report={'status':'completed','n_species':len(data),'n_genera':int(data.genus.nunique()),
            'n_source_earlier_record_species':int(data.source_has_earlier_any_year.sum()),
            'n_added_region_event_species':int(y.sum()),
            'n_repeated_cv':n_repeats,'cv_split':'5 folds held out entire genus, no genus leakage',
            'scores':metrics,'paired_comparisons':paired,
            'historical_source_firstrecord_completeness_groups':by_record.to_dict(orient='records'),
            'repeat_metrics':rows}
    return report,pred

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--source',required=True)
    ap.add_argument('--out',required=True)
    args=ap.parse_args()
    out=Path(args.out)
    guard(not out.exists(),'do not overwrite frozen signed data analysis')
    data,source_meta=input_rows(args.source)
    summary,pred=evaluate(data,5)
    guard(summary['status']=='completed','failed main test')
    out.mkdir(parents=True)
    pred.to_csv(out/'species_historical_source_completeness_OOF_predictions.csv.gz',
                index=False,compression={'method':'gzip','mtime':0})
    pd.DataFrame(summary['repeat_metrics']).to_csv(out/'genus_holdout_repeated_5fold_metrics.csv',index=False)
    pd.DataFrame(summary['historical_source_firstrecord_completeness_groups']).to_csv(
        out/'earlier_any_firstrecord_group_rates.csv',index=False)
    final={'source_recording_opportunity_artifact_SHA256':SOURCE_ZIP_SHA,
           'main':summary,
           'earlier_any_firstrecord_can_include_imputed_casual_or_uncertain':True,
           'source_earliest_strict_not_true_first_record':True,
           'source_recording_opportunity_is_proxy_not_real_effort':True,
           'photo9_unverified_perianth':True,
           'recording_differences_not_causal_biological_dispersal':True,
           'software':{'python':platform.python_version(),'sklearn':sklearn.__version__,
                       'pandas':pd.__version__,'numpy':np.__version__}}
    (out/'summary.json').write_text(json.dumps(final,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    (out/'receipt.json').write_text(json.dumps({
        'source_SHA256':SOURCE_ZIP_SHA,
        'code_SHA256':sha(__file__),
        'output_SHA256':{p.name:sha(p) for p in out.iterdir()}
        },indent=2)+'\n')
    print('MOIS_V40_EARLIER_FIRSTRECORD_PREDICTIVE_AUDIT',json.dumps({
        'samples':summary['historical_source_firstrecord_completeness_groups'],
        'metrics':summary['scores'],
        'deltas':summary['paired_comparisons']},ensure_ascii=False),flush=True)

if __name__=='__main__':main()
