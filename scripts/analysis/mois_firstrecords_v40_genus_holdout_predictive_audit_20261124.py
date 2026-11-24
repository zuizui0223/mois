#!/usr/bin/env python3
"""MOIS first-record history: out-of-genus prediction audit of photo9 hue.

Source: pinned FirstRecords v4.0 x MOIS derived artifact from GitHub Actions
45032061554, SHA256 afa9cbe...; not original records.  This is a sensitivity
and falsification test of added predictive information, not an estimate of
natural dispersal, flower-color evolution, or introduced-species selection.

The exact test question: do photo9 labels improve held-out-genus prediction of
another first-region record within a fully observed 50y follow-up, beyond the
first-observation 25y cohort and recorded starting continent?
"""
import argparse, hashlib, io, json, platform
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss, log_loss, roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import OneHotEncoder
import sklearn

SOURCE_ZIP_SHA256='afa9cbe259ca13d039692e4aa0a0e117045655a3f3fd6ee8e89ddd01b22f41e2'
SOURCE_MEMBER='MOIS_species_recorded_region_accumulation.csv.gz'
REPEATS=5
FOLDS=5
SEEDS=[20261124,20261125,20261126,20261127,20261128]


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def guard(ok,msg):
    if not bool(ok):raise ValueError('MOIS_RECORDING_BIAS_GATE: '+str(msg))


def load_source(zip_path):
    guard(sha(zip_path)==SOURCE_ZIP_SHA256,'source artifact ZIP checksum changed')
    with ZipFile(zip_path) as z:
        receipt=json.load(z.open('receipt.json'))
        for fn,expected in receipt['files_sha256'].items():
            guard(hashlib.sha256(z.read(fn)).hexdigest()==expected,'member checksum changed: '+fn)
        data=pd.read_csv(io.BytesIO(z.read(SOURCE_MEMBER)),compression='gzip',dtype={'key':'string'})
        meta=json.load(z.open('summary.json'))
    guard(len(data)==4072 and data.key.is_unique,'species grain or row count changed')
    guard(meta['source_data_sha']['MOIS_gbif_original_SHA256']=='422f6d6489e61c6e57c5d51db614e47f613039781a901e3fd3481c13b2b085c2','original GBIF provenance changed')
    sampled=data[data.first_year.between(1800,1973)].copy()
    guard(len(sampled)==2813 and int(sampled.any_new_regions_50yr.sum())==793,'fixed 50y sample changed')
    guard(sampled.key.is_unique and sampled.genus.notna().all(),'species and genus not unique')
    guard(sampled.any_new_regions_50yr.isin([0,1]).all(),'invalid binary outcome')
    guard((sampled.new_regions_25yr<=sampled.new_regions_50yr).all() and
          (sampled.new_regions_50yr<=sampled.new_regions_100yr).all(),'non-monotone observations')
    sampled['first_cohort25']=((sampled.first_year.astype(int)-1800)//25).astype(str)
    sampled['cohort_continent']=sampled.first_cohort25+'|'+sampled.first_record_continent
    sampled['white_label']=np.where(sampled.source_hue9.eq('white'),'white','nonwhite')
    return sampled,meta


def make_folds(genus,seed,n_folds=FOLDS):
    rng=np.random.default_rng(seed)
    names=np.unique(np.asarray(genus))
    rng.shuffle(names)
    chunks=np.array_split(names,n_folds)
    allgenus=np.asarray(genus)
    for chunk in chunks:
        test=np.flatnonzero(np.isin(allgenus,chunk))
        train=np.flatnonzero(~np.isin(allgenus,chunk))
        guard(len(test)>0 and len(train)>0,'empty grouped CV fold')
        guard(len(set(allgenus[test])&set(allgenus[train]))==0,'genus leakage')
        yield train,test


def fit_oof(data,seed):
    y=data.any_new_regions_50yr.to_numpy(int)
    base=['first_cohort25','first_record_continent','cohort_continent']
    models={'history_only':base,'history_white':base+['white_label'],
            'history_all9':base+['source_hue9']}
    out={name:np.empty(len(data),dtype=float) for name in models}
    for tr,te in make_folds(data.genus.to_numpy(),seed):
        for name,cols in models.items():
            predictor=make_pipeline(OneHotEncoder(handle_unknown='ignore'),
                                    LogisticRegression(max_iter=1500,C=1.0,solver='lbfgs'))
            predictor.fit(data.iloc[tr][cols],y[tr])
            out[name][te]=predictor.predict_proba(data.iloc[te][cols])[:,1]
    for name,pred in out.items():
        guard(np.isfinite(pred).all() and ((pred>=0)&(pred<=1)).all(),name+' invalid probability')
    return out


def prob_loss(y,p):
    p=np.clip(p,1e-10,1-1e-10)
    return -(y*np.log(p)+(1-y)*np.log(1-p))


def metrics(y,p):
    return {'brier':float(brier_score_loss(y,p)),
            'logloss':float(log_loss(y,p)),
            'AUC':float(roc_auc_score(y,p))}


def bootstrap_diff_by_genus(data,individual_diff,seed=20261124,draws=3000):
    # Paired prediction-loss difference: baseline loss - colour model loss.
    # Reuses fitted OOF predictions.  Does NOT re-fit CV models per bootstrap.
    grouped=pd.DataFrame({'genus':data.genus.to_numpy(),
                          'delta':np.asarray(individual_diff)}).groupby('genus').delta.agg(['sum','count'])
    sums=grouped['sum'].to_numpy(float);counts=grouped['count'].to_numpy(float)
    rng=np.random.default_rng(seed)
    stats=np.empty(draws)
    for j in range(draws):
        idx=rng.integers(0,len(sums),len(sums))
        stats[j]=sums[idx].sum()/counts[idx].sum()
    return {'cluster':'source genus','draws':draws,
            'point':float(individual_diff.mean()),
            'bootstrap_percentile_interval_95pct':np.percentile(stats,[2.5,97.5]).tolist(),
            'conditional_on_fitted_oof_models':True,
            'no_full_training_uncertainty_or_causal_effect_claimed':True}


def analyze(dataset,tag,repeats=REPEATS):
    guard(dataset.genus.nunique()>100 and len(dataset)>300,'subset too small for prespecified genus CV')
    data=dataset.reset_index(drop=True).copy()
    y=data.any_new_regions_50yr.to_numpy(int)
    preds={model:[] for model in ('history_only','history_white','history_all9')}
    rep_metrics=[]
    for seed in SEEDS[:repeats]:
        d=fit_oof(data,seed)
        for model,p in d.items():
            preds[model].append(p)
            rep_metrics.append({'sample':tag,'seed':seed,'model':model,
                                **metrics(y,p)})
    outmeans={k:np.mean(np.stack(p),axis=0) for k,p in preds.items()}
    losses={k:np.mean(np.stack([prob_loss(y,p) for p in v]),axis=0)
            for k,v in preds.items()}
    briers={k:np.mean(np.stack([(p-y)**2 for p in v]),axis=0)
            for k,v in preds.items()}
    comparison={}
    for alt in ('history_white','history_all9'):
        d=losses['history_only']-losses[alt]
        b=briers['history_only']-briers[alt]
        comparison[alt]={'CV_genus_repeats':repeats,
                         'n_species':len(data),'n_genera':int(data.genus.nunique()),
                         'history_minus_photo_model_positive_means_photo_improves':True,
                         'delta_logloss':float(d.mean()),'delta_brier':float(b.mean()),
                         'delta_AUC_on_mean_oof_predictions':float(roc_auc_score(y,outmeans[alt])-roc_auc_score(y,outmeans['history_only'])),
                         'delta_logloss_genus_cluster_bootstrap':bootstrap_diff_by_genus(data,d,draws=3000),
                         'delta_brier_genus_cluster_bootstrap':bootstrap_diff_by_genus(data,b,seed=20261125,draws=3000)}
    per_region=[]
    for c,x in data.groupby('first_record_continent',sort=True):
        ix=x.index.to_numpy()
        record={'sample':tag,'first_continent':c,'n_species':len(ix),
                'n_genera':int(x.genus.nunique()),'mean_first_year':float(x.first_year.mean()),
                'observed_recorded_additional_region_50y_rate':float(y[ix].mean()),
                'white_photo_label_fraction':float(x.source_hue9.eq('white').mean())}
        for model in outmeans:
            record[model+'_oof_logloss']=float(losses[model][ix].mean())
            record[model+'_oof_brier']=float(briers[model][ix].mean())
        per_region.append(record)
    per_era=[]
    for c,x in data.groupby('first_cohort25',sort=True):
        ix=x.index.to_numpy()
        per_era.append({'sample':tag,'first_year_25bin':c,
                        'n_species':len(ix),'n_genera':int(x.genus.nunique()),
                        'recorded_50y_added_region_fraction':float(y[ix].mean()),
                        'white_photo_label_fraction':float(x.source_hue9.eq('white').mean()),
                        'baseline_mean_oof_pred':float(outmeans['history_only'][ix].mean())})
    return {'sample':tag,'n_species':len(data),'n_genera':int(data.genus.nunique()),
            'observed_event_rate':float(y.mean()),'mean_scores_by_model':{k:metrics(y,v) for k,v in outmeans.items()},
            'paired_predictions':comparison,'region_profile':per_region,
            'historical_cohort_profile':per_era,
            'repeat_metric_rows':rep_metrics,
            'source_generation':'Pinned species x region aggregation 2026-11-22; first records not true introductions'},\
           pd.DataFrame({'key':data.key.astype(str),'genus':data.genus,'first_record_continent':data.first_record_continent,
                          'first_year':data.first_year,'outcome_recorded_50y':y,
                          'photo9_label':data.source_hue9,
                          **{k+'_out_of_genus_mean_pred':v for k,v in outmeans.items()},
                          **{k+'_mean_oof_logloss':v for k,v in losses.items()}})


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--source',required=True)
    p.add_argument('--out',required=True)
    a=p.parse_args()
    root=Path(a.out)
    guard(not root.exists(),'refuse overwrite existing published analysis')
    data,source_meta=load_source(a.source)
    all_result,oof=analyze(data,'main_all_50y')
    no_NA=data[data.first_record_continent!='North America'].copy()
    no_AS=data[data.first_record_continent!='Asia'].copy()
    no_AFR=data[data.first_record_continent!='Africa'].copy()
    excluding={}
    for tag,x in [('exclude_NorthAmerica',no_NA),('exclude_Asia',no_AS),('exclude_Africa',no_AFR)]:
        result,_=analyze(x,tag,repeats=3)
        excluding[tag]={'n_species':result['n_species'],'n_genera':result['n_genera'],
                        'observed_event_rate':result['observed_event_rate'],
                        'paired_predictions':result['paired_predictions']}
    root.mkdir(parents=True)
    oof.to_csv(root/'species_genus_holdout_probability_predictions.csv.gz',index=False,
               compression={'method':'gzip','mtime':0})
    pd.DataFrame(all_result['region_profile']).to_csv(root/'historical_first_record_continent_profile.csv',index=False)
    pd.DataFrame(all_result['historical_cohort_profile']).to_csv(root/'historical_first_record_25yr_cohort_profile.csv',index=False)
    pd.DataFrame(all_result['repeat_metric_rows']).to_csv(root/'full_oof_repeat_metrics.csv',index=False)
    summary={'source_zip_SHA256':SOURCE_ZIP_SHA256,
             'status':'completed','primary':all_result,'leave_one_continent_out':excluding,
             'method':'5 random 5-fold genus-group splits for main; 3 repeat splits for LOCO; L2 regularized categorical logistic regression with controls firstrecord 25y cohort, first continent and interaction',
             'time_outcome':'at least one additional original firstrecord location in 1-50y after recorded first',
             'risk':'true first introduction / geography / survey effort / distance not observed',
             'source_photo9_not_independently_verified_petal_colour':True,
             'observational_and_informal_exploratory':True,
             'post_outcome_GBIF_range_area_not_used_as_predictor':True,
             'programming_versions':{'python':platform.python_version(),'pandas':pd.__version__,
                                     'numpy':np.__version__,'sklearn':sklearn.__version__}}
    (root/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n')
    receipt={'source_sha':SOURCE_ZIP_SHA256,'script_sha256':sha(__file__),
             'outputs_sha256':{i.name:sha(i) for i in sorted(root.iterdir()) if i.is_file()},
             'repeated_genus_holdout_prevents_same_genus_train_test':True,
             'uncertainty_bootstrap_reuses_OOF_predictions':True}
    (root/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print('MOIS_V40_GENUS_HOLDOUT_RESULTS',json.dumps({'primary_n':all_result['n_species'],
        'primary_scores':all_result['mean_scores_by_model'],
        'photo9_increment':all_result['paired_predictions'],
        'exclude_continents':{k:v['paired_predictions']['history_white']['delta_logloss'] for k,v in excluding.items()}},ensure_ascii=False))

if __name__=='__main__':main()
