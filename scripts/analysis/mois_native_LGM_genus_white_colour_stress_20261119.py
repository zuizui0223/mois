#!/usr/bin/env python3
"""Audit apparent native LGM climatic-analogue differences by photo9 hue.

Within-genus white vs nonwhite descriptive association, not observed palaeoflora
or independently validated botanical petal colour. Only reads prior signed true
source artifact. Pseudocount was already specified in the source analogue ratio.
"""
from __future__ import annotations
import argparse,hashlib,io,json
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd
import statsmodels.api as sm

SOURCE_ZIP_SHA='955f4c0fdf2ba1be49f3195e065db52b048e28400b4a7938dc8c97222c4e280d'
SOURCE_CSV='WCVP_strict_native_species_LGM_3GCM_climatic_analogue_all_eligible.csv.gz'
MODEL_RESULTS=('LGM_consensus_vs_modern_analog_cell_ratio',
               'CCSM4_LGM_vs_modern_analog_ratio',
               'MIROC_ESM_LGM_vs_modern_analog_ratio',
               'MPI_ESM_P_LGM_vs_modern_analog_ratio')

def checked(ok,message):
    if not bool(ok):raise ValueError('MOIS_NATIVE_LGM_GENUS_COLOUR: '+message)

def sha(data):
    if isinstance(data,(Path,str)):data=Path(data).read_bytes()
    return hashlib.sha256(data).hexdigest()

def read_source(file):
    checked(sha(file)==SOURCE_ZIP_SHA,'Original signed LGM native species ZIP changed')
    with ZipFile(file) as z:
        r=json.loads(z.read('receipt.json'))
        checked(r['n_selected_species']==10882 and r['selection_mode']=='all_eligible'
                and r['strict_WCVP_native_source_species_grid_pairs']==2271852,
                'Source selected species / WCVP native sampling does not match')
        for f,digest in r['outputs_sha256'].items():
            checked(sha(z.read(f))==digest,'Original output member SHA changed: '+f)
        d=pd.read_csv(io.BytesIO(z.read(SOURCE_CSV)),compression='gzip')
    checked(len(d)==10882 and d.source_specieskey.is_unique and
            d.genus.nunique()==3306 and d.source_hue9.nunique()==9,
            'Source native species 9-hue provenance changed')
    checked(d.n_spatial_500km_cv_folds.ge(4).all() and
            d.n_spatial_500km_cv_folds.eq(5).sum()==10877,
            'Source spatial holdout valid-fold denominator changed')
    checked(d.n_modern_native_occupied_cells.ge(50).all() and
            d.n_modern_native_500km_blocks.ge(5).all(),
            'Original native spatial eligibility violated')
    checked(np.isfinite(d[list(MODEL_RESULTS)].to_numpy(float)).all() and
            d[list(MODEL_RESULTS)].gt(0).all().all(),
            'Paleo/current analogue ratio must be positive from 0.5 offset')
    return d,r

def fitted(d,response,weight):
    checked(weight in ('species_equal','genus_equal'),'Unsupported genus weighting')
    d=d.copy();d['white']=d.source_hue9.eq('white').astype(float)
    groups=d.groupby('genus').white.agg(['min','max'])
    d=d.loc[d.genus.isin(groups.index[groups['min']<groups['max']])].copy()
    checked(len(d)>500 and d.genus.nunique()>40,'Too few crossed white/nonwhite genera')
    fields={
      'white':d.white,
      'log_native':np.log1p(d.n_modern_native_occupied_cells),
      'log_candidate_regions':np.log1p(d.n_candidate_cells_in_current_native_regions),
      'log_number_regions':np.log1p(d.n_WCVP_TDWG3_native_regions),
      'width_temperature':d.climate_envelope_T95_C-d.climate_envelope_T5_C,
      'width_log_precip':d.climate_envelope_logP95-d.climate_envelope_logP5,
      'centre_temperature':(d.climate_envelope_T95_C+d.climate_envelope_T5_C)/2,
      'centre_log_precip':(d.climate_envelope_logP95+d.climate_envelope_logP5)/2,
      'log_native_holdout_lift':np.log(np.maximum(
          d.modern_heldout_native_to_region_background_lift,1e-4)),
      'hue_status_sampled_multicolour':d.hue_status.eq('sampled_multicolour').astype(float),
      'hue_status_multicolour_candidate':d.hue_status.eq('multicolour_candidate').astype(float),
      'hue_status_consensus_single':d.hue_status.eq('sampling_consensus_single').astype(float),
    }
    X=pd.DataFrame(fields,index=d.index)
    y=np.log(d[response].to_numpy(float))
    for name in X:
        if name=='white':continue
        sd=X[name].std(ddof=1)
        if sd>1e-10:X[name]=(X[name]-X[name].mean())/sd
    x=X-X.groupby(d.genus).transform('mean')
    yc=y-pd.Series(y,index=d.index).groupby(d.genus).transform('mean').to_numpy()
    x=x.loc[:,x.std(ddof=1).gt(1e-8)]
    checked('white' in x and np.linalg.matrix_rank(x.to_numpy(float))==len(x.columns),
            'Within-genus white effect not identifiable with source controls')
    if weight=='species_equal':
        fit=sm.OLS(yc,x).fit(cov_type='cluster',
            cov_kwds={'groups':d.genus},use_t=False)
    else:
        w=1/d.groupby('genus').genus.transform('size').to_numpy(float)
        fit=sm.WLS(yc,x,weights=w).fit(cov_type='cluster',
            cov_kwds={'groups':d.genus},use_t=False)
    ci=fit.conf_int().loc['white']
    return {'response':response,'weighting':weight,
            'n_source_species':len(d),'n_crossed_genera':int(d.genus.nunique()),
            'n_covariates':len(x.columns),
            'white_minus_nonwhite_log_relative_LGM_climate_analogue_beta':float(fit.params['white']),
            'CI95_cluster_genus':[float(ci.iloc[0]),float(ci.iloc[1])],
            'p_cluster_genus':float(fit.pvalues['white']),
            'not_native_species_colour_causality':True}

def run(args):
    root=Path(args.out)
    checked(not root.exists(),'Refuse overwriting prior signed genus audit')
    d,source=read_source(args.source)
    samples={
      'all':d,
      'exclude_bronze':d.loc[d.source_hue9.ne('bronze')],
      'native_grid100plus':d.loc[d.n_modern_native_occupied_cells.ge(100)],
      'heldout_background_lift_ge1p2':d.loc[
          d.modern_heldout_native_to_region_background_lift.ge(1.2)],
      'provisional_or_consensus_single_photo_status':d.loc[
          d.hue_status.isin(('provisional_single','sampling_consensus_single'))],
      'exclude_20_most_represented_genera':d.loc[
          ~d.genus.isin(d.genus.value_counts().head(20).index)]
    }
    analyses=[]
    for label,subset in samples.items():
        for w in ('species_equal','genus_equal'):
            analyses.append({'cohort':label,**fitted(subset,MODEL_RESULTS[0],w)})
    for gcm in MODEL_RESULTS[1:]:
        analyses.append({'cohort':'all_single_GCM_sensitivity',
                         **fitted(d,gcm,'species_equal')})
    checked(len(analyses)==15 and analyses[0]['n_source_species']==6866
            and analyses[0]['n_crossed_genera']==822,
            'Source multi-colour genus comparison denominator changed')
    root.mkdir(parents=True)
    csv=root/'native_photo9_LGM_genus_white_vs_nonwhite_robustness.csv'
    pd.DataFrame(analyses).to_csv(csv,index=False)
    summary={
       'source_zip_SHA256':SOURCE_ZIP_SHA,
       'full_eligible_wcvp_native_species':len(d),
       'n_total_genera':int(d.genus.nunique()),
       'n_white_nonwhite_crossed_genera':822,
       'n_petal_hue_single_sampling_consensus':int(d.hue_status.eq('sampling_consensus_single').sum()),
       'is_exploratory_not_pre_registered':True,
       'source_image_ROI_not_verified_petal_hue':True,
       'palaeoclimate_analogues_not_fossil_refugia_or_LGM_occupancy':True,
       'only_current_WCVP_native_TDWG3_candidate_regions':True,
       'genus_cluster_SE_does_not_control_all_spatial_and_phylogenetic_correlations':True,
       'all_conditional_native_min50_current_cells_and_min5_500km_blocks':True,
       'source_3GCM_LGM_1960_1990_reconstruction':True,
       'method':'genus FE by within-group demeaning; white indicator and cohort-scaled modern climate breadth/centre, current native richness/region support, heldout niche discrimination and photo-hue status; two weighting regimes',
       'models':analyses}
    js=root/'native_photo9_LGM_genus_white_vs_nonwhite_models.json'
    js.write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    summary['output_sha256']={f.name:sha(f) for f in (csv,js)}
    (root/'receipt.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print('MOIS_NATIVE_LGM_WITHIN_GENUS_WHITE_COLOUR_STRESS_SUCCESS',json.dumps({
        'main_white':analyses[0],'genus_equal':analyses[1],
        'other_sensitivities':[{'cohort':v['cohort'],
          'beta':v['white_minus_nonwhite_log_relative_LGM_climate_analogue_beta'],
          'p':v['p_cluster_genus']} for v in analyses[2:]],
        'not_actual_paleo_refugia':True,
    },allow_nan=False),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--source',required=True)
    p.add_argument('--out',required=True)
    run(p.parse_args())