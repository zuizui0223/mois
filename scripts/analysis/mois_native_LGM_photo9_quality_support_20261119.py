#!/usr/bin/env python3
"""Audit 3-GCM species palaeo-analogue colour interpretation against ROI quality.

Source is immutable verified 10,882 WCVP-native eligible plant species;
photo colours are image-derived modal labels, NOT validated petals.
Species palaeoclimate analogue is NOT observed fossil refugia.
"""
from __future__ import annotations
import argparse,hashlib,io,json
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd

SRC_SHA='955f4c0fdf2ba1be49f3195e065db52b048e28400b4a7938dc8c97222c4e280d'
SRC_NAME='WCVP_strict_native_species_LGM_3GCM_climatic_analogue_all_eligible.csv.gz'
HUE_ORDER=['white','yellow','orange','red','pink','magenta','purple','blue','bronze']
STATES=['sampling_consensus_single','provisional_single','sampled_multicolour','multicolour_candidate']
SEED=20261119

def check(ok,msg):
    if not bool(ok):raise ValueError('MOIS_NATIVE_LGM_PHOTO_QUALITY: '+msg)

def sha(v):
    if isinstance(v,(Path,str)):v=Path(v).read_bytes()
    return hashlib.sha256(v).hexdigest()

def read_source(path):
    check(sha(path)==SRC_SHA,'Original 10882 source zip changed')
    with ZipFile(path) as z:
        meta=json.loads(z.read('receipt.json'))
        check(meta['n_selected_species']==10882 and meta['selection_mode']=='all_eligible',
              'Source cohort no longer all eligible native species')
        for fn,digest in meta['outputs_sha256'].items():
            check(sha(z.read(fn))==digest,'Source member SHA changed: '+fn)
        d=pd.read_csv(io.BytesIO(z.read(SRC_NAME)),compression='gzip')
    check(len(d)==10882 and d.source_specieskey.is_unique and d.genus.nunique()==3306,
          'Source native species identity or genus count changed')
    check(set(d.hue_status)==set(STATES) and set(d.source_hue9)==set(HUE_ORDER),
          'Photo colour hue/status altered')
    check(int(d.hue_status.eq('sampling_consensus_single').sum())==884,
          'Number with independent photo observation consensus changed')
    return d,meta

def assess(d,name):
    d=d.copy()
    d['white']=d.source_hue9.eq('white')
    x=d.groupby('genus').white.agg(['min','max'])
    crossed=x.index[x['min'].ne(x['max'])]
    c=d[d.genus.isin(crossed)].copy()
    cnt=c.source_hue9.value_counts().reindex(HUE_ORDER,fill_value=0).to_dict()
    # The natural weighting of species and genera differs. Display support, not fit a
    # 12-parameter causal model to 25 genera/95 species of consensus photos.
    support={
      'quality_cohort':name,
      'n_species':int(len(d)),
      'n_genera':int(d.genus.nunique()),
      'n_mixed_white_nonwhite_genera':int(len(crossed)),
      'n_species_in_mixed_genera':int(len(c)),
      'mixed_genera_white_species':int(c.white.sum()),
      'mixed_genera_nonwhite_species':int((~c.white).sum()),
      'hue_counts_within_mixed_genera':{k:int(v) for k,v in cnt.items()},
      'source_quality_not_petal_validation':True,
    }
    if len(crossed)>0:
        c['log_3GCM_ratio']=np.log(c.LGM_consensus_vs_modern_analog_cell_ratio.astype(float))
        group=c.groupby('genus')
        diffs=[]
        groups=[]
        for _,chunk in group:
            y=chunk.log_3GCM_ratio.to_numpy(float)
            mask=chunk.white.to_numpy(bool)
            diffs.append(float(y[mask].mean()-y[~mask].mean()))
            groups.append((y,int(mask.sum())))
        differences=np.array(diffs)
        observed=float(differences.mean())
        rng=np.random.default_rng(SEED+len(crossed))
        nperm=1999
        stats=np.zeros(nperm)
        for j in range(nperm):
            pairs=[]
            for values,nwhite in groups:
                perm=rng.permutation(len(values))
                pairs.append(values[perm[:nwhite]].mean()-values[perm[nwhite:]].mean())
            stats[j]=np.mean(pairs)
        bstats=np.mean(differences[rng.integers(0,len(differences),size=(4000,len(differences)))],axis=1)
        p=float((1+(np.abs(stats)>=abs(observed)).sum())/(nperm+1))
        support['unadjusted_genus_equal_difference_log_ratio']=observed
        support['unadjusted_within_genus_colour_exchangeability_permutation_p']=p
        support['genus_resample_95pct_interval']=[float(v) for v in np.percentile(bstats,[2.5,97.5])]
        support['bootstrap_vs_permutation_can_disagree_due_to_different_inference_assumptions']=True
        support['n_colour_label_within_genus_permutations']=nperm
        support['unadjusted_permutation_is_not_climate_covariate_adjusted']=True
        support['genus_resample_not_phylogeny_or_spatial_cluster_model']=True
    return support

def run(args):
    out=Path(args.out)
    check(not out.exists(),'Refuse overwrite of old results')
    d,meta=read_source(args.source)
    cohorts={
      'all_eligible':d,
      'exclude_bronze_photo_class':d[d.source_hue9.ne('bronze')],
      'single_colour_provisional_or_consensus':d[d.hue_status.isin(('provisional_single','sampling_consensus_single'))],
      'consensus_single_only':d[d.hue_status.eq('sampling_consensus_single')]
    }
    results=[assess(v,k) for k,v in cohorts.items()]
    by_hue=pd.crosstab(d.hue_status,d.source_hue9).reindex(index=STATES,columns=HUE_ORDER,fill_value=0)
    check(by_hue.loc['sampling_consensus_single','white']==672 and
          by_hue.loc['sampling_consensus_single','yellow']==198 and
          results[-1]['n_mixed_white_nonwhite_genera']==25 and
          results[-1]['n_species_in_mixed_genera']==95 and
          results[-1]['hue_counts_within_mixed_genera']['yellow']==47,
          'Strict high-quality genus replication denominator unexpectedly changed')
    out.mkdir(parents=True)
    by_hue.to_csv(out/'WCVP_native_photo9_hue_by_quality_counts.csv')
    pd.DataFrame(results).to_csv(out/'WCVP_native_LGM_photo_quality_genus_support.csv',index=False)
    summary={
      'source_zip_sha256':SRC_SHA,
      'source_species':len(d),
      'source_genera':int(d.genus.nunique()),
      'source_original_3GCM_model':True,
      'analysis_intent':'quality denominator and exploratory within-genus hue label randomization, not a refugia hypothesis test',
      'strict_photo_quality_is_not_independently_validated_petal_colour':True,
      'outcome_is_current_climate_analogue_in_LGM_GCMs_not_paleo_species_presence':True,
      '9_hue_comparisons_under_consensus_label_are_not_identifiable':True,
      'results':results,
      'unadjusted_result_not_equivalent_to_prior_covariate_adjusted_genus_FE':True
    }
    report=out/'WCVP_native_LGM_photo_quality_support.json'
    report.write_text(json.dumps(summary,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    summary['outputs_sha256']={p.name:sha(p) for p in out.iterdir()}
    (out/'receipt.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    print('MOIS_NATIVE_LGM_PHOTO9_QUALITY_SUPPORT_VERIFIED',json.dumps({
        'species':len(d),'quality_counts':d.hue_status.value_counts().to_dict(),
        'strict_mixed_support':results[-1],
        'full_mixed_support':results[0]
    },ensure_ascii=False,allow_nan=False),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--source',required=True)
    p.add_argument('--out',required=True)
    run(p.parse_args())