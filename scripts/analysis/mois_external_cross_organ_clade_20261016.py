#!/usr/bin/env python3
"""Audit within-clade organ-colour associations on a curated public flower/fruit CSV.

Dryad 10.5061/dryad.r4xgxd2sc. No original species names are exported.
The main purpose is to distinguish species-pooled associations from clade-
conditional ones; this is not a mechanism, pigment cost or fitness study.
"""
import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

DOI = '10.5061/dryad.r4xgxd2sc'
SOURCE_SHA256 = 'a253308785e4cbd0e361b3ca04cdfdae29c523eefa843375cebf8030ee0874af'
FRUIT_TARGETS = {
    'fruit_red_or_black': frozenset(('red','black')),
    'fruit_white': frozenset(('white',)),
}
DRAW_SEED, BOOT_SEED = 202610162100, 202610162101


def logchoose(n, k):
    return math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1)


def exact_fisher_two_sided(n, a, b, observed):
    """Fisher's probability-ordering two-tailed test given exact margins."""
    lo, hi = max(0,a+b-n), min(a,b)
    logden = logchoose(n,a)
    lobs = logchoose(b,observed)+logchoose(n-b,a-observed)-logden
    return min(1.0, sum(math.exp(logchoose(b,k)+logchoose(n-b,a-k)-logden)
                        for k in range(lo,hi+1)
                        if logchoose(b,k)+logchoose(n-b,a-k)-logden <= lobs + 1e-10))


def read_public_source(path, enforce_digest=True):
    file=Path(path)
    sha=hashlib.sha256(file.read_bytes()).hexdigest()
    if enforce_digest and sha != SOURCE_SHA256:
        raise ValueError(f'Unexpected external source SHA256 {sha}')
    by={}
    with file.open(newline='',encoding='utf-8-sig') as fp:
        reader=csv.DictReader(fp)
        if reader.fieldnames != ['clade','species','flower_color','fruit_color']:
            raise ValueError('Unexpected external CSV schema')
        for i,row in enumerate(reader,1):
            if not row['clade'] or not row['species'] or not row['flower_color'] or not row['fruit_color']:
                raise ValueError(f'Missing clade/species/colour in source row {i}')
            if row['species'] in by:
                raise ValueError('Duplicated original species')
            by[row['species']]=(row['clade'],row['flower_color'].lower(),row['fruit_color'].lower())
    clades=defaultdict(list)
    for clade,flower,fruit in by.values():
        if fruit not in {'red','black','green','blue','yellow','orange','brown','white'}:
            raise ValueError('Unrecognised fruit hue')
        if flower not in {'white','yellow','green','purple','pink','red','orange','black'}:
            raise ValueError('Unrecognised flower hue')
        clades[clade].append((int(flower=='white'),fruit))
    return clades,sha


def conditional_stats(clades, fruit_set, permutations=19999, bootstraps=9999):
    data=[]
    for clade, records in sorted(clades.items()):
        n=len(records)
        a=sum(f for f,fruit in records)
        b=sum(fruit in fruit_set for f,fruit in records)
        overlap=sum(f and fruit in fruit_set for f,fruit in records)
        data.append((n,a,b,overlap))
    n,a,b,v=np.asarray(data,dtype=np.int64).T
    N,A,B,V=map(int,(n.sum(),a.sum(),b.sum(),v.sum()))
    if N != sum(len(v) for v in clades.values()):
        raise AssertionError('Source row mismatch')
    expected=np.asarray(a*b/n,dtype=np.float64)
    numerator=np.asarray(v-expected,dtype=np.float64)
    denominator=np.asarray(a-(a*a/n),dtype=np.float64)
    if denominator.sum() <= 0:
        raise ValueError('No within-clade flower white variation')
    species_slope=numerator.sum()/denominator.sum()
    clade_equal_slope=(numerator/n).sum()/(denominator/n).sum()
    # Exactly equivalent to shuffling fruit colour indicators within clade,
    # conditional on each clade's observed 2x2 row/column margins.
    rng=np.random.default_rng(DRAW_SEED)
    null_o=np.stack([rng.hypergeometric(int(bb),int(nn-bb),int(aa),size=permutations)
                      for nn,aa,bb in zip(n,a,b)],axis=1)
    nullnum=np.asarray(null_o-expected[None,:],dtype=np.float64)
    null_s=nullnum.sum(axis=1)/denominator.sum()
    null_eq=(nullnum/n).sum(axis=1)/(denominator/n).sum()
    p_s=(1+np.count_nonzero(np.abs(null_s)>=abs(species_slope)-1e-14))/(permutations+1)
    p_eq=(1+np.count_nonzero(np.abs(null_eq)>=abs(clade_equal_slope)-1e-14))/(permutations+1)
    bootrng=np.random.default_rng(BOOT_SEED)
    sample=bootrng.integers(0,len(n),size=(bootstraps,len(n)))
    boots=numerator[sample].sum(axis=1)/denominator[sample].sum(axis=1)
    booteq=(numerator[sample]/n[sample]).sum(axis=1)/(denominator[sample]/n[sample]).sum(axis=1)
    matrix=[[N-A-B+V,B-V],[A-V,V]]
    oratio=(matrix[0][0]*matrix[1][1])/(matrix[0][1]*matrix[1][0]) if matrix[0][1]*matrix[1][0]>0 else None
    return {
        'total_species':N,'n_clades':len(n),
        'n_flower_white':A,'n_target_fruit':B,
        'raw_species_2x2_rows_flower_nonwhite_white_columns_fruit_nontarget_target':matrix,
        'raw_species_odds_ratio':oratio,
        'raw_fisher_two_sided_p':exact_fisher_two_sided(N,A,B,V),
        'raw_species_conditional_fruit_prevalence_difference_pp':100*(V/A-(B-V)/(N-A)),
        'overlap_white_flower_target_fruit':V,
        'n_strata_with_both_binary_variation':int(np.count_nonzero((a>0)&(a<n)&(b>0)&(b<n))),
        'clade_conditional_species_weighted_slope_pp':100*species_slope,
        'clade_conditional_species_weighted_two_sided_label_shuffle_p':p_s,
        'clade_conditional_species_weighted_label_shuffle_95_pp':list(map(float,100*np.quantile(null_s,[.025,.975]))),
        'clade_conditional_species_weighted_clade_bootstrap95_pp':list(map(float,100*np.quantile(boots,[.025,.975]))),
        'clade_conditional_equal_clade_slope_pp':100*clade_equal_slope,
        'clade_conditional_equal_clade_two_sided_label_shuffle_p':p_eq,
        'clade_conditional_equal_clade_label_shuffle_95_pp':list(map(float,100*np.quantile(null_eq,[.025,.975]))),
        'clade_conditional_equal_clade_clade_bootstrap95_pp':list(map(float,100*np.quantile(booteq,[.025,.975]))),
    }


def run(path, enforce_digest=True, permutations=19999, bootstraps=9999):
    clades,sha=read_public_source(path,enforce_digest=enforce_digest)
    out={
        'schema':'mois_external_paired_flower_fruit_clade_conditional_v1',
        'source_doi':DOI,'source_sha256':sha,
        'source_rows':sum(map(len,clades.values())), 'source_clades':len(clades),
        'outcomes':{},'label_shuffle_draws':permutations,'clade_bootstrap_draws':bootstraps,
        'label_shuffle_seed':DRAW_SEED,'bootstrap_seed':BOOT_SEED,
        'null':'Fruit target labels randomly reassigned within each original sampled clade with exact fixed group margins (hypergeometric draw); not independent field replication.',
        'limits':'Exploratory fruit-red/black and fruit-white coarsenings of independently curated flower/fruit hue. Not pigment mass, fitness, causal colour allocation, within-species genetics, or independent moisture climate data. Clades are selected and not a random world flora sample.'
    }
    for name,fruit_hues in FRUIT_TARGETS.items():
        out['outcomes'][name]=conditional_stats(clades,fruit_hues,permutations,bootstraps)
    return out


def self_test():
    clades={'A':[(1,'red'),(1,'green'),(0,'red'),(0,'green')],
            'B':[(1,'black'),(1,'green'),(0,'black'),(0,'white')]}
    q=conditional_stats(clades,{'red','black'},permutations=199,bootstraps=399)
    assert q['total_species']==8 and q['n_clades']==2
    assert q['raw_species_2x2_rows_flower_nonwhite_white_columns_fruit_nontarget_target']==[[2,2],[2,2]]
    assert abs(q['clade_conditional_species_weighted_slope_pp']) < 1e-9
    assert 0<=q['clade_conditional_equal_clade_two_sided_label_shuffle_p']<=1
    print('PASS synthetic clade-conditioned organ hue analysis; no original source downloaded')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--csv',type=Path)
    p.add_argument('--out',type=Path)
    p.add_argument('--check',type=Path)
    p.add_argument('--self-test',action='store_true')
    a=p.parse_args()
    if a.self_test: self_test();return
    if a.csv is None: p.error('--csv needed for real external dataset analysis')
    result=run(a.csv)
    text=json.dumps(result,sort_keys=True,indent=2)+'\n'
    if a.out:
        a.out.parent.mkdir(parents=True,exist_ok=True)
        a.out.write_text(text,encoding='utf-8')
    if a.check:
        expected=json.loads(a.check.read_text())
        if expected != result:
            raise AssertionError('Existing source aggregate differs from independent replayer')
        print('PASS exact original public trait CSV and frozen aggregate results')
    else: print(text)


if __name__=='__main__':main()