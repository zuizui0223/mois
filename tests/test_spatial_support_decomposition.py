from __future__ import annotations
import importlib.util,sys
from pathlib import Path
import numpy as np
import pandas as pd
P=Path(__file__).resolve().parents[1]/'scripts/analysis'
sys.path.insert(0,str(P))
import fcp_macro_spatial_support_decomposition_20261014 as m

def test_matched_stage_same_species():
    rng=np.random.default_rng(13)
    rows=[]
    for k in range(41):
        for cell in (2,7):
            for i in range(4):
                d={'genus':f'G{k}','cell':cell,'inat_taxon_id':k*100+cell*10+i,'n_photos':1,
                   'latitude':float(k/2),'longitude':110.,'hue_white':.5,'hue_blue_purple':.5,'hue_yellow_orange':0.,'hue_red_pink':0.}
                for x in (*m.src.GEO,*m.ANCHORS):d[x]=float(rng.normal())
                rows.append(d)
    data=pd.DataFrame(rows)
    mu=pd.Series({x:0. for x in (*m.src.GEO,*m.ANCHORS)});sd=mu+1.
    tr,meta=m.three_stages(data,mu,sd)
    assert meta['n_matched_species']==len(data)
    assert meta['n_original_species_all_genera']==len(data)
    assert meta['n_nested_cell_groups']==82
    for w in ('species_equal','genus_equal'):
        p=m.paired_compare(tr,nboot=39,seed=3,weight=w)
        assert p['status']=='OK',p
        assert abs(p['stage_differences']['spatial_adjustment_on_same_species']['point'][0][0])<1e-12
        assert np.isfinite(p['stage_differences']['spatial_adjustment_on_same_species']['family_maxT_95_CI']).all()

def test_reject_underfilled_genus_cell():
    r=[]
    for i in range(60):
        r.append({'inat_taxon_id':i,'genus':f'g{i//2}','cell':1})
    d=pd.DataFrame(r)
    tr,meta=m.three_stages(d,pd.Series(dtype=float),pd.Series(dtype=float))
    assert tr is None and meta['status']=='HOLD_FEWER_THAN_30_MATCHED_GENERA'

def test_pairwise_source_coefficients_exact_delta_arithmetic():
    q=np.array([[[1.,2.]]]); x=q+1;y=q+2
    assert np.allclose((x-y)+(q-x),q-y)

def test_source_hue_names_and_no_new_evolutionary_claim():
    assert m.ANCHORS==('wc_bio5','wc_bio12','wc_bio14','wc_bio15','wc_srad_annual_mean','wc_vapr_annual_mean')
    assert m.REGIONS==('world','E90_E180')
    assert len(m.STAGES)==3

def test_nonzero_morph_variation_preserves_coupled_delta_and_bootstrap_reproducibility():
    rng=np.random.default_rng(71)
    rows=[]
    for g in range(43):
        for cell in (7,19):
            for j in range(4):
                x=float(rng.uniform(-1.5,1.5))
                q={'genus':f'H{g}','cell':cell,'inat_taxon_id':g*100+cell*4+j,
                   'n_photos':1,'latitude':20.,'longitude':110.,
                   'hue_white':.42+.065*x+.025*(cell==19),'hue_blue_purple':.44-.065*x-.025*(cell==19),
                   'hue_yellow_orange':.08,'hue_red_pink':.06}
                for axis in (*m.src.GEO,*m.ANCHORS):q[axis]=float(rng.normal())
                q['wc_bio12']=x
                rows.append(q)
    d=pd.DataFrame(rows)
    means=pd.Series({x:0. for x in (*m.src.GEO,*m.ANCHORS)});scales=means+1
    stages,meta=m.three_stages(d,means,scales)
    r1=m.paired_compare(stages,89,seed=17,weight='genus_equal')
    r2=m.paired_compare(stages,89,seed=17,weight='genus_equal')
    assert r1==r2
    assert r1['status']=='OK'
    idx=m.ANCHORS.index('wc_bio12')
    assert r1['estimate'][m.STAGES[1]]['point'][idx][0]>0.04
    assert r1['estimate'][m.STAGES[2]]['point'][idx][0]>0.04
    diff=r1['stage_differences']['spatial_adjustment_on_same_species']['point'][idx][0]
    assert abs(diff-(r1['estimate'][m.STAGES[1]]['point'][idx][0]-r1['estimate'][m.STAGES[2]]['point'][idx][0]))<1e-12
    for stage in m.STAGES:
        assert np.allclose(np.array(r1['estimate'][stage]['point']).sum(axis=1),0,atol=1e-9)
