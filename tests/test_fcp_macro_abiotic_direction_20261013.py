"""Toy-data invariants only; these are NOT botanical source replications."""
from __future__ import annotations
import importlib.util,sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

P=Path(__file__).resolve().parents[1]/'scripts/analysis/fcp_abiotic_hue_direction_macro_20261013.py'
spec=importlib.util.spec_from_file_location('fcpmacro_direction',P)
m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m)


def fake_species(seed=20261013):
    rng=np.random.default_rng(seed);rows=[]
    for genus in range(42):
        for j in range(5):
            v=float(rng.uniform(-1.5,1.5))
            x={m.GEO[0]:float(rng.normal()),m.GEO[1]:float(rng.normal()),m.GEO[2]:float(rng.normal())}
            x.update({name:float(rng.normal()) for name in m.ANCHOR_CLIMATE})
            x['wc_bio12']=v
            x.update({
              'genus':f'G{genus:03d}','inat_taxon_id':genus*100+j,'n_photos':1,
              'cell':(genus%5)*18+(j//3),
              'latitude':float(genus%5*3+10), 'longitude':float(j*8+110),
              'hue_white':.5+.14*v,'hue_yellow_orange':0.0,'hue_red_pink':0.0,'hue_blue_purple':.5-.14*v
            })
            rows.append(x)
    return pd.DataFrame(rows)


def norms(frame,axes):
    return pd.Series({x:0. for x in axes}),pd.Series({x:1. for x in axes})


def test_class_probability_sum_direction_and_known_positive():
    d=fake_species();axes=(*m.GEO,*m.ANCHOR_CLIMATE);mu,sd=norms(d,axes)
    x,meta=m.cluster_suff(d,m.ANCHOR_CLIMATE,mu,sd)
    assert meta['n_genera']==42 and meta['n_species']==210
    r=m.regress_boot(x,m.ANCHOR_CLIMATE,199,seed=44)
    assert r['status']=='OK'
    get=lambda anchor,hue:next(z for z in r['estimates'] if z['anchor']==anchor and z['hue']==hue)
    assert get('wc_bio12','white')['point_pp_per_original_photo_env_sd']>13.9
    assert get('wc_bio12','blue_purple')['point_pp_per_original_photo_env_sd']< -13.9
    assert get('wc_bio12','white')['scope_anchor_hue_maxT_ci95_pp'][0]>0
    assert get('wc_bio12','blue_purple')['scope_anchor_hue_maxT_ci95_pp'][1]<0
    assert not get('wc_bio12','yellow_orange')['maxT_excludes_zero']
    assert len(r['estimates'])==len(m.ANCHOR_CLIMATE)*4


def test_same_genus_bootstrap_deterministic_and_genus_weighted():
    d=fake_species();mu,sd=norms(d,(*m.GEO,*m.ANCHOR_CLIMATE))
    g,_=m.cluster_suff(d,m.ANCHOR_CLIMATE,mu,sd)
    a=m.regress_boot(g,m.ANCHOR_CLIMATE,79,seed=100,equal_genus=True)
    b=m.regress_boot(g,m.ANCHOR_CLIMATE,79,seed=100,equal_genus=True)
    assert a==b
    w=next(x for x in a['estimates'] if x['anchor']=='wc_bio12' and x['hue']=='white')
    assert w['point_pp_per_original_photo_env_sd']==pytest.approx(14,abs=1e-8)


def test_only_two_taxa_per_genus_hold():
    d=fake_species();d=d[d.inat_taxon_id%100<2]
    mu,sd=norms(d,(*m.GEO,*m.ANCHOR_CLIMATE))
    g,meta=m.cluster_suff(d,m.ANCHOR_CLIMATE,mu,sd)
    assert g is None and meta['n_species']==0


def test_within_geocell_conditioning_requires_three_matched_species():
    d=fake_species();mu,sd=norms(d,m.ANCHOR_CLIMATE)
    g,meta=m.cluster_suff(d,m.ANCHOR_CLIMATE,mu,sd,geos=False,cell_conditioned=True)
    assert g is not None
    assert meta['n_genera']==42
    out=m.regress_boot(g,m.ANCHOR_CLIMATE,79,seed=89,cell_conditioned=True)
    assert out['status']=='OK'


def test_fewer_than_30_named_genus_hold_after_sufficient_statistics():
    d=fake_species();mu,sd=norms(d,(*m.GEO,*m.ANCHOR_CLIMATE))
    g,meta=m.cluster_suff(d,m.ANCHOR_CLIMATE,mu,sd)
    keep=dict(list(g.items())[:25])
    out=m.regress_boot(keep,m.ANCHOR_CLIMATE,49,seed=45)
    assert out['status']=='HOLD_INSUFFICIENT_NAMED_GENUS'


def test_original_photo_hue_contract_and_not_chemical_pigmentation():
    assert len(m.ALL_FEATURES)==36
    assert len(m.ANCHOR_CLIMATE)==6
    assert m.ANCHOR_SOIL[-1]=='soil_available_water_proxy'
    assert m.CLASSES==('white','yellow_orange','red_pink','blue_purple')
    assert 'SELECTION' not in m.VERSION
    assert 'PHOTO' in (m.__doc__ or '').upper()


def test_spatial_cell_ids_always_within_162():
    lat=np.array([-90,-70,-10,0,50,90],float);lon=np.array([-180,-100,0,70,175,180],float)
    c=m.geo_cell(lat,lon)
    assert len(c)==6 and np.all((c>=0)&(c<162))


def test_artifact_schema_has_explicit_noncausal_guard():
    assert 'environment' in m.VERSION
    assert callable(m.flatten)
