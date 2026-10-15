"""Scientific and arithmetic gates for original-source geography/vapour coefficients.

CI tests aggregate-only receipts and the mathematical omitted-variable identity;
full SHA256-checked original photo and site ZIPs remain owner-held and excluded
from public GitHub Actions.
"""
from __future__ import annotations
import importlib.util
import json
import math
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/'scripts/analysis/mois_geography_vapour_ovb_20261015.py'
REPORT=ROOT/'results/mois_geography_vapour_ovb_20261015.json'
GEOGRAPHY=ROOT/'results/mois_photo_colour_geography_synthesis_20261015.json'
SUPPORT=ROOT/'docs/MOIS_GEB_SUPPORTING_INFORMATION_20261014.md'
MEMO=ROOT/'docs/MOIS_GEOGRAPHY_VAPOUR_SHARED_GRADIENT_20261015.md'
MS=ROOT/'docs/MOIS_GEB_RESEARCH_ARTICLE_DRAFT_20261014.md'

def module():
    spec=importlib.util.spec_from_file_location('mois_ovb',SCRIPT)
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m

def test_weighted_genus_centered_synthetic_omitted_variable_identity():
    m=module()
    rng=np.random.default_rng(951)
    rows=[]
    for g in range(12):
        for i in range(7):
            lat=rng.normal()
            elev=.4*lat+rng.normal()
            rain=.4*lat+.3*rng.normal()
            vap=-.9*lat-.7*elev+.7*rain+.2*rng.normal()
            lon1=rng.normal();lon2=rng.normal()
            white=.2*lat-.1*elev+1.6*vap+.05*rain+.3*lon1+.13*lon2+g*.2+rng.normal(scale=.04)
            rows.append(dict(genus='g'+str(g),abs_lat=lat,lon_sin=lon1,lon_cos=lon2,
                             wc_elevation_m=elev,wc_bio12=rain,
                             wc_vapr_annual_mean=vap,white=white))
    df=pd.DataFrame(rows)
    means={k:{'mean':float(df[k].mean()),'sd':float(df[k].std(ddof=0))} for k in m.VARS}
    G,B=m.genus_moments(df,means)
    A=G.sum(axis=0); rhs=B.sum(axis=0)
    for reduced,full in [('geography_only','geography_plus_vapour'),
                         ('geography_plus_rain','geography_both_moisture')]:
        ix=np.array([m.VARS.index(k) for k in m.MODELS[reduced]])
        fx=np.array([m.VARS.index(k) for k in m.MODELS[full]])
        betar=100*np.linalg.solve(A[np.ix_(ix,ix)],rhs[ix])
        betaf=100*np.linalg.solve(A[np.ix_(fx,fx)],rhs[fx])
        gamma=np.linalg.solve(A[np.ix_(ix,ix)],A[ix,m.VARS.index('wc_vapr_annual_mean')])
        np.testing.assert_allclose(betar-betaf[:-1],betaf[-1]*gamma,rtol=0,atol=1e-10)

def test_original_source_replay_receipt_and_expected_shifts():
    q=json.loads(REPORT.read_text(encoding='utf8'))
    old=json.loads(GEOGRAPHY.read_text(encoding='utf8'))
    assert q['schema']=='mois_within_genus_geography_vapour_omitted_variable_identity_v1'
    assert q['source_photo_sha256']==old['source_original_photo_sha256']
    assert q['source_environment_sha256']==old['source_original_environment_sha256']
    assert (q['source_support']['taxa'],q['source_support']['genera'],
            q['source_support']['original_classified_photos_in_eligible_taxa'])==(20694,1934,38498)
    assert q['paired_genus_bootstrap']['repetitions']==1999
    assert q['paired_genus_bootstrap']['seed']==202610151216
    assert math.isclose(q['original_genus_equal_BIO12_pp_per_sd_reproduced'],2.254545235,abs_tol=1e-8)
    a=q['identity_comparisons']['add_vapour_to_geography']
    b=q['identity_comparisons']['add_vapour_to_geography_and_rain']
    assert math.isclose(a['vapour_variation_explained_by_reduced_geo_and_climate_predictors_R2'],.745324772,abs_tol=1e-8)
    assert math.isclose(b['vapour_variation_explained_by_reduced_geo_and_climate_predictors_R2'],.785732407,abs_tol=1e-8)
    assert math.isclose(a['vapour_collinearity_vif'],3.926569567,abs_tol=1e-8)
    assert b['vapour_collinearity_vif']>a['vapour_collinearity_vif']
    for case,geoA,geoB in ((a,0,2),(b,1,3)):
        for key,oldkey in (('abs_lat','white_vs_10deg_higher_abs_lat_pp'),
                           ('wc_elevation_m','white_vs_1000m_higher_elevation_pp')):
            d=case['geographic_coef_changes'][key]
            q0=old['genus_equal_20694_species_geographically_adjusted'][geoA][oldkey]
            q1=old['genus_equal_20694_species_geographically_adjusted'][geoB][oldkey]
            assert math.isclose(d['reduced_white_coef_pp'],q0,abs_tol=1e-5)
            assert math.isclose(d['with_vapour_white_coef_pp'],q1,abs_tol=1e-5)
            assert math.isclose(d['with_vapour_white_coef_pp']-d['reduced_white_coef_pp'],
                                d['increase_after_adding_vapour_pp'],abs_tol=1e-8)
            assert math.isclose(d['exact_omitted_variable_identity_contribution_pp'],
                                d['increase_after_adding_vapour_pp'],abs_tol=1e-12)
            assert d['source_conditional_paired_genus_bootstrap95_for_coef_increase_pp'][0]>0
    assert math.isclose(a['geographic_coef_changes']['wc_elevation_m']['increase_after_adding_vapour_pp'],2.91359055,abs_tol=1e-8)
    assert math.isclose(a['geographic_coef_changes']['abs_lat']['increase_after_adding_vapour_pp'],2.59978642,abs_tol=1e-8)
    assert q['within_genus_genus_equal_correlations']['absolute_latitude_vs_vapour']<0
    assert q['within_genus_genus_equal_correlations']['elevation_vs_vapour']<0
    assert q['within_genus_genus_equal_correlations']['precipitation_vs_vapour']>0
    assert 'not ecological causal mediation' in q['interpretation']
    assert '"photo_id"' not in REPORT.read_text(encoding='utf8')
    assert '"inat_taxon_id"' not in REPORT.read_text(encoding='utf8')

def test_editorial_limit_on_ecological_interpretation():
    supp=SUPPORT.read_text(encoding='utf8')
    ms=MS.read_text(encoding='utf8')
    memo=MEMO.read_text(encoding='utf8')
    assert '### S2.10 |' in supp
    for v in ('74.53%','78.57%','VIF 3.927','+2.914','+2.600','1,999'):
        assert v in supp
    assert '+2.60 points per 10°' in ms
    assert 'do **not** show causal vapour mediation' in ms
    assert 'coefficient-change algebra' in memo
    assert '*not* a counterfactual intervention' in memo
    assert len(ms[ms.index('## 1 | Introduction'):ms.index('## Data and code availability')].split())<=5000
