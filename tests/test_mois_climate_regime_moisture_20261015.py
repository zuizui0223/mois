"""Original-source aggregate and synthetic checks for climate-regime photo-white slopes.

The licensed FCP source ZIPs are intentionally excluded from GitHub CI; full
original-source --write / --check is done in the authorised execution context.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'scripts/analysis/mois_replay_climate_regime_moisture_20261015.py'
RECEIPT = ROOT / 'results/mois_climate_regime_moisture_20261015.json'
EARLIER = ROOT / 'results/mois_genus_between_within_moisture_20261015.json'


def module():
    spec = importlib.util.spec_from_file_location('climate_regime_white', SCRIPT)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def test_frozen_source_denominators_same_species_and_original_world_coefficients():
    q = json.loads(RECEIPT.read_text(encoding='utf8'))
    original = json.loads(EARLIER.read_text(encoding='utf8'))
    assert q['schema'] == 'mois_moisture_white_genus_climate_regime_v1'
    assert (q['n_original_photo_opportunities'],
            q['n_original_climate26_classifiable_photos'],
            q['n_primary_qualified_species'],
            q['n_primary_qualified_genera'],
            q['n_primary_qualified_classifiable_photos']) == (
                100543, 45836, 20694, 1934, 38498
            )
    for axis in ('wc_bio12', 'wc_vapr_annual_mean'):
        expected = original['estimates']['axes'][axis]['within_nominal_genus']['point_pp_per_original_photo_sd']
        assert abs(q['source_original_genus_equal_within_coefficients_pp'][axis] - expected) < 1e-8
        assert q['source_original_photo_moisture_sd'][axis] > 0
    assert q['bootstrap_replicates'] == 1999
    assert q['bootstrap_seed_base'] == 2026101608
    assert q['frozen_original_species_and_genus_support']
    assert not q['source_photo_colours_relabelled']


def test_synthetic_genus_equal_centred_wls_and_cluster_bootstrap():
    m = module()
    rng = np.random.default_rng(8)
    rows = []
    for gid in range(8):
        for _ in range(12):
            rainfall = rng.uniform(-1, 1)
            geography = rng.normal(size=3)
            rows.append({
                'genus': 'g'+str(gid),
                'white': .5 + .1*rainfall + gid*.005,
                'abs_lat': geography[0],
                'lon_sin': geography[1],
                'lon_cos': geography[2],
                'wc_bio12': rainfall,
            })
    example = pd.DataFrame(rows)
    norm = {k: {'mean':0., 'sd':1.} for k in (*m.GEO,'wc_bio12')}
    genera, aa, bb = m.moments(example, norm, 'wc_bio12')
    assert len(genera) == 8 and aa.shape == (8,4,4) and bb.shape == (8,4)
    point, draws = m.fit(aa, bb, seed=777, nboot=99)
    assert abs(point - 10) < 1e-9
    np.testing.assert_allclose(draws, np.full(99, 10.), atol=1e-8)


def test_all_source_climate_regime_estimates_are_positive_but_differences_unresolved():
    q = json.loads(RECEIPT.read_text(encoding='utf8'))
    groups = q['results_by_genus_mean_environment']
    assert set(groups) == {'wc_srad_annual_mean','wc_bio5','wc_bio12'}
    for grouping, values in groups.items():
        assert values['n_low_mean_genera'] == values['n_high_mean_genera'] == 967
        assert np.isfinite(values['median_genus_mean_in_archived_source_unit'])
        for axis, d in values['axes'].items():
            assert axis in ('wc_bio12','wc_vapr_annual_mean')
            lo = d['low_genus_pp_per_original_photo_environment_sd']
            hi = d['high_genus_pp_per_original_photo_environment_sd']
            assert lo > 0 and hi > 0
            assert d['low_percentile_genus_bootstrap95_pp'][0] > 0
            assert d['high_percentile_genus_bootstrap95_pp'][0] > 0
            diff_ci = d['high_minus_low_percentile_bootstrap95_pp']
            assert diff_ci[0] < 0 < diff_ci[1]
            assert abs(d['high_minus_low_pp']-(hi-lo)) < 1e-6
            assert 0 < d['bootstrap_fraction_high_greater_than_low'] < 1
    assert abs(groups['wc_srad_annual_mean']['axes']['wc_bio12']['low_genus_pp_per_original_photo_environment_sd']-2.076415) < 1e-6
    assert abs(groups['wc_bio5']['axes']['wc_bio12']['high_genus_pp_per_original_photo_environment_sd']-2.485231) < 1e-6
    assert abs(groups['wc_bio12']['axes']['wc_bio12']['low_genus_pp_per_original_photo_environment_sd']-3.380609) < 1e-6


def test_scientific_publication_and_rights_limits():
    receipt = RECEIPT.read_text(encoding='utf8')
    q = json.loads(receipt)
    restrictions = ' '.join(q['limits']).lower()
    for phrase in ('not original 24-hue-axis maxt', 'not causal', 'pollen or seed success',
                   'site coordinates are exported'):
        assert phrase in restrictions
    assert '"photo_id"' not in receipt and '"inat_taxon_id"' not in receipt
    si = (ROOT/'docs/MOIS_GEB_SUPPORTING_INFORMATION_20261014.md').read_text(encoding='utf8')
    memo = (ROOT/'docs/MOIS_CLIMATE_REGIME_WHITE_MOISTURE_20261015.md').read_text(encoding='utf8')
    assert '### S2.8 |' in si
    for term in ('20,694','1,934','+2.076','+2.481','+1.781','+2.485','+3.381','+2.021'):
        assert term in si and term in memo
    assert 'source photo' in memo.lower() or 'photo-site' in memo.lower()
    assert 'not a causal' in si.lower()
