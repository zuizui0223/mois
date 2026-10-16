"""Audit genus-heldout pairwise congeneric white-photo difference forecasting.

Both source ZIPs remain owner-held; check aggregate receipt without exporting
private taxon names, photo IDs, sites or original photo media.
"""
from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/"results/mois_pairwise_congeneric_photo_white_transfer_20261015.json"
SCRIPT=ROOT/"scripts/analysis/mois_pairwise_congeneric_transfer_20261015.py"
PRIOR=ROOT/"results/mois_env_shapley_and_genus_residual_transfer_20261015.json"
SI=ROOT/"docs/MOIS_GEB_SUPPORTING_INFORMATION_20261014.md"
MEMO=ROOT/"docs/MOIS_PAIRWISE_CONGENERIC_TRANSFER_INDEPENDENT_VALIDATION_20261015.md"


def module():
    spec=importlib.util.spec_from_file_location("mois_pairwise",SCRIPT)
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_analytically_exact_all_species_pairs_without_test_colour_mean_as_input():
    m=module()
    rng=np.random.default_rng(2048)
    for n in (3,4,7,15,59):
        truth=rng.random(n)
        predicted=rng.normal(size=n)
        a=m.pairwise_group_squared_error(truth,predicted)
        b=m.direct_pairwise_squared_error(truth,predicted)
        assert math.isclose(a,b,rel_tol=0,abs_tol=1e-12)
        assert math.isclose(
            m.pairwise_group_squared_error(truth+7,predicted-2),a,
            rel_tol=0,abs_tol=1e-12
        )
    try:
        m.pairwise_group_squared_error([0.3],[0.1])
    except ValueError:
        pass
    else:
        raise AssertionError("A single species is not a congeneric pair")
    assert len(m.MODELS)==6
    assert len(m.SALTS)==4


def test_original_source_frozen_species_genera_and_distinct_pairwise_target():
    a=json.loads(RECEIPT.read_text(encoding="utf8"))
    b=json.loads(PRIOR.read_text(encoding="utf8"))
    assert a["schema"]=="mois_original_photo_white_congeneric_pairwise_transfer_v1"
    assert a["source_sha256"]["photo"]==b["original_photo_sha256"]
    assert a["source_sha256"]["environment"]==b["original_environment_sha256"]
    assert (a["population"]["qualified_photo_white_taxa"],
            a["population"]["qualifying_genera"],
            a["population"]["classified_source_photo_rows"])==(16801,1665,26163)
    assert a["ridge_penalty"]==.01
    assert a["bootstrap_draws"]==1999
    assert a["n_source_genus_assignments"]==4
    assert a["folds"]==5
    assert "NO phenotypic mean of the test genus" in a["target"]
    assert "after prediction" in a["limits"][2].lower() or "EVALUATE prediction error" in a["limits"][2]


def test_pairwise_original_genus_heldout_scores_and_robustness():
    a=json.loads(RECEIPT.read_text(encoding="utf8"))
    r=a["comparisons"]
    assert set(r)=={"original","salt202610151","salt202610152","salt202610153"}
    assert r["original"]["heldout_fold_species_counts"]==[3212,3505,3803,3306,2975]
    assert all(sum(x["heldout_fold_species_counts"])==16801 for x in r.values())
    scores=r["original"]["pairwise_white_difference_squared_error"]
    assert math.isclose(scores["geography"],.316284948033,abs_tol=1e-12)
    assert math.isclose(scores["geography_moisture"],.315620364594,abs_tol=1e-12)
    assert math.isclose(scores["null"],.31675040559,abs_tol=1e-12)
    historical=b=json.loads(PRIOR.read_text(encoding="utf8"))
    old=b["within_genus_residual_out_of_genus"]["model_results"]
    assert old["geography_moisture"]["gain_over_geography"]>0
    for key,case in r.items():
        pt=case["moisture_increment_over_geography"]
        assert pt["point"]>0
        assert pt["whole_genus_fixed_prediction_bootstrap95"][0]<0<pt["whole_genus_fixed_prediction_bootstrap95"][1]
        assert 0<pt["fraction_bootstrap_improved"]<1
        soil=case["soil5_increment_over_fullclimate"]
        assert soil["whole_genus_fixed_prediction_bootstrap95"][0]<0<soil["whole_genus_fixed_prediction_bootstrap95"][1]
        mse=case["pairwise_white_difference_squared_error"]
        for name in ("null","geography","geography_moisture","geography_fullclimate","geography_fullclimate_soil5","geography_soil5"):
            assert name in mse and .28<mse[name]<.35
        assert math.isclose(pt["point"],mse["geography"]-mse["geography_moisture"],abs_tol=2e-12)
        assert math.isclose(soil["point"],mse["geography_fullclimate"]-mse["geography_fullclimate_soil5"],abs_tol=2e-12)
    assert r["original"]["soil5_increment_over_fullclimate"]["point"]<0
    assert all(r[k]["soil5_increment_over_fullclimate"]["point"]>0 for k in ("salt202610151","salt202610152","salt202610153"))
    text=RECEIPT.read_text(encoding="utf8")
    for sensitive in ('"inat_taxon_id"','"photo_id"','"latitude"','"longitude"'):
        assert sensitive not in text


def test_supplement_makes_correct_distinction_between_absolute_residual_and_pairs():
    si=SI.read_text(encoding="utf8")
    memo=MEMO.read_text(encoding="utf8")
    assert "### S2.17 |" in si
    for term in ("16,801","1,665","26,163","+0.000665","−0.000227","+0.001506","−0.000205","4"):
        assert term in si and term in memo,term
    assert "no test-genus" in memo.lower() or "without test-genus" in memo.lower()
    assert "not" in memo.lower() and "causal" in memo.lower()
    assert "not" in si.lower()
