"""Aggregate scientific gates for exploratory nominal-genus contribution heterogeneity.

Original licensed photo/environment archives are NOT in public CI. CI runs
synthetic invariants and checks the aggregate frozen original-source receipt.
"""
from __future__ import annotations

import csv
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/analysis/replay_mois_genus_sign_mosaic_20261015.py"
RESULT = ROOT / "results/mois_congeneric_moisture_genus_sign_mosaic_20261015.json"
LEDGER = ROOT / "results/fcp_macro_direction_20261013/direction_ledger.csv"


def source_module():
    spec = importlib.util.spec_from_file_location("mois_genus_sign_mosaic", SCRIPT)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def test_synthetic_positive_negative_zero_genus_contributions():
    model = source_module()
    rng = np.random.default_rng(57)
    samples = []
    for genus in range(20):
        for k, x in enumerate([-2, -1, 0, 1, 2]):
            y = float((k + 1) / 6)
            if genus == 0:
                y = 0.0
            elif genus == 1:
                y = 1.0
            elif genus % 2 == 0:
                y = 1.0 - y
            samples.append({
                "genus": "genus" + str(genus),
                "hue_white": y,
                "wc_bio12": 1000 + 200 * x + genus * 30,
                "wc_vapr_annual_mean": 1.1 + 0.2 * x + genus * 0.03,
                "abs_lat": float(rng.uniform(5, 65)),
                "lon_sin": float(rng.uniform(-1, 1)),
                "lon_cos": float(rng.uniform(-1, 1)),
            })
    data = pd.DataFrame(samples)
    norm = {k: {"mean": float(data[k].mean()), "sd": float(data[k].std(ddof=0))}
            for k in (*model.GEO, *model.AXES)}
    r = model.model(data, norm, minimum=3, seed=2026101603, nperm=59)
    assert r["n_genera"] == 20 and r["n_taxa"] == 100
    assert r["n_photo_white_variable_genera"] == 18
    assert r["n_photo_white_invariant_genera"] == 2
    assert r["permutations"] == 59
    for axis in model.AXES:
        x = r["axes"][axis]
        assert x["n_no_within_genus_photo_white_variation"] == 2
        assert x["n_positive_contributing_genera"] + x["n_negative_contributing_genera"] == 18
        assert 0 < x["fraction_positive_among_variable_genera"] < 1
        assert x["permuted_positive_fraction_upper_tail_p"] >= 1 / 60
        assert 0 < x["share_of_positive_covariance_from_top_tenth_of_genera"] <= 1


def test_source_receipt_matches_original_global_slopes_and_support():
    q = json.loads(RESULT.read_text(encoding="utf-8"))
    assert q["schema"] == "mois_congeneric_moisture_genus_sign_mosaic_v1"
    assert q["source_support"] == {
        "n_original_photos": 100543,
        "n_climate26_photos": 45836,
        "n_primary_species": 20694,
        "n_primary_genera": 1934,
        "n_primary_eligible_classified_photos": 38498,
    }
    assert set(q["cases"]) == {"3", "5", "10", "20"}
    assert "not separate genus regression significance tests" in q["estimand"]
    with LEDGER.open(newline="", encoding="utf-8") as stream:
        old = list(csv.DictReader(stream))
    for ax in ("wc_bio12", "wc_vapr_annual_mean"):
        qmain = q["cases"]["3"]["axes"][ax]
        found = [r for r in old if
                 r["mode"] == "climate" and r["quality"] == "all"
                 and r["region"] == "world" and r["weight"] == "genus_equal"
                 and r["adjustment"] == "genus_geo"
                 and r["anchor"] == ax and r["hue"] == "white"]
        assert len(found) == 1
        assert math.isclose(
            qmain["source_wls_slope_pp_per_environment_sd"],
            float(found[0]["direction_pp_per_env_global_sd"]),
            abs_tol=1e-6,
        )
        assert float(found[0]["maxT_ci_lo_pp"]) > 0

    for key, row in q["cases"].items():
        assert int(key) == row["minimum_photo_classifiable_species_per_genus"]
        assert row["n_genera"] >= row["n_photo_white_variable_genera"] > 0
        assert row["n_genera"] - row["n_photo_white_variable_genera"] == row["n_photo_white_invariant_genera"]
        assert row["permutations"] == 1999
        assert row["permutation_seed"] == 2026101600 + int(key)
        for axis, a in row["axes"].items():
            assert a["n_positive_contributing_genera"] + a["n_negative_contributing_genera"] == row["n_photo_white_variable_genera"]
            assert a["n_no_within_genus_photo_white_variation"] == row["n_photo_white_invariant_genera"]
            assert math.isclose(a["fraction_positive_among_all_genera"],
                                a["n_positive_contributing_genera"] / row["n_genera"], abs_tol=1e-6)
            assert math.isclose(a["fraction_positive_among_variable_genera"],
                                a["n_positive_contributing_genera"] / row["n_photo_white_variable_genera"], abs_tol=1e-6)
            assert len(a["permuted_positive_fraction_interval95"]) == 2
            assert a["permuted_positive_fraction_upper_tail_p"] >= 1 / 2000
            assert 0 < a["share_of_positive_covariance_from_top_tenth_of_genera"] < 1
            assert a["source_wls_slope_pp_per_environment_sd"] > 0


def test_paper_reports_heterogeneity_without_evolutionary_overreach():
    s = (ROOT / "docs/MOIS_GEB_RESEARCH_ARTICLE_DRAFT_20261014.md").read_text(encoding="utf8")
    si = (ROOT / "docs/MOIS_GEB_SUPPORTING_INFORMATION_20261014.md").read_text(encoding="utf8")
    memo = (ROOT / "docs/MOIS_CONGENERIC_MOISTURE_GENUS_SIGN_MOSAIC_20261015.md").read_text(encoding="utf8")
    assert "830 positive, 719 negative and 385 photo-white-invariant" in s
    assert "### S2.7 |" in si
    for phrase in ("1,549", "830", "719", "385", "53.58%", "54.29%", "60.28%", "0.0040"):
        assert phrase in si and phrase in memo
    assert "not a universal genus response" in s
    assert "not" in memo.lower() and "fitness" in memo.lower()
    assert "No original photo IDs" in si or "no raw per-genus" in si
