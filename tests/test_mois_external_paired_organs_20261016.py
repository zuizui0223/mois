"""Clade-controlled reproductive organ-colour validation on independent Dryad source.

Hosted GitHub CI has only safe aggregate receipts and synthetic tests.  The
2,960 original labelled species remain in the licensed/source-controlled file.
"""
from pathlib import Path
import json
import math

ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/"results/mois_external_paired_flower_fruit_20261016.json"
SCRIPT=ROOT/"scripts/analysis/mois_external_cross_organ_clade_20261016.py"
SI=ROOT/"docs/MOIS_GEB_SUPPORTING_INFORMATION_20261014.md"


def test_raw_vs_clade_adjustment_and_sample_identity():
    r=json.loads(RECEIPT.read_text(encoding="utf-8"))
    assert r["source_doi"]=="10.5061/dryad.r4xgxd2sc"
    assert r["source_sha256"]=="a253308785e4cbd0e361b3ca04cdfdae29c523eefa843375cebf8030ee0874af"
    assert (r["source_rows"],r["source_clades"])==(2960,51)
    assert r["label_shuffle_draws"]==19999 and r["clade_bootstrap_draws"]==9999
    assert set(r["outcomes"])=={"fruit_red_or_black","fruit_white"}
    red=r["outcomes"]["fruit_red_or_black"]
    assert red["raw_species_2x2_rows_flower_nonwhite_white_columns_fruit_nontarget_target"]==[[800,711],[630,819]]
    assert red["n_target_fruit"]==1530 and red["n_flower_white"]==1449
    assert red["n_strata_with_both_binary_variation"]==43
    assert math.isclose(red["raw_species_odds_ratio"],1.4627285513361463,abs_tol=1e-12)
    assert red["raw_fisher_two_sided_p"]<1e-6
    assert math.isclose(red["raw_species_conditional_fruit_prevalence_difference_pp"],9.466808620838485,abs_tol=1e-12)
    assert math.isclose(red["clade_conditional_species_weighted_slope_pp"],2.8281341470691967,abs_tol=1e-12)
    assert red["clade_conditional_species_weighted_two_sided_label_shuffle_p"]>.05
    assert red["clade_conditional_species_weighted_clade_bootstrap95_pp"][0]<0<red["clade_conditional_species_weighted_clade_bootstrap95_pp"][1]
    assert math.isclose(red["clade_conditional_equal_clade_slope_pp"],-1.452377766881527,abs_tol=1e-12)
    assert red["clade_conditional_equal_clade_two_sided_label_shuffle_p"]>.05
    white=r["outcomes"]["fruit_white"]
    assert white["raw_species_2x2_rows_flower_nonwhite_white_columns_fruit_nontarget_target"]==[[1469,42],[1405,44]]
    assert white["n_target_fruit"]==86
    assert white["raw_fisher_two_sided_p"]>.05
    assert white["clade_conditional_species_weighted_two_sided_label_shuffle_p"]>.05
    assert "Not pigment mass" in r["limits"]


def test_method_and_manuscript_keep_prediction_and_causal_boundaries():
    s=SCRIPT.read_text(encoding="utf-8")
    note=SI.read_text(encoding="utf-8")
    for token in ("rng.hypergeometric","fruit_red_or_black","fruit_white","SOURCE_SHA256","--self-test"):
        assert token in s
    for token in ("S2.19","2,960","1,530","9.47","2.83","0.125"):
        assert token in note
    assert "not" in note.lower()
