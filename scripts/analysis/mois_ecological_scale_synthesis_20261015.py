#!/usr/bin/env python3
"""Source-receipt-only ecological scale synthesis for the mois GEB manuscript.

This script recomputes arithmetic across EXISTING, independently source-replayed
aggregate receipts. It does not fit a new model, supply missing flower-photo
labels, reproduce source archives, or turn any association into selection.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"


def load(name: str) -> dict:
    return json.loads((RESULTS / name).read_text(encoding="utf-8"))


def diagnostic() -> dict:
    between = load("mois_genus_between_within_moisture_20261015.json")
    sign = load("mois_congeneric_moisture_genus_sign_mosaic_20261015.json")
    cluster = load("mois_original_genus_clustering_20261015.json")
    shapley = load("mois_env_shapley_and_genus_residual_transfer_20261015.json")
    pair = load("mois_pairwise_congeneric_photo_white_transfer_20261015.json")
    source = between["source_cohort"]
    within_rain = between["estimates"]["axes"]["wc_bio12"]["within_nominal_genus"]
    genus_rain = between["estimates"]["axes"]["wc_bio12"]["between_genus_means"]
    paired = between["estimates"]["axes"]["wc_bio12"]["between_minus_within"]
    signs = sign["cases"]["3"]["axes"]["wc_bio12"]
    clustered = cluster["comparisons"]["union_global"]
    rsq = shapley["shapley_explanatory_in_sample"]["in_sample_total_R2"]

    assert source["n_primary_taxa"] == sign["cases"]["3"]["n_taxa"] == 20694
    assert source["n_primary_genera"] == sign["cases"]["3"]["n_genera"] == 1934
    assert (clustered["n_species"], clustered["n_genera"]) == (20725, 1935)
    assert (shapley["eligible_species"], shapley["eligible_genera"]) == (16801, 1665)
    assert (pair["population"]["qualified_photo_white_taxa"], pair["population"]["qualifying_genera"]) == (16801, 1665)
    assert signs["n_positive_contributing_genera"] + signs["n_negative_contributing_genera"] + signs["n_no_within_genus_photo_white_variation"] == 1934
    assert math.isclose(genus_rain["point_pp_per_original_photo_sd"] - within_rain["point_pp_per_original_photo_sd"], paired["point_pp_per_original_photo_sd"], abs_tol=1e-9)
    assert paired["whole_genus_percentile_CI95_pp"][0] > 0
    assert math.isclose(rsq, sum(x["shapley_fraction_of_total_within_genus_white_variance"] for x in shapley["shapley_explanatory_in_sample"]["block_attribution"].values()), abs_tol=1e-12)

    fold_summaries = []
    soil_increments = []
    for fold, rec in pair["comparisons"].items():
        losses = rec["pairwise_white_difference_squared_error"]
        gain = losses["geography"] - losses["geography_moisture"]
        archived = rec["moisture_increment_over_geography"]
        assert math.isclose(gain, archived["point"], abs_tol=2e-9), fold
        ci = archived["whole_genus_fixed_prediction_bootstrap95"]
        assert ci[0] <= 0 <= ci[1], fold
        soil = rec["soil5_increment_over_fullclimate"]["point"]
        assert math.isclose(soil, losses["geography_fullclimate"] - losses["geography_fullclimate_soil5"], abs_tol=2e-9)
        soil_increments.append(soil)
        fold_summaries.append({
            "fold": fold,
            "null_pairwise_mse": losses["null"],
            "geography_pairwise_mse": losses["geography"],
            "geography_moisture_pairwise_mse": losses["geography_moisture"],
            "incremental_relative_mse_reduction_pct": 100 * gain / losses["geography"],
            "total_relative_mse_reduction_vs_null_pct": 100 * (losses["null"] - losses["geography_moisture"]) / losses["null"],
            "incremental_gain_95_interval_contains_zero": True,
            "soil5_increment_over_climate": soil,
        })
    assert len(fold_summaries) == 4
    assert min(soil_increments) < 0 < max(soil_increments)
    assert all(f["incremental_relative_mse_reduction_pct"] < 0.5 for f in fold_summaries)

    return {
        "scope": "Existing photo-classifier-conditional ecological results; arithmetic consistency only, no new raw-photo refit",
        "taxonomic_pattern": {
            "clustering_cohort_taxa": clustered["n_species"],
            "clustering_cohort_genera": clustered["n_genera"],
            "observed_genus_composition_r2_pct": 100 * clustered["metrics"]["genus_r2_species_white_fraction"]["observed"],
            "genus_size_shuffle_null_r2_pct": 100 * clustered["metrics"]["genus_r2_species_white_fraction"]["null_mean"],
        },
        "same_20694_taxa_moisture_slopes": {
            "within_genus_rain_pp_per_sd": within_rain["point_pp_per_original_photo_sd"],
            "between_genus_rain_pp_per_sd": genus_rain["point_pp_per_original_photo_sd"],
            "between_minus_within_rain_pp_per_sd": paired["point_pp_per_original_photo_sd"],
            "paired_difference_ordinary_95": paired["whole_genus_percentile_CI95_pp"],
            "ratio_between_to_within_descriptive_only": genus_rain["point_pp_per_original_photo_sd"] / within_rain["point_pp_per_original_photo_sd"],
            "positive_negative_invariant_contributing_genera": [
                signs["n_positive_contributing_genera"],
                signs["n_negative_contributing_genera"],
                signs["n_no_within_genus_photo_white_variation"],
            ],
        },
        "different_soil36_cohort_within_genus_environment": {
            "taxa": shapley["eligible_species"],
            "genera": shapley["eligible_genera"],
            "descriptive_insample_13_predictor_r2_pct": 100 * rsq,
            "unseen_genus_pairwise_four_nondependent_folds": fold_summaries,
        },
        "noncomparison": "Clustering R2, congeneric 13-predictor R2 and heldout pairwise MSE reductions have different source populations/estimands; do NOT subtract or divide them as a common variance partition or selection measure.",
        "not_measured": ["true global white-morph incidence", "pigment cost", "pollinator fitness", "natural selection", "source-linked DEM terrain"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Verify arithmetic, ranges and manuscript claim without accessing original-source ZIPs")
    args = parser.parse_args()
    outcome = diagnostic()
    if args.check:
        doc = (ROOT / "docs/MOIS_ECOLOGICAL_SCALE_SEPARATION_SYNTHESIS_20261015.md").read_text(encoding="utf-8")
        paper = (ROOT / "docs/MOIS_GEB_RESEARCH_ARTICLE_DRAFT_20261014.md").read_text(encoding="utf-8")
        for token in ("34.38%", "9.33%", "+8.555", "+2.255", "0.6202%", "0.21–0.26%", "different cohorts"):
            assert token in doc, token
        for token in ("0.316285", "0.315620", "0.21%", "different soil-complete cohort"):
            assert token in paper, token
        print("PASS source-receipt ecological scale consistency, 4 reused genus-fold intervals, and manuscript grain boundaries")
    else:
        print(json.dumps(outcome, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
