#!/usr/bin/env python3
"""Source-replayed signed genus contributions to global congeneric moisture.

Source: two SHA256-verified original FCP measured-photo/environment ZIP files.
Reuses the verified original species loader; releases aggregate-only results.
No photo IDs, genus names, species rows, source sites or photographer attributes
are written to the JSON. No pigment costs or biological fitness are estimated.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np

AXES = ("wc_bio12", "wc_vapr_annual_mean")
GEO = ("abs_lat", "lon_sin", "lon_cos")
MIN_G = (3, 5, 10, 20)
N_PERM = 1999
SEED_BASE = 2026101600


def load_original(photo_zip: Path, env_zip: Path):
    # All original-photo deduplication, environmental completeness, SHA256
    # checks and 20694/1934 source cohort assertions reside in this loader.
    source_path = Path(__file__).with_name("replay_mois_genus_between_within_moisture_20261015.py")
    spec = importlib.util.spec_from_file_location("mois_source_taxa", source_path)
    source_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(source_module)
    return source_module.load_species(photo_zip, env_zip)


def model(species, norm, minimum: int, seed: int, nperm: int = N_PERM):
    count = species.genus.value_counts()
    data = species.loc[species.genus.isin(count.index[count >= minimum])].reset_index(drop=True)
    genus, gi = np.unique(data.genus.to_numpy(), return_inverse=True)
    k = len(genus)
    size = np.bincount(gi).astype(float)
    weights = 1.0 / size[gi]
    n = len(data)

    y = data.hue_white.to_numpy(dtype=float)
    y_mean = np.bincount(gi, weights=y) / size
    yc = y - y_mean[gi]
    variable = (np.bincount(gi, weights=yc ** 2) / size) > 1e-14
    n_variable = int(variable.sum())

    G = np.stack([(data[t].to_numpy(float) - norm[t]["mean"]) / norm[t]["sd"]
                  for t in GEO], axis=1)
    g_means = np.stack([np.bincount(gi, weights=G[:, j]) / size
                        for j in range(len(GEO))], axis=1)
    Gc = G - g_means[gi]
    geo_inv = np.linalg.inv((Gc.T * weights) @ Gc)
    xr = []
    for axis in AXES:
        x = (data[axis].to_numpy(float) - norm[axis]["mean"]) / norm[axis]["sd"]
        xc = x - (np.bincount(gi, weights=x) / size)[gi]
        xr.append(xc - Gc @ (geo_inv @ (Gc.T @ (weights * xc))))
    X = np.stack(xr, axis=1)

    denominator = np.stack([
        np.bincount(gi, weights=weights * X[:, j] ** 2, minlength=k)
        for j in range(len(AXES))
    ], axis=1)
    # Frisch-Waugh-Lovell contribution: residualize focal moisture globally
    # against within-genus centered geography; KEEP y only genus-centered.
    # This ensures a genus with constant photo-white has exactly 0 numerator.
    numerator = np.stack([
        np.bincount(gi, weights=weights * X[:, j] * yc, minlength=k)
        for j in range(len(AXES))
    ], axis=1)
    fitted = 100 * numerator.sum(axis=0) / denominator.sum(axis=0)
    positive = (numerator > 1e-14) & variable[:, None]
    negative = (numerator < -1e-14) & variable[:, None]
    share = positive.sum(axis=0) / n_variable

    # One common within-genus source-label shuffle per permutation, shared
    # between both environmental axes; genus size and colour fractions fixed.
    order = np.argsort(gi, kind="stable")
    group_order = gi[order]
    yc_order = yc[order]
    w_order = weights[order]
    X_order = X[order]
    rng = np.random.default_rng(seed)
    null_share = np.empty((nperm, len(AXES)))
    null_fitted = np.empty((nperm, len(AXES)))
    for t in range(nperm):
        perm = yc_order[np.lexsort((rng.random(n), group_order))]
        group_num = np.stack([
            np.bincount(group_order, weights=w_order * X_order[:, j] * perm, minlength=k)
            for j in range(len(AXES))
        ], axis=1)
        null_share[t] = ((group_num > 1e-14) & variable[:, None]).sum(axis=0) / n_variable
        null_fitted[t] = 100 * group_num.sum(axis=0) / denominator.sum(axis=0)

    results = {}
    for j, axis in enumerate(AXES):
        n_top = max(1, int(np.ceil(.1 * k)))
        positive_cov = np.maximum(numerator[:, j], 0)
        results[axis] = {
            "source_wls_slope_pp_per_environment_sd": round(float(fitted[j]), 6),
            "n_positive_contributing_genera": int(positive[:, j].sum()),
            "n_negative_contributing_genera": int(negative[:, j].sum()),
            "n_no_within_genus_photo_white_variation": k - n_variable,
            "fraction_positive_among_all_genera": round(float(positive[:, j].sum() / k), 6),
            "fraction_positive_among_variable_genera": round(float(share[j]), 6),
            "permuted_positive_fraction_mean": round(float(null_share[:, j].mean()), 6),
            "permuted_positive_fraction_interval95": np.round(
                np.quantile(null_share[:, j], [.025, .975]), 6).tolist(),
            "permuted_positive_fraction_upper_tail_p": round(float(
                (1 + (null_share[:, j] >= share[j] - 1e-15).sum()) / (nperm + 1)), 6),
            "permuted_slope_interval95_pp": np.round(
                np.quantile(null_fitted[:, j], [.025, .975]), 6).tolist(),
            "permuted_global_slope_upper_tail_p": round(float(
                (1 + (null_fitted[:, j] >= fitted[j] - 1e-12).sum()) / (nperm + 1)), 6),
            "share_of_positive_covariance_from_top_tenth_of_genera": round(
                float(np.sort(positive_cov)[::-1][:n_top].sum() / positive_cov.sum()), 6),
        }
    return {
        "minimum_photo_classifiable_species_per_genus": minimum,
        "n_taxa": n,
        "n_genera": k,
        "n_photo_white_variable_genera": n_variable,
        "n_photo_white_invariant_genera": k - n_variable,
        "permutations": nperm,
        "permutation_seed": seed,
        "axes": results,
    }


def calculate(photo_zip: Path, env_zip: Path):
    species, norm, support = load_original(photo_zip, env_zip)
    assert (support["n_primary_taxa"], support["n_primary_genera"]) == (20694, 1934)
    cases = {str(g): model(species, norm, g, SEED_BASE + g) for g in MIN_G}
    assert abs(cases["3"]["axes"]["wc_bio12"]["source_wls_slope_pp_per_environment_sd"] - 2.254545) < 2e-6
    assert abs(cases["3"]["axes"]["wc_vapr_annual_mean"]["source_wls_slope_pp_per_environment_sd"] - 3.951850) < 2e-6
    return {
        "schema": "mois_congeneric_moisture_genus_sign_mosaic_v1",
        "date": "2026-10-15",
        "original_photo_sha256": "2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0",
        "original_environment_sha256": "6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca",
        "source_support": {"n_original_photos": 100543, "n_climate26_photos": 45836,
                           "n_primary_species": 20694, "n_primary_genera": 1934,
                           "n_primary_eligible_classified_photos": 38498},
        "estimand": "Signed contribution of each nominal genus to one worldwide genus-centred species-level photographic white/moisture WLS coefficient after the same global pooled within-genus geographic adjustment; individual signs are not separate genus regression significance tests.",
        "randomization": "1999 within-genus permutations of source species photographic-white fractions with genus sample sizes and white composition, original photo depth and source species environment/geography held fixed. Both environmental axes share each permutation. Conditional retrospective tests, not programme-wide multiple-testing adjustment.",
        "cases": cases,
        "limits": "Results are classifier-conditional and opportunistically sampled. Genus names are not a molecular phylogeny; positive covariance does not demonstrate all genera have positive ecological slopes or any water-driven selection. Top-tenth contribution is the share of POSITIVE covariance mass, not causal explained variance. No original photo IDs, named taxa, sites, genetic morphs, pollinator fitness or pigment costs are deposited."
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--photo-zip", required=True, type=Path)
    ap.add_argument("--env-zip", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    modes = ap.add_mutually_exclusive_group(required=True)
    modes.add_argument("--write", action="store_true")
    modes.add_argument("--check", action="store_true")
    args = ap.parse_args()
    report = calculate(args.photo_zip, args.env_zip)
    if args.write:
        args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        print("Wrote aggregate-only genus contribution sign mosaic", args.output)
    else:
        assert json.loads(args.output.read_text(encoding="utf-8")) == report
        print("PASS original-source matched genus contribution signs and shuffled baselines")


if __name__ == "__main__":
    main()
