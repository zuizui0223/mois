#!/usr/bin/env python3
"""External, independently scored flower-colour clade benchmark for Mois.

Requires the researcher to supply the exact published final_dataset.csv from:
  Sinnott-Armstrong et al. (2025) Dryad DOI: 10.5061/dryad.r4xgxd2sc

This script NEVER downloads from the website, guesses blocked source data,
modifies Mois historical hue labels, joins to the Mois photos, infers climate
coefficients, or exposes species/clade row-level source data in its output.
A passing --self-test is NOT an external empirical replication.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

DOI = "10.5061/dryad.r4xgxd2sc"
# The archived 2025 CSV uses 'black' and 'purple', not the website's
# descriptive 'black/dark' or 'blue/purple'. All are treated as nonwhite.
# No unrecorded synonyms or inferred whiteness are silently accepted.
COLOURS = frozenset({"black", "black/dark", "blue/purple", "purple", "green", "orange", "pink", "red", "white", "yellow"})
MISSING = frozenset({"", "na", "n/a", "none", "null", "unknown"})


def _fraction(values):
    return sum(values) / len(values)


def percentile(values, probability):
    values = sorted(values)
    index = probability * (len(values) - 1)
    low = int(math.floor(index))
    high = int(math.ceil(index))
    return values[low] + (values[high] - values[low]) * (index - low)


def load_csv(path):
    by_clade = defaultdict(list)
    encountered = {}
    n_rows = 0
    n_missing = 0
    hue_counts = Counter()
    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError("Input has no CSV header")
        columns = {col.strip().casefold(): col for col in reader.fieldnames}
        required = {"clade", "species", "flower_color"}
        if not required.issubset(columns):
            raise ValueError(f"Expected CSV columns {sorted(required)}, got {reader.fieldnames!r}")
        for row in reader:
            n_rows += 1
            clade = (row[columns["clade"]] or "").strip()
            species = (row[columns["species"]] or "").strip()
            colour = (row[columns["flower_color"]] or "").strip().casefold()
            if not clade or not species:
                raise ValueError(f"Missing clade or species at CSV data row {n_rows}")
            if species in encountered:
                raise ValueError(f"Species appears in more than one source row: {species}; use unambiguous input first")
            encountered[species] = clade
            hue_counts[colour] += 1
            if colour in MISSING:
                n_missing += 1
                continue
            if colour not in COLOURS:
                raise ValueError(f"Unexpected flower_color={colour!r} at row {n_rows}; no implicit merging of pale/white hues")
            by_clade[clade].append(1 if colour == "white" else 0)
    if not by_clade:
        raise ValueError("No classifiable flowers")
    return by_clade, n_rows, n_missing, hue_counts


def _stats(grouped):
    groups = [observations for observations in grouped.values()]
    sizes = [len(group) for group in groups]
    n = sum(sizes)
    positive = sum(sum(group) for group in groups)
    white = positive / n
    if n < 2 or white in (0.0, 1.0):
        raise ValueError("Binary-white variance is zero; no between-clade variance decomposition")
    means = [_fraction(group) for group in groups]
    between_ss = sum(size * (mean - white) ** 2 for size, mean in zip(sizes, means))
    total_ss = n * white * (1.0 - white)
    eta = between_ss / total_ss
    return {"n_species": n, "n_clades": len(groups),
            "white_species_equal": white, "white_clade_equal": _fraction(means),
            "clade_mean_eta_squared": eta,
            "near_single_hue_fraction_ge3": sum(x <= 0.1 or x >= 0.9 for x in means) / len(means)}


def run(path, *, draws=9999, seed=202610152328, min_clade_species=3):
    if draws < 99:
        raise ValueError("At least 99 permutations required for meaningful null diagnostic")
    by_clade, n_rows, n_missing, hue_counts = load_csv(path)
    eligible = {name: observations for name, observations in by_clade.items() if len(observations) >= min_clade_species}
    if len(eligible) < 2:
        raise ValueError("At least two qualifying independent sampled clades required")
    estimate = _stats(eligible)
    all_white = sum(sum(v) for v in by_clade.values()) / sum(len(v) for v in by_clade.values())
    ge5 = {c: x for c, x in eligible.items() if len(x) >= 5}
    ge5_nearpure = None
    if ge5:
        ge5_nearpure = sum((sum(v)/len(v) <= .1 or sum(v)/len(v) >= .9) for v in ge5.values()) / len(ge5)
    labels = [z for group in eligible.values() for z in group]
    sizes = [len(v) for v in eligible.values()]
    rng = random.Random(seed)
    eta_null = []
    near_null = []
    near5_null = []
    for _ in range(draws):
        rng.shuffle(labels)
        start = 0
        simulated = {}
        for i, size in enumerate(sizes):
            simulated[i] = labels[start:start+size]
            start += size
        s = _stats(simulated)
        eta_null.append(s["clade_mean_eta_squared"])
        near_null.append(s["near_single_hue_fraction_ge3"])
        if ge5:
            ge5_means = [sum(v) / len(v) for v in simulated.values() if len(v) >= 5]
            near5_null.append(sum(z <= .1 or z >= .9 for z in ge5_means) / len(ge5_means))

    def perm_result(observed, null):
        return {"observed": observed,
                "global_species_colour_shuffle_mean": sum(null) / len(null),
                "global_species_colour_shuffle_95": [percentile(null, .025), percentile(null, .975)],
                "one_sided_upper_tail_permutation_p": (1 + sum(x >= observed - 1e-14 for x in null)) / (draws + 1)}

    # Nonparametric resampling of the 51 observed SOURCE clades. This describes
    # variation across selected clades, not worldwide angiosperm uncertainty.
    bootstrap_seed, bootstrap_draws = 202610152329, 9999
    clade_rng = random.Random(bootstrap_seed)
    observed_groups = list(eligible.values())
    eta_boot, white_boot, clade_white_boot = [], [], []
    for _ in range(bootstrap_draws):
        replicate = {i: observed_groups[clade_rng.randrange(len(observed_groups))]
                     for i in range(len(observed_groups))}
        # A tiny synthetic clade bootstrap can sample only all-white or
        # only nonwhite source clades: white prevalence is defined, but
        # between-clade eta² is undefined when binary variance is zero.
        groups = list(replicate.values())
        sample_p = sum(sum(group) for group in groups) / sum(len(group) for group in groups)
        white_boot.append(sample_p)
        clade_white_boot.append(sum(sum(group) / len(group) for group in groups) / len(groups))
        if 0 < sample_p < 1:
            eta_boot.append(_stats(replicate)["clade_mean_eta_squared"])

    return {
        "schema": "mois_external_independent_clade_flower_white_replication_v1",
        "external_data_doi": DOI,
        "external_csv_sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest(),
        "external_csv_size_bytes": Path(path).stat().st_size,
        "original_external_label_values": sorted(x for x in hue_counts if x not in MISSING),
        "external_eight_hue_counts": dict(sorted((k, v) for k, v in hue_counts.items() if k not in MISSING)),
        "source_row_count": n_rows,
        "source_rows_without_classifiable_white_category": n_missing,
        "classifiable_species_in_all_clades": sum(len(x) for x in by_clade.values()),
        "n_original_named_clades": len(by_clade),
        "qualified_clade_min_species": min_clade_species,
        "qualified_classifiable_species": estimate["n_species"],
        "qualified_clades": estimate["n_clades"],
        "n_white_species_in_all_classifiable": sum(sum(v) for v in by_clade.values()),
        "white_species_equal_all_classifiable": all_white,
        "white_species_equal_qualified_clades": estimate["white_species_equal"],
        "white_equal_qualified_clade_means": estimate["white_clade_equal"],
        "clade_mean_binary_white_variance_fraction": perm_result(
            estimate["clade_mean_eta_squared"], eta_null),
        "nearly_uniform_white_or_nonwhite_clade_share_ge3": perm_result(
            estimate["near_single_hue_fraction_ge3"], near_null),
        "n_qualified_clades_ge5": len(ge5),
        "nearly_uniform_white_or_nonwhite_clade_share_ge5": ge5_nearpure,
        "nearly_uniform_white_or_nonwhite_clade_share_ge5_null": (perm_result(ge5_nearpure, near5_null) if ge5 else None),
        "clade_resampling_diagnostic": {
            "seed": bootstrap_seed, "draws": bootstrap_draws,
            "selected_clade_bootstrap95_white_species_equal": [percentile(white_boot, .025), percentile(white_boot, .975)],
            "selected_clade_bootstrap95_white_clade_equal": [percentile(clade_white_boot, .025), percentile(clade_white_boot, .975)],
            "selected_clade_bootstrap95_eta_squared": [percentile(eta_boot, .025), percentile(eta_boot, .975)],
            "warning": "Clades were not drawn at random from the world flora; intervals are conditional descriptive resampling, not population-universal confidence limits.",
        },
        "permutations": draws, "random_seed": seed,
        "comparison_warning": ("External phylogenetic clades are not Mois nominal genera. "
                               "This tests taxonomy-conditional colour clustering in a distinct "
                               "animal-pollinated/fleshy-fruited source sample; no independent "
                               "moisture gradient, BIO12 slope or reproductive fitness is estimated."),
        "scope": "Independent external taxon colour labels, not validation of Mois photo-classification/UV or a causal phylogenetic effect",
    }


def self_test():
    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp) / "synthetic.csv"
        lines = ["clade,species,flower_color,fruit_color",
                 "A,A_1,white,red", "A,A_2,white,blue/purple", "A,A_3,white,yellow",
                 "B,B_1,pink,red", "B,B_2,pink,blue/purple", "B,B_3,red,yellow"]
        p.write_text("\n".join(lines) + "\n", encoding="utf-8")
        s = run(p, draws=199, seed=1)
        assert s["source_row_count"] == 6
        assert s["qualified_clades"] == 2 and s["qualified_classifiable_species"] == 6
        assert math.isclose(s["white_species_equal_all_classifiable"], .5, abs_tol=1e-12)
        assert math.isclose(s["white_equal_qualified_clade_means"], .5, abs_tol=1e-12)
        assert math.isclose(s["clade_mean_binary_white_variance_fraction"]["observed"], 1.0, abs_tol=1e-12)
        assert math.isclose(s["nearly_uniform_white_or_nonwhite_clade_share_ge3"]["observed"], 1.0, abs_tol=1e-12)
        assert s["clade_mean_binary_white_variance_fraction"]["one_sided_upper_tail_permutation_p"] > 0
        p.write_text("clade,species,flower_color\nA,A_1,blue\nB,B_1,white\n", encoding="utf-8")
        try:
            run(p, draws=199, min_clade_species=1)
            raise AssertionError("Unknown colour was silently accepted")
        except ValueError as e:
            assert "Unexpected flower_color" in str(e)
        p.write_text("clade,species,flower_color\nA,A_1,white\nB,A_1,red\n", encoding="utf-8")
        try:
            run(p, draws=199, min_clade_species=1)
            raise AssertionError("Duplicate species counted as independent")
        except ValueError as e:
            assert "appears in more than one source row" in str(e)
        # Published source uses the exact labels 'purple' and 'black', not
        # only the website's explanatory umbrella labels.
        p.write_text("clade,species,flower_color\nA,A_1,white\nA,A_2,purple\nA,A_3,black\nB,B_1,white\nB,B_2,green\nB,B_3,pink\n", encoding="utf-8")
        z = run(p, draws=199, seed=1)
        assert z["n_white_species_in_all_classifiable"] == 2
        assert len(z["clade_resampling_diagnostic"]["selected_clade_bootstrap95_white_species_equal"]) == 2
    print("PASS synthetic external colour counts, published-category aliases, clade bootstrap, permutation and ambiguity gates (NOT external replication)")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, help="Researcher-supplied official Dryad final_dataset.csv; never scraped")
    parser.add_argument("--aggregate-output", type=Path,
                        help="Optional output JSON with only population-level statistics; no taxon/clade identities")
    parser.add_argument("--draws", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=202610152328)
    parser.add_argument("--expected-sha256", help="Optional independently recorded expected SHA256 of locally provided CSV")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if args.csv is None:
        parser.error("--csv is required unless --self-test is used")
    digest = hashlib.sha256(args.csv.read_bytes()).hexdigest()
    if args.expected_sha256 and digest.lower() != args.expected_sha256.lower():
        raise ValueError(f"External CSV SHA256 mismatch: expected {args.expected_sha256}, got {digest}")
    summary = run(args.csv, draws=args.draws, seed=args.seed)
    report = json.dumps(summary, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
    if args.aggregate_output:
        args.aggregate_output.parent.mkdir(parents=True, exist_ok=True)
        args.aggregate_output.write_text(report, encoding="utf-8")
    print(report)


if __name__ == "__main__":
    main()
