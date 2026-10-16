#!/usr/bin/env python3
"""Exploratory sensitivity: visible flower white clustering by nominal genus token.

The original Dryad 'clade' can span more than one genus. This separate
same-colour-label audit groups the public final_dataset.csv by the first token
of its species column, preserving the original data and never uploading
individual source rows to GitHub. Name token != phylogenetically verified genus.
The four-column dataset contains no independent climate variable.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

from mois_external_clade_flower_colour_benchmark_20261015 import _stats, run

ORIGINAL_DOI = "10.5061/dryad.r4xgxd2sc"


def nominal_genus(s: str) -> str:
    tokens = str(s).strip().split("_", 1)
    if len(tokens) != 2 or not tokens[0].strip() or not tokens[1].strip():
        raise ValueError("Invalid species genus_species token; do not invent a genus")
    return tokens[0].strip()


def analyse(csv_path: Path, *, permutations=9999, seed=202610152330) -> dict:
    content_hash = hashlib.sha256(csv_path.read_bytes()).hexdigest()
    with csv_path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        expected = {"clade", "species", "flower_color", "fruit_color"}
        if not reader.fieldnames or not expected.issubset(reader.fieldnames):
            raise ValueError("Expected 2025 Dryad clade,species,flower_color,fruit_color")
        rows = list(reader)
    if not rows:
        raise ValueError("Empty external table")
    counts = Counter(nominal_genus(row["species"]) for row in rows)
    kept = [row for row in rows if counts[nominal_genus(row["species"])] >= 3]
    if not kept:
        raise ValueError("No external genus has >=3 source species")
    unmatched = sum(nominal_genus(row["species"]) != row["clade"].strip() for row in rows)
    by_clade = defaultdict(list)
    for row in kept:
        by_clade[row["clade"]].append(int(row["flower_color"].strip().casefold() == "white"))
    matched_clade_eta = _stats(by_clade)["clade_mean_eta_squared"]

    # Explicitly tagged temporary derivative: the only data transformation is
    # changing source clade to nominal genus and requiring >=3 genus species.
    # The temporary file is deleted before printing or returning any result.
    with tempfile.TemporaryDirectory(prefix="mois_genus_source_") as tmp:
        derived = Path(tmp) / "genus_grouped.csv"
        with derived.open("w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["clade","species","flower_color","fruit_color"])
            w.writeheader()
            for row in kept:
                w.writerow({**row, "clade": nominal_genus(row["species"])})
        outcome = run(derived, draws=permutations, seed=seed, min_clade_species=3)
    return {
        "schema": "mois_external_nominal_genus_colour_sensitivity_v1",
        "external_data_doi": ORIGINAL_DOI,
        "original_csv_sha256": content_hash,
        "original_source_rows": len(rows),
        "source_clades": len(set(row["clade"] for row in rows)),
        "all_unique_nominal_genus_tokens": len(counts),
        "source_rows_genus_token_differs_from_study_clade": unmatched,
        "minimum_species_per_nominal_genus": 3,
        "qualifying_nominal_genera": outcome["qualified_clades"],
        "qualifying_source_species": outcome["qualified_classifiable_species"],
        "white_species_in_eligible_genera": outcome["n_white_species_in_all_classifiable"],
        "white_species_equal_eligible_genera": outcome["white_species_equal_all_classifiable"],
        "white_equal_nominal_genus_means": outcome["white_equal_qualified_clade_means"],
        "nominal_genus_eta2": outcome["clade_mean_binary_white_variance_fraction"],
        "nominal_genus_near_pure_fraction_ge3": outcome["nearly_uniform_white_or_nonwhite_clade_share_ge3"],
        "n_eligible_nominal_genera_ge5": outcome["n_qualified_clades_ge5"],
        "nominal_genus_near_pure_fraction_ge5": outcome["nearly_uniform_white_or_nonwhite_clade_share_ge5"],
        "nominal_genus_near_pure_fraction_ge5_null": outcome["nearly_uniform_white_or_nonwhite_clade_share_ge5_null"],
        "source_clade_eta2_on_identical_species": matched_clade_eta,
        "group_bootstrap": outcome["clade_resampling_diagnostic"],
        "fixed_genus_size_permutations": permutations,
        "permutation_seed": seed,
        "limits": "Nominal genus is only the first species-name token, not a vetted taxonomic synonym or phylogenetic genus. This source does not independently test rainfall, water vapour, UV pigment costs, pollinators, selection, or fitness.",
    }


def self_test():
    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp) / "test.csv"
        p.write_text(
            "clade,species,flower_color,fruit_color\n"
            "C,A_1,white,red\nC,A_2,white,red\nC,A_3,red,red\n"
            "D,B_1,green,red\nD,B_2,green,red\nD,B_3,white,red\n",
            encoding="utf-8",
        )
        a = analyse(p, permutations=199, seed=11)
        assert a["qualifying_nominal_genera"] == 2
        assert a["qualifying_source_species"] == 6
        assert a["source_rows_genus_token_differs_from_study_clade"] == 6
        assert a["white_species_in_eligible_genera"] == 3
        try:
            nominal_genus("InvalidBinomial")
            raise AssertionError("Bad species name accepted")
        except ValueError:
            pass
    print("PASS synthetic nominal-genus grouping and colour clustering; not a real source replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path)
    parser.add_argument("--expected-sha256")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--permutations", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=202610152330)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if args.csv is None:
        parser.error("--csv required for actual original trait table")
    d = analyse(args.csv, permutations=args.permutations, seed=args.seed)
    if args.expected_sha256 and d["original_csv_sha256"].lower() != args.expected_sha256.lower():
        raise ValueError("Original external source SHA256 mismatch")
    output = json.dumps(d, ensure_ascii=False, sort_keys=True, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(output, encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
