#!/usr/bin/env python3
"""Compare native minus introduced photo-hue Shannon in the same fixed cells.

Each randomization assigns one hue to each original source species, shared by
all cells, both status groups, and both strict/moderate populations. The primary
comparison uses a common global genus x occupancy x six-sector-affinity null.
Positive residual difference means introduced assemblages are relatively more
concentrated. A separate native-only donor analysis has a different estimand.

The two primary comparisons are strict 9-hue and white-fixed 8-hue, each with
both assemblages having at least 30 eligible species, equal cell weighting.
Holm2 applies to those two comparisons only. Other weights, richness thresholds,
regions, the simpler null, and moderate status are sensitivity analyses.
"""
from __future__ import annotations

import argparse
import io
import json
import platform
import sys
from importlib import metadata
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd

from scripts.analysis.mois_global_hue9_community_assembly_genus_null_20261118 import (
    HUES, SEED, baseline, check, empirical_entropy, sha, signed_original,
)

SCHEMA = "mois_native_intro_samecell_global_georegion_null_20261122_v1"
NATIVE999_SHA = "c031cf8fe3f80d6cf868ee20477fdb4c550b836435b79ae16d75e771c564f951"
PRIMARY = "both_richness_ge30_equal_cells"
POPULATIONS = ("strict", "moderate")
COHORTS = tuple(p + "_" + status for p in POPULATIONS for status in ("native", "introduced"))
MODELS = {
    "photo9": ("genus_breadth", "genus_breadth_sixzones"),
    "whitefixed8": ("original_genus_breadth_whitefixed8", "sixzone_genus_breadth_whitefixed8"),
}
ATOL = 1e-10
TAIL_TOLERANCE = 1e-12
SECTOR_NAMES = ("Americas_N", "Americas_S", "Africa_Eurasia_W_N",
                "Africa_Eurasia_W_S", "Asia_Oceania_E_N", "Asia_Oceania_E_S")


def counts_for_cohort(c, labels):
    """Count fixed membership, with coherent original species identities."""
    labels = np.asarray(labels)
    check(labels.shape == c.original_labels.shape and
          np.issubdtype(labels.dtype, np.integer) and
          ((labels >= 0) & (labels < 9)).all(), "Invalid shared species hue vector")
    check(c.n_hues in (8, 9), "Invalid hue dimension")
    if c.n_hues == 8:
        check(np.array_equal(labels == 0, c.original_labels == 0),
              "White species identity changed")
        check((labels[c.species] > 0).all(), "White species entered fixed nonwhite membership")
    offset = 1 if c.n_hues == 8 else 0
    codes = c.cell.astype(np.int64) * c.n_hues + labels[c.species] - offset
    return np.bincount(codes, minlength=c.n_cells * c.n_hues).reshape(c.n_cells, c.n_hues)


def permutation_stats(observed, draws):
    draws = np.asarray(draws, dtype=float)
    check(draws.ndim == 1 and len(draws) >= 2 and np.isfinite(draws).all()
          and np.isfinite(observed), "Invalid observed or permutation statistic")
    lower_count = int(np.count_nonzero(draws <= observed + TAIL_TOLERANCE))
    upper_count = int(np.count_nonzero(draws >= observed - TAIL_TOLERANCE))
    lower, upper = (1 + lower_count) / (len(draws) + 1), (1 + upper_count) / (len(draws) + 1)
    sd = float(draws.std(ddof=1))
    return {
        "n_permutations": len(draws), "observed": float(observed),
        "null_mean": float(draws.mean()), "observed_minus_null": float(observed - draws.mean()),
        "null_sd": sd, "null_mean_monte_carlo_se": sd / np.sqrt(len(draws)),
        "null_reference_95pct_interval": np.percentile(draws, [2.5, 97.5]).tolist(),
        "interval_interpretation": "null reference interval, not an effect confidence interval",
        "n_draws_at_or_below_observed": lower_count, "n_draws_at_or_above_observed": upper_count,
        "p_rank_lower": float(lower), "p_rank_upper": float(upper),
        "p_rank_two_sided": float(min(1., 2 * min(lower, upper))),
        "inference_status": "conditional_common_global_species_label_permutation",
    }


def holm_adjust(pvalues):
    p = np.asarray(pvalues, dtype=float)
    check(p.ndim == 1 and len(p) > 0 and np.isfinite(p).all()
          and ((p >= 0) & (p <= 1)).all(), "Invalid Holm family")
    order = np.argsort(p, kind="stable")
    adjusted = np.minimum(1, np.maximum.accumulate(p[order] * np.arange(len(p), 0, -1)))
    result = np.empty(len(p))
    result[order] = adjusted
    return result


def make_pair_specs(native_richness, introduced_richness, zone, common8=None):
    """Create one common fixed mask and common weights for both status groups."""
    nr, ir, zone = map(np.asarray, (native_richness, introduced_richness, zone))
    check(nr.ndim == 1 and nr.shape == ir.shape == zone.shape and
          np.isfinite(nr).all() and np.isfinite(ir).all() and
          (nr >= 0).all() and (ir >= 0).all() and np.isin(zone, np.arange(6)).all(),
          "Invalid paired richness or geographic coordinates")
    common = (nr >= 30) & (ir >= 30)
    check(common.any(), "No jointly eligible cells for native/introduced comparison")
    specs = {}

    def add(name, mask, raw_weights=None, wording="equal weight for each same eligible cell"):
        ids = np.flatnonzero(mask)
        if not len(ids):
            return
        weights = np.ones(len(ids)) if raw_weights is None else np.asarray(raw_weights)[ids].astype(float)
        check(np.isfinite(weights).all() and (weights > 0).all(), "Invalid common positive weights")
        specs[name] = {"indices": ids, "weights": weights / weights.sum(), "weighting": wording}

    add(PRIMARY, common)
    add("both_richness_ge30_min_richness_weight", common, np.minimum(nr, ir),
        "same-cell weight proportional to min(native richness, introduced richness)")
    add("both_richness_ge100_equal_cells", (nr >= 100) & (ir >= 100))
    for k, name in enumerate(SECTOR_NAMES):
        mask = common & (zone == k)
        if mask.sum() >= 30:
            add(name + "_both_richness_ge30", mask)
    for key, parity in (("northern_hemisphere", 0), ("southern_hemisphere", 1)):
        mask = common & (zone % 2 == parity)
        if mask.sum() >= 30:
            add(key + "_both_richness_ge30", mask)
    sector_sizes = np.bincount(zone[common], minlength=6)
    if (sector_sizes >= 30).all():
        weights = np.zeros(len(zone))
        weights[common] = 1 / sector_sizes[zone[common]]
        add("six_coordinate_sectors_equal_weight", common, weights,
            "equal weight for each of six sector-specific same-cell means")
    if common8 is not None:
        common8 = np.asarray(common8)
        check(common8.shape == common.shape and common8.dtype == bool
              and not (common8 & ~common).any(), "Eight-hue mask is not a subset of nine-hue eligibility")
        add("both_nonwhite_richness_ge30_equal_cells", common8)
    return specs


def weighted_statistic(values, spec):
    values = np.asarray(values)
    ids, weights = spec["indices"], spec["weights"]
    check(len(ids) == len(weights) and len(ids) > 0 and
          np.isfinite(values[ids]).all() and abs(weights.sum() - 1) < ATOL,
          "Invalid fixed paired summary")
    return float(np.dot(values[ids], weights))


def simulate_paired(cohorts, groups, specs_by_population, permutations, seed,
                    model="", log_every=100):
    check(isinstance(permutations, (int, np.integer)) and permutations >= 2,
          "At least two randomizations are required")
    source = next(iter(cohorts.values()))
    original = source.original_labels
    mobile = np.concatenate(groups) if groups else np.array([], dtype=int)
    check(np.issubdtype(mobile.dtype, np.integer) and
          ((mobile >= 0) & (mobile < len(original))).all() and
          len(np.unique(mobile)) == len(mobile), "Overlapping or invalid donor groups")
    if source.n_hues == 8:
        check((original[mobile] > 0).all(), "White species in a white-fixed donor group")
    for c in cohorts.values():
        check(c.n_cells == source.n_cells and c.n_hues == source.n_hues
              and np.array_equal(c.original_labels, original),
              "Status populations do not share the source species catalogue")
    for population in specs_by_population:
        check(population + "_native" in cohorts and population + "_introduced" in cohorts,
              "Paired status population is missing")
    sums = {name: np.zeros(c.n_cells) for name, c in cohorts.items()}
    first199 = {name: np.zeros(c.n_cells) for name, c in cohorts.items()}
    first999 = {name: np.zeros(c.n_cells) for name, c in cohorts.items()}
    draws = {population: {stat: np.empty((permutations, 3)) for stat in specs}
             for population, specs in specs_by_population.items()}
    rng = np.random.default_rng(seed)
    for draw_id in range(permutations):
        labels = original.copy()
        for group in groups:
            labels[group] = original[group][rng.permutation(len(group))]
        entropy = {}
        for name, c in cohorts.items():
            counts = counts_for_cohort(c, labels)
            check(np.array_equal(counts.sum(axis=1), c.richness),
                  "Randomization altered fixed status-group richness")
            h = empirical_entropy(counts)
            entropy[name] = h
            sums[name] += h
            if draw_id < 199:
                first199[name] += h
            if draw_id < 999:
                first999[name] += h
        for population, specs in specs_by_population.items():
            for stat, spec in specs.items():
                native = weighted_statistic(entropy[population + "_native"], spec)
                introduced = weighted_statistic(entropy[population + "_introduced"], spec)
                draws[population][stat][draw_id] = (native, introduced, native - introduced)
        if log_every and (draw_id + 1) % log_every == 0:
            print("MOIS_NATIVE_INTRO_SAMECELL_PROGRESS", source.n_hues, model,
                  draw_id + 1, "/", permutations, flush=True)
    return {"expected": {k: v / permutations for k, v in sums.items()},
            "first199": {k: v / min(199, permutations) for k, v in first199.items()},
            "first999": {k: v / min(999, permutations) for k, v in first999.items()},
            "draws": draws}


def support_for_cohort(groups, c, genus_codes):
    mobile = np.concatenate(groups) if groups else np.array([], dtype=int)
    occupancy = np.bincount(c.species, minlength=len(c.original_labels))
    present = occupancy > 0
    affected = int(occupancy[mobile].sum())
    return {
        "global_variable_donor_strata": len(groups),
        "global_movable_species": int(len(mobile)),
        "global_genera_with_variable_strata": int(len(np.unique(genus_codes[mobile]))),
        "status_unique_species": int(present.sum()),
        "status_unique_movable_species": int(present[mobile].sum()),
        "status_species_cell_pairs": int(len(c.species)),
        "status_pairs_eligible_for_label_change": affected,
        "status_pair_fraction_eligible_for_label_change": affected / max(1, len(c.species)),
        "mobility_definition": "eligible for label change; actual hue may remain unchanged in a draw",
    }


def status_overlap(cohorts):
    result = {}
    for population in POPULATIONS:
        native, introduced = (cohorts[population + "_" + s] for s in ("native", "introduced"))
        n = np.bincount(native.species, minlength=len(native.original_labels))
        i = np.bincount(introduced.species, minlength=len(introduced.original_labels))
        intro_without_native = (i > 0) & (n == 0)
        encoded_n = native.cell.astype(np.int64) * len(n) + native.species
        encoded_i = introduced.cell.astype(np.int64) * len(n) + introduced.species
        check(not np.intersect1d(encoded_n, encoded_i).size,
              "The same source species-cell pair has both native and introduced status")
        result[population] = {
            "native_unique_species": int((n > 0).sum()),
            "introduced_unique_species": int((i > 0).sum()),
            "species_native_somewhere_and_introduced_elsewhere": int(((n > 0) & (i > 0)).sum()),
            "introduced_species_without_native_footprint_in_this_population": int(intro_without_native.sum()),
            "introduced_pairs_without_native_footprint_in_this_population": int(i[intro_without_native].sum()),
            "native_species_cell_pairs": int(n.sum()), "introduced_species_cell_pairs": int(i.sum()),
            "native_and_introduced_same_cell_species_overlap": 0,
            "introduced_species_without_native_footprint_dropped": False,
        }
    return result


def load_native999(path):
    check(sha(path) == NATIVE999_SHA, "Pinned native999 ZIP SHA changed")
    with ZipFile(path) as z:
        receipt = json.loads(z.read("receipt.json"))
        check(receipt["permutations_per_model"] == 999 and
              receipt["all_reproducibility_gates_passed"], "Native999 reference run invalid")
        for member, expected in receipt["output_SHA256"].items():
            check(sha(z.read(member)) == expected, "Native999 reference member SHA changed: " + member)
        frames = {outcome: pd.read_csv(io.BytesIO(z.read(member)), compression="gzip",
                                      float_precision="round_trip", low_memory=False)
                  for outcome, member in (("photo9", "native_photo9_georegion_cells.csv.gz"),
                                          ("whitefixed8", "native_whitefixed8_georegion_cells.csv.gz"))}
        summary = json.loads(z.read("native_georegion_summary.json"))
    return frames, summary, receipt


def code_provenance():
    root = Path(__file__).resolve().parents[2]
    paths = {Path(__file__).resolve()}
    for name, module in tuple(sys.modules.items()):
        if name.startswith("scripts.analysis.") and getattr(module, "__file__", None):
            path = Path(module.__file__).resolve()
            if path.suffix == ".py" and path.is_relative_to(root):
                paths.add(path)
    test = root / "tests/test_native_intro_samecell_georegion_null_20261122.py"
    if test.exists():
        paths.add(test)
    return {str(p.relative_to(root)): sha(p) for p in sorted(paths)}


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def run(args):
    # Geographic source dependencies are only needed for the real source run;
    # the invariant tests use NumPy/Pandas and the Python standard library.
    from scripts.analysis.mois_native_hue9_whitefixed8_georegion_null_20261121 import (
        compare_gate, load_reference, prepare_assemblage, WCVP_RECEIPT_SHA,
    )
    from scripts.analysis.mois_global_hue9_regional_pool_conditional_null_20261121 import (
        geographic_regions, genus_breadth_codes, nulls, species_affinity, THRESH,
    )
    from scripts.analysis.mois_global_photo9_whitefixed8_georegion_null_20261121 import calculate_groups
    from scripts.analysis.mois_wcvp_native_hue9_geographical_null_20261118 import (
        source_guild_rows, MANIFEST, NAME_FILE, DIST_FILE, TDWG_GEO_SHA,
    )

    check(isinstance(args.permutations, int) and args.permutations >= 999,
          "At least 999 draws are required for the canonical native999 gate")
    out = Path(args.out)
    check(not out.exists(), "Refuse overwriting a frozen same-cell analysis directory")
    old_ref, old_summary, old_receipt = load_reference(args.native_reference, "native")
    latest_ref, latest_summary, latest_receipt = load_native999(args.native999_reference)
    wcvp = Path(args.wcvp)
    check(sha(wcvp / MANIFEST) == WCVP_RECEIPT_SHA, "Pinned WCVP receipt changed")
    check(sha(args.tdwg) == TDWG_GEO_SHA, "Pinned TDWG3 geometry changed")
    pairs, species = signed_original(args.gbif, args.env)
    x = baseline(pairs, species, null_model="genus_breadth")
    coords = pd.DataFrame(np.asarray(x["cell_key"].tolist(), dtype=np.int64),
                          columns=["grid_x", "grid_y"])
    gates = [compare_gate("old_native_exact_grid_order", coords.to_numpy(),
                          old_ref[["grid_x", "grid_y"]].to_numpy(), atol=0)]
    for outcome, ref in latest_ref.items():
        gates.append(compare_gate(outcome + "_native999_exact_grid_order", coords.to_numpy(),
                                  ref[["grid_x", "grid_y"]].to_numpy(), atol=0))
    zone, longitude, latitude = geographic_regions(x)
    base, occupied, genus_codes, _ = genus_breadth_codes(x, species)
    affinity, share, _ = species_affinity(x, zone)
    groups9 = {model: nulls(x, species, model, zone, affinity, base, occupied,
                           genus_codes, share) for model in MODELS["photo9"]}
    groups8, genus_codes8, affinity8 = calculate_groups(x, species, zone)
    check(np.array_equal(genus_codes8, genus_codes) and np.array_equal(affinity8, affinity),
          "White-fixed donor model geography differs")
    source_cohorts, source_audit = source_guild_rows(pairs, species, coords, wcvp, args.tdwg, coords)
    check(source_audit == old_summary["source_WCVP_exact_match_and_TDWG3_status"]
          and source_audit == latest_summary["native_source_audit"],
          "Frozen WCVP status source audit changed")
    prepared = {
        n: {name: prepare_assemblage(name, source_cohorts[name]["cell_index"],
                                     source_cohorts[name]["species_index"], x["hue"], x["ngrids"], n)
            for name in COHORTS} for n in (9, 8)
    }
    for n, outcome in ((9, "photo9"), (8, "whitefixed8")):
        for name, c in prepared[n].items():
            check(np.array_equal(counts_for_cohort(c, x["hue"]), c.counts),
                  "Independent counting does not reproduce frozen source counter")
            if n == 9:
                gates += [compare_gate(name + "_observed_H9", c.observed,
                                       old_ref[name + "_Shannon_observed"].to_numpy(float)),
                          compare_gate(name + "_richness9", c.richness,
                                       old_ref[name + "_n_photo9_species"].to_numpy(int), atol=0)]
            if name.endswith("_native"):
                gates += [compare_gate(name + "_native999_observed_" + outcome, c.observed,
                                       latest_ref[outcome][name + "_observed_Shannon"].to_numpy(float)),
                          compare_gate(name + "_native999_richness_" + outcome, c.richness,
                                       latest_ref[outcome][name + "_richness"].to_numpy(int), atol=0)]
    masks8 = {p: (prepared[8][p + "_native"].richness >= 30) &
                 (prepared[8][p + "_introduced"].richness >= 30) for p in POPULATIONS}
    specs = {n: {p: make_pair_specs(prepared[n][p + "_native"].richness,
                                   prepared[n][p + "_introduced"].richness, zone,
                                   common8=masks8[p] if n == 9 else None)
                 for p in POPULATIONS} for n in (9, 8)}
    for p, expected in (("strict", 3106), ("moderate", 3602)):
        gates.append(compare_gate(p + "_samecell9_richness30_denominator",
                                  [len(specs[9][p][PRIMARY]["indices"])], [expected], atol=0))
    gates.append(compare_gate("strict_samecell9_richness100_denominator",
                              [len(specs[9]["strict"]["both_richness_ge100_equal_cells"]["indices"])],
                              [1299], atol=0))
    overlap = {"photo9": status_overlap(prepared[9]), "whitefixed8": status_overlap(prepared[8])}
    code_sha = code_provenance()
    results, cells_out, draw_frames, statistic_rows, composition_rows, masks_out = {}, {}, [], [], [], []
    print("MOIS_NATIVE_INTRO_SAMECELL_SOURCES_VERIFIED", json.dumps({
        "source_audit": source_audit,
        "paired_cells": {str(n): {p: len(s[PRIMARY]["indices"]) for p, s in specs[n].items()}
                         for n in (9, 8)},
        "overlap": overlap["photo9"],
    }, sort_keys=True), flush=True)
    for n, outcome, all_groups in ((9, "photo9", groups9), (8, "whitefixed8", groups8)):
        cohorts = prepared[n]
        seed = SEED if n == 9 else SEED + 803
        frame = coords.copy()
        frame["center_longitude"], frame["center_latitude"] = longitude, latitude
        frame["coordinate_sector"] = np.asarray(SECTOR_NAMES)[zone]
        hue_names = HUES if n == 9 else HUES[1:]
        for name, c in cohorts.items():
            frame[name + "_richness"] = c.richness
            frame[name + "_observed_Shannon"] = c.observed
            for j, hue in enumerate(hue_names):
                frame[name + "_observed_" + hue + "_species"] = c.counts[:, j]
        for p, pspecs in specs[n].items():
            native, introduced = cohorts[p + "_native"], cohorts[p + "_introduced"]
            for stat, spec in pspecs.items():
                ids = spec["indices"]
                masks_out.append(pd.DataFrame({
                    "outcome": outcome, "population": p, "statistic_name": stat,
                    "grid_x": coords.grid_x.iloc[ids].to_numpy(),
                    "grid_y": coords.grid_y.iloc[ids].to_numpy(),
                    "cell_index": ids, "normalized_common_weight": spec["weights"],
                    "native_richness": native.richness[ids],
                    "introduced_richness": introduced.richness[ids],
                }))
                for j, hue in enumerate(hue_names):
                    pn = native.counts[:, j] / np.maximum(native.richness, 1)
                    pi = introduced.counts[:, j] / np.maximum(introduced.richness, 1)
                    mean_n, mean_i = weighted_statistic(pn, spec), weighted_statistic(pi, spec)
                    composition_rows.append({
                        "outcome": outcome, "population": p, "statistic_name": stat,
                        "hue": hue, "n_cells": len(ids), "weighting": spec["weighting"],
                        "native_mean_cell_species_fraction": mean_n,
                        "introduced_mean_cell_species_fraction": mean_i,
                        "native_minus_introduced_percentage_points": 100 * (mean_n - mean_i),
                        "inference_status": "observed_composition_descriptive_no_hue_specific_p",
                    })
        model_results = {}
        for model in MODELS[outcome]:
            groups = all_groups[model]
            simulation = simulate_paired(cohorts, groups, specs[n], args.permutations, seed, model)
            for p in POPULATIONS:
                name = p + "_native"
                gates.append(compare_gate(name + "_native999_expected_" + model,
                                          simulation["first999"][name],
                                          latest_ref[outcome][name + "_expected_Shannon_" + model].to_numpy(float)))
            if n == 9 and model == "genus_breadth":
                for name in COHORTS:
                    gates.append(compare_gate(name + "_baseline_first199_expected",
                                              simulation["first199"][name],
                                              old_ref[name + "_Shannon_null"].to_numpy(float)))
            population_results = {}
            for name, c in cohorts.items():
                frame[name + "_expected_Shannon_" + model] = simulation["expected"][name]
                frame[name + "_observed_minus_expected_" + model] = c.observed - simulation["expected"][name]
            for p, pspecs in specs[n].items():
                native, introduced = cohorts[p + "_native"], cohorts[p + "_introduced"]
                summaries = {}
                for stat, spec in pspecs.items():
                    values = simulation["draws"][p][stat]
                    observed_n = weighted_statistic(native.observed, spec)
                    observed_i = weighted_statistic(introduced.observed, spec)
                    components = {"native": permutation_stats(observed_n, values[:, 0]),
                                  "introduced": permutation_stats(observed_i, values[:, 1])}
                    direct = permutation_stats(observed_n - observed_i, values[:, 2])
                    direct.update(observed_contrast=direct["observed"],
                                  null_mean_contrast=direct["null_mean"],
                                  residual_difference=direct["observed_minus_null"],
                                  contrast_sign="native minus introduced",
                                  positive_residual_interpretation="introduced relatively more concentrated",
                                  synchronous_species_labels=True,
                                  n_cells=len(spec["indices"]), weighting=spec["weighting"],
                                  mean_native_richness=weighted_statistic(native.richness, spec),
                                  mean_introduced_richness=weighted_statistic(introduced.richness, spec),
                                  primary_test=(p == "strict" and model == MODELS[outcome][1] and stat == PRIMARY),
                                  population_role="primary" if p == "strict" else "overlapping sensitivity")
                    expected_n = weighted_statistic(simulation["expected"][p + "_native"], spec)
                    expected_i = weighted_statistic(simulation["expected"][p + "_introduced"], spec)
                    check(abs(expected_n - components["native"]["null_mean"]) <= ATOL
                          and abs(expected_i - components["introduced"]["null_mean"]) <= ATOL
                          and abs((components["native"]["observed_minus_null"] -
                                   components["introduced"]["observed_minus_null"]) -
                                  direct["residual_difference"]) <= ATOL,
                          "Same-cell component/contrast/cell expectation identity failed")
                    summaries[stat] = {"native": components["native"],
                                       "introduced": components["introduced"],
                                       "native_minus_introduced": direct}
                    draw_frames.append(pd.DataFrame({
                        "outcome": outcome, "model": model, "population": p,
                        "statistic_name": stat, "seed": seed,
                        "draw_id": np.arange(1, args.permutations + 1),
                        "native_Shannon": values[:, 0], "introduced_Shannon": values[:, 1],
                        "native_minus_introduced_Shannon": values[:, 2],
                    }))
                population_results[p] = {
                    "exchangeability_support": {s: support_for_cohort(groups, cohorts[p + "_" + s], genus_codes)
                                                for s in ("native", "introduced")},
                    "statistics": summaries,
                }
            model_results[model] = {"seed": seed, "populations": population_results}
            print("MOIS_NATIVE_INTRO_SAMECELL_MODEL_COMPLETED", outcome, model,
                  json.dumps(population_results["strict"]["statistics"][PRIMARY], sort_keys=True), flush=True)
        results[outcome] = {"n_hues": n, "primary_model": MODELS[outcome][1],
                            "models": model_results}
        cells_out[outcome] = frame
    primary_stats = [results[o]["models"][MODELS[o][1]]["populations"]["strict"]
                     ["statistics"][PRIMARY]["native_minus_introduced"] for o in ("photo9", "whitefixed8")]
    adjusted = holm_adjust([s["p_rank_two_sided"] for s in primary_stats])
    family = []
    for outcome, stats, corrected in zip(("photo9", "whitefixed8"), primary_stats, adjusted):
        stats["p_rank_two_sided_holm2"] = float(corrected)
        family.append({"outcome": outcome, "population": "strict", "model": MODELS[outcome][1],
                       "statistic_name": PRIMARY, "p_rank_two_sided": stats["p_rank_two_sided"],
                       "p_rank_two_sided_holm2": float(corrected)})
    for outcome, result in results.items():
        for model, model_result in result["models"].items():
            for population, pop_result in model_result["populations"].items():
                for statistic, stats in pop_result["statistics"].items():
                    for component, record in stats.items():
                        row = {"outcome": outcome, "model": model, "population": population,
                               "statistic_name": statistic, "component": component, **record}
                        lo, hi = row.pop("null_reference_95pct_interval")
                        row["null_reference_2p5_percentile"] = lo
                        row["null_reference_97p5_percentile"] = hi
                        statistic_rows.append(row)
    check(code_provenance() == code_sha, "Analysis code changed during execution")
    sources = {"GBIF_zip": Path(args.gbif), "ENV_zip": Path(args.env),
               "native199_reference": Path(args.native_reference),
               "native999_reference": Path(args.native999_reference),
               "WCVP_receipt": wcvp / MANIFEST, "WCVP_names": wcvp / NAME_FILE,
               "WCVP_distribution": wcvp / DIST_FILE, "TDWG3_geojson": Path(args.tdwg)}
    source_sha = {k: {"filename": v.name, "SHA256": sha(v)} for k, v in sources.items()}
    runtime = {"python": platform.python_version(), "platform": platform.platform()}
    for package in ("numpy", "pandas", "pyarrow", "geopandas", "shapely", "pyproj"):
        runtime[package] = metadata.version(package)
    summary = {
        "schema": SCHEMA, "design_status": "choices fixed before this run; not formal preregistration",
        "source_species_unit": "original accepted GBIF specieskey; WCVP supplies regional status only",
        "source_species": len(species), "source_species_cell_pairs": len(pairs), "source_cells": x["ngrids"],
        "permutations_per_model": args.permutations,
        "primary_population": "strict", "sensitivity_population": "moderate; overlapping, not independent",
        "primary_statistic": "equal-cell mean(native Shannon minus introduced Shannon), both richness>=30",
        "positive_residual_difference": "introduced relatively more concentrated under the same null",
        "donor_pool": "full GBIF occupied species shared by every population and status group",
        "breadth_and_affinity_source": "full GBIF footprints, including introduced occurrences",
        "occupancy_bins": ["1", "2-5", "6-20", "21-100", ">=101"],
        "affinity_threshold": THRESH, "sectors_are_coordinate_partitions_not_biological_realms": True,
        "sector_names": list(SECTOR_NAMES),
        "region_sensitivity_minimum_paired_cells": 30,
        "same_species_has_one_shared_hue_per_draw": True,
        "same_cell_mask_and_weights_applied_to_both_statuses": True,
        "whitefixed8": "original white species IDs and original nonwhite memberships fixed",
        "nine_vs_eight_mask_control": "additional 9-hue result on the identical 8-hue-eligible cells",
        "weight_sensitivity": "min(native richness, introduced richness); weights fixed in every draw",
        "rank_p": "(1 + tail count)/(B+1), two-sided = min(1, 2*min(tails)), ties tolerance1e-12",
        "uncertainty": "null reference intervals are not effect confidence intervals",
        "primary_multiplicity": {"method": "Holm2", "family": family,
                                 "scope": "two same-cell outcome comparisons; not all project analyses"},
        "other_statistics": "sensitivity analyses; raw p-values do not receive confirmatory status",
        "colour_compositions": "species fractions within each cell, then same fixed common cell weights",
        "status_overlap": overlap, "source_status_audit": source_audit,
        "source_photo_hues_are_not_verified_petal_colours": True,
        "status_is_TDWG3_region_attribute_not_individual_wildness": True,
        "residence_time_dispersal_rate_speciation_or_adaptation_not_measured": True,
        "source_files": source_sha, "code_SHA256": code_sha, "runtime": runtime,
        "references": {"native199": old_receipt, "native999_ZIP_SHA256": NATIVE999_SHA,
                       "native999_signed_members": latest_receipt["output_SHA256"]},
        "outcomes": results,
    }
    out.mkdir(parents=True, exist_ok=False)
    gzip = {"method": "gzip", "mtime": 0}
    for outcome, frame in cells_out.items():
        frame.to_csv(out / ("native_intro_" + outcome + "_samecell_cells.csv.gz"),
                     index=False, compression=gzip)
    pd.concat(draw_frames, ignore_index=True).to_csv(out / "native_intro_samecell_permutation_draws.csv.gz",
                                                   index=False, compression=gzip)
    pd.concat(masks_out, ignore_index=True).to_csv(out / "native_intro_samecell_masks_weights.csv.gz",
                                                 index=False, compression=gzip)
    pd.DataFrame(statistic_rows).to_csv(out / "native_intro_samecell_statistics.csv", index=False)
    pd.DataFrame(composition_rows).to_csv(out / "native_intro_samecell_colour_composition.csv", index=False)
    write_json(out / "native_intro_samecell_summary.json", summary)
    write_json(out / "reproducibility_gates.json", {"all_passed": all(g["passed"] for g in gates),
               "native999_full_cell_vector_gates": 8, "native_introduced_baseline199_vector_gates": 4,
               "gates": gates})
    write_json(out / "receipt.json", {"schema": SCHEMA, "status": "complete",
               "source_files": source_sha, "code_SHA256": code_sha, "runtime": runtime,
               "output_SHA256": {p.name: sha(p) for p in sorted(out.iterdir()) if p.is_file()},
               "command_arguments": vars(args), "all_reproducibility_gates_passed": True,
               "draw_id_interpretation": "1-based; shared across status and population within an outcome/model only"})
    print("MOIS_NATIVE_INTRO_SAMECELL_DONE", json.dumps({
        "output": str(out), "primary_Holm2_family": family, "source_gates_passed": len(gates),
        "paired_statistic_draw_rows": sum(len(f) for f in draw_frames),
    }, sort_keys=True), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("gbif", "env", "wcvp", "tdwg", "native-reference", "native999-reference", "out"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--permutations", type=int, default=999)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
