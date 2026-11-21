#!/usr/bin/env python3
"""Source-matched native photo-hue nulls with global geographic donor pools.

The estimand is the mean photographic hue Shannon entropy of recorded
WCVP-TDWG3-native species in fixed 100-km cells. Strict native status is the
primary population; moderate native status is an overlapping sensitivity
population. These are source photo-hue classes, not verified petal colours.

Reuse the original global genus/breadth and six-coordinate-sector pool
definitions without changing their donor species, grouping order, or seeds.
Every draw assigns one species label shared by all cells and both cohorts.
Eight-hue analyses fix original white IDs and original nonwhite memberships.
The matched-mobility models condition on observed hue heterogeneity and are
reported descriptively, without an exact randomization-test interpretation.

The first 199 draws reproduce three full-global 9-hue and three full-global
8-hue references, and the prior strict/moderate native 9-hue baseline. The
main run defaults to 999 draws. No frozen output directory is overwritten.
"""
from __future__ import annotations

import argparse
import io
import json
import platform
import sys
from dataclasses import dataclass, field
from importlib import metadata
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd

from scripts.analysis.mois_global_hue9_community_assembly_genus_null_20261118 import (
    ENV_SHA, GBIF_SHA, HUES, SEED, baseline, check, counts_by_cell,
    empirical_entropy, sha, signed_original,
)
from scripts.analysis.mois_global_hue9_regional_pool_conditional_null_20261121 import (
    THRESH, ZONE_NAMES, geographic_regions, genus_breadth_codes, nulls,
    species_affinity,
)
from scripts.analysis.mois_global_photo9_whitefixed8_georegion_null_20261121 import (
    calculate_groups,
)
from scripts.analysis.mois_nonwhite8_assembly_white_fixed_null_20261118 import (
    count_chromatic,
)
from scripts.analysis.mois_wcvp_native_hue9_geographical_null_20261118 import (
    ARCTIC_LAT, DIST_FILE, MANIFEST, NAME_FILE, TDWG_GEO_SHA, source_guild_rows,
)


SCHEMA = "mois_native_hue9_whitefixed8_global_georegion_null_20261121_v1"
REFERENCE_PERMUTATIONS = 199
TAIL_TOLERANCE = 1e-12
REPRODUCIBILITY_ATOL = 1e-10
NATIVE_SHA = "cbe1eef381e3e1a0062dbf5134b047f947338b7c499d2ecc0cd22ee8855909fd"
GLOBAL9_SHA = "3516b3d232f2bbe8d0358de4d9cca677d3eba6c255b80a87f60c43c954dc6d33"
GLOBAL8_SHA = "b1243a8505efb7837781341383ebf41cd0f72c0a0b94c6a2052b1b2f85087683"
WCVP_RECEIPT_SHA = "a3af69d99e389c20eac4d8bc2c5ab4a586a58c401d65743d9dbac1189baaaf9c"
NATIVE_CELL = "global_WCVP_TDWG3_native_only_photo9_island_alpine_polar_south_cells.csv.gz"
GLOBAL9_CELL = "MOIS_global_photo9_100km_conditional_geography_18577_cells.csv.gz"
GLOBAL8_CELL = "global_100km_nonwhite8_whitefixed_geography_conditioned_cells.csv.gz"
MODELS9 = ("genus_breadth", "matched_exchangeable_genus_breadth", "genus_breadth_sixzones")
MODELS8 = ("original_genus_breadth_whitefixed8", "matched_exchangeable_whitefixed8",
           "sixzone_genus_breadth_whitefixed8")
COHORTS = ("strict_native", "moderate_native")
SECTOR_GROUPS = tuple(name + "_richness_ge30" for name in ZONE_NAMES)
EQUAL_SECTORS = "six_coordinate_sectors_equal_weight_mean"
NORTH = "northern_hemisphere_richness_ge30"
SOUTH = "southern_hemisphere_richness_ge30"
ARCTIC = "Arctic_lat_ge66p562_richness_ge30"
BOREAL = "north_boreal_lat50to66p562_richness_ge30"
CONTRASTS = {
    "south_minus_north": (SOUTH, NORTH),
    "Arctic_minus_boreal": (ARCTIC, BOREAL),
}


@dataclass
class Assemblage:
    name: str
    n_hues: int
    n_cells: int
    original_labels: np.ndarray
    cell: np.ndarray
    species: np.ndarray
    source_pair_count_before_nonwhite_filter: int
    counts: np.ndarray
    richness: np.ndarray
    observed: np.ndarray
    indices: dict[str, np.ndarray] = field(default_factory=dict)


def count_assemblage(c: Assemblage, labels: np.ndarray) -> np.ndarray:
    """Count a fixed species/cell membership using one catalogue label vector."""
    labels = np.asarray(labels)
    check(labels.shape == c.original_labels.shape and
          np.issubdtype(labels.dtype, np.integer) and
          ((labels >= 0) & (labels < 9)).all(), "Invalid species hue labels")
    x = {"cell": c.cell, "species": c.species, "ngrids": c.n_cells}
    if c.n_hues == 9:
        return counts_by_cell(x, labels)
    check(np.array_equal(labels == 0, c.original_labels == 0),
          "White species identities changed in the white-fixed null")
    # The unchanged legacy counter requires min/max of a nonempty vector.
    # Empty native/chromatic cells (or a completely empty toy population)
    # have entropy zero and cannot pass any richness inclusion threshold.
    if not len(c.species):
        return np.zeros((c.n_cells, 8), dtype=np.int64)
    return count_chromatic(x, labels, c.cell, c.species)


def prepare_assemblage(name, cell, species, labels, n_cells, n_hues):
    labels = np.asarray(labels)
    cell, species = np.asarray(cell), np.asarray(species)
    check(isinstance(n_cells, (int, np.integer)) and n_cells > 0 and
          n_hues in (8, 9), "Invalid grid count or hue dimension")
    check(labels.ndim == 1 and len(labels) > 0 and
          np.issubdtype(labels.dtype, np.integer) and
          ((labels >= 0) & (labels < 9)).all(), "Invalid original hue labels")
    check(cell.ndim == species.ndim == 1 and cell.shape == species.shape and
          np.issubdtype(cell.dtype, np.integer) and
          np.issubdtype(species.dtype, np.integer), "Invalid species-cell index vectors")
    check(((cell >= 0) & (cell < n_cells)).all() and
          ((species >= 0) & (species < len(labels))).all(),
          "Species-cell index outside source catalogue")
    encoded = cell.astype(np.int64) * len(labels) + species
    check(len(np.unique(encoded)) == len(encoded), "Duplicate species-cell pair")
    source_n = len(cell)
    if n_hues == 8:
        keep = labels[species] > 0
        cell, species = cell[keep], species[keep]
    c = Assemblage(name, n_hues, int(n_cells), labels.copy(),
                   cell.astype(np.int32), species.astype(np.int32), source_n,
                   np.empty((0, n_hues)), np.empty(0), np.empty(0))
    c.counts = count_assemblage(c, labels)
    c.richness = c.counts.sum(axis=1)
    c.observed = empirical_entropy(c.counts)
    check(int(c.richness.sum()) == len(c.cell), "Fixed source richness changed")
    return c


def permutation_statistics(observed, draws, descriptive_only=False):
    draws = np.asarray(draws, dtype=float)
    check(draws.ndim == 1 and len(draws) >= 2 and np.isfinite(draws).all()
          and np.isfinite(observed), "Invalid finite permutation statistic")
    mean = float(draws.mean())
    sd = float(draws.std(ddof=1))
    lower = float((1 + np.count_nonzero(draws <= observed + TAIL_TOLERANCE)) /
                  (len(draws) + 1))
    upper = float((1 + np.count_nonzero(draws >= observed - TAIL_TOLERANCE)) /
                  (len(draws) + 1))
    legacy = float((1 + np.count_nonzero(
        np.abs(draws - mean) >= abs(observed - mean) - TAIL_TOLERANCE
    )) / (len(draws) + 1))
    return {
        "n_permutations": len(draws),
        "observed": float(observed),
        "null_mean": mean,
        "observed_minus_null": float(observed - mean),
        "null_sd": sd,
        "null_mean_monte_carlo_se": sd / float(np.sqrt(len(draws))),
        "null_reference_95pct_interval": np.percentile(draws, [2.5, 97.5]).tolist(),
        "interval_interpretation": "null reference interval, not an effect confidence interval",
        "p_rank_lower": None if descriptive_only else lower,
        "p_rank_upper": None if descriptive_only else upper,
        "p_rank_two_sided": None if descriptive_only else min(1.0, 2 * min(lower, upper)),
        "legacy_simulated_mean_centered_p_two_sided_diagnostic_only":
            None if descriptive_only else legacy,
        "inference_status": "descriptive_only_outcome_selected_mobility" if descriptive_only
            else "conditional_global_species_label_permutation",
    }


def paired_contrast(observed_a, observed_b, draws_a, draws_b, descriptive_only=False):
    a, b = np.asarray(draws_a, dtype=float), np.asarray(draws_b, dtype=float)
    check(a.shape == b.shape and a.ndim == 1, "Paired draws are not aligned")
    draws = a - b
    stats = permutation_statistics(observed_a - observed_b, draws, descriptive_only)
    stats.update(observed_contrast=stats["observed"],
                 null_mean_contrast=stats["null_mean"],
                 residual_difference=stats["observed_minus_null"],
                 synchronous_within_model_draws=True)
    return stats, draws


def adjust_pvalues(pvalues):
    p = np.asarray(pvalues, dtype=float)
    check(p.ndim == 1 and len(p) > 0 and np.isfinite(p).all() and
          ((p >= 0) & (p <= 1)).all(), "Invalid multiplicity family p-values")
    order = np.argsort(p, kind="stable")
    ordered = p[order]
    n = len(p)
    holm = np.minimum(1, np.maximum.accumulate(ordered * np.arange(n, 0, -1)))
    bh = np.minimum(1, np.minimum.accumulate(
        (ordered * n / np.arange(1, n + 1))[::-1])[::-1])
    h, b = np.empty(n), np.empty(n)
    h[order], b[order] = holm, bh
    return {"holm": h, "bh_exploratory": b}


def compare_gate(name, actual, reference, atol=REPRODUCIBILITY_ATOL):
    a, b = np.asarray(actual), np.asarray(reference)
    check(a.shape == b.shape and np.isfinite(a).all() and np.isfinite(b).all(),
          "Source gate shape/finite failure: " + name)
    difference = float(np.max(np.abs(a.astype(float) - b.astype(float)), initial=0))
    check(difference <= atol, f"Source gate {name} differs by {difference:.17g}")
    return {"name": name, "passed": True, "n_values": int(a.size),
            "maximum_absolute_difference": difference, "absolute_tolerance": atol,
            "relative_tolerance": 0}


def group_means(c, values):
    result = {name: float(values[ii].mean()) for name, ii in c.indices.items()}
    if all(name in result for name in SECTOR_GROUPS):
        result[EQUAL_SECTORS] = float(np.mean([result[name] for name in SECTOR_GROUPS]))
    return result


def simulate_model(global_source, cohorts, shuffle_groups, permutations, seed, model,
                   reference_global_expected=None, reference_native_expected=None):
    """Draw once per species, then evaluate every population with that draw."""
    check(isinstance(permutations, (int, np.integer)) and permutations >= 2,
          "At least two draws are needed for Monte Carlo SD")
    original = global_source.original_labels
    mobile = np.concatenate(shuffle_groups) if shuffle_groups else np.array([], dtype=int)
    check(np.issubdtype(mobile.dtype, np.integer) and
          ((mobile >= 0) & (mobile < len(original))).all() and
          len(np.unique(mobile)) == len(mobile), "Overlapping or invalid shuffle groups")
    if global_source.n_hues == 8:
        check((original[mobile] > 0).all(), "White IDs included in chromatic shuffle groups")
    for c in cohorts.values():
        check(c.n_cells == global_source.n_cells and c.n_hues == global_source.n_hues
              and np.array_equal(c.original_labels, original),
              "Cohorts do not share the same source species label catalogue")
        check(bool(c.indices) and all(len(ii) > 0 for ii in c.indices.values()),
              "A summary group has no eligible native cells")
    if reference_global_expected is not None or reference_native_expected is not None:
        check(permutations >= REFERENCE_PERMUTATIONS,
              "The signed reference gate requires at least 199 permutations")
    rng = np.random.default_rng(seed)
    sums = {name: np.zeros(c.n_cells) for name, c in cohorts.items()}
    first = {name: np.zeros(c.n_cells) for name, c in cohorts.items()}
    draws = {name: {key: np.empty(permutations) for key in group_means(c, c.observed)}
             for name, c in cohorts.items()}
    full_sum = np.zeros(global_source.n_cells)
    gates = []
    for i in range(permutations):
        labels = original.copy()
        # Preserve group order and RNG calls exactly as in both source scripts.
        for group in shuffle_groups:
            labels[group] = original[group][rng.permutation(len(group))]
        if i < REFERENCE_PERMUTATIONS:
            full_sum += empirical_entropy(count_assemblage(global_source, labels))
        for name, c in cohorts.items():
            counts = count_assemblage(c, labels)
            check(np.array_equal(counts.sum(axis=1), c.richness),
                  "A shuffle changed fixed native or nonwhite species membership")
            entropy = empirical_entropy(counts)
            sums[name] += entropy
            if i < REFERENCE_PERMUTATIONS:
                first[name] += entropy
            for key, value in group_means(c, entropy).items():
                draws[name][key][i] = value
        if i + 1 == REFERENCE_PERMUTATIONS:
            if reference_global_expected is not None:
                gates.append(compare_gate(
                    f"full_global_{global_source.n_hues}_{model}_first199_expected",
                    full_sum / REFERENCE_PERMUTATIONS, reference_global_expected))
            for name, expected in (reference_native_expected or {}).items():
                gates.append(compare_gate(f"{name}_{model}_native9_first199_expected",
                                          first[name] / REFERENCE_PERMUTATIONS, expected))
            if reference_global_expected is not None:
                print("MOIS_NATIVE_NULL_FIRST199_GATES_PASSED", model,
                      json.dumps(gates, sort_keys=True), flush=True)
        if (i + 1) % 100 == 0:
            print("MOIS_NATIVE_GEO_NULL_PROGRESS", global_source.n_hues,
                  model, i + 1, "/", permutations, flush=True)
    n_first = min(permutations, REFERENCE_PERMUTATIONS)
    return {
        "expected": {name: value / permutations for name, value in sums.items()},
        "expected_first199": {name: value / n_first for name, value in first.items()},
        "draws": draws,
        "global_expected_first199": full_sum / n_first,
        "n_global_entropy_evaluations": n_first,
        "gates": gates,
    }


def summary_indices(richness, zone, latitude):
    richness, zone, latitude = map(np.asarray, (richness, zone, latitude))
    check(richness.ndim == 1 and richness.shape == zone.shape == latitude.shape
          and np.isfinite(richness).all() and (richness >= 0).all()
          and np.isin(zone, np.arange(6)).all(), "Invalid fixed geographic summary inputs")
    eligible = richness >= 30
    masks = {"richness_ge30": eligible, "richness_ge100": richness >= 100}
    masks.update({name: eligible & (zone == k) for k, name in enumerate(SECTOR_GROUPS)})
    # geographic_regions uses signed projected northing for hemisphere. This
    # remains defined for all 27 polar edge cells with undefined latitude.
    masks[NORTH] = eligible & (zone % 2 == 0)
    masks[SOUTH] = eligible & (zone % 2 == 1)
    masks[ARCTIC] = eligible & np.isfinite(latitude) & (latitude >= ARCTIC_LAT)
    masks[BOREAL] = eligible & np.isfinite(latitude) & (latitude >= 50) & (latitude < ARCTIC_LAT)
    return {name: np.flatnonzero(mask) for name, mask in masks.items()}


def exchangeability_support(groups, global_source, native, genus_codes):
    mobile = np.unique(np.concatenate(groups)) if groups else np.array([], dtype=int)
    native_occ = np.bincount(native.species, minlength=len(native.original_labels))
    global_occ = np.bincount(global_source.species, minlength=len(native.original_labels))
    native_mobile = mobile[native_occ[mobile] > 0]
    affected = int(native_occ[mobile].sum())
    n_pairs = len(native.species)
    return {
        "global_variable_donor_strata": len(groups),
        "global_movable_species": int(len(mobile)),
        "global_genera_with_variable_donor_strata": int(len(np.unique(genus_codes[mobile]))),
        "global_species_cell_pairs_eligible_for_label_change": int(global_occ[mobile].sum()),
        "global_species_cell_pair_fraction_eligible_for_label_change":
            float(global_occ[mobile].sum() / len(global_source.species)),
        "native_unique_species": int(np.count_nonzero(native_occ)),
        "native_unique_movable_species": int(len(native_mobile)),
        "native_unique_movable_species_fraction":
            float(len(native_mobile) / max(1, np.count_nonzero(native_occ))),
        "native_species_cell_pairs": n_pairs,
        "native_species_cell_pairs_eligible_for_label_change": affected,
        "native_species_cell_pair_fraction_eligible_for_label_change":
            float(affected / max(1, n_pairs)),
        "mobility_definition": "eligible for label change; a draw can retain the original label",
    }


def load_reference(path, kind):
    specs = {
        "native": (NATIVE_SHA, NATIVE_CELL,
                   "WCVP_native_only_geographical_gradients_and_nulls.json"),
        "global9": (GLOBAL9_SHA, GLOBAL9_CELL,
                    "ecological_species_affinity_conditional_photo9_permutation_summary.json"),
        "global8": (GLOBAL8_SHA, GLOBAL8_CELL,
                    "global_nonwhite8_species_region_whitefixed_null_summary.json"),
    }
    expected_sha, cell_name, summary_name = specs[kind]
    check(sha(path) == expected_sha, "Pinned reference ZIP SHA changed: " + kind)
    with ZipFile(path) as archive:
        receipt = json.loads(archive.read("receipt.json"))
        if kind == "native":
            signed = receipt["output_sha256"]
        elif kind == "global9":
            signed = {**receipt["outputs_SHA256"], summary_name: receipt["summary_SHA256"]}
        else:
            signed = {cell_name: receipt["source_global_8chromatic_cells_SHA256"],
                      summary_name: receipt["analysis_summary_SHA256"]}
        check(cell_name in signed and summary_name in signed,
              "Reference receipt does not sign required files: " + kind)
        for member, digest in signed.items():
            check(sha(archive.read(member)) == digest,
                  "Reference member SHA changed: " + kind + "/" + member)
        cells = pd.read_csv(io.BytesIO(archive.read(cell_name)), compression="gzip",
                            float_precision="round_trip", low_memory=False)
        summary = json.loads(archive.read(summary_name))
    check(len(cells) == 18577 and not cells[["grid_x", "grid_y"]].duplicated().any(),
          "Reference no longer contains 18,577 unique coordinates: " + kind)
    if kind == "native":
        check(summary["source_37104_accepted_photo_hue_taxa"] == 37104 and
              summary["source_full_GBIF_species_cell_pairs"] == 5493931 and
              summary["source_null_199_same_global_genus_by_width_colour_label_permutations"]
              is True, "Native reference source or permutation design changed")
    else:
        check(summary["null_permutations"] == REFERENCE_PERMUTATIONS and
              summary["random_seed" if kind == "global9" else "seed"] ==
              (SEED if kind == "global9" else SEED + 803),
              "Global reference permutation design changed: " + kind)
    return cells, summary, {"ZIP_SHA256": expected_sha, "members_SHA256": signed}


def summarize_model(cohorts, result, groups, global_source, genus_codes, model, outcome, seed):
    descriptive = model.startswith("matched_exchangeable")
    model_result = {
        "inference_status": "descriptive_only_outcome_selected_mobility" if descriptive
            else "conditional_global_species_label_permutation",
        "seed": seed, "global_entropy_evaluations_for_reference_only":
            result["n_global_entropy_evaluations"],
        "cohorts": {},
    }
    draw_frames = []
    for name, c in cohorts.items():
        observed = group_means(c, c.observed)
        statistics = {}
        for key, value in observed.items():
            stats = permutation_statistics(value, result["draws"][name][key], descriptive)
            ids = c.indices["richness_ge30"] if key == EQUAL_SECTORS else c.indices[key]
            stats.update(
                n_cells=int(len(ids)),
                median_native_richness_of_included_cells=float(np.median(c.richness[ids])),
                mean_native_richness_of_included_cells=float(np.mean(c.richness[ids])),
                weighting="equal weight for each of six sector cell means" if key == EQUAL_SECTORS
                    else "equal weight for each eligible cell",
                richness_denominator="original nonwhite native species" if c.n_hues == 8
                    else "all native source species",
            )
            expected_from_cells = group_means(c, result["expected"][name])[key]
            check(abs(stats["null_mean"] - expected_from_cells) <= REPRODUCIBILITY_ATOL,
                  "Per-cell and per-draw mean expectations disagree")
            statistics[key] = stats
        contrasts = {}
        for key, (a, b) in CONTRASTS.items():
            stats, values = paired_contrast(observed[a], observed[b],
                                            result["draws"][name][a],
                                            result["draws"][name][b], descriptive)
            stats.update(group_A=a, group_B=b, n_cells_A=len(c.indices[a]),
                         n_cells_B=len(c.indices[b]))
            check(abs(stats["residual_difference"] -
                      (statistics[a]["observed_minus_null"] -
                       statistics[b]["observed_minus_null"])) <= REPRODUCIBILITY_ATOL,
                  "Synchronous contrast residual identity failed")
            contrasts[key] = stats
            result["draws"][name][key] = values
        for key, values in result["draws"][name].items():
            draw_frames.append(pd.DataFrame({
                "outcome": outcome, "model": model, "cohort": name, "seed": seed,
                "statistic_kind": "contrast" if key in CONTRASTS else "group_mean",
                "statistic_name": key, "draw_id": np.arange(1, len(values) + 1),
                "draw_value": values,
            }))
        model_result["cohorts"][name] = {
            "population_role": "primary" if name == "strict_native"
                else "overlapping sensitivity, not independent replication",
            "exchangeability_support": exchangeability_support(
                groups, global_source, c, genus_codes),
            "group_statistics": statistics, "paired_contrasts": contrasts,
        }
    return model_result, draw_frames


def apply_regional_multiplicity(models, primary_model):
    strict = models[primary_model]["cohorts"]["strict_native"]
    members = [("group_statistics", key) for key in SECTOR_GROUPS]
    members += [("paired_contrasts", key) for key in CONTRASTS]
    check(len(members) == 8, "The strict regional multiplicity family must have eight tests")
    p = [strict[kind][key]["p_rank_two_sided"] for kind, key in members]
    corrected = adjust_pvalues(p)
    family = []
    for i, (kind, key) in enumerate(members):
        strict[kind][key]["p_rank_two_sided_holm8"] = float(corrected["holm"][i])
        strict[kind][key]["q_BH8_exploratory_no_dependence_guarantee"] = float(
            corrected["bh_exploratory"][i])
        family.append({"statistic_kind": kind, "statistic_name": key,
                       "p_rank_two_sided": p[i],
                       "p_holm8": float(corrected["holm"][i]),
                       "q_BH8_exploratory": float(corrected["bh_exploratory"][i])})
    return {
        "cohort": "strict_native", "model": primary_model, "family_size": 8,
        "main_adjustment": "Holm familywise error control under arbitrary dependence",
        "input_p": "rank two-sided, twice the smaller finite Monte Carlo tail",
        "BH_status": "exploratory diagnostic only; PRDS dependence not established",
        "members": family,
    }


def flatten_statistics(outcomes, kind):
    records = []
    for outcome, outcome_result in outcomes.items():
        for model, model_result in outcome_result["models"].items():
            for cohort, c in model_result["cohorts"].items():
                for name, stats in c[kind].items():
                    row = {"outcome": outcome, "model": model, "cohort": cohort,
                           "statistic_name": name, **stats}
                    lo, hi = row.pop("null_reference_95pct_interval")
                    row["null_reference_2p5_percentile"] = lo
                    row["null_reference_97p5_percentile"] = hi
                    records.append(row)
    return pd.DataFrame(records)


def code_provenance():
    root = Path(__file__).resolve().parents[2]
    paths = {Path(__file__).resolve()}
    for module_name, module in tuple(sys.modules.items()):
        if module_name.startswith("scripts.analysis.") and getattr(module, "__file__", None):
            path = Path(module.__file__).resolve()
            if path.suffix == ".py" and path.is_relative_to(root):
                paths.add(path)
    test = root / "tests/test_native_georegion_null_20261121.py"
    if test.exists():
        paths.add(test)
    return {str(path.relative_to(root)): sha(path) for path in sorted(paths)}


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def run(args):
    check(isinstance(args.permutations, int) and args.permutations >= REFERENCE_PERMUTATIONS,
          "At least 199 permutations are required for the original-source gates")
    out = Path(args.out)
    check(not out.exists(), "Refuse overwriting a frozen native geographic null directory")
    native_ref, native_summary, native_receipt = load_reference(args.native_reference, "native")
    global9_ref, _, global9_receipt = load_reference(args.global9_reference, "global9")
    global8_ref, _, global8_receipt = load_reference(args.global8_reference, "global8")
    wcvp = Path(args.wcvp)
    check(sha(wcvp / MANIFEST) == WCVP_RECEIPT_SHA,
          "Original WCVP receipt differs from the native-status source")
    check(sha(args.tdwg) == TDWG_GEO_SHA, "Pinned original TDWG3 polygons changed")
    pairs, species = signed_original(args.gbif, args.env)
    x = baseline(pairs, species, null_model="genus_breadth")
    coords = pd.DataFrame(np.asarray(x["cell_key"].tolist(), dtype=np.int64),
                          columns=["grid_x", "grid_y"])
    gates = []
    for name, reference in (("native", native_ref), ("global9", global9_ref),
                            ("global8", global8_ref)):
        gates.append(compare_gate(name + "_exact_18577_coordinate_order",
                                  coords.to_numpy(), reference[["grid_x", "grid_y"]].to_numpy(),
                                  atol=0))
    zone, longitude, latitude = geographic_regions(x)
    check(np.array_equal(np.asarray(ZONE_NAMES)[zone], global9_ref.macrozone.to_numpy())
          and np.array_equal(np.asarray(ZONE_NAMES)[zone], global8_ref.macrozone.to_numpy()),
          "Original coordinate sector assignments changed")
    gates.append({"name": "six_coordinate_sectors_and_27_undefined_polar_latitudes",
                  "passed": True, "n_undefined_latitudes": int((~np.isfinite(latitude)).sum())})
    base, occupied, genus_codes, _ = genus_breadth_codes(x, species)
    affinity, share, _ = species_affinity(x, zone)
    groups9 = {model: nulls(x, species, model, zone, affinity, base, occupied,
                           genus_codes, share) for model in MODELS9}
    groups8, genus_codes8, affinity8 = calculate_groups(x, species, zone)
    check(tuple(groups8) == MODELS8 and np.array_equal(genus_codes8, genus_codes)
          and np.array_equal(affinity8, affinity), "Original white-fixed model definitions changed")
    for models, matched, regional in ((groups9, MODELS9[1], MODELS9[2]),
                                      (groups8, MODELS8[1], MODELS8[2])):
        check(np.array_equal(np.sort(np.concatenate(models[matched])),
                             np.sort(np.concatenate(models[regional]))),
              "Matched and regional model movable species identities differ")
    global_sources = {
        n: prepare_assemblage("global", x["cell"], x["species"], x["hue"], x["ngrids"], n)
        for n in (9, 8)
    }
    gates.extend([
        compare_gate("full_global9_observed_entropy", global_sources[9].observed,
                     global9_ref.observed_hue9_Shannon.to_numpy(float)),
        compare_gate("full_global9_richness", global_sources[9].richness,
                     global9_ref.n_species.to_numpy(int), atol=0),
        compare_gate("full_global8_observed_entropy", global_sources[8].observed,
                     global8_ref.observed_nonwhite8_Shannon.to_numpy(float)),
        compare_gate("full_global8_nonwhite_richness", global_sources[8].richness,
                     global8_ref.n_nonwhite_source_photo9_species.to_numpy(int), atol=0),
    ])
    # This source loader only consumes grid coordinates from its historical
    # `lgm` argument. Reuse the exact source coordinates without adding any
    # palaeoclimate dependency or replacing the original 25-point assignment.
    all_cohorts, native_audit = source_guild_rows(
        pairs, species, coords, wcvp, args.tdwg, coords)
    check(native_audit["n_WCVP_accepted_unique_exact_matched_species"] == 36762
          and all_cohorts["strict_native"]["n_source_species_grid_pairs"] == 2271852
          and all_cohorts["moderate_native"]["n_source_species_grid_pairs"] == 2604870,
          "Original exact WCVP mapping or strict/moderate native pair counts changed")
    check(native_audit == native_summary["source_WCVP_exact_match_and_TDWG3_status"],
          "Reconstructed source WCVP lineage differs from the native reference")
    cohort_sources = {name: all_cohorts[name] for name in COHORTS}
    prepared = {}
    for n_hues in (9, 8):
        prepared[n_hues] = {}
        for name, source in cohort_sources.items():
            c = prepare_assemblage(name, source["cell_index"], source["species_index"],
                                    x["hue"], x["ngrids"], n_hues)
            c.indices = summary_indices(c.richness, zone, latitude)
            check(all(len(ii) > 0 for ii in c.indices.values()),
                  f"No eligible native cells in a required {n_hues}-hue geographic group")
            prepared[n_hues][name] = c
            if n_hues == 9:
                gates.extend([
                    compare_gate(name + "_native9_observed_entropy", c.observed,
                                 native_ref[name + "_Shannon_observed"].to_numpy(float)),
                    compare_gate(name + "_native9_richness", c.richness,
                                 native_ref[name + "_n_photo9_species"].to_numpy(int), atol=0),
                ])
    for name, counts in {"strict_native": (7231, 5389, 5799, 3051),
                         "moderate_native": (8155, 6076, 6527, 3440)}.items():
        observed_counts = [len(prepared[n][name].indices[k]) for n in (9, 8)
                           for k in ("richness_ge30", "richness_ge100")]
        gates.append(compare_gate(name + "_richness_threshold_denominators",
                                  observed_counts, counts, atol=0))
    gates.append(compare_gate("strict_native9_primary_observed_mean",
                              [group_means(prepared[9]["strict_native"],
                                           prepared[9]["strict_native"].observed)["richness_ge30"]],
                              [1.2152685262485332]))
    gates.append({"name": "source_WCVP_accepted_and_native_pair_lineage", "passed": True,
                  "accepted_exact_matches": 36762, "strict_native_pairs": 2271852,
                  "moderate_native_pairs": 2604870})
    print("MOIS_NATIVE_GEOGRAPHIC_SOURCES_VERIFIED", json.dumps({
        "source_audit": native_audit,
        "richness_denominators": {str(n): {
            name: {key: len(c.indices[key]) for key in ("richness_ge30", "richness_ge100")}
            for name, c in prepared[n].items()} for n in (9, 8)},
    }, sort_keys=True), flush=True)
    code_sha = code_provenance()
    outcomes = {}
    tables = {}
    draw_frames = []
    full_validation = coords.copy()
    for n_hues, outcome, groups, models, reference, prefix in (
        (9, "photo9", groups9, MODELS9, global9_ref, "expected_Shannon_"),
        (8, "whitefixed8", groups8, MODELS8, global8_ref, "null_expected_nonwhite8_H_"),
    ):
        cohorts = prepared[n_hues]
        seed = SEED if n_hues == 9 else SEED + 803
        cells = coords.copy()
        cells["center_longitude"] = longitude
        cells["center_latitude"] = latitude
        cells["coordinate_sector"] = np.asarray(ZONE_NAMES)[zone]
        for name, c in cohorts.items():
            cells[name + "_richness"] = c.richness
            cells[name + "_observed_Shannon"] = c.observed
            hue_names = HUES if n_hues == 9 else HUES[1:]
            for j, hue_name in enumerate(hue_names):
                cells[name + "_observed_" + hue_name + "_species"] = c.counts[:, j]
            cells[name + "_richness_ge30"] = c.richness >= 30
            cells[name + "_richness_ge100"] = c.richness >= 100
        model_results = {}
        for model in models:
            prior_native = {name: native_ref[name + "_Shannon_null"].to_numpy(float)
                            for name in COHORTS} if n_hues == 9 and model == "genus_breadth" else None
            simulation = simulate_model(
                global_sources[n_hues], cohorts, groups[model], args.permutations, seed, model,
                reference_global_expected=reference[prefix + model].to_numpy(float),
                reference_native_expected=prior_native)
            gates.extend(simulation["gates"])
            model_results[model], frames = summarize_model(
                cohorts, simulation, groups[model], global_sources[n_hues], genus_codes,
                model, outcome, seed)
            draw_frames.extend(frames)
            for name, c in cohorts.items():
                cells[name + "_expected_Shannon_" + model] = simulation["expected"][name]
                cells[name + "_observed_minus_expected_" + model] = (
                    c.observed - simulation["expected"][name])
                cells[name + "_expected_first199_" + model] = simulation["expected_first199"][name]
            full_validation[outcome + "_expected_first199_" + model] = simulation["global_expected_first199"]
            print("MOIS_NATIVE_GEOGRAPHIC_MODEL_COMPLETED", outcome, model,
                  json.dumps({name: model_results[model]["cohorts"][name]["group_statistics"]
                              ["richness_ge30"] for name in COHORTS}, sort_keys=True), flush=True)
        outcomes[outcome] = {
            "n_hues": n_hues, "seed": seed, "permutations": args.permutations,
            "primary_regional_model": models[2], "models": model_results,
            "strict_regional_eight_test_family": apply_regional_multiplicity(model_results, models[2]),
        }
        tables[outcome] = cells
    check(code_sha == code_provenance(), "Analysis code changed during the run")
    source_paths = {
        "GBIF_zip": Path(args.gbif), "ENV_zip": Path(args.env),
        "native_reference_zip": Path(args.native_reference),
        "global9_reference_zip": Path(args.global9_reference),
        "global8_reference_zip": Path(args.global8_reference),
        "WCVP_receipt": wcvp / MANIFEST, "WCVP_exact_names": wcvp / NAME_FILE,
        "WCVP_native_distribution": wcvp / DIST_FILE, "TDWG3_geojson": Path(args.tdwg),
    }
    sources_sha = {name: {"filename": path.name, "SHA256": sha(path)}
                   for name, path in source_paths.items()}
    runtime = {"python": platform.python_version(), "platform": platform.platform()}
    for distribution in ("numpy", "pandas", "pyarrow", "geopandas", "shapely", "pyproj"):
        runtime[distribution] = metadata.version(distribution)
    summary = {
        "schema": SCHEMA,
        "design_status": "analysis choices recorded for this run; no formal preregistration claim",
        "source_species_unit": "original accepted GBIF specieskey; WCVP crosswalk supplies status only",
        "source_species": len(species), "source_species_cell_pairs": len(pairs),
        "source_cells": x["ngrids"], "permutations_per_model": args.permutations,
        "reference_prefix_permutations": REFERENCE_PERMUTATIONS,
        "primary_population": "strict_native", "sensitivity_population": "moderate_native",
        "populations_overlap_and_are_not_independent_replication": True,
        "primary_statistic": "photo9 strict_native richness>=30 cell-equal mean Shannon",
        "primary_regional_model": MODELS9[2],
        "primary_direction": "lower-tail concentration; rank two-sided also reported",
        "donor_pool": "original full-present-GBIF occupied species; never restricted to native cohorts",
        "range_breadth_and_affinity_source": "full present GBIF footprints, including introduced occurrences",
        "absent_GBIF_species_excluded_from_donor_groups": 346,
        "occupancy_breadth_bins": ["1", "2-5", "6-20", "21-100", ">=101"],
        "species_affinity_threshold": THRESH,
        "species_affinity_rule": ">=65% of full-present occupied cells in one sector; otherwise multizone",
        "sectors_are_coordinate_partitions_not_biological_realms": True,
        "sector_definition": "longitude <-30, -30 to <60, >=60 crossed with signed-northing hemisphere",
        "sector_names": list(ZONE_NAMES),
        "polar_handling": "27 source edge-cell latitudes remain missing; hemisphere uses signed northing; missing latitudes excluded from Arctic/boreal cuts",
        "polar_edge_cells_with_undefined_latitude": int((~np.isfinite(latitude)).sum()),
        "whitefixed8": "white IDs and original nonwhite species-cell memberships fixed in every draw",
        "photo9": "original nine-label shuffle, including white; native cohort memberships fixed",
        "all_cells_and_both_native_cohorts_share_one_species_label_per_model_draw": True,
        "shared_seed_across_different_models_does_not_define_paired_model_contrasts": True,
        "matched_mobility": "eligible IDs selected from observed hue heterogeneity; descriptive only, p-values suppressed; equal IDs do not equate palettes or null variance",
        "rank_p_formulas": {
            "lower": "(1 + count(draw <= observed + 1e-12)) / (B + 1)",
            "upper": "(1 + count(draw >= observed - 1e-12)) / (B + 1)",
            "two_sided": "min(1, 2 * min(lower, upper))",
            "legacy": "simulated-mean-centered tail is a separately named diagnostic only",
        },
        "uncertainty": "null mean Monte Carlo SE = sample SD / sqrt(B); 95% null reference interval is not an effect CI",
        "regional_contrasts": "synchronous within-model draw differences; never differences between p-values",
        "multiplicity": "strict regional six sectors plus two contrasts, Holm8 separately for photo9 and whitefixed8; BH8 exploratory only",
        "source_photo_hues_are_not_independently_verified_petal_colours": True,
        "native_status_is_broad_TDWG3_not_individual_wildness_or_endemicity": True,
        "causal_selection_pollination_and_phylogenetic_inference_not_identified": True,
        "native_source_audit": native_audit,
        "source_files": sources_sha,
        "reference_integrity": {"native": native_receipt, "global9": global9_receipt,
                                "global8": global8_receipt},
        "runtime": runtime, "code_SHA256": code_sha, "outcomes": outcomes,
    }
    # Outputs are written only after all source gates and models succeed.
    out.mkdir(parents=True, exist_ok=False)
    gzip = {"method": "gzip", "mtime": 0}
    tables["photo9"].to_csv(out / "native_photo9_georegion_cells.csv.gz", index=False, compression=gzip)
    tables["whitefixed8"].to_csv(out / "native_whitefixed8_georegion_cells.csv.gz", index=False, compression=gzip)
    full_validation.to_csv(out / "full_global_first199_validation_cells.csv.gz", index=False, compression=gzip)
    pd.concat(draw_frames, ignore_index=True).to_csv(
        out / "native_georegion_permutation_draws.csv.gz", index=False, compression=gzip)
    flatten_statistics(outcomes, "group_statistics").to_csv(
        out / "native_georegion_group_statistics.csv", index=False)
    flatten_statistics(outcomes, "paired_contrasts").to_csv(
        out / "native_georegion_contrast_statistics.csv", index=False)
    write_json(out / "native_georegion_summary.json", summary)
    write_json(out / "reproducibility_gates.json", {
        "schema": SCHEMA, "all_passed": all(g["passed"] for g in gates),
        "first199_reference_gate_count": 8,
        "reference_prefix_permutations": REFERENCE_PERMUTATIONS, "gates": gates,
    })
    receipt = {
        "schema": SCHEMA, "status": "complete", "permutations_per_model": args.permutations,
        "source_files": sources_sha, "code_SHA256": code_sha, "runtime": runtime,
        "output_SHA256": {path.name: sha(path) for path in sorted(out.iterdir()) if path.is_file()},
        "command_arguments": vars(args),
        "draw_id_interpretation": "1-based; paired only within outcome/model/seed across groups and cohorts",
        "all_reproducibility_gates_passed": True,
    }
    write_json(out / "receipt.json", receipt)
    print("MOIS_NATIVE_HUE9_WHITEFIXED8_GEOREGION_DONE", json.dumps({
        "output": str(out), "permutations_per_model": args.permutations,
        "source_gates_passed": len(gates),
        "preserved_individual_statistics": sum(len(frame) for frame in draw_frames),
    }, sort_keys=True), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("gbif", "env", "native-reference", "global9-reference",
                 "global8-reference", "wcvp", "tdwg", "out"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--permutations", type=int, default=999)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
