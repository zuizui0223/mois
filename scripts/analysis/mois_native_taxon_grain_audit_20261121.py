#!/usr/bin/env python3
"""Read-only source-taxonomy audit; writes only its JSON report.

Run from the repository:
  python -m scripts.analysis.mois_native_taxon_grain_audit_20261121 \\\n    --source source --out results/wcvp_taxon_grain_audit.json
No species identities, source files, or permutation results are modified.
"""
from pathlib import Path
import argparse
import hashlib
import itertools
import json
import numpy as np
import pandas as pd

from scripts.analysis.mois_global_hue9_community_assembly_genus_null_20261118 import (
    HUES, empirical_entropy, signed_original,
)
from scripts.analysis.mois_wcvp_native_hue9_geographical_null_20261118 import (
    get_wcvp_source, match_gbif_to_wcvp, source_guild_rows,
)

SOURCE_FILES = (
    "GBIF.zip", "ENV.zip", "TDWG3.geojson",
    "WCVP/mois350399_wcvp_TDWG3_source_only_taxon_native_receipt.json",
    "WCVP/mois350399_wcvp_source_taxon_exact_name_resolution.csv.gz",
    "WCVP/mois350399_wcvp_accepted_native_distribution_TDWG3.csv.gz",
)
INTERPRETATION = {
    "primary_unit": "Frozen original COLXR/GBIF source species key, once per grid cell",
    "wcvp_role": "Exact-name crosswalk supplies regional native status; it does not redefine the source taxonomic unit",
    "many_to_one_crosswalk_is_not_duplicate_source_rows": True,
    "frozen_source_species_keys_unchanged": True,
    "conflicting_hues_make_automatic_taxon_collapse_ambiguous": True,
    "native_status_membership_reconstructed_with_original_25_subpoint_rule": True,
    "no_permutations_run_by_this_audit": True,
    "observed_entropy_bound_scope": (
        "Fixed original >=30 source-species cell mask; retain any one observed hue "
        "per co-occurring WCVP accepted ID. Per-cell choices may vary, giving a "
        "conservative bound on observed mean Shannon only. This is not a new "
        "taxonomic analysis and does not bound changes to null means, p values, "
        "or a changed eligibility mask."
    ),
    "current_action": "Preserve frozen primary analysis and disclose the source-unit/status-crosswalk distinction",
}


def native_summary(c, wid, hue, collision_ids, ngrids):
    ci, si = c["cell_index"], c["species_index"]
    richness = np.bincount(ci, minlength=ngrids)
    eligible = richness >= 30
    rows = pd.DataFrame({"cell": ci, "source": si, "wid": wid[si], "hue": hue[si]})
    collision = rows.loc[rows.wid.isin(collision_ids)]
    grouped = collision.groupby(["cell", "wid"]).agg(
        n_source=("source", "nunique"), n_hue=("hue", "nunique")
    )
    duplicated = grouped.loc[grouped.n_source > 1]
    surplus = duplicated.n_source - 1
    removed = surplus.groupby("cell").sum().reindex(np.arange(ngrids), fill_value=0).to_numpy()
    duplicate30 = duplicated.loc[
        duplicated.index.get_level_values("cell").isin(np.flatnonzero(eligible))
    ]
    counts = np.bincount(ci.astype(np.int64) * 9 + hue[si], minlength=ngrids * 9).reshape(-1, 9)
    old_h = empirical_entropy(counts)
    min_h, max_h = old_h.copy(), old_h.copy()
    options = {}
    for (cell, _), group in collision.groupby(["cell", "wid"]):
        if len(group) < 2:
            continue
        hist = np.bincount(group.hue.to_numpy(), minlength=9)
        choices = []
        for retained_hue in np.flatnonzero(hist):
            removed_hues = hist.copy()
            removed_hues[retained_hue] -= 1
            choices.append(removed_hues)
        options.setdefault(int(cell), []).append(choices)
    max_choices = 0
    for cell, groups in options.items():
        combinations = list(itertools.product(*groups))
        max_choices = max(max_choices, len(combinations))
        retained = np.array([counts[cell] - np.sum(choice, axis=0) for choice in combinations])
        entropy = empirical_entropy(retained)
        min_h[cell], max_h[cell] = entropy.min(), entropy.max()
    return {
        "source_pairs": len(ci), "source_species": len(np.unique(si)),
        "cooccurring_duplicate_taxon_cells": len(duplicated),
        "excess_pairs": int(surplus.sum()),
        "excess_pair_fraction": float(surplus.sum() / len(ci)),
        "affected_cells_all": int(duplicated.index.get_level_values("cell").nunique()),
        "affected_cells_ge30": int(duplicate30.index.get_level_values("cell").nunique()),
        "conflicting_hue_duplicate_groups_all": int((duplicated.n_hue > 1).sum()),
        "conflicting_hue_duplicate_groups_ge30": int((duplicate30.n_hue > 1).sum()),
        "n_ge30_source": int(eligible.sum()),
        "n_ge30_after_acceptedspecies_dedup": int(((richness - removed) >= 30).sum()),
        "same_source_mask_observed_mean_H": float(old_h[eligible].mean()),
        "same_source_mask_any_per_cell_dedup_observed_mean_H_min": float(min_h[eligible].mean()),
        "same_source_mask_any_per_cell_dedup_observed_mean_H_max": float(max_h[eligible].mean()),
        "same_source_mask_max_absolute_change_in_observed_mean_H": float(max(
            abs((min_h - old_h)[eligible].mean()), abs((max_h - old_h)[eligible].mean())
        )),
        "maximum_hue_retention_combinations_per_cell": max_choices,
    }


def main(args):
    root = Path(args.source)
    output = Path(args.out)
    if output.exists():
        raise FileExistsError("Refusing to overwrite source taxonomy audit")
    pairs, species = signed_original(root / "GBIF.zip", root / "ENV.zip")
    names, _, _ = get_wcvp_source(root / "WCVP")
    crosswalk = match_gbif_to_wcvp(species, names)
    matched = crosswalk.loc[crosswalk.wcvp_match_exact_binomial]
    grouped = matched.groupby("wcvp_accepted_plant_name_id").agg(
        n_source=("key", "nunique"), n_hue=("source_hue9", "nunique")
    )
    collision_ids = set(grouped.index[grouped.n_source > 1])
    collisions = matched.loc[matched.wcvp_accepted_plant_name_id.isin(collision_ids)]
    catalog = {
        "catalog_species": len(species), "mapped_species": len(matched),
        "unique_wcvp_species": len(grouped),
        "wcvp_ids_receiving_multiple_source_ids": len(collision_ids),
        "source_ids_in_many_to_one_wcvp_mappings": len(collisions),
        "many_to_one_ids_with_different_hues": int(((grouped.n_source > 1) & (grouped.n_hue > 1)).sum()),
    }
    subset = pairs.merge(
        collisions[["key", "wcvp_accepted_plant_name_id", "source_hue9"]],
        left_on="specieskey", right_on="key", how="inner", validate="many_to_one",
    )
    grouped_pairs = subset.groupby(["grid_x", "grid_y", "wcvp_accepted_plant_name_id"]).agg(
        n_source=("specieskey", "nunique"), n_hue=("source_hue9", "nunique")
    )
    duplicate_pairs = grouped_pairs.loc[grouped_pairs.n_source > 1]
    global_counts = {
        "all_source_pairs": len(pairs), "collision_source_pairs": len(subset),
        "duplicate_accepted_taxon_grid_groups": len(duplicate_pairs),
        "source_pairs_in_duplicate_groups": int(duplicate_pairs.n_source.sum()),
        "excess_source_pairs_vs_wcvp_unit": int((duplicate_pairs.n_source - 1).sum()),
        "affected_cells": len(duplicate_pairs.reset_index()[["grid_x", "grid_y"]].drop_duplicates()),
        "conflicting_hue_groups": int((duplicate_pairs.n_hue > 1).sum()),
        "maximum_duplicate_ids_per_taxon_cell": int(duplicate_pairs.n_source.max()) if len(duplicate_pairs) else 0,
    }
    coords = pairs[["grid_x", "grid_y"]].drop_duplicates().sort_values(["grid_x", "grid_y"]).reset_index(drop=True)
    cohorts, _ = source_guild_rows(pairs, species, coords, root / "WCVP", root / "TDWG3.geojson", coords)
    wid = crosswalk.wcvp_accepted_plant_name_id.fillna("").to_numpy(str)
    hue = np.array([HUES.index(k) for k in species.source_hue9], dtype=int)
    result = {
        "schema": "mois_wcvp_source_taxon_grain_audit_20261121_v1",
        "source_sha256": {name: hashlib.sha256((root / name).read_bytes()).hexdigest() for name in SOURCE_FILES},
        "catalog_crosswalk": catalog,
        "global_cooccurrence": global_counts,
        "native_cohorts": {name: native_summary(cohorts[name], wid, hue, collision_ids, len(coords))
                           for name in ("strict_native", "moderate_native")},
        "interpretation": INTERPRETATION,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True)
    parser.add_argument("--out", required=True)
    main(parser.parse_args())

