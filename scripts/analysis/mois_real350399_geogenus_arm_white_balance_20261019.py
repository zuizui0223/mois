#!/usr/bin/env python3
"""350399 REAL photographic ROI outcomes: does source-arm whiteness survive spatial/taxonomic matching?

Compares the original deeply photographed core with independently recruited,
mostly one-photo-per-taxon broadened sources. These recruitment arms have
UNKNOWN cross-species inclusion probabilities; balancing named genera and
geography changes a DESCRIPTIVE estimand, not a population bias correction.

Every original pixel outcome is re-used from SHA-verified Releases. Nothing
labels unmeasured metadata or technically unknown ROIs as coloured flowers.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis.mois_actual_photo_cohort_season_audit_20261019 import (
    HUES, assemble, require, sha256, summary,
)

OUT_JSON = "mois350399_real_photo_arm_geogenus_balance.json"
OUT_CSV = "mois350399_matched_geogenus_photo_hue_contrasts.csv"
OUT_GEO = "mois350399_geographic_only_matched_taxon_photo_hues.csv"
OUT_ARM = "mois350399_photo_arm_hue_yield.csv"
N_PHOTOS = 350399
N_CORE = 317631
N_NEW = 32768
TILES = (5, 10, 20)
GATES = (1, 2)
GENUS_PATTERN = re.compile(r"^[A-Z][A-Za-z-]+$")


def arm_stats(d):
    s = summary(d)
    return {k: v for k, v in s.items() if k != "original_unclassified_morph_label_counts"}


def verified(d):
    """Actual photo classifications only; old FCP and new one-taxon-arm stay separate."""
    require(len(d) == N_PHOTOS and d.photo_id.is_unique, "Not full 350399 original photo pixels")
    require(int(d.arm.eq("core").sum()) == N_CORE and
            int(d.arm.ne("core").sum()) == N_NEW, "Unexpected old/new source cohort sizes")
    x = d[d.classified].copy()
    x["cohort"] = np.where(x.arm.eq("core"), "core", "new_species_first")
    x["genus"] = x.species.str.strip().str.split().str[0]
    x["valid_genus"] = x.genus.map(lambda z: isinstance(z, str) and bool(GENUS_PATTERN.fullmatch(z)))
    x["valid_geography"] = x.latitude.between(-90, 90) & x.longitude.between(-180, 180)
    x["valid_taxon"] = x.inat_taxon_id.ne("")
    return x


def group_hue(d, by):
    """One equal-weight named taxon per stratum; photos within taxon not independent taxa."""
    x = d.copy()
    for hue in HUES:
        x["p_" + hue] = x.morph.eq(hue).astype(float)
    cols = ["p_" + h for h in HUES]
    taxon = x.groupby([*by, "inat_taxon_id"], observed=True, sort=False).agg(
        **{h: (h, "mean") for h in cols},
        photo_n=("photo_id", "size")
    ).reset_index()
    aggs = {h: (h, "mean") for h in cols}
    return taxon.groupby(by, observed=True, sort=False).agg(
        **aggs, taxa=("inat_taxon_id", "size"),
        photos=("photo_n", "sum")
    ).reset_index(), len(taxon)


def matched_contrasts(x, *, tile, genus, min_taxa, draws=999, seed=20261019):
    """Compute paired geographic or genus×geographic source strata.

    The strata are paired by NAME/location, *not* by an identical species
    or biological individual, hence any residual gap is not phenotype change.
    """
    require(tile in TILES and min_taxa in GATES and draws >= 1, "Invalid prespecified matching grid")
    z = x.loc[x.valid_geography & x.valid_taxon &
              (x.valid_genus if genus else True)].copy()
    z["ilat"] = np.floor((z.latitude + 90) / tile).astype(int)
    z["ilon"] = np.floor((z.longitude + 180) / tile).astype(int)
    keys = (["genus"] if genus else []) + ["ilat", "ilon"]
    # Each original taxon gets equal weight within its recruitment arm and cell;
    # the groups share genus/geography only when genus is requested.
    cells, ntaxon_cells = group_hue(z, ["cohort", *keys])
    parent = cells.loc[cells.cohort.eq("core")].set_index(keys)
    newer = cells.loc[cells.cohort.eq("new_species_first")].set_index(keys)
    m = parent.join(newer, how="inner", lsuffix="_core", rsuffix="_new")
    m = m.loc[m.taxa_core.ge(min_taxa) & m.taxa_new.ge(min_taxa)].copy()
    if len(m):
        m["white_diff_pp"] = 100 * (m.p_white_new - m.p_white_core)
    n_geocells = int(z.groupby(["ilat", "ilon"], observed=True).ngroups)
    taxa_total = {str(k): int(v) for k, v in
                  cells.groupby("cohort", observed=True).taxa.sum().items()}
    n = len(m)
    selected = {
        "core": int(m.taxa_core.sum()) if n else 0,
        "new_species_first": int(m.taxa_new.sum()) if n else 0,
    }
    support = dict(
        matching=("nominal_genus_and_geographic_cell" if genus else "geographic_cell_only"),
        tile_degrees=tile, minimum_distinct_named_taxa_per_arm_per_cell=min_taxa,
        matched_strata=n, matched_nominal_genera=int(m.reset_index().genus.nunique())
        if n and genus else None,
        matched_geographic_cells=int(m.reset_index()[["ilat", "ilon"]].drop_duplicates().shape[0])
        if n else 0,
        photographed_geographic_cells=n_geocells,
        source_taxon_cells_all_cohorts=ntaxon_cells,
        matched_strata_named_taxon_cells_by_arm=selected,
        total_classifiable_taxon_cells_by_arm=taxa_total,
        retained_taxon_cell_share_by_arm={
            k: (selected[k] / taxa_total[k] if taxa_total.get(k, 0) else None)
            for k in ("core", "new_species_first")
        },
        matched_original_photo_ROIs_by_arm={
            "core": int(m.photos_core.sum()) if n else 0,
            "new_species_first": int(m.photos_new.sum()) if n else 0,
        },
        exact_same_species_required=False,
        true_inclusion_probabilities_known=False,
    )
    if not n:
        return {**support, "matched_hue_differences_pp": None,
                "genus_equal_bootstrap_white_95_pp": None,
                "core_taxon_cell_weighted_white_pp": None,
                "new_taxon_cell_weighted_white_pp": None}
    hue_differences = {
        hue: float(100 * (m["p_" + hue + "_new"] - m["p_" + hue + "_core"]).mean())
        for hue in HUES
    }
    require(abs(sum(hue_differences.values())) < 1e-7,
            "Actual photographed original four-hue stratum composition does not sum to 1")
    if genus:
        g = m.reset_index().groupby("genus", observed=True).white_diff_pp.mean()
        values = g.to_numpy(float)
        # Whole nominal genera, not photo counts; exploratory, without phylogenetic independence.
        rng = np.random.default_rng(seed + 100 * tile + min_taxa)
        b = np.array([values[rng.integers(0, len(values), len(values))].mean()
                      for _ in range(draws)])
        genus_mean = float(values.mean())
        bounds = [float(a) for a in np.quantile(b, [0.025, 0.975])]
    else:
        genus_mean, bounds = None, None
    return {
        **support,
        "matched_hue_differences_pp": hue_differences,
        "genus_equal_white_difference_pp": genus_mean,
        "genus_equal_bootstrap_white_95_pp": bounds,
        "core_taxon_cell_weighted_white_pp": float(np.average(
            m.white_diff_pp, weights=m.taxa_core)),
        "new_taxon_cell_weighted_white_pp": float(np.average(
            m.white_diff_pp, weights=m.taxa_new)),
        "matched_stratum_equal_white_difference_pp": float(m.white_diff_pp.mean()),
    }


def source_composition(x, d):
    core = x.loc[x.cohort.eq("core")]
    new = x.loc[x.cohort.eq("new_species_first")]
    raw_gap = 100 * (float(new.morph.eq("white").mean()) -
                     float(core.morph.eq("white").mean()))
    # Unique nominal taxa are not necessarily taxonomically validated species.
    old_species = set(core.inat_taxon_id)
    new_species = set(new.inat_taxon_id)
    old_genus = set(core.loc[core.valid_genus, "genus"])
    new_genus = set(new.loc[new.valid_genus, "genus"])
    return {
        "source_photo_attempts": len(d),
        "source_classifiable_four_hue_photos": len(x),
        "source_arm_photo_attempts": {
            "core": int(d.arm.eq("core").sum()),
            "new_species_first": int(d.arm.ne("core").sum())
        },
        "source_arm_classifiable_photo_ROIs": {
            "core": len(core),
            "new_species_first": len(new)
        },
        "raw_new_minus_core_white_percentage_points": raw_gap,
        "raw_core_white_fraction": float(core.morph.eq("white").mean()),
        "raw_new_white_fraction": float(new.morph.eq("white").mean()),
        "overlapping_classifiable_named_taxon_ids": len(old_species & new_species),
        "overlapping_classifiable_nominal_genera": len(old_genus & new_genus),
        "classifiable_new_arm_taxa_not_in_core": len(new_species - old_species),
        "classifiable_new_arm_genera_not_in_core": len(new_genus - old_genus),
        "valid_genus_and_coordinate_photo_count": {
            group: int((s.valid_genus & s.valid_geography & s.valid_taxon).sum())
            for group, s in [("core", core), ("new_species_first", new)]
        },
    }


def run(core, broad, wave2, cached8192, tail16384, output):
    require(not output.exists() or not any(output.iterdir()),
            "Refusing to overwrite prior finished photographic-source audit")
    d, arms = assemble(core, broad, wave2, cached8192, tail16384, strict=True)
    x = verified(d)
    report = {
        "schema": "mois_actual350399_geography_genus_matched_source_arm_four_hue_v1",
        "immutable_original_release_photo_sha256": {
            label: sha256(f) for label, f in [
                ("primary317631", core), ("first_broad_cumulative", broad),
                ("second4096", wave2), ("cached8192", cached8192),
                ("finished_tail16384", tail16384)]
        },
        "source": source_composition(x, d),
        "all_photos": arm_stats(d),
        "by_source_arm": {
            name: arm_stats(frame) for name, frame in zip(
                ("core", "broad4096", "wave24096", "cached8192", "tail16384"), arms)
        },
        "matched_geogenus": [], "matched_geography_only": [],
        "limitations": [
            "All outcomes are original automatically inferred image ROI hue classes, not independently verified open-flower petals or calibrated anthocyanin/UV colour.",
            "Historical deeply sampled and new deliberately genus/species-first arms have unknown global sampling inclusion probabilities; these are descriptive source comparisons.",
            "Matched named genus is not matched phylogenetic lineage, species, site population, or paired biological sample.",
            "Equal taxon in genus x geographic cells does not remove observation/photo selection, species turnover, invasive/non-native taxa, season or unknown hue missingness.",
            "Genus bootstrap percentiles are exploratory and cannot establish independent phylogenetic lineages or multiple-test-adjusted effect.",
            "Unmatched geographic/genus cells are outside the matched-support estimand; contrasting only shared support is not a worldwide population white estimate."
        ],
        "no_source_photo_IDs_without_measured_pixels_assigned_hues": True,
        "all_source_original_photo_IDs_unique": True,
    }
    for t in TILES:
        for gate in GATES:
            report["matched_geogenus"].append(matched_contrasts(
                x, tile=t, genus=True, min_taxa=gate))
            report["matched_geography_only"].append(matched_contrasts(
                x, tile=t, genus=False, min_taxa=gate))
    output.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(report["matched_geogenus"]).to_csv(output/OUT_CSV, index=False)
    pd.DataFrame(report["matched_geography_only"]).to_csv(output/OUT_GEO, index=False)
    pd.DataFrame([
        {"arm": k, **{kk: vv for kk, vv in value.items()
                      if not isinstance(vv, dict)}}
        for k, value in report["by_source_arm"].items()
    ]).to_csv(output/OUT_ARM, index=False)
    report["csv_sha256"] = {f: sha256(output/f) for f in
                            (OUT_CSV, OUT_GEO, OUT_ARM)}
    (output/OUT_JSON).write_text(json.dumps(
        report, ensure_ascii=False, indent=2, sort_keys=True,
        allow_nan=False) + "\n", encoding="utf-8")
    print("MOIS_REAL350399_MATCHED_GENUS_GEOGRAPHY_SOURCE_BIAS", json.dumps({
        "raw_new_minus_core_white_pp": report["source"]["raw_new_minus_core_white_percentage_points"],
        "shared_source_named_taxa": report["source"]["overlapping_classifiable_named_taxon_ids"],
        "shared_nominal_genera": report["source"]["overlapping_classifiable_nominal_genera"],
        "geogenus": [
            {"tile": r["tile_degrees"], "min_taxa": r["minimum_distinct_named_taxa_per_arm_per_cell"],
             "matched_strata": r["matched_strata"],
             "matched_genera": r["matched_nominal_genera"],
             "equal_genus_white_difference_pp": r.get("genus_equal_white_difference_pp"),
             "bootstrap_95": r["genus_equal_bootstrap_white_95_pp"],
             "core_coverage": r["retained_taxon_cell_share_by_arm"]["core"],
             "new_coverage": r["retained_taxon_cell_share_by_arm"]["new_species_first"]}
            for r in report["matched_geogenus"]
        ],
    }, ensure_ascii=False, sort_keys=True), flush=True)
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("core", "broad", "wave2", "cached8192", "tail16384", "output"):
        p.add_argument("--" + name, type=Path, required=True)
    a = p.parse_args()
    run(a.core, a.broad, a.wave2, a.cached8192, a.tail16384, a.output)


if __name__ == "__main__":
    main()
