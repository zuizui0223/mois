#!/usr/bin/env python3
"""MOIS: what part of the summer photographic-white contrast is taxon mix?

Use ONLY authenticated, genuinely classified original photo pixels, not photo
opportunities. A symmetric Kitagawa partition is exact for two seasons among
identical taxa but is observational photo-weighted, not an ecological causal
attribution or an effect within a true population.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis.mois_seasonal_267k_climate_normals import (
    HUES, MOIS_HUES, MOIS_ROWS, PHOTO_DATA, PHOTO_RECEIPT,
    labelled_calendar, read_signed, require,
)


def symmetric_photo_decomposition(counts: np.ndarray) -> np.ndarray:
    """[within-taxon, source-taxon-composition, pooled delta, taxon-equal delta].

    Input columns: summer number, winter number, summer white, winter white;
    only taxa observed in BOTH seasons are admitted. All outputs are p.p.
    """
    a = np.asarray(counts, dtype=float)
    require(a.ndim == 2 and a.shape[1] == 4 and len(a) > 0,
            "Taxon × season 4-column counts are required")
    summer_n, winter_n, summer_w, winter_w = a.T
    require(np.all(summer_n > 0) and np.all(winter_n > 0)
            and np.all(summer_w >= 0) and np.all(winter_w >= 0)
            and np.all(summer_w <= summer_n) and np.all(winter_w <= winter_n),
            "Invalid photographic numerator/denominator or one-season-only species")
    ps, pw = summer_w / summer_n, winter_w / winter_n
    ws, ww = summer_n / summer_n.sum(), winter_n / winter_n.sum()
    within = np.sum((ps - pw) * (ws + ww) / 2)
    composition = np.sum((ws - ww) * (ps + pw) / 2)
    raw = np.sum(ps * ws) - np.sum(pw * ww)
    require(np.isclose(within + composition, raw, rtol=0, atol=1e-12),
            "Symmetric taxon-mixture decomposition is not additive")
    equal = np.mean(ps - pw)
    return np.array([within, composition, raw, equal], dtype=float) * 100


def bootstrap_taxa(a: np.ndarray, *, n_draws: int = 999,
                   seed: int = 20261017) -> tuple[np.ndarray, np.ndarray]:
    point = symmetric_photo_decomposition(a)
    rng = np.random.default_rng(seed)
    distribution = np.vstack([
        symmetric_photo_decomposition(a[rng.integers(0, len(a), size=len(a))])
        for _ in range(n_draws)
    ])
    limits = np.quantile(distribution, [0.025, 0.975], axis=0).T
    return point, limits


def summarize_zone(d: pd.DataFrame, *, min_per_season: int,
                   seed: int = 20261017) -> dict:
    # d is restricted to ONE extratropical hemisphere's classified, dated,
    # genuinely sourced photo ROIs in summer or winter.
    require(d.observation_id.notna().all() and d.observation_id.is_unique,
            "Duplicate or missing original observation IDs require explicit observation-level review")
    tab = d.groupby(["inat_taxon_id", "season"]).white.agg(["sum", "size"]).unstack(fill_value=0)
    for s in ("summer", "winter"):
        require(s in tab["size"], "Source lacks a meteorological season")
    same_taxon = tab[(tab["size"]["summer"] >= min_per_season)
                     & (tab["size"]["winter"] >= min_per_season)]
    a = np.stack((same_taxon["size"]["summer"].to_numpy(dtype=float),
                  same_taxon["size"]["winter"].to_numpy(dtype=float),
                  same_taxon["sum"]["summer"].to_numpy(dtype=float),
                  same_taxon["sum"]["winter"].to_numpy(dtype=float)), axis=1)
    require(len(a) >= 2, "At least 2 species must be represented in BOTH seasons")
    point, ci = bootstrap_taxa(a, n_draws=999, seed=seed)
    names = ("within_taxon_hue_rate_component", "taxon_photo_frequency_composition_component",
             "source_photo_weighted_common_taxa_difference", "taxon_equal_mean_difference")
    result = {name + "_pp": round(float(value), 6) for name, value in zip(names, point)}
    result.update({name + "_taxon_bootstrap_percentile95_pp":
                   [round(float(v), 6) for v in ci[i]] for i, name in enumerate(names)})
    result.update({
        "n_common_taxa": int(len(a)),
        "photographs_summer_common_taxa": int(a[:, 0].sum()),
        "photographs_winter_common_taxa": int(a[:, 1].sum()),
        "n_distinct_observations_in_common_taxa": int(a[:, 0].sum() + a[:, 1].sum()),
        "minimum_original_observations_per_taxon_per_season": min_per_season,
        "n_taxon_bootstrap_draws": 999,
        "identical_named_taxon_support": True,
        "exact_symmetric_Kitagawa_decomposition": True,
        "causal_or_population_within_species_response": False,
        "time_and_location_not_conditioned": True,
    })
    return result


def run(photo_zip: Path, outdir: Path) -> dict:
    require(not outdir.exists() or not any(outdir.iterdir()),
            "Never overwrite source-locked output")
    d, rec, source_sha = read_signed(photo_zip, PHOTO_DATA, PHOTO_RECEIPT,
        "source_photo_level_csv_sha256", ["photo_id", "observation_id", "inat_taxon_id",
        "morph", "observed_on", "latitude", "longitude"])
    require(len(d) == MOIS_ROWS and int(d.morph.isin(HUES).sum()) == MOIS_HUES,
            "Frozen original 267425 actual photograph/colour denominator changed")
    require(rec["combined_unique_photo_IDs"] == MOIS_ROWS,
            "Source signature does not refer to actual MOIS 267425 photographed pixels")
    d = labelled_calendar(d)
    require(int((d.classified & d.date.notna()).sum()) == 81867,
            "Date and colour denominator changed")
    d = d[d.classified & d.date.notna() & d.season.isin(["summer", "winter"])
          & d.zone.isin(["north_extratropical", "south_extratropical"])].copy()
    zones = {}
    for name in ("north_extratropical", "south_extratropical"):
        q = d[d.zone.eq(name)]
        season = q.groupby("season").white.agg(["sum", "size"])
        observed_photo_pp = 100 * (season.loc["summer", "sum"] / season.loc["summer", "size"]
                                   - season.loc["winter", "sum"] / season.loc["winter", "size"])
        contrasts = {str(m): summarize_zone(q, min_per_season=m, seed=20261017 + m)
                     for m in (1, 2, 4)}
        zones[name] = {
            "all_original_photo_summer_vs_winter_difference_pp": round(float(observed_photo_pp), 6),
            "all_photo_summer_n": int(season.loc["summer", "size"]),
            "all_photo_winter_n": int(season.loc["winter", "size"]),
            "common_taxa_thresholds": contrasts,
        }
    outdir.mkdir(parents=True, exist_ok=True)
    result = {
        "schema": "mois_267425_real_photo_same_taxon_summer_winter_Kitagawa_v1",
        "source_photo_csv_sha256": source_sha,
        "original_actual_photo_ids": MOIS_ROWS,
        "classifiable_four_hue_original_photo_ids": MOIS_HUES,
        "dated_classified_original_photo_ids": 81867,
        "calculations": zones,
        "method": "Symmetric photographic composition vs within-taxon white-label-rate Kitagawa decomposition, reweighted only to taxa observed in both seasons, 999 whole-taxon bootstrap percentile intervals",
        "important_limits": [
            "Photographic white label only, not botanically verified open flower organs",
            "Photo source sampling and observer effort are not independent population counts",
            "Common taxa may occupy different latitudes/locations and observation years across seasons",
            "Only within matching-taxon source support; effect is not exact attribution of full-atlas all-taxon seasonal difference",
            "Bootstrap percentile intervals not multiple-testing familywise intervals",
            "No FCP D/W polymorphism or within-species genetic selection estimates imported",
        ],
    }
    out = outdir / "mois_267425_real_photo_summer_winter_taxon_frequency_decomposition.json"
    out.write_text(json.dumps(result, ensure_ascii=False, sort_keys=True, indent=2) + "\n")
    print("MOIS_REAL_PHOTO_TAXON_MIX_VS_WHITE_RATE", json.dumps({
        k: {"raw_pp": v["all_original_photo_summer_vs_winter_difference_pp"],
            "common_taxa_at2": v["common_taxa_thresholds"]["2"]["n_common_taxa"],
            "composition_pp_at2": v["common_taxa_thresholds"]["2"]["taxon_photo_frequency_composition_component_pp"]}
        for k, v in zones.items()}, sort_keys=True), flush=True)
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--photo-zip", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    run(args.photo_zip, args.output)
