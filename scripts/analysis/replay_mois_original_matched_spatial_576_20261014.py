#!/usr/bin/env python3
"""Independent real-ZIP replay of all 576 original FCP geographic-support values.

Reads original SHA-verified photograph-colour and site-environment TABLE ZIPs.
Recreates exact same-species broad-versus-genus-cell model comparisons in
world/East and both genus/species equal weights; resamples original genera 999
times with the original seeds. It exports only aggregate numeric diagnostics;
photo IDs, species/genus names and plant coordinates NEVER leave memory.
Retrospective observational model check, NOT selection or independent biology.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

PHOTO_SHA = "2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0"
ENV_SHA = "6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca"
KEYS = ["inat_taxon_id", "observation_id", "photo_id"]
HUES = ["white", "yellow_orange", "red_pink", "blue_purple"]
GEO = ["abs_lat", "lon_sin", "lon_cos"]
AXES = ["wc_bio5", "wc_bio12", "wc_bio14", "wc_bio15", "wc_srad_annual_mean", "wc_vapr_annual_mean"]
CLIMATE = ["wc_elevation_m"] + [f"wc_bio{i}" for i in range(1, 20)] + [
    "wc_srad_annual_mean", "wc_srad_monthly_cv", "wc_wind_annual_mean",
    "wc_wind_monthly_cv", "wc_vapr_annual_mean", "wc_vapr_monthly_cv"]
STAGES = ["full_photo_in_matched_genera", "same_species_with_broad_geography",
          "same_species_with_genus_cell"]
DIFFS = ["opportunity_within_matched_genera",
         "spatial_adjustment_on_same_species",
         "combined_coverage_and_spatial_adjustment"]


def read_member(z, member, columns):
    return pd.read_csv(io.BytesIO(z.read(member)), compression="gzip",
                       usecols=columns, low_memory=False)


def original_source(photo_path, env_path):
    for path, sha in ((photo_path, PHOTO_SHA), (env_path, ENV_SHA)):
        with path.open("rb") as stream:
            if hashlib.file_digest(stream, "sha256").hexdigest() != sha:
                raise ValueError("Original full-source ZIP hash mismatch: " + str(path))
    with zipfile.ZipFile(photo_path) as archive:
        cols = KEYS + ["species", "morph", "measurement_status"]
        a = read_member(archive, "fcp42111_real_site_abiotic/breadth_flower_colour_climate_soil.csv.gz", cols)
        b = read_member(archive, "fcp42111_real_site_abiotic/taxon_cell_flower_colour_climate_soil.csv.gz", cols)
    assert (len(a), len(b)) == (42111, 85337)
    photographs = pd.concat((a, b), ignore_index=True).drop_duplicates("photo_id")
    assert len(photographs) == 100543 and photographs.photo_id.is_unique
    with zipfile.ZipFile(env_path) as archive:
        environment = read_member(
            archive, "fcp42111_original_photo_all_climate_soil/all_original_100543_photo_sites_expanded_environment.csv.gz",
            KEYS + ["latitude", "longitude", *CLIMATE])
    assert len(environment) == 100543 and environment.photo_id.is_unique
    data = photographs.merge(environment, on=KEYS, validate="one_to_one", indicator=True)
    assert data._merge.eq("both").all()
    data["abs_lat"] = data.latitude.abs()
    data["lon_sin"] = np.sin(np.deg2rad(data.longitude))
    data["lon_cos"] = np.cos(np.deg2rad(data.longitude))
    data["sector"] = np.clip(np.floor((data.longitude + 180) / 90), 0, 3).astype("Int64")
    classified = data.measurement_status.eq("classified_four_state_morph") & data.morph.isin(HUES)
    geo_ok = data.latitude.between(-90, 90) & data.longitude.between(-180, 180)
    climate_ok = classified & geo_ok & data[GEO + CLIMATE].notna().all(axis=1)
    assert int(classified.sum()) == 45953 and int(climate_ok.sum()) == 45836
    data = data.loc[climate_ok].copy()
    data["genus"] = data.species.str.split().str[0]
    for hue in HUES:
        data["hue_" + hue] = data.morph.eq(hue).astype(float)
    means = data[GEO + CLIMATE].mean()
    sds = data[GEO + CLIMATE].std(ddof=0).where(lambda x: x >= 1e-12, 1.0)
    return data, means, sds


def species_frame(photos):
    grouped = photos.groupby("inat_taxon_id", sort=True)
    z = grouped[GEO + CLIMATE + ["hue_" + h for h in HUES] + ["latitude", "longitude"]].mean()
    z["genus"] = grouped.genus.first()
    z["n_photos"] = grouped.photo_id.size()
    latcell = np.clip(np.floor((np.sin(np.deg2rad(z.latitude)) + 1) * 4.5), 0, 8).astype(int)
    loncell = np.clip(np.floor((z.longitude + 180) / 20), 0, 17).astype(int)
    z["cell"] = 18 * latcell + loncell
    return z


def genus_moments(frame, genera, means, sds, cell=False):
    xtx = np.zeros((len(genera), 6, 4, 4))
    xty = np.zeros_like(xtx)
    for i, genus in enumerate(genera):
        src = frame.loc[frame.genus.eq(genus)]
        n = len(src)
        assert n >= 3
        y = src[["hue_" + h for h in HUES]].to_numpy(float)
        if cell:
            x = ((src[AXES] - means[AXES]) / sds[AXES]).to_numpy(float)
            _, idx = np.unique(src.cell.to_numpy(int), return_inverse=True)
            count = np.bincount(idx)
            assert count.min() >= 3
            x -= np.stack([np.bincount(idx, weights=x[:, k]) / count
                           for k in range(6)], axis=1)[idx]
            y -= np.stack([np.bincount(idx, weights=y[:, k]) / count
                           for k in range(4)], axis=1)[idx]
            for j in range(6):
                xtx[i, j, 3, 3] = (x[:, j] @ x[:, j]) / n
                xty[i, j, 3, :] = (x[:, j] @ y) / n
        else:
            x = ((src[GEO + AXES] - means[GEO + AXES]) / sds[GEO + AXES]).to_numpy(float)
            x -= x.mean(axis=0)
            y -= y.mean(axis=0)
            for j in range(6):
                v = x[:, [0, 1, 2, j + 3]]
                xtx[i, j] = v.T @ v / n
                xty[i, j] = v.T @ y / n
    return xtx, xty


def solve_bootstrap(draws, xtx, xty, cell=False):
    x = xtx.sum(axis=0)
    y = xty.sum(axis=0)
    bx = np.einsum("bg,gjkl->bjkl", draws, xtx, optimize=True)
    by = np.einsum("bg,gjkl->bjkl", draws, xty, optimize=True)
    if cell:
        point = y[:, 3, :] / x[:, 3, 3, None]
        bootstrap = by[:, :, 3, :] / bx[:, :, 3, 3, None]
    else:
        point = np.linalg.solve(x, y)[:, 3, :]
        bootstrap = np.linalg.solve(bx, by)[:, :, 3, :]
    assert np.allclose(point.sum(axis=1), 0, atol=1e-9)
    return point, bootstrap


def replay(photo_path: Path, env_path: Path, ledger_path: Path):
    photos, mean, sd = original_source(photo_path, env_path)
    ledger = pd.read_csv(ledger_path)
    assert len(ledger) == 576
    checked = {"points": 0, "percentile_interval_endpoints": 0,
               "maxT_interval_endpoints": 0}
    maximum_error = 0.0
    estimates = {}
    for region, weight in (("world", "genus_equal"), ("E90_E180", "genus_equal"),
                           ("world", "species_equal"), ("E90_E180", "species_equal")):
        regional = photos if region == "world" else photos.loc[photos.sector.eq(3)]
        full = species_frame(regional)
        key = full.genus.astype(str) + "|" + full.cell.astype(str)
        counts = key.value_counts()
        chosen = full.loc[key.map(counts).ge(3)].copy()
        genera = sorted(chosen.genus.unique())
        wide = full.loc[full.genus.isin(genera)].copy()
        target = (9370, 845, 14268) if region == "world" else (2886, 280, 4149)
        assert (len(chosen), len(genera), int(chosen.n_photos.sum())) == target
        count_a = wide.genus.value_counts().reindex(genera).to_numpy(int)
        count_b = chosen.genus.value_counts().reindex(genera).to_numpy(int)
        seed = 2026101422 + len(region) * 11
        draws = np.random.default_rng(seed).multinomial(
            len(genera), np.repeat(1 / len(genera), len(genera)), size=999)
        fits = {}
        for name, frame, iscell, count in (
            ("A", wide, False, count_a), ("B", chosen, False, count_b),
            ("C", chosen, True, count_b)):
            x, y = genus_moments(frame, genera, mean, sd, cell=iscell)
            if weight == "species_equal":
                x *= count[:, None, None, None]
                y *= count[:, None, None, None]
            fits[name] = solve_bootstrap(draws, x, y, cell=iscell)
        entries = [("stage", STAGES[i], *fits[name]) for i, name in enumerate(("A", "B", "C"))]
        entries += [
            ("paired_difference", DIFFS[0], fits["A"][0] - fits["B"][0], fits["A"][1] - fits["B"][1]),
            ("paired_difference", DIFFS[1], fits["B"][0] - fits["C"][0], fits["B"][1] - fits["C"][1]),
            ("paired_difference", DIFFS[2], fits["A"][0] - fits["C"][0], fits["A"][1] - fits["C"][1])
        ]
        for kind, label, point, bs in entries:
            low, high = np.quantile(bs, [0.025, 0.975], axis=0)
            if kind == "paired_difference":
                dev = bs - point[None, :, :]
                devsd = dev.std(axis=0, ddof=1)
                keep = devsd > 1e-12
                crit = np.quantile(np.max(np.abs(dev[:, keep] / devsd[keep]), axis=1), 0.95)
                simultaneous_low = point - crit * devsd
                simultaneous_high = point + crit * devsd
            for j, axis in enumerate(AXES):
                for k, hue in enumerate(HUES):
                    match = ledger.loc[(ledger.region == region) & (ledger.weight == weight) &
                        (ledger.type == kind) & (ledger.comparison == label) &
                        (ledger.anchor == axis) & (ledger.hue == hue)]
                    assert len(match) == 1
                    orig = match.iloc[0]
                    components = [(point[j, k] * 100, orig.point_pp, "points"),
                                  (low[j, k] * 100, orig.CI95_lo_pp, "percentile_interval_endpoints"),
                                  (high[j, k] * 100, orig.CI95_hi_pp, "percentile_interval_endpoints")]
                    if kind == "paired_difference":
                        components += [
                            (simultaneous_low[j, k] * 100, orig.maxT95_lo_pp, "maxT_interval_endpoints"),
                            (simultaneous_high[j, k] * 100, orig.maxT95_hi_pp, "maxT_interval_endpoints")]
                    for now, then, field in components:
                        delta = abs(now - then)
                        maximum_error = max(maximum_error, delta)
                        checked[field] += 1
                        if delta > 1e-8:
                            raise AssertionError((region, weight, label, axis, hue, field, delta))
            if kind == "paired_difference" and label == DIFFS[1]:
                estimates[region + "|" + weight] = {
                    "n_matched_species": len(chosen),
                    "n_matched_genera": len(genera),
                    "n_matched_photos": int(chosen.n_photos.sum()),
                    "B_white_BIO12_pp": float(fits["B"][0][1, 0] * 100),
                    "C_white_BIO12_pp": float(fits["C"][0][1, 0] * 100),
                    "B_minus_C_white_BIO12_pp": float(point[1, 0] * 100),
                    "B_minus_C_white_BIO12_maxT95_pp": [
                        float(simultaneous_low[1, 0] * 100),
                        float(simultaneous_high[1, 0] * 100)]}
    assert checked == {"points": 576, "percentile_interval_endpoints": 1152,
                       "maxT_interval_endpoints": 576}
    return {"schema": "mois_source_zip_independent_spatial_576_ci_20261014_v1",
        "status": "PASS_ORIGINAL_576_POINTS_PERCENTILE_AND_MAXT_INTERVALS",
        "photo_sha256": PHOTO_SHA, "environment_sha256": ENV_SHA,
        "n_original_source_rows": 576, "n_bootstrap": 999,
        "comparison_fields": checked, "max_abs_difference_pp": maximum_error,
        "matched_results": estimates,
        "not_independent_biology": True, "not_causal_geography_effect": True,
        "programme_wide_selection_not_adjusted": True,
        "no_photo_or_taxon_locations_exported": True,
        "source_licence_review_required": True}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--photo", type=Path, required=True)
    p.add_argument("--environment", type=Path, required=True)
    p.add_argument("--ledger", type=Path, required=True)
    p.add_argument("--out", type=Path)
    args = p.parse_args()
    result = replay(args.photo, args.environment, args.ledger)
    if args.out:
        args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf8")
    print(json.dumps({k:v for k,v in result.items() if k != "matched_results"}, indent=2))


if __name__ == "__main__":
    main()
