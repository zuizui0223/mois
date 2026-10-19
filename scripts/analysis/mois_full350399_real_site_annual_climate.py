#!/usr/bin/env python3
"""Source-authentic WorldClim 1970-2000 extension to 350399 genuine MOIS ROI photos.

Never resample the existing 317631 photographed site pixels/environment,
never assign hue to 93m candidate photo IDs, and never fill unmeasured soil.
Original added 32768 pixels were source-measured across four independent
taxon-prioritised image acquisition cohorts. All 33k *new* photo sites
receive raster lookups only after the five input source SHA checks.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis.mois_actual_photo_cohort_season_audit_20261019 import (
    HUES, assemble, require, sha256, summary
)
from scripts.analysis import mois_original313535_full19_bio_vapr_terrain_lake as old
from scripts.analysis import mois_global313535_water_light_heat_habitat_partition as climate
from scripts.analysis.mois_sample_248312_worldclim_static_and_monthnormals import (
    sample_tiff, raster_terrain
)

N_CORE = 317631
N_NEW = 32768
N_TOTAL = N_CORE + N_NEW
N_CLASSIFIED = 164297
N_WHITE = 78026
PHOTO = "mois350399_original_rgb_photo_hues_sourcearms.parquet"
ENVIRONMENT = "mois350399_original_site_worldclim_bio19_vapr_solar_terrain.parquet"
RECEIPT = "mois350399_original_rgb_photo_full_environment_receipt.json"
COVERAGE = "mois350399_photo_environment_coverage_by_arm.csv"
RASTER_FILES = {
    "bio": "wc2.1_10m_bio.zip",
    "srad": "wc2.1_10m_srad.zip",
    "vapr": "wc2.1_10m_vapr.zip",
    "elev": "wc2.1_2.5m_elev.zip",
}


def raster_site_variables(photos: pd.DataFrame, raster: Path):
    """Sample immutable official annual raster normals at NEW actual-photo sites."""
    ids = photos.photo_id.astype(str).to_numpy()
    require(len(ids) == N_NEW and len(set(ids)) == N_NEW,
            "Expected precisely 32768 new source photo sites, not metadata-only candidates")
    lon = pd.to_numeric(photos.longitude, errors="coerce").to_numpy(float)
    lat = pd.to_numeric(photos.latitude, errors="coerce").to_numpy(float)
    data = {"photo_id": ids}
    for i, name in enumerate(old.BIO, 1):
        path = raster / "bio" / f"wc2.1_10m_bio_{i}.tif"
        require(path.is_file(), "Missing pinned original WorldClim BIO raster: " + str(path))
        data[name] = sample_tiff(path, lon, lat)
    sums = {k: np.zeros(N_NEW, dtype=np.float64) for k in ("vapr", "solar")}
    seen = {k: np.zeros(N_NEW, dtype=np.int8) for k in ("vapr", "solar")}
    for m in range(1, 13):
        for folder, kind, factor in (("vapr", "vapr", 1), ("srad", "solar", 86.4)):
            path = raster / folder / f"wc2.1_10m_{folder}_{m:02d}.tif"
            require(path.is_file(), "Missing original 1970-2000 monthly archive raster " + str(path))
            a = sample_tiff(path, lon, lat).astype(float) / factor
            ok = np.isfinite(a)
            sums[kind][ok] += a[ok]
            seen[kind][ok] += 1
    for kind, out in (("vapr", old.VAPR), ("solar", climate.SOLAR)):
        result = np.full(N_NEW, np.nan, dtype=np.float32)
        ok = seen[kind] == 12
        result[ok] = (sums[kind][ok] / 12).astype(np.float32)
        data[out] = result
    dem = raster / "elev" / "wc2.1_2.5m_elev.tif"
    require(dem.is_file(), "Original site terrain raster not restored")
    terrain = raster_terrain(dem, lon, lat)
    for field in old.TERRAIN:
        data[field] = terrain[field]
    result = pd.DataFrame(data)
    require(len(result) == N_NEW and result.photo_id.is_unique,
            "One-to-one original 32768 site raster lookups failed")
    return result


def assemble_lake(photo_files, core_environment: Path, core_receipt: Path,
                  raster: Path, output: Path, sampler=raster_site_variables):
    require(not output.exists() or not any(output.iterdir()),
            "Cannot overwrite a completed source-verified original-photo result")
    all_pixels, arms = assemble(*photo_files, strict=True)
    require(len(all_pixels) == N_TOTAL and all_pixels.photo_id.is_unique
            and int(all_pixels.classified.sum()) == N_CLASSIFIED
            and int(all_pixels.white.sum()) == N_WHITE,
            "Actual original RGB photographs/hues differ from signed 350399 source")
    core, *extras = arms
    new_pixels = pd.concat(extras, ignore_index=True)
    require(len(new_pixels) == N_NEW and
            not new_pixels.photo_id.isin(core.photo_id).any(),
            "Original 32768 separately measured photo identity not verified")
    require(core_environment.is_file() and core_receipt.is_file(),
            "Missing original immutable 317631 photo-environment receipt")
    receipt = json.loads(core_receipt.read_text())
    require(receipt["schema"] == "mois_real317631_fifth4096_incremental_RGB_annual_bio19_terrain_v1"
            and receipt["all_original_photo_attempts"] == N_CORE
            and receipt["all_classified_four_hue_photo_ROIs"] == 152827
            and receipt["all_white_photo_ROIs"] == 72175,
            "Previous 317631 real original RGB/photo population changed")
    expected_core_sha = receipt["photo_ID_aligned_Zstd_files_sha256"].get(core_environment.name)
    require(expected_core_sha == sha256(core_environment),
            "Previous 317631 photo environment Parquet modified")
    previous_env = pd.read_parquet(core_environment)
    require(len(previous_env) == N_CORE and previous_env.photo_id.is_unique
            and previous_env.photo_id.astype(str).tolist() == core.photo_id.tolist(),
            "Previous source photo IDs and annual climate table were misaligned")
    needed = ("photo_id", *old.BIO, old.VAPR, climate.SOLAR, *old.TERRAIN)
    require(set(needed) == set(previous_env.columns),
            "Unexpected previously verified all19/VAPR/solar/terrain schema")
    expected_source = receipt["source_original_raster_sha256"]
    sources = {}
    for category, filename in RASTER_FILES.items():
        src = raster / category / filename
        require(src.is_file(), "Missing physically extracted source WorldClim ZIP")
        actual = sha256(src)
        require(actual == expected_source[category],
                "New raster snapshot differs from source used for earlier photos: " + category)
        sources[category] = actual
    require(sources["bio"] == old.PIN_BIO
            and sources["srad"] == old.PIN_SRAD,
            "Exact original annual climate/solar source or units changed")
    newly_sampled = sampler(new_pixels, raster)
    newly_sampled = newly_sampled[list(previous_env.columns)]
    combined_env = pd.concat([previous_env, newly_sampled], ignore_index=True)
    combined_env["photo_id"] = combined_env.photo_id.astype(str)
    require(len(combined_env) == N_TOTAL and combined_env.photo_id.is_unique
            and combined_env.photo_id.tolist() == all_pixels.photo_id.tolist(),
            "Original photo/annual climate source is not one-to-one 350399")
    require(combined_env.iloc[:N_CORE].drop(columns=["photo_id"]).reset_index(drop=True).equals(
                previous_env.drop(columns=["photo_id"]).reset_index(drop=True)),
            "Frozen 317631 previous climate source values were recomputed")
    photo_data = all_pixels[
        ["photo_id", "inat_taxon_id", "species", "morph", "observation_id",
         "observed_on", "latitude", "longitude", "classified", "white", "arm"]
    ].copy()
    require(photo_data.photo_id.tolist() == combined_env.photo_id.tolist()
            and int(photo_data.classified.sum()) == N_CLASSIFIED,
            "Frozen ROI hues/source arm and climate sites no longer align")
    output.mkdir(parents=True, exist_ok=True)
    photo_data.to_parquet(output / PHOTO, index=False, compression="zstd")
    combined_env.to_parquet(output / ENVIRONMENT, index=False, compression="zstd")
    env_cols = [*old.BIO, old.VAPR, climate.SOLAR, *old.TERRAIN]
    rows = []
    for arm, x in photo_data.groupby("arm", sort=False):
        ids = x.index
        s = combined_env.iloc[ids]
        q = x.classified.to_numpy()
        rows.append({
            "arm": arm,
            "original_rgb_attempts": len(x),
            "four_hue_classified": int(q.sum()),
            "white_four_hue": int(x.white.sum()),
            "annual_rain_available_photo_sites": int(s[climate.RAIN].notna().sum()),
            "all21_climate_available_photo_sites": int(s[[*old.BIO, old.VAPR, climate.SOLAR]].notna().all(axis=1).sum()),
            "all_7_terrain_available_photo_sites": int(s[list(old.TERRAIN)].notna().all(axis=1).sum()),
        })
    pd.DataFrame(rows).to_csv(output / COVERAGE, index=False)
    source_hashes = {
        key: sha256(path) for key, path in zip(
            ("core317631", "broad4096_cumulative", "wave24096", "cached8192", "tail16384"), photo_files)
    }
    s = summary(all_pixels)
    rec = {
        "schema": "mois350399_real_original_photo_bio19_vapr_solar_dem_full_lake_v1",
        "actually_processed_original_photo_pixels": N_TOTAL,
        "genuinely_new_photo_sites_sampled_now": N_NEW,
        "earlier_original_photo_sites_climate_unchanged": N_CORE,
        "actually_four_hue_classifiable_source_photos": N_CLASSIFIED,
        "white_original_photo_ROIs": N_WHITE,
        "taxa_with_classified_hues": s["nominal_taxa_with_classifiable_hues"],
        "genera_with_classified_hues": s["nominal_genera_with_classifiable_hues"],
        "photo_source_parquet_sha256": source_hashes,
        "previous317631_annual_environment_original_sha256": expected_core_sha,
        "worldclim1970_2000_source_original_raster_sha256": sources,
        "annual_solar_kJ_m2_day_converted_to_Wm2": True,
        "annual_climate_not_actual_observation_year_weather": True,
        "only_new_photo_sites_raster_sampled": True,
        "old_real_photo_hues_and_climate_unchanged": True,
        "soil_unmeasured_for_new_taxa_not_imputed": True,
        "source_arms_separate_not_globally_probability_weighted": True,
        "no_original_licensed_photo_bytes_redistributed": True,
        "not_botanically_confirmed_flower_petal_anthesis": True,
        "output_sha256": {name: sha256(output/name)
                          for name in (PHOTO, ENVIRONMENT, COVERAGE)},
    }
    (output / RECEIPT).write_text(
        json.dumps(rec, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False) + "\n")
    print("MOIS_FULL_350399_ACTUAL_PIXEL_SITE_CLIMATE_COMPLETE", json.dumps({
        "photo_attempts": N_TOTAL,
        "new_real_site_climate": N_NEW,
        "four_hue": N_CLASSIFIED,
        "white": N_WHITE,
        "raster_source_sha256": sources,
        "coverage": rows
    }, sort_keys=True), flush=True)
    return rec


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("core", "broad", "wave2", "cached8192", "tail16384",
                 "core_environment", "core_receipt", "raster", "output"):
        p.add_argument("--" + name.replace("_", "-"), type=Path, required=True)
    a = p.parse_args()
    assemble_lake(
        (a.core, a.broad, a.wave2, a.cached8192, a.tail16384),
        a.core_environment, a.core_receipt, a.raster, a.output)


if __name__ == "__main__":
    main()
