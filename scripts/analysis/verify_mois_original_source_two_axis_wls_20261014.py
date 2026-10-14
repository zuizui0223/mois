#!/usr/bin/env python3
"""Independent original-source WLS audit for the global mois moisture paper.

Source-only check, separate from migrated FCP estimators; no bootstrap refit.
Use two local full historical GitHub Actions ZIPs; NEVER publish plant-site/photo rows.
Only aggregate cell-free diagnostics are output.
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

ROOT = Path(__file__).resolve().parents[2]
PHOTO_SHA = "2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0"
ENV_SHA = "6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca"
CLASSES = ("white", "yellow_orange", "red_pink", "blue_purple")
CLIMATE = ["wc_elevation_m", *(f"wc_bio{i}" for i in range(1, 12)),
           *(f"wc_bio{i}" for i in range(12, 20)), "wc_srad_annual_mean",
           "wc_srad_monthly_cv", "wc_wind_annual_mean", "wc_wind_monthly_cv",
           "wc_vapr_annual_mean", "wc_vapr_monthly_cv"]
SOIL = ["soil_pH", "soil_SOC", "soil_N", "soil_clay", "soil_available_water_proxy",
        "soil_cec_0_30cm_source_raw", "soil_sand_0_30cm_source_raw",
        "soil_silt_0_30cm_source_raw", "soil_bdod_0_30cm_source_raw",
        "soil_cfvo_0_30cm_source_raw"]
GEO = ("abs_lat", "lon_sin", "lon_cos")
TARGETS = ("wc_bio12", "wc_vapr_annual_mean")


def verify_sha(path: Path, sha: str) -> None:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(chunk)
    assert h.hexdigest() == sha, f"Source ZIP mismatch: {path}"


def read_member(archive: zipfile.ZipFile, name: str, columns: list[str]) -> pd.DataFrame:
    return pd.read_csv(io.BytesIO(archive.read(name)), compression="gzip",
                       usecols=columns, low_memory=False)


def audit(photo_path: Path, environment_path: Path) -> dict:
    verify_sha(photo_path, PHOTO_SHA)
    verify_sha(environment_path, ENV_SHA)
    cols = ["inat_taxon_id", "observation_id", "photo_id", "species", "morph",
            "measurement_status", "photo_license", "site_geo_status"]
    with zipfile.ZipFile(photo_path) as archive:
        assert archive.testzip() is None
        breadth = read_member(archive, "fcp42111_real_site_abiotic/breadth_flower_colour_climate_soil.csv.gz", cols)
        regional = read_member(archive, "fcp42111_real_site_abiotic/taxon_cell_flower_colour_climate_soil.csv.gz", cols)
    assert (len(breadth), len(regional)) == (42111, 85337)
    original = pd.concat([breadth, regional], ignore_index=True)
    assert len(original) - original.photo_id.nunique() == 26905
    original = original.drop_duplicates("photo_id")
    assert len(original) == 100543 and original.inat_taxon_id.nunique() == 42111

    features = CLIMATE + SOIL
    assert len(CLASSES) == 4 and len(CLIMATE) == 26 and len(features) == 36
    envcols = ["inat_taxon_id", "observation_id", "photo_id",
               "latitude", "longitude", *features]
    with zipfile.ZipFile(environment_path) as archive:
        assert archive.testzip() is None
        env = read_member(archive, "fcp42111_original_photo_all_climate_soil/all_original_100543_photo_sites_expanded_environment.csv.gz", envcols)
    assert len(env) == 100543 and env.photo_id.is_unique
    data = original.merge(env, on=["inat_taxon_id", "observation_id", "photo_id"],
                          how="left", validate="one_to_one", indicator=True)
    assert data["_merge"].eq("both").all()
    data["abs_lat"] = data.latitude.abs()
    data["lon_sin"] = np.sin(np.deg2rad(data.longitude))
    data["lon_cos"] = np.cos(np.deg2rad(data.longitude))
    classifiable = data.measurement_status.eq("classified_four_state_morph") & data.morph.isin(CLASSES)
    valid_geo = data.latitude.between(-90, 90) & data.longitude.between(-180, 180)
    climate = classifiable & valid_geo & data[list(GEO) + CLIMATE].notna().all(axis=1)
    soil = climate & data[SOIL].notna().all(axis=1)
    assert (int(classifiable.sum()), int(climate.sum()), int(soil.sum())) == (45953, 45836, 31925)
    assert (data.loc[climate, "inat_taxon_id"].nunique(),
            data.loc[soil, "inat_taxon_id"].nunique()) == (25028, 20786)

    sample = data.loc[climate].copy()
    sample["genus"] = sample.species.str.split().str[0]
    for hue in CLASSES:
        sample["hue_" + hue] = sample.morph.eq(hue).astype(float)
    nuisance = list(GEO) + list(TARGETS)
    mean = sample[list(GEO) + CLIMATE].mean()
    sd = sample[list(GEO) + CLIMATE].std(ddof=0).where(lambda x: x >= 1e-12, 1)
    cols = nuisance + ["hue_" + h for h in CLASSES]
    group = sample.groupby("inat_taxon_id", sort=True)
    species = group[cols].mean()
    species["genus"] = group.genus.first()
    counts = species.genus.value_counts()
    species = species.loc[species.genus.isin(counts[counts >= 3].index)]
    assert len(species) == 20694 and species.genus.nunique() == 1934

    table = pd.read_csv(ROOT / "results/fcp_macro_direction_20261013/direction_ledger.csv")
    assert len(table) == 1280
    estimated = {}
    max_diff = 0.0
    for feature in TARGETS:
        x = ((species[list(GEO) + [feature]] - mean[list(GEO) + [feature]]) / sd[list(GEO) + [feature]]).to_numpy(float)
        y = species[["hue_" + h for h in CLASSES]].to_numpy(float)
        _, ids = np.unique(species.genus.to_numpy(str), return_inverse=True)
        n = np.bincount(ids)
        xbar = np.stack([np.bincount(ids, weights=x[:,i])/n for i in range(x.shape[1])], axis=1)
        ybar = np.stack([np.bincount(ids, weights=y[:,i])/n for i in range(y.shape[1])], axis=1)
        cx,cy = x-xbar[ids],y-ybar[ids]
        w = 1.0 / n[ids]
        beta = np.linalg.solve(cx.T@(w[:,None]*cx), cx.T@(w[:,None]*cy))[-1]
        estimated[feature] = {}
        for hue, val in zip(CLASSES, beta):
            actual = float(val * 100)
            q = table.loc[(table["mode"]=="climate") & (table["quality"]=="all") &
                          (table["region"]=="world") & (table["adjustment"]=="genus_geo") &
                          (table["weight"]=="genus_equal") & (table["anchor"]==feature) &
                          (table["hue"]==hue)]
            assert len(q) == 1
            diff = abs(actual - float(q.iloc[0]["direction_pp_per_env_global_sd"]))
            max_diff = max(max_diff, diff)
            assert diff < 1e-9, (feature, hue, diff)
            estimated[feature][hue] = actual
        assert abs(sum(estimated[feature].values())) < 1e-9

    return {
        "schema": "mois_original_source_independent_two_axis_four_hue_wls_20261014_v1",
        "status": "PASS_INDEPENDENT_ORIGINAL_ZIP_JOIN_AND_EIGHT_WLS_COEFFICIENTS",
        "photo_source_sha256": PHOTO_SHA, "environment_source_sha256": ENV_SHA,
        "source_photo_opportunities": 100543, "named_taxa": 42111,
        "classifiable_photo_hues": 45953, "climate26_classifiable_photo_rows": 45836,
        "climate26_source_species": 25028, "soil36_classifiable_photo_rows": 31925,
        "soil36_source_species": 20786, "genus_eligible_species": 20694,
        "genus_eligible_genera": 1934,
        "eight_point_slopes_pp_per_source_sd": estimated,
        "max_abs_difference_from_original_1280_ledger_pp": max_diff,
        "licence_counts": {str(k): int(v) for k,v in original.photo_license.fillna("UNSPECIFIED").value_counts().items()},
        "source_geo_status_counts": {str(k): int(v) for k,v in original.site_geo_status.fillna("UNSPECIFIED").value_counts().items()},
        "point_estimates_only": True, "bootstrap_CIs_recomputed": False,
        "all_original_photo_pixels_present": False,
        "per_photo_ids_coordinates_and_attribution_released": False,
        "source_rights_and_geoprivacy_human_review_required": True,
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--photo", required=True, type=Path)
    p.add_argument("--environment", required=True, type=Path)
    p.add_argument("--out", type=Path)
    args = p.parse_args()
    result = audit(args.photo, args.environment)
    if args.out:
        args.out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf8")
    print(json.dumps({"status": result["status"], "source_species": result["genus_eligible_species"],
                      "max_abs_coefficient_error_pp": result["max_abs_difference_from_original_1280_ledger_pp"]}))
