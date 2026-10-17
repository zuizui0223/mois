#!/usr/bin/env python3
"""MOIS: honest seasonal photo-hue diagnostics on signed 267425 / 248312 frames.

Outcomes are PHOTO ROI-v4 labels, not botanically verified open flowers.
WorldClim 1970-2000 monthly normals are NOT observation-year weather.
The existing frozen 248312 raster artifact has an unconverted WorldClim srad
column falsely named Wm2; its kJ m-2 day-1 is divided by 86.4 here.
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

HUES = ("white", "yellow_orange", "red_pink", "blue_purple")
MOIS_ROWS, MOIS_HUES = 267425, 129795
ENV_ROWS = 248312
SITE_DEGREES = 0.25
ENERGY_PER_WATT_DAY_KJ = 86.4
PHOTO_DATA = "mois_fcp_wave4_PARTIAL_completed_pixel_original_photo_union.csv.gz"
PHOTO_RECEIPT = "mois_fcp_wave4_PARTIAL_real_completed_original_ROI_receipt.json"
ENV_DATA = "mois_248312_original_photo_worldclim_static_monthly_normal_terrain.csv.gz"
ENV_RECEIPT = "mois_248312_worldclim_static_month_normal_receipt.json"
CLIMATE = {
    "temperature_5C": ("wc2_month_tavg_normal_10m_C", 5.0, 1.0),
    "precipitation_50mm": ("wc2_month_prec_normal_10m_mm", 50.0, 1.0),
    "solar_50Wm2": ("wc2_month_srad_normal_10m_Wm2", 50.0, ENERGY_PER_WATT_DAY_KJ),
}


def require(ok: bool, reason: str) -> None:
    if not ok:
        raise ValueError(reason)


def read_signed(zpath: Path, data_name: str, receipt_name: str,
                sha_key: str, columns: list[str]) -> tuple[pd.DataFrame, dict, str]:
    with zipfile.ZipFile(zpath) as z:
        receipt = json.loads(z.read(receipt_name))
        raw = z.read(data_name)
    sha = hashlib.sha256(raw).hexdigest()
    require(sha == receipt[sha_key], "Input original photo/raster SHA mismatch")
    d = pd.read_csv(io.BytesIO(raw), compression="gzip", usecols=columns,
                    dtype={"photo_id": str, "observation_id": str, "inat_taxon_id": str},
                    low_memory=False)
    require(d.photo_id.notna().all() and d.photo_id.is_unique,
            "Original photo IDs must be unique")
    return d, receipt, sha


def labelled_calendar(d: pd.DataFrame) -> pd.DataFrame:
    """Meteorological seasons local to latitude, tropical months left unbinned.

    December belongs to the NEXT seasonal year so winter/summer spanning the
    New Year is not silently split into two different seasonal years.
    """
    d = d.copy()
    d["date"] = pd.to_datetime(d.observed_on, errors="coerce")
    d["month"] = d.date.dt.month
    d["season_year"] = d.date.dt.year + d.month.eq(12).fillna(False).astype(int)
    d["classified"] = d.morph.isin(HUES)
    d["white"] = d.morph.eq("white")
    d["zone"] = np.select([d.latitude.ge(23.5), d.latitude.le(-23.5)],
                          ["north_extratropical", "south_extratropical"],
                          default="tropical_or_unknown")
    season = ((d.month % 12) // 3 + np.where(d.latitude.lt(0), 2, 0)) % 4
    d["season"] = pd.Series(season, index=d.index).map(
        {0: "winter", 1: "spring", 2: "summer", 3: "autumn"})
    d.loc[d.zone.eq("tropical_or_unknown"), "season"] = "month_only"
    return d


def photo_season_table(d: pd.DataFrame) -> pd.DataFrame:
    d = d[d.date.notna() & d.latitude.between(-90, 90) &
          d.longitude.between(-180, 180)].copy()
    g = d.groupby(["zone", "season"], observed=True, dropna=False)
    out = g.agg(photo_attempts=("photo_id", "size"),
                classified_hues=("classified", "sum"),
                white_hues=("white", "sum"),
                photographed_taxa=("inat_taxon_id", "nunique")).reset_index()
    out["white_share_of_classified"] = out.white_hues / out.classified_hues
    out["classified_fraction_of_attempts"] = out.classified_hues / out.photo_attempts
    return out.sort_values(["zone", "season"]).reset_index(drop=True)


def matched_climate_panel(d: pd.DataFrame, min_per_season: int,
                          cell_width: float = SITE_DEGREES) -> pd.DataFrame:
    """Eligible repeated photos within taxon, approximate cell, and season-year.

    The same named source observation must not be counted twice within a season.
    This is a sampling-support gate, NOT observer or botanical verification.
    """
    q = d[d.classified & d.date.notna() &
          d.zone.isin(["north_extratropical", "south_extratropical"]) &
          d.latitude.between(-90, 90) & d.longitude.between(-180, 180)].copy()
    q["latcell"] = np.floor((q.latitude + 90) / cell_width).astype(int)
    q["loncell"] = np.floor((q.longitude + 180) / cell_width).astype(int)
    q["group"] = (q.inat_taxon_id.astype(str) + ":" + q.latcell.astype(str)
                  + ":" + q.loncell.astype(str) + ":" + q.season_year.astype(int).astype(str))
    season_observations = q.groupby(["group", "season"])["observation_id"].nunique().unstack(fill_value=0)
    admissible = season_observations.index[(season_observations >= min_per_season).sum(axis=1) >= 2]
    q = q[q.group.isin(admissible)].copy()
    # All three physical monthly normals must exist on a matched original photo.
    q = q.dropna(subset=[spec[0] for spec in CLIMATE.values()])
    complete = q.groupby("group").season.nunique()
    return q[q.group.isin(complete.index[complete.ge(2)])]


def equal_group_within_estimate(d: pd.DataFrame, field: str,
                                photo_unit: float, raw_divisor: float,
                                *, bootstrap: int = 999, seed: int = 20261017) -> dict:
    q = d.copy()
    q["x"] = q[field].astype(float) / raw_divisor
    q["y"] = q.white.astype(float)
    g = q.groupby("group")
    q["xx"] = q.x - g.x.transform("mean")
    q["yy"] = q.y - g.y.transform("mean")
    q["weight"] = 1 / g.x.transform("size")
    q["num"] = q.weight * q.xx * q.yy
    q["den"] = q.weight * q.xx * q.xx
    taxa = q.groupby("inat_taxon_id")[["num", "den"]].sum()
    nums = taxa.num.to_numpy()
    dens = taxa.den.to_numpy()
    require(len(taxa) >= 2 and np.sum(dens) > 0,
            "Too little within-site seasonal climatic contrast")
    rng = np.random.default_rng(seed)
    sampled = rng.integers(0, len(taxa), size=(bootstrap, len(taxa)))
    den_boot = dens[sampled].sum(axis=1)
    require(np.all(den_boot > 0), "Bootstrap has no remaining temporal exposure contrast")
    sample = 100 * photo_unit * (nums[sampled].sum(axis=1) / den_boot)
    point = 100 * photo_unit * nums.sum() / dens.sum()
    return {
        "point_white_pp_per_unit": round(float(point), 6),
        "taxon_cluster_bootstrap_percentile95_pp": [round(float(t), 6) for t in np.quantile(sample, [0.025, 0.975])],
        "photo_count": int(len(q)), "taxa": int(len(taxa)),
        "site_taxon_season_year_groups": int(q.group.nunique()),
        "site_taxon_season_year_groups_with_both_hues": int(q.groupby("group").white.nunique().gt(1).sum()),
        "bootstrap_draws": bootstrap, "equal_weight_to_each_taxon_site_season_year_group": True,
        "all_three_predictors_separate_univariate_descriptive_models": True,
        "observed_year_weather_estimated": False,
    }


def analyze(photo_zip: Path, worldclim_zip: Path, outdir: Path) -> dict:
    require(not outdir.exists() or not any(outdir.iterdir()), "Do not overwrite source-locked results")
    a, p, sha_a = read_signed(photo_zip, PHOTO_DATA, PHOTO_RECEIPT,
        "source_photo_level_csv_sha256", ["photo_id", "observation_id", "inat_taxon_id", "morph",
        "observed_on", "latitude", "longitude"])
    require(len(a) == MOIS_ROWS and int(a.morph.isin(HUES).sum()) == MOIS_HUES and
            p["combined_unique_photo_IDs"] == MOIS_ROWS and
            p["combined_four_hue_classified_photo_IDs"] == MOIS_HUES,
            "Measured 267425 original image-photo colour cohort changed")
    b, e, sha_b = read_signed(worldclim_zip, ENV_DATA, ENV_RECEIPT,
        "full_photo_static_monthly_normal_csv_sha256", ["photo_id", "inat_taxon_id", "morph",
        "observed_on", "latitude", "longitude", *[k[0] for k in CLIMATE.values()]])
    require(len(b) == ENV_ROWS and e["n_original_photo_ids"] == ENV_ROWS,
            "Historical WorldClim photo-site cohort changed")
    require(e["schema"] == "mois_248312_real_original_photo_worldclim_static_and_monthly_climate_normals_v1",
            "This reanalysis requires frozen original v1 *unconverted* radiation provenance")
    merged = a.merge(b, on="photo_id", suffixes=("", "_env"), how="inner", validate="one_to_one")
    require(len(merged) == ENV_ROWS, "WorldClim row IDs not subset of expanded photo source")
    for col in ("morph", "inat_taxon_id", "observed_on"):
        require((merged[col].fillna("") == merged[col + "_env"].fillna("")).all(),
                "Incompatible photo identity, taxonomy, hue or observed date across sources")
    a = labelled_calendar(a)
    b = labelled_calendar(merged)
    season = photo_season_table(a)
    models = {}
    for minimum in (1, 2, 3):
        sample = matched_climate_panel(b, minimum)
        rec = {"minimum_distinct_observations_per_season": minimum,
               "site_taxon_season_year_groups": int(sample.group.nunique()),
               "eligible_source_photos": int(len(sample)),
               "eligible_taxa": int(sample.inat_taxon_id.nunique()),
               "n_groups_with_both_hue_states": int(sample.groupby("group").white.nunique().gt(1).sum())}
        if sample.inat_taxon_id.nunique() >= 2:
            rec["univariate_photo_white_vs_month_climate_normal"] = {
                name: equal_group_within_estimate(sample, field, unit, divisor)
                for name, (field, unit, divisor) in CLIMATE.items()}
        models[str(minimum)] = rec
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / "mois_real_267425_dated_photo_season_white_rates.csv"
    season.to_csv(path, index=False, lineterminator="\n", float_format="%.8f")
    result = {
        "schema": "mois_real_original_photo_season_normals_QC_20261017_v1",
        "photo_source_sha256": sha_a, "worldclim_v1_source_sha256": sha_b,
        "date_present_267425": int(a.date.notna().sum()),
        "classified_and_dated_267425": int((a.date.notna() & a.classified).sum()),
        "classified_and_dated_248312_worldclim": int((b.date.notna() & b.classified).sum()),
        "total_real_measured_photo_ids": MOIS_ROWS,
        "classified_real_measured_photo_ids": MOIS_HUES,
        "observed_year_month_climate": "NO; these are WorldClim 1970-2000 MONTH normals",
        "v1_worldclim_srad_column_existing_wrong_unit": "Stored as kJ m-2 day-1, despite misleading _Wm2 name",
        "corrected_srad_unit": "W m-2 = (kJ m-2 day-1) / 86.4",
        "four_season_northern_and_southern_hemisphere_reversed": True,
        "tropical_sites_reported_as_month_only": True,
        "december_assigned_to_following_local_season_year": True,
        "season_table_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "equal_site_year_photo_models": models,
        "inference_limit": ["Only automated photographic ROIs, not botanically independently verified open petals",
            "Original observational photographic effort, taxa and classifiability vary by season",
            "Exact observer IDs not available here; distinct observation does not establish photographer independence",
            "0.25-degree geographic cells are coarse, not identified biological populations",
            "Univariate associations with monthly CLIMATE NORMALS, not actual weather, are noncausal",
            "Exploratory taxon-cluster bootstrap percentile intervals are not multiple-testing maxT intervals"],
    }
    output = outdir / "mois_real_267425_season_climate_normal_decision.json"
    output.write_text(json.dumps(result, ensure_ascii=False, sort_keys=True, indent=2) + "\n")
    print("MOIS_AUTHENTIC_PHOTO_SEASON_QC", json.dumps({
        "dated_classifiable_267k": result["classified_and_dated_267425"],
        "paired_2plus_groups": models["2"]["site_taxon_season_year_groups"],
        "n_taxa": models["2"]["eligible_taxa"]}, sort_keys=True), flush=True)
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--photo-zip", type=Path, required=True)
    p.add_argument("--worldclim-zip", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    analyze(a.photo_zip, a.worldclim_zip, a.output)
