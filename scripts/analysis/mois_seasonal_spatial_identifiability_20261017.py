#!/usr/bin/env python3
"""Source-locked MOIS seasonal photographic taxon turnover and geographic positivity.

Photographic ROI hues are NOT confirmed petals or biological flower polymorphism.
Do not mistake all-site taxon seasonal overlap for same-place/year resampling.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scripts.analysis.mois_seasonal_267k_climate_normals import (
    PHOTO_DATA, PHOTO_RECEIPT, MOIS_ROWS, MOIS_HUES,
    HUES, read_signed, labelled_calendar, require,
)
from scripts.analysis.mois_taxon_season_composition_decomposition import symmetric_photo_decomposition

ZONES = ("north_extratropical", "south_extratropical")
GRIDS = (0.25, 0.5, 1.0, 2.0, 5.0)


def full_photo_turnover(q: pd.DataFrame) -> dict:
    """Additive exact photographic source decomposition INCLUDING season-only taxa."""
    s = q[q.season.eq("summer")].copy()
    w = q[q.season.eq("winter")].copy()
    summer_taxa, winter_taxa = set(s.inat_taxon_id), set(w.inat_taxon_id)
    common = summer_taxa & winter_taxa
    sc = s[s.inat_taxon_id.isin(common)]
    wc = w[w.inat_taxon_id.isin(common)]
    sx = s[~s.inat_taxon_id.isin(common)]
    wx = w[~w.inat_taxon_id.isin(common)]
    tab = q[q.inat_taxon_id.isin(common)].groupby(["inat_taxon_id", "season"]).white.agg(["size", "sum"]).unstack(fill_value=0)
    a = np.column_stack((tab["size"]["summer"], tab["size"]["winter"],
                         tab["sum"]["summer"], tab["sum"]["winter"])).astype(float)
    within, mix, common_delta, _ = symmetric_photo_decomposition(a)
    summer_only = 100 * len(sx) / len(s) * (sx.white.mean() - sc.white.mean())
    winter_only = -100 * len(wx) / len(w) * (wx.white.mean() - wc.white.mean())
    full_delta = 100 * (s.white.mean() - w.white.mean())
    require(np.isclose(full_delta, within + mix + summer_only + winter_only, atol=1e-9),
            "Complete measured-photo original seasonal hue partition is not additive")
    return {
        "all_photo_summer_minus_winter_white_pp": round(float(full_delta), 6),
        "n_taxa_summer": len(summer_taxa), "n_taxa_winter": len(winter_taxa),
        "n_taxa_common_both_seasons": len(common),
        "n_taxa_observed_only_summer": len(summer_taxa - winter_taxa),
        "n_taxa_observed_only_winter": len(winter_taxa - summer_taxa),
        "photographic_taxa_jaccard": round(len(common) / len(summer_taxa | winter_taxa), 6),
        "summer_photos": len(s), "winter_photos": len(w),
        "summer_photo_ids_only_season_taxa": len(sx),
        "winter_photo_ids_only_season_taxa": len(wx),
        "photographic_only_season_summer_fraction": round(len(sx)/len(s), 6),
        "photographic_only_season_winter_fraction": round(len(wx)/len(w), 6),
        "common_taxa_within_taxon_white_photo_rate_pp": round(float(within), 6),
        "common_taxa_photo_frequency_mix_pp": round(float(mix), 6),
        "common_taxa_photo_weighted_total_pp": round(float(common_delta), 6),
        "observed_only_summer_taxa_support_component_pp": round(float(summer_only), 6),
        "observed_only_winter_taxa_support_component_pp": round(float(winter_only), 6),
        "is_exact_full_source_photo_algebra_not_phenological_causation": True,
    }


def eligible_site_photo_support(q: pd.DataFrame, degrees: float, same_year: bool,
                                min_observations_per_season: int) -> dict:
    """Minimum distinct original observations for BOTH local winter and summer."""
    d = q.copy()
    d["latcell"] = np.floor((d.latitude + 90) / degrees).astype("int32")
    d["loncell"] = np.floor((d.longitude + 180) / degrees).astype("int32")
    fields = ["inat_taxon_id", "latcell", "loncell"] + (["season_year"] if same_year else [])
    tab = d.groupby(fields + ["season"], observed=True).agg(
        obs=("observation_id", "nunique"), n=("photo_id", "size"), white=("white", "sum")).unstack(fill_value=0)
    common = tab[(tab["obs"]["summer"] >= min_observations_per_season) &
                 (tab["obs"]["winter"] >= min_observations_per_season)]
    if common.empty:
        return {"matched_groups": 0, "matched_taxa": 0, "matched_photos": 0,
                "matched_original_observations": 0}
    return {"matched_groups": int(len(common)),
            "matched_taxa": int(common.index.get_level_values("inat_taxon_id").nunique()),
            "matched_photos": int(common["n"]["summer"].sum() + common["n"]["winter"].sum()),
            "matched_original_observations": int(common["obs"]["summer"].sum() + common["obs"]["winter"].sum())}


def run(archive: Path, output: Path) -> dict:
    require(not output.exists() or not any(output.iterdir()), "Do not overwrite frozen ecological results")
    d, rec, source_sha = read_signed(archive, PHOTO_DATA, PHOTO_RECEIPT,
        "source_photo_level_csv_sha256",
        ["photo_id", "observation_id", "inat_taxon_id", "morph", "observed_on",
         "latitude", "longitude", "old_mois_or_new_fcp"])
    require(len(d) == MOIS_ROWS and int(d.morph.isin(HUES).sum()) == MOIS_HUES
            and rec["combined_unique_photo_IDs"] == MOIS_ROWS,
            "Signed original full photo and classified hue denominator changed")
    d = labelled_calendar(d)
    require(int((d.classified & d.date.notna()).sum()) == 81867, "Measured original date+hue denominator changed")
    d = d[d.classified & d.date.notna() & d.season.isin(("summer","winter")) &
          d.zone.isin(ZONES) & d.latitude.between(-90,90) &
          d.longitude.between(-180,180)].copy()
    regions, rows = {}, []
    for name in ZONES:
        q = d[d.zone.eq(name)].copy()
        require(q.observation_id.notna().all() and q.observation_id.is_unique,
                "Repeated original observation IDs; naive distinct photographer/photo counts prohibited")
        raw = full_photo_turnover(q)
        by_source = {}
        for source, p in q.groupby("old_mois_or_new_fcp"):
            seasons = p.groupby("season").white.agg(["sum","size"])
            if {"summer","winter"} <= set(seasons.index):
                ss = seasons.loc["summer"]
                ww = seasons.loc["winter"]
                by_source[str(source)] = {
                    "original_summer_photos": int(ss["size"]), "original_winter_photos": int(ww["size"]),
                    "summer_white_label_share": round(float(ss["sum"]/ss["size"]),6),
                    "winter_white_label_share": round(float(ww["sum"]/ww["size"]),6),
                    "summer_minus_winter_white_pp": round(float(100*(ss["sum"]/ss["size"] - ww["sum"]/ww["size"])),6),
                }
        regions[name] = {"full_source_photo_turnover": raw, "source_stratified": by_source}
        for grid in GRIDS:
            for year in (False, True):
                for minimum in (1,2):
                    rows.append({"hemisphere":name, "cell_degrees":grid,
                        "same_local_meteorological_season_year":year,
                        "min_distinct_observations_each_season":minimum,
                        **eligible_site_photo_support(q, grid, year, minimum)})
    output.mkdir(parents=True, exist_ok=True)
    csv = output / "mois_seasonal_original_photo_spatial_resolution_support.csv"
    pd.DataFrame(rows).to_csv(csv, index=False, lineterminator="\n")
    result = {
        "schema": "mois_photo_white_original267425_temporal_spatial_support_20261017_v1",
        "source_photo_union_SHA256": source_sha, "total_original_photo_ids": MOIS_ROWS,
        "four_hue_original_photo_ids": MOIS_HUES,
        "dated_and_classified_original_photo_ids": 81867,
        "regions": regions,
        "spatial_support_sha256": hashlib.sha256(csv.read_bytes()).hexdigest(),
        "spatial_result_count": len(rows),
        "site_resolution_degrees": list(GRIDS),
        "same_year_means_local_meteorological_season_year_not_calendar_photo_year": True,
        "interpretation": "Measured-photo taxon turnover / photo-frequency composition, not confirmed flowering phenology or pigment selection",
        "flower_stage_independently_verified": False,
        "historical_100543_dated_photos_in_this_source": 0,
        "actual_observed_month_weather_or_anomaly_modelled": False,
        "same_observation_does_not_guarantee_distinct_photographer_or_population": True,
        "coarse_cell_is_not_biological_population": True,
        "no_FCP_intraspecific_W_D_or_50km_spatial_colour_estimator": True,
    }
    (output/"mois_seasonal_original_photo_spatial_identifiability_receipt.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True)+"\n")
    print("MOIS_SIGNED_ORIGINAL_SPATIAL_POSITIVITY", json.dumps({
        z:{"taxa_both_summer_winter":regions[z]["full_source_photo_turnover"]["n_taxa_common_both_seasons"],
           "site_year_cell_0p25_2photo_groups":next(r["matched_groups"] for r in rows
                if r["hemisphere"]==z and r["cell_degrees"]==0.25 and
                   r["same_local_meteorological_season_year"] and
                   r["min_distinct_observations_each_season"]==2)}
         for z in ZONES}, sort_keys=True))
    return result

if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--photo-archive",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    run(args.photo_archive,args.output)
