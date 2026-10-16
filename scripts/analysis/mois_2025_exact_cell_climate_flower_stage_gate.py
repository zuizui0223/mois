#!/usr/bin/env python3
"""MOIS: exact-photo observed-month climate join and seasonal inference-feasibility audit.

This is a deterministic source-validated *opportunity/measurement* diagnostic,
NOT a flower-colour frequency, colour-environment relationship or causal test.
The photo & 40-grid-month inputs originate in two authenticated 2026-10-16
GitHub Actions artifacts and retain their original immutable SHA256 checks.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

VARS = (
    "tmean_C", "ppt_mm", "srad_W_m2", "def_mm", "vpd_kPa",
    "tmean_C_anomaly_1991_2020_same_month",
    "ppt_mm_anomaly_1991_2020_same_month",
    "srad_W_m2_anomaly_1991_2020_same_month",
    "def_mm_anomaly_1991_2020_same_month",
    "vpd_kPa_anomaly_1991_2020_same_month",
)
KNOWN = {"white", "yellow_orange", "red_pink", "blue_purple"}
MONTH_TO_SEASON = {f"2025-{n:02d}": ("summer" if n <= 8 else "autumn") for n in range(6, 12)}
MIN_BOTANICAL_CONFIRMED_EACH_SEASON = 20


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def require(test: bool, message: str) -> None:
    if not test:
        raise ValueError(message)


def average(rows: list[dict], key: str) -> float:
    require(bool(rows), "Empty observed-photo group")
    return round(sum(float(r[key]) for r in rows) / len(rows), 6)


def seasonal_support(photos: list[dict], climate: list[dict]) -> dict:
    """Join strictly by actual photographed calendar month and 1/24-degree cell."""
    require(len(photos) == 64 and len(climate) == 40, "Original 64/40 source denominators drift")
    cm = {}
    for r in climate:
        key = (r["year_month"], r["tc_lat_idx"], r["tc_lon_idx"])
        require(key not in cm and key[0] in MONTH_TO_SEASON, "Duplicate/invalid actual climate key")
        require(r["climate_is_actual_2025_not_a_historical_normal"] == "True",
                "Calendar normals mistaken for original observation-year climate")
        require(r["flowering_organ_is_botanically_unreviewed"] == "True",
                "Botanical anthesis cannot be inferred from RGB")
        for field in VARS:
            float(r[field])
        cm[key] = r

    joined = []
    photo_ids, obs_ids, observer_ids = set(), set(), set()
    counts_by_key = Counter()
    for r in photos:
        ym = r["observed_year_month"]
        season = r["local_season"]
        require(ym in MONTH_TO_SEASON and season == MONTH_TO_SEASON[ym],
                "Calendar month/season mismatch")
        require(r["observed_on"].startswith(ym), "Observation date does not match month")
        require(r["botanical_open_flower_review"] == "UNREVIEWED",
                "No external botanical review supplied in source")
        require(r["no_white_rate_until_open_flower_confirmation"] == "True",
                "Upstream source prematurely promoted automated flower colour")
        for col, seen in (("photo_id", photo_ids), ("observation_id", obs_ids),
                          ("observer_id", observer_ids)):
            require(r[col] and r[col] not in seen, "Photo/observation/observer not independent: " + col)
            seen.add(r[col])
        key = (ym, r["tc_lat_idx"], r["tc_lon_idx"])
        require(key in cm, "Original photo lacks exact year-month climate cell")
        counts_by_key[key] += 1
        label = r["roi_morph_unverified"]
        classified = r["roi_measurement_status"] == "classified_four_state_morph"
        require(not classified or label in KNOWN, "Unknown classifiable ROI label")
        mixed = r["same_roi_organ_mix"] == "True"
        require(not mixed or classified, "Mixed-ROI flag on unclassified source")
        joined.append({**r, **{v: cm[key][v] for v in VARS},
                       "is_automatic_classified": classified,
                       "is_unreviewed_clean_photo_colour": classified and not mixed,
                       "is_mixed_roi": mixed})

    require(len(photo_ids) == len(obs_ids) == len(observer_ids) == 64,
            "All 64 source records must be independent")
    require(set(counts_by_key) == set(cm), "40 climate month-keys are not all observed")
    require(all(counts_by_key[k] == int(cm[k]["n_photo_observations"]) for k in cm),
            "Month×cell source counts do not agree with original signed climate requests")
    cells = {(r["tc_lat_idx"], r["tc_lon_idx"]) for r in photos}
    require(len(cells) == 21, "Original geography changed")

    groups = {s: [r for r in joined if r["local_season"] == s]
              for s in ("summer", "autumn")}
    require(all(len(group) == 32 for group in groups.values()),
            "Seasonal source counts are not 32+32")
    support = {}
    for season, rows in groups.items():
        auto = sum(r["is_automatic_classified"] for r in rows)
        clean = sum(r["is_unreviewed_clean_photo_colour"] for r in rows)
        mix = sum(r["is_mixed_roi"] for r in rows)
        support[season] = {
            "n_unique_original_observations": len(rows),
            "n_independent_original_observers": len({r["observer_id"] for r in rows}),
            "n_automatic_four_colour_roi": auto,
            "n_automatic_four_colour_roi_without_mixed_organs": clean,
            "n_mixed_organ_roi": mix,
            "n_unclassifiable_roi": len(rows) - auto,
            "n_botanically_verified_open_flower_colour": 0,
            "max_current_primary_colour_photos_even_if_all_clean_roi_validated": clean,
            "n_distinct_photo_months": len({r["observed_year_month"] for r in rows}),
            "n_distinct_1over24_degree_grid_cells": len({
                (r["tc_lat_idx"], r["tc_lon_idx"]) for r in rows}),
            "observation_weighted_actual_2025_month_climate": {
                var: average(rows, var) for var in VARS},
        }
    require((support["summer"]["n_automatic_four_colour_roi"],
             support["autumn"]["n_automatic_four_colour_roi"]) == (12, 6),
            "Source ROI eligibility not preserved")
    require((support["summer"]["n_automatic_four_colour_roi_without_mixed_organs"],
             support["autumn"]["n_automatic_four_colour_roi_without_mixed_organs"]) == (9, 6),
            "Original same-ROI mixed-colour exclusions changed")

    cross = {
        s: {(r["tc_lat_idx"], r["tc_lon_idx"]) for r in group}
        for s, group in groups.items()
    }
    common = cross["summer"] & cross["autumn"]
    require(len(common) == 7, "Matched source geography drift")
    matched = {s: [r for r in group if (r["tc_lat_idx"], r["tc_lon_idx"]) in common]
               for s, group in groups.items()}
    require((len(matched["summer"]), len(matched["autumn"])) == (15, 24),
            "Matched source site-season observation counts drift")
    paired_cell_diffs = {}
    for var in VARS:
        delta = []
        for cell in sorted(common):
            a = [r for r in matched["autumn"] if (r["tc_lat_idx"], r["tc_lon_idx"]) == cell]
            s = [r for r in matched["summer"] if (r["tc_lat_idx"], r["tc_lon_idx"]) == cell]
            delta.append(average(a, var) - average(s, var))
        paired_cell_diffs[var] = {
            "autumn_minus_summer_cell_equal_mean": round(sum(delta)/len(delta), 6),
            "n_paired_cells": len(delta),
            "n_positive_cells": sum(x > 0 for x in delta),
            "n_negative_cells": sum(x < 0 for x in delta),
        }
    months = {}
    for ym in sorted(MONTH_TO_SEASON):
        rows = [r for r in joined if r["observed_year_month"] == ym]
        months[ym] = {
            "original_photo_attempts": len(rows),
            "automatically_classified_roi": sum(r["is_automatic_classified"] for r in rows),
            "open_flower_botanically_confirmed": 0,
        }
    return {
        "schema": "mois_original_photo_2025_month_climate_join_and_design_gate_v1",
        "status": "SOURCE_CLIMATE_JOIN_VERIFIED_BIOLOGICAL_FLOWERING_OUTCOME_NOT_ESTIMABLE",
        "n_original_photo_attempts": 64,
        "n_distinct_climate_month_cells": len(cm),
        "n_distinct_climate_grid_cells": len(cells),
        "n_unique_observers": len(observer_ids),
        "seasonal_photo_and_climate": support,
        "monthly_original_photo_attempts": months,
        "matched_1over24_degree_cells": {
            "n_cells_with_both_seasons": len(common),
            "n_original_photos_summer_in_matched_cells": len(matched["summer"]),
            "n_original_photos_autumn_in_matched_cells": len(matched["autumn"]),
            "within_cell_actual_climate_contrasts_not_flower_effects": paired_cell_diffs,
        },
        "scientific_gates": {
            "minimum_botanically_verified_primary_photos_per_season": MIN_BOTANICAL_CONFIRMED_EACH_SEASON,
            "current_max_below_threshold_both_seasons": all(
                r["max_current_primary_colour_photos_even_if_all_clean_roi_validated"]
                < MIN_BOTANICAL_CONFIRMED_EACH_SEASON for r in support.values()),
            "all_actual_anthesis_reviews_pending": True,
            "no_repeated_source_observer_to_enable_observer_paired_seasonal_test": True,
            "single_named_species_one_year": True,
            "global_white_flower_seasonality_estimable": False,
            "flower_colour_vs_monthly_climate_coefficient_estimable": False,
            "no_new_colour_selection_or_replacement": True,
        },
        "interpretation": "Real 2025 climate joined to 64 photographed opportunities at the original climate cells. The seven matched cells establish only physical seasonal conditions, not plant flowering hue responses. All source photo flower organs remain unreviewed, classification coverage differs by season, and no white-frequency or climate-effect model was fitted.",
    }


def run(photos: Path, source_keys: Path, source_receipt: Path,
        climate: Path, climate_receipt: Path, output: Path) -> dict:
    s = json.loads(source_receipt.read_text(encoding="utf-8"))
    c = json.loads(climate_receipt.read_text(encoding="utf-8"))
    require(s["status"] == "REAL_DATE_AND_ROI_ATTEMPTS_NO_BOTANICAL_WHITE_CLIMATE_EFFECT",
            "Unrecognised upstream original source state")
    require(s["original_photo_attempts_sha256"] == sha(photos) and
            s["original_month_climate_grid_csv_sha256"] == sha(source_keys),
            "Authentic observed-on photo or month-key input checksum mismatch")
    require(c["status"] == "AUTHENTIC_2025_MONTHLY_CLIMATE_VALUES_NOT_WHITE_FLOWER_ASSOCIATION"
            and c["input_original_2025_photo_climate_request_sha256"] == sha(source_keys)
            and c["input_original_2025_roi_source_receipt_sha256"] == sha(source_receipt)
            and c["decoded_month_grid_climate_sha256"] == sha(climate)
            and c["genuine_original_monthly_photo_observations"] == 64
            and c["provider_request_count"] == 126
            and c["no_botanically_confirmed_open_flower_colour_data"] is True,
            "Original CF-corrected 2025 TerraClimate receipt mismatch")
    photos_rows = read_csv(photos)
    key_rows = read_csv(source_keys)
    climate_rows = read_csv(climate)
    require(len(key_rows) == 40, "Original photo-climate source key count changed")
    climate_by_key = {(r["year_month"], r["tc_lat_idx"], r["tc_lon_idx"]): r
                      for r in climate_rows}
    require(len(climate_by_key) == 40, "Original 40 climate keys are not unique")
    for row in key_rows:
        key = (row["year_month"], row["tc_lat_idx"], row["tc_lon_idx"])
        require(key in climate_by_key and
                int(row["n_photo_observations"]) ==
                int(climate_by_key[key]["n_photo_observations"]),
                "Climate data do not match photo-only original grid-month keys")
    result = seasonal_support(photos_rows, climate_rows)
    result["provenance"] = {
        "photo_observation_table_sha256": sha(photos),
        "observed_month_source_grid_sha256": sha(source_keys),
        "original_post_roi_source_receipt_sha256": sha(source_receipt),
        "decoded_original_2025_climate_sha256": sha(climate),
        "physical_climate_provider_release_receipt_sha256": sha(climate_receipt),
        "original_provider_monthly_requests": c["provider_request_count"],
        "original_provider_CF_version": "TerraClimate V1.1",
        "FCP_photo_input_run": 38909164068,
        "MOIS_photo_month_join_run": 38914197937,
        "MOIS_40gridmonth_climate_run": 38914456992,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2,
                                 sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original-photo-table", type=Path, required=True)
    p.add_argument("--original-grid-month-keys", type=Path, required=True)
    p.add_argument("--original-photo-source-receipt", type=Path, required=True)
    p.add_argument("--decoded-climate-month-table", type=Path, required=True)
    p.add_argument("--physical-climate-source-receipt", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    r = run(a.original_photo_table, a.original_grid_month_keys,
            a.original_photo_source_receipt, a.decoded_climate_month_table,
            a.physical_climate_source_receipt, a.output)
    print("MOIS_2025_REAL_SITE_MONTH_SOURCE_QA_PASS",
          json.dumps({"photos": r["n_original_photo_attempts"],
                      "cells": r["n_distinct_climate_grid_cells"],
                      "common_cells": r["matched_1over24_degree_cells"]["n_cells_with_both_seasons"],
                      "white_flower_effect_estimable": r["scientific_gates"]["flower_colour_vs_monthly_climate_coefficient_estimable"]},
                     sort_keys=True))
