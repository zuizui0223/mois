#!/usr/bin/env python3
"""Audit genuinely measured MOIS photographic colour cohorts and source seasons.

Never pool intentionally taxon-enriched cohorts as a probability sample.
Never turn unsuccessful original photo ROI attempts into non-white flowers.
Only the original-photo pixel outcomes, not metadata-only candidates, enter
the hue statistics. Same-taxon/local/year seasonal overlap is a diagnostic
of observational support, not a test of colour plasticity or climate causality.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

HUES = ("white", "yellow_orange", "red_pink", "blue_purple")
COUNTS = {"core": 317631, "broad": 4096, "wave2": 4096, "cached8192": 8192,
          "tail16384": 16384}
EXPECTED_BASE_BROAD_OVERLAP = 313535


def require(value, message):
    if not value:
        raise ValueError(message)


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def load_pixels(path, arm):
    raw = pd.read_parquet(path)
    required = {"photo_id", "inat_taxon_id", "species", "morph",
                "observation_id", "latitude", "longitude"}
    require(required <= set(raw), f"{arm}: missing actual photo fields {required - set(raw)}")
    d = raw.copy()
    for key in ("photo_id", "inat_taxon_id", "species", "morph", "observation_id"):
        d[key] = d[key].fillna("").astype(str)
    require(d.photo_id.str.fullmatch(r"\d+").all() and d.photo_id.is_unique,
            f"{arm}: blank/duplicate/non-numeric original photo IDs")
    # Historical core releases use multiple legitimate technical unknown codes.
    # Only these four exact hue labels count as classified; retain all other
    # source status strings verbatim for review, never assume they are chromatic.
    if "classified" in d:
        reported = d["classified"].astype(str).str.lower().eq("true")
        require(reported.eq(d.morph.isin(HUES)).all(),
                f"{arm}: original four-hue classified flag contradicts source morph")
    if "classified_four_hue" in d:
        actual = d.classified_four_hue.astype(str).str.lower().eq("true")
        require(actual.eq(d.morph.isin(HUES)).all(),
                f"{arm}: classifier success flag contradicts original hue")
    if "measurement_status" in d:
        known = d.morph.isin(HUES)
        require(known.eq(d.measurement_status.eq("classified_four_state_morph")).all(),
                f"{arm}: original photo technical status contradicts hue")
    d["observed_on"] = (d["qualified_date"] if "qualified_date" in d
                        else d["observed_on"] if "observed_on" in d
                        else "")
    d["date"] = pd.to_datetime(d.observed_on, errors="coerce")
    d["latitude"] = pd.to_numeric(d.latitude, errors="coerce")
    d["longitude"] = pd.to_numeric(d.longitude, errors="coerce")
    d["arm"] = arm
    d["white"] = d.morph.eq("white")
    d["classified"] = d.morph.isin(HUES)
    d["technical_failure"] = (d.measurement_status.eq("image_acquisition_failed")
                              if "measurement_status" in d else pd.NA)
    return d


def assemble(core, broad, wave2, cached8192, tail16384=None, strict=True):
    c = load_pixels(core, "core")
    b = load_pixels(broad, "broad")
    require(not strict or len(c) == COUNTS["core"], "Core is not 317631 original pixels")
    ids = set(c.photo_id)
    overlap = b.loc[b.photo_id.isin(ids)]
    require(not strict or len(overlap) == EXPECTED_BASE_BROAD_OVERLAP,
            "Expected frozen 313535 shared baseline between independent ROI arms")
    same = overlap.merge(
        c[["photo_id", "inat_taxon_id", "species", "morph"]],
        on="photo_id", suffixes=("_broad", "_core"), validate="one_to_one")
    for col in ("inat_taxon_id", "species", "morph"):
        require(same[f"{col}_broad"].equals(same[f"{col}_core"]),
                f"Shared originally photographed {col} was altered")
    extra = b.loc[~b.photo_id.isin(ids)].copy()
    d2 = load_pixels(wave2, "wave2")
    d8 = load_pixels(cached8192, "cached8192")
    waves = [c, extra, d2, d8]
    if tail16384 is not None:
        waves.append(load_pixels(tail16384, "tail16384"))
    for name, frame in zip(("core", "broad", "wave2", "cached8192", "tail16384"), waves):
        if strict:
            require(len(frame) == COUNTS[name], f"{name}: incomplete authentic pixel source")
    joined = pd.concat(waves, ignore_index=True)
    require(joined.photo_id.is_unique,
            "Two acquired cohorts contain the same original photo ID")
    require(joined.classified.eq(joined.morph.isin(HUES)).all(),
            "An unknown ROI was treated as a four-hue flower")
    return joined, waves


def summary(d):
    classified = d.loc[d.classified]
    n = len(d)
    nc = len(classified)
    require(n > 0, "No actually processed photographic pixels")
    taxa = classified.loc[classified.inat_taxon_id.ne("")].copy()
    taxa["genus"] = taxa.species.str.split().str[0]
    named = taxa.groupby("inat_taxon_id", sort=False).agg(
        fraction=("white", "mean"), genus=("genus", "first"))
    gen = named.loc[named.genus.fillna("").ne("")].groupby("genus").fraction.mean()
    technical = (None if d.technical_failure.isna().all()
                 else int(d.technical_failure.fillna(False).sum()))
    white = int(classified.white.sum())
    return {
        "attempted_original_photos": n,
        "four_hue_classified": nc,
        "classifiable_fraction": nc / n,
        "unclassifiable_or_failed": n - nc,
        "acquisition_failures_when_status_recorded": technical,
        "hue_counts": {h: int(classified.morph.eq(h).sum()) for h in HUES},
        "original_unclassified_morph_label_counts": {
            str(k): int(v) for k, v in d.loc[~d.classified, "morph"].value_counts(
                dropna=False).items()},
        "white_fraction_among_classifiable": white / nc if nc else None,
        "taxon_equal_white_fraction": float(named.fraction.mean()) if len(named) else None,
        "genus_equal_white_fraction": float(gen.mean()) if len(gen) else None,
        "nominal_taxa_with_classifiable_hues": len(named),
        "nominal_genera_with_classifiable_hues": len(gen),
        "dated_classifiable_photo_rows": int(classified.date.notna().sum()),
        "dated_classifiable_fraction": (float(classified.date.notna().mean()) if nc else None),
        "white_lower_bound_if_all_unknown_not_white": white / n,
        "white_upper_bound_if_all_unknown_white": (white + n - nc) / n,
    }


def season_support(d, *, min_abs_lat=23.5):
    """North and south analysed separately, including observation-level repeats."""
    sub = d.loc[d.classified & d.date.notna() &
                d.latitude.between(-90, 90) &
                d.longitude.between(-180, 180) &
                d.latitude.abs().ge(min_abs_lat) &
                d.inat_taxon_id.ne("") &
                d.observation_id.ne("")].copy()
    if sub.empty:
        return {"north": {"photo_count": 0}, "south": {"photo_count": 0}}
    month = sub.date.dt.month
    north = sub.latitude.ge(0)
    summer = (north & month.isin((6, 7, 8))) | (~north & month.isin((12, 1, 2)))
    winter = (north & month.isin((12, 1, 2))) | (~north & month.isin((6, 7, 8)))
    sub = sub.loc[summer | winter].copy()
    if sub.empty:
        return {"north": {"photo_count": 0}, "south": {"photo_count": 0}}
    sub["hemisphere"] = np.where(sub.latitude.ge(0), "north", "south")
    month = sub.date.dt.month
    sub["season"] = np.where(
        (sub.latitude.ge(0) & month.isin((6, 7, 8))) |
        (sub.latitude.lt(0) & month.isin((12, 1, 2))), "summer", "winter")
    # Local DJF spans December of the preceding calendar year.
    sub["season_year"] = sub.date.dt.year + month.eq(12).astype(int)
    sub["lat_cell025"] = np.floor((sub.latitude + 90) / .25).astype(int)
    sub["lon_cell025"] = np.floor((sub.longitude + 180) / .25).astype(int)
    # A documented photo set can contain multiple pictures of ONE observation;
    # do not count these as distinct local seasonal replication.
    obs = sub.groupby(["hemisphere", "observation_id", "inat_taxon_id",
                       "season", "season_year", "lat_cell025", "lon_cell025"],
                      as_index=False).agg(white=("white", "mean"),
                                          photo_count=("photo_id", "size"))
    result = {}
    for hemi in ("north", "south"):
        p = sub.loc[sub.hemisphere.eq(hemi)]
        o = obs.loc[obs.hemisphere.eq(hemi)]
        out = {"photo_count": len(p), "distinct_observations": len(o)}
        season_by_photo = {}
        for ss in ("summer", "winter"):
            x = p.loc[p.season.eq(ss)]
            season_by_photo[ss] = {
                "photos": len(x), "white_photo_fraction": (
                    float(x.white.mean()) if len(x) else None),
                "taxa": int(x.inat_taxon_id.nunique()),
                "observations": int(o.season.eq(ss).sum()),
            }
        out["seasons"] = season_by_photo
        if season_by_photo["summer"]["photos"] and season_by_photo["winter"]["photos"]:
            out["raw_summer_minus_winter_white_pp"] = 100 * (
                season_by_photo["summer"]["white_photo_fraction"] -
                season_by_photo["winter"]["white_photo_fraction"])
        else:
            out["raw_summer_minus_winter_white_pp"] = None
        summer_taxa = set(p.loc[p.season.eq("summer"), "inat_taxon_id"])
        winter_taxa = set(p.loc[p.season.eq("winter"), "inat_taxon_id"])
        out["taxa_in_both_seasons"] = len(summer_taxa & winter_taxa)
        out["fraction_summer_photos_from_taxa_absent_in_winter"] = (
            float((~p.loc[p.season.eq("summer"), "inat_taxon_id"]
                   .isin(winter_taxa)).mean()) if summer_taxa else None)
        out["fraction_winter_photos_from_taxa_absent_in_summer"] = (
            float((~p.loc[p.season.eq("winter"), "inat_taxon_id"]
                   .isin(summer_taxa)).mean()) if winter_taxa else None)
        matches = {}
        for needed in (1, 2):
            g = o.groupby(["inat_taxon_id", "lat_cell025", "lon_cell025",
                           "season_year", "season"]).white.agg(["mean", "size"])
            pairs = g.unstack("season")
            if pairs.empty or not {"summer", "winter"} <= set(pairs["size"]):
                matches[str(needed)] = {"taxon_cell_year_pairs": 0, "observations": 0,
                                        "average_pair_summer_minus_winter_white_pp": None}
                continue
            valid = pairs.loc[pairs[("size", "summer")].ge(needed) &
                              pairs[("size", "winter")].ge(needed)].copy()
            matches[str(needed)] = {
                "taxon_cell_year_pairs": len(valid),
                "observations": int((valid[("size", "summer")] +
                                     valid[("size", "winter")]).sum()),
                "average_pair_summer_minus_winter_white_pp": (
                    float(100 * (valid[("mean", "summer")] -
                                 valid[("mean", "winter")]).mean()) if len(valid) else None),
            }
        out["same_taxon_0p25deg_same_season_year_support"] = matches
        result[hemi] = out
    return result


def run(core, broad, wave2, cached8192, output, tail16384=None, strict=True):
    d, waves = assemble(core, broad, wave2, cached8192, tail16384, strict)
    output = Path(output)
    require(not output.exists() or not any(output.iterdir()),
            "Never overwrite an already published source audit")
    output.mkdir(parents=True, exist_ok=True)
    report = {
        "schema": "mois_actual_original_photo_cohort_composition_season_support_v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            name: {"sha256": sha256(path), "file": Path(path).name}
            for name, path in [("core", core), ("broad_cumulative", broad),
                               ("wave2", wave2), ("cached8192", cached8192)] +
             ([("tail16384", tail16384)] if tail16384 is not None else [])
        },
        "actual_original_pixel_photo_union": summary(d),
        "cohort_measured_photo_statistics": {
            arm: summary(frame) for arm, frame in
            zip(("core", "broad", "wave2", "cached8192", "tail16384"), waves)
        },
        "seasonal_source_photo_support_core_only": season_support(waves[0]),
        "seasonal_source_photo_support_all_cohorts": season_support(d),
        "recruitment_note": "Core and taxon/species-first broad/cached arms have different unknown acquisition probabilities; pooled fractions are NOT globally representative.",
        "season_note": "Phenological flower anthesis, observed-year temperature/solar/water, identical individuals and pigment costs are NOT measured here.",
        "unknown_note": "Unknown/failure original-photo ROIs never counted as chromatic.",
        "all_original_source_photo_IDs_unique_after_old_baseline_reconciliation": True,
        "broad_cumulative_313535_original_photo_ROIs_NOT_double_counted": True,
        "tail_included": tail16384 is not None,
    }
    (output / "actual_photo_colour_and_season_support.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True,
                   allow_nan=False) + "\n")
    rows = []
    for key, value in [("all_real_cohorts", report["actual_original_pixel_photo_union"])] + list(
            report["cohort_measured_photo_statistics"].items()):
        rows.append({"arm": key, **{k: v for k, v in value.items()
                                    if not isinstance(v, dict)}})
    pd.DataFrame(rows).to_csv(output / "cohort_white_balance.csv", index=False)
    # This human-readable receipt is extracted from the SAME SHA-verified
    # 334k real-photo calculation and permits full auditing without copying
    # the 10KB released JSON. The full lossless report remains the Release.
    print("MOIS_ACTUAL_RGB_NEW_COHORT_AUDIT", json.dumps({
        "photo_union": len(d),
        "classifiable": report["actual_original_pixel_photo_union"]["four_hue_classified"],
        "tail_included": tail16384 is not None,
        "full_photo_hue_stats": report["actual_original_pixel_photo_union"],
        "by_actual_acquisition_cohort": {
            arm: {"attempted": r["attempted_original_photos"],
                  "four_hue": r["four_hue_classified"],
                  "white_photo_fraction": r["white_fraction_among_classifiable"],
                  "taxa": r["nominal_taxa_with_classifiable_hues"],
                  "genera": r["nominal_genera_with_classifiable_hues"]}
            for arm, r in report["cohort_measured_photo_statistics"].items()},
        "season_by_hemisphere": {
            hemi: {"season_photos": r.get("seasons", {}),
                   "taxa_in_both_seasons": r.get("taxa_in_both_seasons"),
                   "summer_minus_winter_white_pp":
                       r.get("raw_summer_minus_winter_white_pp"),
                   "taxa_season_nonoverlap": {
                       "fraction_summer_photos_from_taxa_absent_in_winter":
                           r.get("fraction_summer_photos_from_taxa_absent_in_winter"),
                       "fraction_winter_photos_from_taxa_absent_in_summer":
                           r.get("fraction_winter_photos_from_taxa_absent_in_summer")},
                   "same_taxon_local_year_observed_support":
                       r.get("same_taxon_0p25deg_same_season_year_support", {})}
            for hemi, r in report["seasonal_source_photo_support_all_cohorts"].items()},
    }, sort_keys=True, allow_nan=False), flush=True)
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ("core", "broad", "wave2", "cached8192", "output"):
        p.add_argument("--" + arg, type=Path, required=True)
    p.add_argument("--tail16384", type=Path)
    args = p.parse_args()
    run(args.core, args.broad, args.wave2, args.cached8192,
        args.output, args.tail16384)


if __name__ == "__main__":
    main()
