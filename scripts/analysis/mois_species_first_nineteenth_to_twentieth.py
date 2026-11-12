#!/usr/bin/env python3
"""Append genuine nineteenth original image ROI9 and freeze a twentieth 4096-photo source.

The nineteenth consists of 1182 independent original photos for 1182 still-unresolved prior-attempt taxa.
A retry can have existing strict9 but unresolved mixed hue: retain all previous
independent observations and update the mean; never call a retry a new species.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis.mois_species_first_palette9 import PALETTE, COLUMNS, class_labels, sha256
from scripts.analysis.mois_species_first_incremental_first256 import require, parse_boolean
from scripts.analysis.mois_species_first_real1280_cumulative import choose_third_wave

SCHEMA = "mois_species_first_nineteenth4096_append_twentieth4096_v1"
AUDIT_SCHEMA = "mois_species_first_nineteenth_authentic_originalROI9_audit_v1"
SOURCE_ARM = "MOIS_SPECIES_FIRST_NINETEENTH_PRIORITIZED4096"
PRIOR_ATLAS_SHA = "b4b4af6b1597f8e15bb60b818f4c484a9d6988b0ce3ccd711f965e39ef877182"
PRIOR_PHOTOS = 1224828
PRIOR_NUMERIC_PHOTOS = 465729
PRIOR_SUPPORTED_TAXA = 39637
PRIOR_LABELLED_TAXA = 38102
HUES = {"white", "yellow_orange", "red_pink", "blue_purple"}
PIX = ["flower_fraction_" + x for x in PALETTE]


def read_nineteenth(source, workers, audited):
    """Independently recompute authentic nine colour fractions from signed workers."""
    source, workers, audited = map(Path, (source, workers, audited))
    manifest = json.loads((source / "manifest.json").read_text())
    require(
        manifest.get("schema") == "mois_full_source_observation_sampling_cohort_v1"
        and manifest.get("source_arm") == SOURCE_ARM
        and manifest.get("colour_measured") is False
        and manifest.get("source_previous_attempts") == PRIOR_PHOTOS
        and manifest.get("source_prior_no_reused_original_photo_or_observation") is True
        and manifest.get("work_units") == 5
        and manifest.get("worker_ids") == list(range(5))
        and manifest.get("new_photo_attempts_planned") == 1182,
        "Wrong nineteenth 1182-photo frozen manifest",
    )
    signed = json.loads((audited / "receipt.json").read_text())
    require(
        signed.get("schema") == AUDIT_SCHEMA
        and signed.get("actual_independent_original_photo_attempts") == 1182
        and signed.get("photo_observation_source_taxon_ID_collisions") == 0
        and signed.get("frozen_original_manifest_sha256") == sha256(source / "manifest.json")
        and signed.get("audited_derivative_sha256") ==
        sha256(audited / "mois_nineteenth_actual_original_nine_roi.csv.gz")
        and len(signed.get("original_photo_worker_receipt_sha256", {})) == 5,
        "Nineteenth actual original nineROI audit is missing or invalid",
    )
    rows = []
    for i in range(5):
        src = source / "shards" / f"shard-{i:03d}.csv"
        folder = workers / f"mois-nineteenth18-real-original-worker-{i//4}" / f"shard_{i:03d}"
        terminal, receipt = folder / "photos.csv", folder / "receipt.json"
        require(src.is_file() and terminal.is_file() and receipt.is_file(),
                f"Missing actual nineteenth worker source/photo receipt {i}")
        require(sha256(src) == manifest["photo_shard_sha256"][str(i)]
                and sha256(receipt) ==
                signed["original_photo_worker_receipt_sha256"][str(i)],
                f"Nineteenth frozen source or worker receipt SHA drift {i}")
        r = json.loads(receipt.read_text())
        require(
            r.get("schema") == "mois_full_original93m_real_RGB_photo_unit_v1"
            and r.get("shard") == i and r.get("attempted") == len(pd.read_csv(src))
            and r.get("photo_colours_measured") is True
            and r.get("source_sha256") == sha256(src)
            and r.get("photos_sha256") == sha256(terminal),
            f"Nineteenth original ROI numerical worker signature mismatch {i}",
        )
        a = pd.read_csv(src, dtype=str, keep_default_na=False)
        b = pd.read_csv(terminal, dtype=str, keep_default_na=False)
        require(1 <= len(a) == len(b) <= 256 and a.photo_id.is_unique and
                b.photo_id.is_unique and set(a.photo_id) == set(b.photo_id),
                f"Nineteenth source/photo terminal ID mismatch {i}")
        for key in ("observation_id", "inat_taxon_id"):
            require(a.set_index("photo_id").sort_index()[key].equals(
                b.set_index("photo_id").sort_index()[key]),
                "Nineteenth source observation/taxon ID drift: " + key)
        rows.append(b)
    raw = pd.concat(rows, ignore_index=True)
    require(len(raw) == 1182 and raw.photo_id.is_unique and
            raw.observation_id.is_unique and raw.inat_taxon_id.is_unique,
            "Nineteenth original photo/observation/taxon collision")
    audit = pd.read_csv(audited / "mois_nineteenth_actual_original_nine_roi.csv.gz",
                        dtype=str, keep_default_na=False)
    require(len(audit) == 1182 and audit.photo_id.is_unique and
            set(audit.photo_id) == set(raw.photo_id),
            "Nineteenth audit derivative source photo IDs differ")
    x, y = (t.set_index("photo_id").sort_index() for t in (raw, audit))
    for field in ("inat_taxon_id", "observation_id", "morph", "measurement_status", "image_sha256"):
        if field in x.columns and field in y.columns:
            require(x[field].equals(y[field]),
                    "Audited original photo identity/provenance drift " + field)
    rgb = x[PIX].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    strict = (x.morph.isin(HUES).to_numpy(bool)
              & np.isfinite(rgb).all(axis=1)
              & (rgb >= 0).all(axis=1)
              & (rgb <= 1.00001).all(axis=1)
              & (np.abs(np.nan_to_num(rgb).sum(axis=1) - 1) <= 0.001))
    flags = parse_boolean(y.strict9_original_photo_numeric).to_numpy(bool)
    labels = np.full(len(x), "", dtype=object)
    labels[strict] = class_labels(rgb[strict])
    require(np.array_equal(flags, strict)
            and np.array_equal(labels, y.dominant_hue9_label.to_numpy())
            and int(strict.sum()) == signed["actual_genuine_strict9_numeric_success"]
            and int((labels != "").sum()) == signed["actual_clear_nine_palette_hues"]
            and int((~strict).sum()) == signed["actual_technical_unknowns"],
            "Nineteenth audit numeric fractions/hue/technical unknown counters drifted")
    require(y.loc[~flags, PIX].eq("").all().all(),
            "Technical unknown nineteenth photos must carry no fabricated nine-fraction values")
    measured = y.loc[flags, PIX].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    require(np.allclose(rgb[flags], measured, atol=1e-6, rtol=0)
            and np.allclose(measured.sum(axis=1), 1, atol=0.001),
            "Nineteenth original numeric nine fractions do not sum to one / changed in audit")
    out = pd.DataFrame({"photo_id": x.index.to_numpy(),
                        "observation_id": x.observation_id.to_numpy(),
                        "inat_taxon_id": x.inat_taxon_id.to_numpy(),
                        "strict9": strict, "dominant_label9": labels})
    for index, name in enumerate(COLUMNS):
        out[name] = np.where(strict, rgb[:, index], np.nan)
    proof = {
        "original_manifest_sha256": sha256(source / "manifest.json"),
        "original_audit_receipt_sha256": sha256(audited / "receipt.json"),
        "audited_derivative_sha256": sha256(audited / "mois_nineteenth_actual_original_nine_roi.csv.gz"),
        "five_original_worker_receipt_sha256":
            signed["original_photo_worker_receipt_sha256"],
    }
    return out, proof


def append_with_retries(atlas, nineteenth, previous, options, *, quota=4096,
                        expected_previous=PRIOR_PHOTOS):
    """Update first ROIs and independent unresolved-image retries without erasing history."""
    require(len(atlas) == 42111 and atlas.inat_taxon_id.is_unique and
            int(atlas.after_strict9.gt(0).sum()) == PRIOR_SUPPORTED_TAXA and
            int(atlas.dominant_hue9.fillna("").ne("").sum()) == PRIOR_LABELLED_TAXA,
            "Previous ninth 42,111 photo-source atlas totals drifted")
    require(len(previous) == expected_previous and previous.photo_id.is_unique,
            "Prior actual original-photo ID registry incomplete")
    require(len(nineteenth) == 1182 and nineteenth.photo_id.is_unique and
            nineteenth.observation_id.is_unique and nineteenth.inat_taxon_id.is_unique,
            "Nineteenth wave actual original photo/observation/taxon duplications")
    used_photo = set(previous.photo_id.fillna("").astype(str))
    used_obs = set(previous.observation_id.fillna("").astype(str)) - {""}
    require(not used_photo.intersection(set(nineteenth.photo_id.astype(str))) and
            not used_obs.intersection(set(nineteenth.observation_id.astype(str))),
            "Nineteenth real original photo or observation reused an old measurement")
    prev = atlas.copy()
    for field in ("dominant_hue9", "observed_hue9_set", "hue_status"):
        prev[field] = prev[field].fillna("").astype(str)
    prev = prev.set_index("inat_taxon_id")
    ids = nineteenth.inat_taxon_id.astype(str)
    require(ids.isin(prev.index).all(), "Nineteenth real original species missing from fixed source census")
    old = prev.loc[ids].copy()
    require(old.dominant_hue9.eq("").all() and
            old.hue_status.isin(["no_strict9", "mixed_image_unresolved"]).all(),
            "Nineteenth run remeasured a species with an assigned sampled photo hue")
    prev.loc[ids, "after_photo_attempts"] = (
        pd.to_numeric(old.after_photo_attempts).to_numpy(int) + 1)
    good = nineteenth.loc[nineteenth.strict9].copy()
    if not good.empty:
        gids = good.inat_taxon_id.astype(str)
        state = prev.loc[gids].copy()
        old_counts = pd.to_numeric(state.independent_roi9_observations).to_numpy(int)
        require(((old_counts > 0) == (pd.to_numeric(state.after_strict9).to_numpy(int) > 0)).all(),
                "Old strict9 ROI observation counts inconsistent")
        for hue in PALETTE:
            key = "observed_mean_" + hue
            values = pd.to_numeric(good["strict9_" + hue], errors="coerce").to_numpy(float)
            oldmean = pd.to_numeric(state[key], errors="coerce").to_numpy(float)
            require(np.isfinite(values).all() and
                    np.isfinite(oldmean[old_counts > 0]).all(),
                    "Nineteenth authentic nine fractions or prior mean invalid")
            mean = np.where(old_counts > 0,
                            (np.nan_to_num(oldmean) * old_counts + values) / (old_counts + 1),
                            values)
            prev.loc[gids, key] = mean
        prev.loc[gids, "after_strict9"] = 1
        prev.loc[gids, "independent_roi9_observations"] = old_counts + 1
        labels = good.dominant_label9.fillna("").astype(str).to_numpy()
        labelled = labels != ""
        prev.loc[gids, "dominant_hue9"] = labels
        prev.loc[gids, "observed_hue9_set"] = labels
        prev.loc[gids, "hue_status"] = np.where(labelled, "provisional_single",
                                                   "mixed_image_unresolved")
        prev.loc[gids, "distinct_palette_votes"] = (
            pd.to_numeric(state.distinct_palette_votes).to_numpy(int) + labelled.astype(int))
        prev.loc[gids, "ambiguous_observation_rois"] = (
            pd.to_numeric(state.ambiguous_observation_rois).to_numpy(int) +
            (~labelled).astype(int))
        prev.loc[gids, "next_measure_priority"] = np.where(labelled, 2, 1)
    updated = prev.reset_index()
    updated["routine_hue_assigned"] = updated.dominant_hue9.fillna("").ne("")
    require(len(updated) == 42111 and updated.inat_taxon_id.is_unique,
            "42,111 original source taxon atlas conservation lost")
    increased = int((nineteenth.strict9.to_numpy(bool) &
                     old.after_strict9.eq(0).to_numpy(bool)).sum())
    # Old is ordered the same as nineteenth; strict photo success can come from
    # an already numerically supported taxon whose prior hue was inconclusive.
    require(int(updated.after_strict9.gt(0).sum()) == PRIOR_SUPPORTED_TAXA + increased and
            int(updated.dominant_hue9.ne("").sum()) ==
            PRIOR_LABELLED_TAXA + int(nineteenth.dominant_label9.ne("").sum()),
            "Numerically supported / clearly labelled source species totals drifted")
    used = pd.concat([previous[["photo_id", "observation_id", "inat_taxon_id"]],
                      nineteenth[["photo_id", "observation_id", "inat_taxon_id"]]],
                     ignore_index=True)
    require(len(used) == expected_previous + 1182 and used.photo_id.is_unique,
            "Nineteenth appended source inventory not unique")
    # The frozen candidate metadata may be exhausted before 4,096 distinct
    # eligible taxa remain. Never reuse photo/observation IDs or invent rows:
    # acquire the largest complete 256-photo multiple actually available.
    known_photos = set(used.photo_id.fillna("").astype(str))
    known_obs = set(used.observation_id.fillna("").astype(str)) - {""}
    pool = options.loc[~options.photo_id.astype(str).isin(known_photos) &
                       options.observation_id.astype(str).ne("") &
                       ~options.observation_id.astype(str).isin(known_obs)].copy()
    pool = pool.merge(updated[["inat_taxon_id", "hue_status", "dominant_hue9"]],
                      on="inat_taxon_id", how="inner", validate="many_to_one")
    pool = pool.loc[pool.dominant_hue9.eq("") &
                    pool.hue_status.isin(["no_strict9", "mixed_image_unresolved"])]
    available = int(pool.inat_taxon_id.nunique())
    eligible_quota = min(quota, available)
    print("MOIS_TWENTIETH_ELIGIBLE_FROZEN_ORIGINAL_PHOTO_TAXA",
          {"eligible_distinct_taxa": available, "requested": quota,
           "safe_photo_quota": eligible_quota}, flush=True)
    require(eligible_quota > 0,
            "Signed original-photo fallback inventory exhausted; new authentic photo metadata is necessary")
    choices = choose_third_wave(updated, used, options, quota=eligible_quota)
    require(len(choices) == eligible_quota and choices.photo_id.is_unique and
            choices.observation_id.is_unique and choices.inat_taxon_id.is_unique and
            choices.acquisition_priority.isin(["first_unmeasured", "retry_unclear"]).all(),
            "Twentieth original photo source queue is not complete/unique")
    selected = choices.merge(
        updated[["inat_taxon_id", "dominant_hue9"]],
        on="inat_taxon_id", validate="one_to_one", suffixes=("", "_signed"))
    require(selected.dominant_hue9_signed.fillna("").eq("").all() and
            not set(choices.photo_id).intersection(set(used.photo_id)),
            "Twentieth attempted original photo or previously hued taxon reselected")
    choices["measurement_id"] = "MOIS-SPECIES-FIRST-W20-" + choices.photo_id.astype(str)
    choices["image_filename"] = choices.measurement_id + ".jpg"
    choices["shard"] = (np.arange(len(choices)) // 256).astype(str)
    return updated, choices, increased


def run(atlas, source, workers, audited, previous, global_options, output, *, quota=4096):
    output = Path(output)
    require(not output.exists(), "Refuse overwrite nineteenth genuine source outputs")
    require(sha256(atlas) == PRIOR_ATLAS_SHA, "Tenth genuine signed photo atlas SHA drifted")
    oldatlas = pd.read_csv(atlas, compression="infer", dtype={"inat_taxon_id": str})
    photos, proof = read_nineteenth(source, workers, audited)
    oldids = pd.concat([pd.read_parquet(p, columns=[
        "photo_id", "observation_id", "inat_taxon_id"]) for p in previous.values()],
        ignore_index=True)
    options = pd.read_csv(global_options, compression="infer", dtype=str)
    updated, next_queue, fresh = append_with_retries(
        oldatlas, photos, oldids, options, quota=quota)
    output.mkdir(parents=True)
    paths = {"atlas": output / "mois_species_first_42111_after_nineteenth.csv.gz",
             "delta": output / "mois_nineteenth4096_real_original9_photo_delta.parquet",
             "next": output / "mois_twentieth4096_prioritized_original_photo_queue.csv"}
    updated.to_csv(paths["atlas"], index=False, compression="gzip")
    photos.to_parquet(paths["delta"], index=False)
    next_queue.to_csv(paths["next"], index=False)
    report = {
        "schema": SCHEMA,
        "source_taxa_total": 42111,
        "prior_attempts": PRIOR_PHOTOS,
        "prior_cumulative_original_numeric_photo_roi9": PRIOR_NUMERIC_PHOTOS,
        "nineteenth_actual_independent_original_photo_attempts": 1182,
        "nineteenth_actual_strict9_photo_successes": int(photos.strict9.sum()),
        "nineteenth_new_source_taxa_with_first_numeric9": fresh,
        "nineteenth_strict9_on_previous_mixed_unresolved": int(photos.strict9.sum()) - fresh,
        "nineteenth_new_source_taxa_with_clear_dominant_hue9": int(photos.dominant_label9.ne("").sum()),
        "nineteenth_actual_technical_unknowns": int((~photos.strict9).sum()),
        "cumulative_original_photo_attempts": PRIOR_PHOTOS + 1182,
        "cumulative_original_numeric_photo_roi9": PRIOR_NUMERIC_PHOTOS + int(photos.strict9.sum()),
        "source_taxa_with_numerical_nineROI": int(updated.after_strict9.gt(0).sum()),
        "source_taxa_with_clear_sampled_nine_palette": int(updated.dominant_hue9.ne("").sum()),
        "source_taxa_without_nineROI": int(updated.after_strict9.eq(0).sum()),
        "twentieth_source_distinct_new_photo_observation_taxon": len(next_queue),
        "twentieth_source_photo_quota": len(next_queue),
        "twentieth_source_requested_photo_quota": quota,
        "twentieth_source_first_unmeasured": int(next_queue.acquisition_priority.eq("first_unmeasured").sum()),
        "twentieth_source_retry_unclear": int(next_queue.acquisition_priority.eq("retry_unclear").sum()),
        "prior_photo_or_observation_ID_reused": False,
        "nineteenth_original_photo_provenance": proof,
        "output_sha256": {k: sha256(p) for k, p in paths.items()},
        "GBIF_occurrences_downloaded": 0,
        "species_wide_flower_colour_fixed_proven": False,
    }
    (output / "receipt.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print("MOIS_REAL_NINETEENTH4096_APPEND_TWENTIETH4096", json.dumps(report, sort_keys=True), flush=True)
    return report


if __name__ == "__main__":
    a = argparse.ArgumentParser(description=__doc__)
    for k in ("atlas", "source", "workers", "audited", "baseline",
              "first", "second", "third", "fourth", "fifth", "sixth", "seventh", "eighth", "ninth", "tenth", "eleventh", "twelfth", "thirteenth", "fourteenth", "fifteenth", "sixteenth", "seventeenth", "eighteenth",
              "global-options", "output"):
        a.add_argument("--" + k, required=True)
    a.add_argument("--quota", type=int, default=4096)
    v = a.parse_args()
    run(v.atlas, v.source, v.workers, v.audited,
        {k: getattr(v, k) for k in ("baseline", "first", "second",
                                   "third", "fourth", "fifth", "sixth", "seventh", "eighth", "ninth", "tenth", "eleventh", "twelfth", "thirteenth", "fourteenth", "fifteenth", "sixteenth", "seventeenth", "eighteenth")},
        v.global_options, v.output, quota=v.quota)
