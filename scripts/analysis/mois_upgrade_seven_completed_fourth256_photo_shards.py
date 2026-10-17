#!/usr/bin/env python3
"""Seven completed original RGB/ROI worker shards supersede the prior 512 subset.

Targeted 256 site-season photos remain a separate, unrepresentative
ascertainment arm. Counts are original photographed pixels, NOT flower census.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import tempfile
from pathlib import Path
import pandas as pd
from scripts.analysis import mois_checkpoint_real_fourth_photo_pixel_shards as fourth
from scripts.analysis import mois_reconcile_real_disjoint_fourth512_targeted256 as prev
from scripts.analysis.mois_admit_targeted_256_real_pixels_and_siteyear_support_20261017 import support

SEVEN = tuple(range(7))
NEW_PHOTOS, NEW_FOUR_HUE = 1792, 674
TOTAL = (311487, 150407, 150398)
SOURCE = "mois_311487_actual_original_rgb_photo_hues_plus_targeted256.csv.gz"
RECEIPT = "mois_311487_real_seven_worker_and_targeted256_source_qc.json"
COVERAGE = "mois_311487_season_siteyear_source_support_by_sampling_arm.csv"
ARMS = "mois_311487_original_photo_separate_sampling_arms.csv"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def validate_supersession(old512, all1792):
    """Return only the 1,280 *new* original image pixels from seven workers."""
    require(len(old512) == 512 and len(all1792) == NEW_PHOTOS and
            all1792.photo_id.is_unique, "Incomplete seven-shard image pixels")
    require(all1792.iloc[:512].reset_index(drop=True).equals(
            old512.reset_index(drop=True)), "Earlier512 real pixel rows changed")
    difference = all1792.iloc[512:].reset_index(drop=True)
    require(len(difference) == 1280 and
            int(difference.morph.isin(prev.HUES).sum()) == 471,
            "Seven-worker additional photo classification numbers changed")
    return difference


def source_coverage(old, ordinary, targeted):
    stages = {
        "original309439": old,
        "original_plus_7_fourth_shards": pd.concat((old, ordinary), ignore_index=True),
        "original_plus_targeted256_only": pd.concat((old, targeted), ignore_index=True),
        "all_completed_disjoint_photo_pixels_UNREPRESENTATIVE": prev.join_disjoint(
            old, ordinary, targeted),
    }
    min1 = "same_original_taxon_025degree_cell_same_seasonyear_both_summer_winter_min1_obs"
    min2 = "same_original_taxon_025degree_cell_same_seasonyear_both_summer_winter_min2_obs"
    actual = {k: support(v)[0] for k, v in stages.items()}
    for hemi in ("north_extratropical", "south_extratropical"):
        for label in (min1, min2):
            b = actual["original309439"][hemi][label]
            f = actual["original_plus_7_fourth_shards"][hemi][label]
            t = actual["original_plus_targeted256_only"][hemi][label]
            a = actual["all_completed_disjoint_photo_pixels_UNREPRESENTATIVE"][hemi][label]
            require(b <= f <= a and b <= t <= a,
                    "Adding actually classified photos cannot destroy paired source support")
    rows = []
    for arm, hemispheres in actual.items():
        for hemi, counts in hemispheres.items():
            rows.append({
                "photo_source_arm": arm,
                "hemisphere": hemi,
                "summer_winter_same_taxon_cell_seasonyear_min1_obs": counts[min1],
                "summer_winter_same_taxon_cell_seasonyear_min2_obs": counts[min2],
            })
    return actual, rows


def run(prior, plan_zip, workers, earlier512, targeted256, earlier310207, output):
    require(not output.exists() or not any(output.iterdir()),
            "Never overwrite prior immutable source-authenticated image results")
    old, first512, new256, receipt512, receipt256 = prev.load_real_sources(
        prior, earlier512, targeted256)
    old310207 = earlier310207 / prev.OUT_SOURCE
    oldreceipt = json.loads((earlier310207 / prev.OUT_RECEIPT).read_text())
    require(oldreceipt["actual_combined_photo_source_gzip_SHA256"] == sha(old310207)
            and oldreceipt["individual_original_photo_IDs_disjoint_across_all_arms"] is True,
            "Earlier 310207 genuine pixel/photo-source receipt not verified")
    oldjoined = prev.source_table(old310207)
    require(oldjoined.equals(prev.join_disjoint(old, first512, new256)),
            "Earlier 310207 source pixel rows have changed")
    with tempfile.TemporaryDirectory(prefix="mois_seven_original_rgb_") as td:
        tempout = Path(td) / "actual_seven_shards"
        observed = fourth.run(
            prior / "original309439-source-archive.zip",
            plan_zip, workers, tempout)
        shards = observed["actual_finished_shard_source_receipts"]
        require(tuple(r["shard"] for r in shards) == SEVEN and
                observed["completed_authentic_ROI_shards"] == 7 and
                observed["original_actual_photo_IDs_with_real_pixel_measurement_added"] == NEW_PHOTOS and
                observed["original_four_hue_ROI_source_photo_labels_added"] == NEW_FOUR_HUE and
                observed["original_photo_hue_plus_observation_day_added"] == NEW_FOUR_HUE,
                "Seven completed original photo/model worker receipts are not source-valid")
        completed = prev.source_table(tempout / fourth.RESULT_PHOTOS)
        require(completed.iloc[:len(old)].reset_index(drop=True).equals(old),
                "Old original309439 photo/day/colour prefix changed")
        new1792 = completed.iloc[len(old):].reset_index(drop=True)
        extra1280 = validate_supersession(first512, new1792)
        require(not new1792.photo_id.isin(new256.photo_id).any() and
                new1792.observation_id.is_unique,
                "Targeted source photo/observation overlaps the fourth image wave")
        combined = prev.join_disjoint(old, new1792, new256)
        tally = prev.counted(combined)
        require(tuple(tally[k] for k in (
            "measured_source_photo_IDs", "automatic_four_hue_ROIs",
            "dated_and_automatic_four_hue_ROIs")) == TOTAL,
                "Genuine photo pixel/source date/ROI label count mismatch")
        prior_ids = set(oldjoined.photo_id)
        require(prior_ids.issubset(set(combined.photo_id)) and
                set(combined.photo_id)-prior_ids == set(extra1280.photo_id),
                "Previous 512 checkpoint double-counted or source-only photos included")
        cov, covrows = source_coverage(old, new1792, new256)
        require(cov["original_plus_targeted256_only"] ==
                receipt256["after_named_taxon_cell_seasonyear_support"],
                "Prior targeted seasonal source groups changed")
        output.mkdir(parents=True, exist_ok=True)
        sourcefile = output / SOURCE
        combined.to_csv(sourcefile, index=False, lineterminator="\n",
                        compression={"method": "gzip", "mtime": 0})
        csvcoverage = output / COVERAGE
        pd.DataFrame(covrows).to_csv(csvcoverage, index=False, lineterminator="\n")
        csvarms = output / ARMS
        armrows = []
        for arm, data in (
            ("original_309439_frozen_global_climate_model", old),
            ("fourth4096_FIRST_SEVEN_workers_partial_supersedes_512", new1792),
            ("targeted256_unequal_recruitment_exclude_global_prevalence", new256)):
            armrows.append({"photographic_sampling_arm": arm, **prev.counted(data)})
        pd.DataFrame(armrows).to_csv(csvarms, index=False, lineterminator="\n")
        report = {
            "schema": "mois_actual311487_seven_fourth_worker_source_safely_superseded_v1",
            "genuine_image_workers": list(SEVEN),
            "genuinely_completed_new_original_photo_pixels_in_fourth_wave": NEW_PHOTOS,
            "new_four_hue_orig_photo_ROIs_in_seven_workers": NEW_FOUR_HUE,
            "old512_is_subset_of_seven_workers_not_separate_cohort": True,
            "only_new_1280_photo_pixels_beyond_existing_310207": 1280,
            "only_new_471_four_hue_date_ROIs_beyond_existing_310207": 471,
            "existing_targeted_256_genuine_RGB_photo_ids": 256,
            "existing_targeted_71_original_date_plus_four_hues": 71,
            "overall_actual_photo_pixels_with_disjoint_unique_IDs": TOTAL[0],
            "overall_automatic_four_hue_original_photo_ROIs": TOTAL[1],
            "overall_date_plus_automatic_four_hue_photo_ROIs": TOTAL[2],
            "base309439_original_photo_source_sha256": prev.BASE_SHA,
            "verified_previous310207_photo_source_sha256": sha(old310207),
            "verified_previous512_source_sha256":
                receipt512["new_original_photo_source_full_GZIP_SHA256"],
            "verified_targeted256_source_sha256":
                receipt256["expanded_actual_original_photo_source_SHA256"],
            "complete_seven_worker_original_photo_source_sha256":
                observed["new_original_photo_source_full_GZIP_SHA256"],
            "verified_worker_terminal_receipts": shards,
            "expanded_original_photo_source_sha256": sha(sourcefile),
            "season_siteyear_source_coverage_sha256": sha(csvcoverage),
            "sampling_arm_counts_sha256": sha(csvarms),
            "named_taxon_025degree_same_season_year_photo_support": cov,
            "no_source_only_candidate_counted_as_processed_photo_RGB": True,
            "original_prior_photos_hues_dates_preserved": True,
            "season_enriched_targeted_pixels_not_globally_representative": True,
            "primary_global_six_climate_axes_999genus_maxT_remains_309439_photos": True,
            "not_botanically_verified_true_open_petals_or_individuals": True,
            "no_observed_year_month_weather_or_reproductive_fitness_claim": True,
            "FCP_within_species_photographic_white_darkness_not_reused": True,
            "full_fourth4096_on_success_must_supersede_partial1792": True,
        }
        (output / RECEIPT).write_text(json.dumps(
            report, indent=2, ensure_ascii=False, sort_keys=True) + "\n")
    print("MOIS_REAL_ORIGINAL_311487_PHOTO_SOURCE",
          json.dumps({"photos": TOTAL[0], "four_hues": TOTAL[1],
                      "dated_hues": TOTAL[2],
                      "photo_season_support": cov}, sort_keys=True), flush=True)
    return report


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--prior309439", type=Path, required=True)
    p.add_argument("--photo4096-source-zip", type=Path, required=True)
    p.add_argument("--seven-completed-workers", type=Path, required=True)
    p.add_argument("--earlier512", type=Path, required=True)
    p.add_argument("--targeted256", type=Path, required=True)
    p.add_argument("--earlier310207", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    run(a.prior309439, a.photo4096_source_zip, a.seven_completed_workers,
        a.earlier512, a.targeted256, a.earlier310207, a.output)
