#!/usr/bin/env python3
"""Reconcile two genuinely completed, disjoint, non-representative MOIS ROI arms.

The 512 fourth-wave images are a checkpoint subset of the ongoing 4096 wave;
they MUST be superseded, not added again, when the full wave finishes.
The 256 targeted images were selected to enrich source photo-season support;
their white labels MUST NOT be pooled into the primary global climatic or
unweighted seasonal-prevalence estimands. Open petals are not validated.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from scripts.analysis.mois_admit_targeted_256_real_pixels_and_siteyear_support_20261017 import support
from scripts.analysis.mois_reassemble_original_new4096_measurements import HUES, BASE_COLS

BASE_SHA = "546266e2bbaa340a2523864c28c5a877842eb4b0904a7e0dff90870c6519a47e"
BASE_FILE = "mois_original309439_fullFCP256_actual_photo_hues_and_source_dates.csv.gz"
FOURTH_FILE = "mois_fourth_partial_finished_real_RGB_photo_hues_plus309439.csv.gz"
FOURTH_RECEIPT = "mois_fourth_partial_finished_pixel_source_checkpoint_receipt.json"
TARGET_FILE = "mois_targeted256_actual_ROI_309695_original_photo_and_dates.csv.gz"
TARGET_RECEIPT = "mois_targeted256_original_photo_repeated_season_siteyear_QC_receipt.json"
OUT_SOURCE = "mois_310207_source_disjoint_fourth512_targeted256_real_original_photo_hues.csv.gz"
OUT_RECEIPT = "mois_310207_source_disjoint_768_real_photo_pixel_reconciliation.json"
OUT_SUPPORT = "mois_310207_disjoint_real_photo_season_support.csv"
OUT_ARMS = "mois_310207_real_photo_sampling_arms_NOT_global_prevalence.csv"

BASE_N, BASE_HUES, BASE_DATED = 309439, 149662, 149653
FOURTH_N, FOURTH_HUES = 512, 203
TARGET_N, TARGET_HUES = 256, 71
TOTAL_N, TOTAL_HUES, TOTAL_DATED = 310207, 149936, 149927


def require(condition, explanation):
    if not condition:
        raise ValueError(explanation)


def sha256(path: Path) -> str:
    dig = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            dig.update(chunk)
    return dig.hexdigest()


def source_table(path: Path) -> pd.DataFrame:
    require(path.is_file(), f"Missing actually-measured photo source: {path}")
    d = pd.read_csv(path, compression="gzip", dtype=str, keep_default_na=False,
                    low_memory=False)
    require(tuple(d.columns) ==
            (*BASE_COLS, "qualified_date", "original_source_date_provenance"),
            "Source original photo/schema/date columns differ")
    require(d.photo_id.is_unique and d.photo_id.str.fullmatch(r"[0-9]+").all(),
            "Duplicated or invalid original photographic photo_id")
    require(d.morph.isin((*HUES, "mixed_uncertain")).all(),
            "Source photographic four-hue/unknown model outcomes drift")
    return d


def counted(d: pd.DataFrame) -> dict:
    yes = d.morph.isin(HUES)
    return {
        "measured_source_photo_IDs": len(d),
        "automatic_four_hue_ROIs": int(yes.sum()),
        "dated_and_automatic_four_hue_ROIs":
            int((yes & d.qualified_date.str.fullmatch(r"\d{4}-\d{2}-\d{2}")).sum()),
        "white_hue_ROIs": int(d.morph.eq("white").sum()),
        "technical_hue_unknown_ROIs": int((~yes).sum()),
    }


def original_prefix_and_new_photos(baseline: pd.DataFrame,
                                   source_with_extra: pd.DataFrame,
                                   expected_new: int, expected_classified: int,
                                   arm: str) -> pd.DataFrame:
    n = len(baseline)
    require(len(source_with_extra) == n + expected_new and
            baseline.columns.equals(source_with_extra.columns) and
            source_with_extra.iloc[:n].reset_index(drop=True).equals(baseline),
            f"{arm}: historical photographed photo IDs, hue or day changed")
    new = source_with_extra.iloc[n:].reset_index(drop=True)
    require(new.photo_id.is_unique and
            not new.photo_id.isin(baseline.photo_id).any() and
            new.qualified_date.str.fullmatch(r"\d{4}-\d{2}-\d{2}").all() and
            int(new.morph.isin(HUES).sum()) == expected_classified,
            f"{arm}: unmeasured/duplicate pixels, dates, or classified ROIs")
    return new


def join_disjoint(baseline: pd.DataFrame, fourth_new: pd.DataFrame,
                  targeted_new: pd.DataFrame) -> pd.DataFrame:
    for label, extra in (("ordinary fourth wave", fourth_new),
                         ("season-targeted", targeted_new)):
        require(extra.photo_id.is_unique and
                not extra.photo_id.isin(baseline.photo_id).any() and
                baseline.columns.equals(extra.columns),
                f"{label}: source photo identity, order or schema differs")
    require(not fourth_new.photo_id.isin(targeted_new.photo_id).any(),
            "Fourth and targeted photo arms overlap: cannot sum pixel outcomes")
    merged = pd.concat((baseline, fourth_new, targeted_new), ignore_index=True)
    require(merged.photo_id.is_unique and
            merged.iloc[:len(baseline)].reset_index(drop=True).equals(baseline),
            "Original measured source overwritten during disjoint merge")
    return merged


def load_real_sources(prior: Path, fourth: Path, targeted: Path):
    old_file, fourth_file, target_file = (
        prior / BASE_FILE, fourth / FOURTH_FILE, targeted / TARGET_FILE)
    fourth_proof = json.loads((fourth / FOURTH_RECEIPT).read_text())
    target_proof = json.loads((targeted / TARGET_RECEIPT).read_text())
    require(sha256(old_file) == BASE_SHA,
            "Original 309439 immutable real-photo baseline byte SHA mismatch")
    require(
        fourth_proof["schema"] ==
          "mois_fourth_original_actual_finished_partial_photo_pixel_checkpoint_v1" and
        fourth_proof["source_original_real_309439_photo_SHA256"] == BASE_SHA and
        fourth_proof["new_original_photo_source_full_GZIP_SHA256"] ==
          sha256(fourth_file) and
        fourth_proof["completed_authentic_ROI_shards"] == 2 and
        fourth_proof["original_actual_photo_IDs_with_real_pixel_measurement_added"] ==
          FOURTH_N and
        fourth_proof["original_four_hue_ROI_source_photo_labels_added"] ==
          FOURTH_HUES and
        fourth_proof["original_photo_hue_plus_observation_day_added"] ==
          FOURTH_HUES and
        fourth_proof["previous_photo_ID_hue_taxon_and_source_day_UNCHANGED"] is True and
        fourth_proof["NOT_AN_INDEPENDENT_REPLICATION_AND_WILL_BE_SUPERSEDED_BY_FULL_FOURTH4096"] is True,
        "Fourth512 has no immutable completed real-photo ROI provenance")
    require(
        target_proof["schema"] ==
          "mois_targeted256_actual_original_ROI_photo_season_support_v1" and
        target_proof["source_prior_309439_real_photo_SHA256"] == BASE_SHA and
        target_proof["expanded_actual_original_photo_source_SHA256"] ==
          sha256(target_file) and
        target_proof["number_original_targeted_photo_pixels_actually_measured"] ==
          TARGET_N and
        target_proof["additional_four_hue_photo_ROIs"] == TARGET_HUES and
        target_proof["additional_four_hue_photo_ROIs_with_original_observation_day"] ==
          TARGET_HUES and
        target_proof["target_source_selection_prior_hue_label_not_read"] is True and
        target_proof["original_photo_colour_model_source_attested"] is True and
        target_proof["strict_targeted_sampling_not_a_globally_representative_colour_census"] is True and
        target_proof["NO_original_global_moisture_genus_water_maxT_refit_or_causal_claim"] is True,
        "Targeted256 has no immutable complete real-photo and sampling provenance")
    old, f, t = (source_table(old_file), source_table(fourth_file),
                 source_table(target_file))
    require(counted(old)["measured_source_photo_IDs"] == BASE_N and
            counted(old)["automatic_four_hue_ROIs"] == BASE_HUES and
            counted(old)["dated_and_automatic_four_hue_ROIs"] == BASE_DATED,
            "Previously frozen full measured real-photo source drifted")
    fourth_new = original_prefix_and_new_photos(
        old, f, FOURTH_N, FOURTH_HUES, "fourth512")
    targeted_new = original_prefix_and_new_photos(
        old, t, TARGET_N, TARGET_HUES, "targeted256")
    require(fourth_new.old_mois_or_new_fcp.eq("mois_new_original_ROI").all() and
            targeted_new.old_mois_or_new_fcp.eq(
                "mois_targeted_opposite_photo_season").all(),
            "New original photo sampling-arm identity lost")
    return old, fourth_new, targeted_new, fourth_proof, target_proof


def support_rows(before, after) -> list[dict]:
    rows = []
    first = "same_original_taxon_025degree_cell_same_seasonyear_both_summer_winter_min1_obs"
    second = "same_original_taxon_025degree_cell_same_seasonyear_both_summer_winter_min2_obs"
    for hemi in ("north_extratropical", "south_extratropical"):
        b, a = before[hemi], after[hemi]
        require(a[first] >= b[first] and a[second] >= b[second],
                "Adding genuine photos cannot erase prior matching source observations")
        rows.append({
            "hemisphere": hemi,
            "before_taxon_site_seasonyear_min1": b[first],
            "after_taxon_site_seasonyear_min1": a[first],
            "gain_min1_source_groups": a[first] - b[first],
            "before_taxon_site_seasonyear_min2": b[second],
            "after_taxon_site_seasonyear_min2": a[second],
            "gain_min2_source_groups": a[second] - b[second],
        })
    return rows


def run(prior: Path, fourth: Path, targeted: Path, output: Path):
    require(not output.exists() or not any(output.iterdir()),
            "Do not replace archived source-verified original real-photo result")
    original, f, t, fproof, tproof = load_real_sources(prior, fourth, targeted)
    combined = join_disjoint(original, f, t)
    stats = counted(combined)
    require((stats["measured_source_photo_IDs"],
             stats["automatic_four_hue_ROIs"],
             stats["dated_and_automatic_four_hue_ROIs"]) ==
            (TOTAL_N, TOTAL_HUES, TOTAL_DATED),
            "768 genuinely measured extra photo pixels/date+hue count mismatch")
    old_support, _ = support(original)
    with_fourth_support, _ = support(
        pd.concat((original, f), ignore_index=True))
    with_target_support, _ = support(
        pd.concat((original, t), ignore_index=True))
    full_support, _ = support(combined)
    require(with_fourth_support == old_support,
            "Ordinary first512 checkpoint site/year support unexpectedly changed")
    require(with_target_support == tproof["after_named_taxon_cell_seasonyear_support"],
            "Actual targeted256 precomputed source support cannot be reproduced")
    require(old_support == tproof["before_named_taxon_cell_seasonyear_support"],
            "Signed prior taxon-site/year summer/winter source support differs")
    require(full_support == with_target_support,
            "Disjoint fourth512 unexpectedly affects targeted paired source support")
    rows = support_rows(old_support, full_support)
    require([row["after_taxon_site_seasonyear_min1"] for row in rows] ==
            [87, 54] and
            [row["after_taxon_site_seasonyear_min2"] for row in rows] ==
            [2, 1],
            "Real independent seasonal source support differs from signed run")
    output.mkdir(parents=True, exist_ok=True)
    out = output / OUT_SOURCE
    combined.to_csv(out, index=False, lineterminator="\n",
                    compression={"method": "gzip", "mtime": 0})
    support_file = output / OUT_SUPPORT
    pd.DataFrame(rows).to_csv(support_file, index=False, lineterminator="\n")
    arms = []
    for arm, group in (
        ("historical_and_fullFCP256_unbiased_estimand_frozen", original),
        ("ordinary_fourth512_PARTIAL_do_not_count_again", f),
        ("targeted256_ENRICHED_exclude_from_global_prevalence", t)):
        counts = counted(group)
        arms.append({"source_arm": arm, **counts,
                     "photographic_white_share_among_classified_in_arm_ONLY":
                       counts["white_hue_ROIs"] / counts["automatic_four_hue_ROIs"]})
    arms_file = output / OUT_ARMS
    pd.DataFrame(arms).to_csv(arms_file, index=False, lineterminator="\n")
    receipt = {
        "schema": "mois_310207_completed_disjoint_photo_source_sampling_arm_reconciliation_v1",
        "source_prior309439_real_photo_SHA256": BASE_SHA,
        "source_fourth512_completed_image_pixel_GZIP_SHA256":
            fproof["new_original_photo_source_full_GZIP_SHA256"],
        "source_targeted256_completed_image_pixel_GZIP_SHA256":
            tproof["expanded_actual_original_photo_source_SHA256"],
        "prior_global_climate_and_unpaired_season_photo_denominator_stays_309439": True,
        "only_original_photo_IDs_with_actual_completed_image_pixel_ROI_counted": True,
        "individual_original_photo_IDs_disjoint_across_all_arms": True,
        "latest_disjoint_photo_pixel_source_counts_FOR_PROVENANCE_NOT_GLOBAL_PREVALENCE":
            stats,
        "new_ordinary_completed_photo_pixel_attempts": FOURTH_N,
        "new_ordinary_four_hue": FOURTH_HUES,
        "new_targeted_completed_photo_pixel_attempts": TARGET_N,
        "new_targeted_four_hue": TARGET_HUES,
        "prior_same_taxon_siteyear_support": old_support,
        "with_both_nonoverlapping_source_arms_same_taxon_siteyear_support": full_support,
        "added_matched_taxon_025degree_site_seasonyear_min1":
            sum(row["gain_min1_source_groups"] for row in rows),
        "added_matched_taxon_025degree_site_seasonyear_min2":
            sum(row["gain_min2_source_groups"] for row in rows),
        "actual_combined_photo_source_gzip_SHA256": sha256(out),
        "matched_source_support_figure_csv_SHA256": sha256(support_file),
        "sampling_arm_transparency_csv_SHA256": sha256(arms_file),
        "exclusions_and_noncausal_limits": [
            "Targeted pixels selected via previous season/taxon/site photo support, not independent prevalence sampling",
            "Fourth512 is contained in the running fourth4096 real-pixel study: supersede, never add together",
            "No independently botanically confirmed open flowers or pigment chemistry",
            "Named taxon by 0.25-degree site/year min1 usually has one photographed observation per season",
            "Photos and observation IDs do not imply independent individuals or pollinators",
            "Existing 309439 original WorldClim 1970-2000 climate model is NOT refitted on enriched data",
            "Observed-on season is not measured observed-year monthly weather or fitness",
            "No FCP within-species W/D/50km results are imported",
        ],
    }
    (output / OUT_RECEIPT).write_text(json.dumps(
        receipt, indent=2, ensure_ascii=False, sort_keys=True) + "\n")
    print("MOIS_DISJOINT_ORIGINAL_REAL_PHOTO_SOURCE_RECONCILIATION",
          json.dumps({
              "photo_pixels_measured": stats["measured_source_photo_IDs"],
              "real_four_hue_ROIs": stats["automatic_four_hue_ROIs"],
              "original_source_day_plus_four_hue": stats["dated_and_automatic_four_hue_ROIs"],
              "gain_matched_source_groups_min1":
                  receipt["added_matched_taxon_025degree_site_seasonyear_min1"],
              "gain_matched_source_groups_min2":
                  receipt["added_matched_taxon_025degree_site_seasonyear_min2"],
          }, sort_keys=True), flush=True)
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prior309439", type=Path, required=True)
    parser.add_argument("--fourth512", type=Path, required=True)
    parser.add_argument("--targeted256", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.prior309439, args.fourth512, args.targeted256, args.output)
