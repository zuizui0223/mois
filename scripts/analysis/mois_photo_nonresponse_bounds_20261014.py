#!/usr/bin/env python3
"""Sharp partial-identification of raw wet/dry photographed-white differences.

Uses ONLY archived aggregate counts, no individual photo IDs, no original image
pixels, no new fitted ecological coefficients, and no missing-at-random assumption.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "results/source_receipts/photo_colour_missingness_mnar.json"
OUTPUT = ROOT / "results/mois_photo_nonresponse_bounds_20261014.json"


def build_receipt() -> dict:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    assert source["all_original_unique_photo_ids"] == 100543
    assert (source["source_all_classifiable"], source["source_all_unclassifiable"]) == (45953, 54590)
    bins = source["rain_bin_all_photo_source_counts"]
    assert len(bins) == 10 and sum(row["n"] for row in bins) == 100255

    def decile(index: int) -> dict:
        row = next(b for b in bins if b["rain_decile"] == index)
        classified_share = source["original_photo_wet_dry_unweighted"][str(index)]["white_fraction"]
        classified_white = round(row["n_class"] * classified_share)
        assert math.isclose(classified_white, row["n_class"] * classified_share, abs_tol=1e-7)
        n = row["n"]
        unknown = n - row["n_class"]
        assert n > 0 and unknown >= 0
        return {
            "n_photo_opportunities": n,
            "n_classified": row["n_class"],
            "n_unclassified": unknown,
            "n_classified_white": classified_white,
            "classified_white_fraction": round(classified_share, 8),
            "hypothetical_all_photo_white_share_percent_bounds": [
                round(100 * classified_white / n, 6),
                round(100 * (classified_white + unknown) / n, 6),
            ],
        }

    dry = decile(0)
    wet = decile(9)
    classified_diff = (
        source["original_photo_wet_dry_unweighted"]["9"]["white_fraction"]
        - source["original_photo_wet_dry_unweighted"]["0"]["white_fraction"]
    )
    assert math.isclose(classified_diff, source["original_photo_wet_minus_dry_unweighted"], abs_tol=1e-12)
    lower = (
        wet["hypothetical_all_photo_white_share_percent_bounds"][0]
        - dry["hypothetical_all_photo_white_share_percent_bounds"][1]
    )
    upper = (
        wet["hypothetical_all_photo_white_share_percent_bounds"][1]
        - dry["hypothetical_all_photo_white_share_percent_bounds"][0]
    )
    dry_all_share = (
        dry["n_classified_white"] + 0.5 * dry["n_unclassified"]
    ) / dry["n_photo_opportunities"]
    wet_tipping = (
        dry_all_share * wet["n_photo_opportunities"] - wet["n_classified_white"]
    ) / wet["n_unclassified"]
    original_tipping = source["latent_hue_missing_not_at_random_tipping_curve"]["0.5"]["required_unclassified_wet_white_fraction_for_EQUAL_true_dry_wet_white"]
    assert math.isclose(wet_tipping, original_tipping, abs_tol=1e-12)

    return {
        "schema": "mois_photo_nonresponse_sharp_bounds_v1",
        "role": "Reporting-only hypothetical latent-white sharp bounds; no reclassification, no model refit, no new biological observations",
        "source_receipt": "results/source_receipts/photo_colour_missingness_mnar.json",
        "original_photo_source_sha256": source["source_original_archive_sha256"],
        "n_all_original_photo_opportunities": source["all_original_unique_photo_ids"],
        "n_georeferenced_rain_photo_opportunities": sum(b["n"] for b in bins),
        "comparison": "unweighted original photo opportunity rainfall deciles 0 and 9; classified-only observed white fractions",
        "dry_decile": dry,
        "wet_decile": wet,
        "observed_classified_only_wet_minus_dry_pp": round(100 * classified_diff, 6),
        "hypothetical_all_photo_wet_minus_dry_sharp_bounds_pp": [round(lower, 6), round(upper, 6)],
        "original_mnar_tipping_wet_unknown_white_if_dry_unknown_white_is_0_5": round(wet_tipping, 6),
        "identification_assumption": "Every unclassifiable photo opportunity is hypothetically assigned a latent binary white/nonwhite status; photo without a diagnosable floral organ may have no such well-defined latent trait.",
        "scope_restriction": "Extreme logical bounds, not bootstrap confidence intervals or sensitivity bounds for genus-adjusted BIO12 regression; no adjusted causal or ecological population inference.",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    actions = parser.add_mutually_exclusive_group(required=True)
    actions.add_argument("--write", action="store_true")
    actions.add_argument("--check", action="store_true")
    args = parser.parse_args()
    receipt = build_receipt()
    if args.write:
        OUTPUT.write_text(json.dumps(receipt, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print("Wrote", OUTPUT)
    else:
        assert OUTPUT.exists(), "missing source-derived bounds receipt"
        assert json.loads(OUTPUT.read_text(encoding="utf-8")) == receipt, "partial-identification receipt differs from original source"
        print("PASS original-photo nonresponse bounds: wet-dry sign unidentified")


if __name__ == "__main__":
    main()
