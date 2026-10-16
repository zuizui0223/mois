#!/usr/bin/env python3
"""Fail closed when FCP and mois primary analysis levels are accidentally mixed.

Runs offline on committed frozen receipts. It does not refit models or inspect
unavailable historical photo ZIPs. The original receipt files are read-only.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load(relative: str) -> dict:
    return json.loads((ROOT / relative).read_text(encoding="utf-8"))


def require(ok: bool, message: str) -> None:
    if not ok:
        raise AssertionError(message)


def main() -> None:
    scope = load("study_scope.json")
    study = scope["study"]
    require(study in ("fcp", "mois"), "Unknown scientific study")
    require(scope["repository"] == f"zuizui0223/{study}", "Wrong repository identity")
    require(scope["shared_upstream"]["photo_id_intersection_status"] == "not_measured",
            "Claimed cross-paper photo-ID overlap requires an actual ID-level audit")
    require(scope["shared_upstream"]["independent_source_replication"] is False,
            "Shared iNaturalist photos are not independent source replication")
    require(scope["shared_upstream"]["can_pool_model_estimates"] is False,
            "Within- and between-species effects cannot be pooled as one estimand")
    require((ROOT / scope["active_manuscript"]).is_file(), "Active manuscript missing")
    for path in scope["primary_result_paths"]:
        require((ROOT / path).is_file(), f"Frozen result missing: {path}")

    if study == "fcp":
        require(scope["estimand_level"] == "within_named_species_photographic_colour_structure",
                "FCP must remain a within-species inference")
        frame = scope["fixed_analysis_frame"]
        require((frame["photographs"], frame["named_species"]) == (149900, 1499),
                "FCP primary high-depth frame changed without a new study protocol")
        require(frame["cohorts"] == {"discovery": 500, "validation": 500, "prospective": 499},
                "FCP cohort composition changed")
        require(frame["is_global_polymorphism_prevalence_denominator"] is False,
                "An opportunity sample is not global polymorphism prevalence")

        h2 = load("results/polymorphism_h2_third_cohort_prospective_white_axis_20260917/result.json")
        require(h2["stage_receipt"]["rows"] == 49900, "Third-cohort photo rows changed")
        p = h2["thresholds"]["primary_0_10"]
        require(p["species"] == 158 and p["structured_null_upper_p"] == 0.001,
                "Prospective W claim must remain 158-species pre-frozen test")
        require(h2["decision"]["verdict"] == "H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED",
                "Original H2 frozen verdict changed")

        geo = load("results/polymorphism_distributed_polymorphism_posthoc_20261007/result.json")
        require(geo["primary_radius_km"] == 50, "Primary geographic scale changed")
        require(geo["confirmatory_decisions_changed"] is False, "Post hoc is not confirmatory")
        require([geo["primary"][c]["n_species"] for c in ("discovery", "validation", "third")]
                == [166, 181, 204], "FCP within-species geographic support changed")
        require(all(geo["primary"][c]["mean_local_depletion"] > 0
                    for c in ("discovery", "validation", "third")),
                "Original species-conditioned local depletion signs changed")

        foreign = [
            "results/fcp_macro_direction_20261013",
            "results/fcp_m1_matched_spatial_support_20261014",
            "scripts/analysis/fcp_abiotic_hue_direction_macro_20261013.py",
            "docs/MOIS_GLOBAL_MOISTURE_MANUSCRIPT_20261014.md",
        ]
        print("PASS fcp: 1,499-species within-species D/local structure/W; not mois rainfall")
    else:
        require(scope["estimand_level"] == "between_named_species_global_photographed_hue_environment",
                "mois must remain a between-species global environment analysis")
        frame = scope["fixed_analysis_frame"]
        require((frame["unique_original_photos"], frame["named_taxa"],
                 frame["colour_classified_photos"], frame["unclassified_photos"])
                == (100543, 42111, 45953, 54590), "mois original photo frame changed")
        require((frame["global_climate26_comparison_species"], frame["global_climate26_genera"])
                == (20694, 1934), "mois congeneric primary comparison changed")

        macro = load("results/fcp_macro_direction_20261013/decision_receipt.json")
        require((macro["original_photo_opportunities"], macro["source_photo_classified"],
                 macro["directional_estimates_and_sensitivities"]) == (100543, 45953, 1280),
                "mois primary original-source receipt changed")
        key = ("climate", "all", "world", "genus_geo", "wc_bio12", "white")
        matched = [r for r in macro["key_white_directional_observations"]
                   if tuple(r[k] for k in ("mode", "quality", "region",
                                          "adjustment", "axis", "hue")) == key]
        require(len(matched) == 1, "Global genus-centred white/BIO12 slope not unique")
        rain = matched[0]
        require((rain["n_photo_species"], rain["n_genera"]) == (20694, 1934),
                "mois global model must be between named species")
        require(abs(rain["point_pp"] - 2.254545235123587) < 1e-10,
                "Original congeneric rain coefficient changed")

        strict_key = ("climate", "strict", "world", "genus_geo", "wc_bio12", "white")
        strict = [r for r in macro["key_white_directional_observations"]
                  if tuple(r[k] for k in ("mode", "quality", "region",
                                         "adjustment", "axis", "hue")) == strict_key]
        require(len(strict) == 1, "Strict-photo sensitivity missing")
        require(strict[0]["family_95_interval_pp"][0] < 0 < strict[0]["family_95_interval_pp"][1],
                "Strict-photo negative uncertainty gate lost")
        spatial = load("results/fcp_m1_matched_spatial_support_20261014/decision_receipt.json")
        east = spatial["BIO12_white_by_region_and_weight"]["E90_E180|genus_equal"]
        require(east["difference_excludes_zero_after_24_maxT"] is False,
                "Spatially matched nonconfirmation cannot be omitted")
        require(east["paired_24_hue_axis_maxT_95_interval_pp"][0] < 0
                < east["paired_24_hue_axis_maxT_95_interval_pp"][1],
                "Spatial matched simultaneous interval must span zero")
        foreign = [
            "results/polymorphism_h2_third_cohort_prospective_white_axis_20260917",
            "results/polymorphism_distributed_polymorphism_posthoc_20261007",
            "docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md",
        ]
        print("PASS mois: original 20,694-species congeneric moisture slopes and negative gates; not FCP W")

    for path in foreign:
        require(not (ROOT / path).exists(), f"Cross-study primary output leaked: {path}")
    print("Boundary checked from frozen receipts only; original photographs were not reanalysed.")


if __name__ == "__main__":
    main()
