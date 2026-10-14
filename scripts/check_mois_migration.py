#!/usr/bin/env python3
"""Verify that mois preserves the original FCP *interspecific* moisture source outputs.

This fast offline check validates original Git blob bytes and receipt arithmetic.
It does NOT re-download photo/environment ZIPs, refit models, recompute bootstrap,
or serve as an independent source or prospective biological confirmation.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SOURCE_SHA={
"docs/FCP_GLOBAL_COLOUR_ENVIRONMENT_DIRECTION_AND_CAUSAL_GATES_20261013.md":"fd03cd7582addc08186e26cfe39730202edb6246",
"docs/FCP_MACRO_HUE_MATCHED_SPATIAL_SUPPORT_20261014.md":"9ea26aca218882ef594144f4098f85e8a4977da7",
"scripts/analysis/fcp_abiotic_hue_direction_macro_20261013.py":"e973abc64b0a834f1ad64412bffd998f0d3bf5e6",
"scripts/analysis/build_fcp_macro_direction_summary_20261013.py":"0e96d6c3543ce6db93d92b60f0d0239095097037",
"scripts/analysis/verify_fcp_macro_direction_20261013.py":"1025ceb8a5fac62519e5e9590ddc107aebb0f2e0",
"scripts/analysis/fcp_macro_spatial_support_decomposition_20261014.py":"5fdbb11a2b19eea2e0f070d848255a842007c50d",
"scripts/analysis/build_fcp_macro_spatial_support_decision_20261014.py":"bae06c182753bb947745901da1b5f236ec24336e",
"scripts/analysis/verify_fcp_macro_spatial_support_20261014.py":"f03f3940f35fb68d56faea6c4c32042243b3fd2f",
"tests/test_fcp_macro_abiotic_direction_20261013.py":"a49140ee65cab44f53de3bd77fedf64d43989017",
"tests/test_spatial_support_decomposition.py":"981af47d8f49011c22955c06c555c6bc87246180",
"results/fcp_macro_direction_20261013/decision_receipt.json":"89f0fd9cc914edcc3d490b1fd629ab08553745f6",
"results/fcp_macro_direction_20261013/independent_verified.json":"72d89d9dd0ea186095c8dd314bb02aae6fb2e594",
"results/fcp_macro_direction_20261013/direction_ledger.csv":"868a8994bfa8db4ef0a9d945bcfe08720975c3ed",
"results/fcp_m1_matched_spatial_support_20261014/decision_receipt.json":"fce7ea6f3a52fa76d08c914d9c53a9a8457ddf6d",
"results/fcp_m1_matched_spatial_support_20261014/independent_verified.json":"0d73a7e580ea5ed598ea31286d7f0bc6e22a8892",
"results/fcp_m1_matched_spatial_support_20261014/matched_spatial_support_differences.csv":"f372fa669097cac75c2b7ada9d68e9f518ef9391",
"results/fcp_m1_matched_spatial_support_20261014/result.json":"60020db92af1c8c64ac61e56f1ecddbaf6140a73",
"results/figure_source_data.json":"2f9c3869ce417e144a894afc94d06c8d7b0bbd6f",
"figures/global_moisture_colour.svg":"88afb9ffab77b3cd78421110b080d0c2d2d9cd4d",
"results/source_receipts/nhm_independent_botanical_colour_summary.json":"fac2c427bb9f652cacc1c4f57ae8610efdf22f72",
"results/source_receipts/photo_colour_missingness_mnar.json":"4f1ccddbff6a43ecc0196d8849d41c4922d95975",
"results/source_receipts/matched_congeneric_vs_withinspecies_rainfall.json":"3a0f7fa2172a931a62318c7c0285459c5284295b",
"results/source_receipts/continuous_visible_white_roi_rainfall.json":"4e7b680ac37301c5b93afcf0d9262c3418660cd1",
"results/source_receipts/global_classified_white_fraction_and_compensation.json":"d1037003f7d586fbf6b6fe58eb511851ba561ebb",
"results/source_receipts/photo_white_latitudinal_null.json":"58c12575f0bbdf318b26c09c0ce56f5c6ca5a909",
"results/source_receipts/matched_photo_reporting_summary.json":"e450b0a4d6ba1ede14214ea160493108f9289d16",
"results/source_receipts/soil_addition_same_species.json":"7af6978ff476dac5a2dceb5f1c9336026dda5634",
"results/source_receipts/unseen_genus_other_geography_500km.json":"a9920c05628842e479abae1a40efcc866c1eb4ba",
"results/source_receipts/nonlinear_geographic_transport.json":"d99c4bee8f51b303af297f56039b3e176dedec45",
"docs/supporting/FOUR_HUE_ENVIRONMENTAL_PREDICTION_DECOMPOSITION.md":"d701856259b5ad9e4731140d086162ad65072f4c",
"docs/supporting/NONLINEAR_GEOGRAPHIC_TRANSPORT_LIMIT.md":"a17dbc0e5623bf2ae290fd9bb40bb0e5f1b56c83",
}

PHOTO_SHA="2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0"
ENV_SHA="6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca"

def source(path):
    return json.loads((ROOT/path).read_text(encoding="utf8"))

def git_blob_sha(b: bytes) -> str:
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()

def main():
    for path,expected in SOURCE_SHA.items():
        file=ROOT/path
        assert file.is_file(),f"Migration incomplete: {path}"
        actual=git_blob_sha(file.read_bytes())
        assert actual==expected,f"Source bytes differ: {path}, {actual} != {expected}"

    macro=source("results/fcp_macro_direction_20261013/decision_receipt.json")
    spatial=source("results/fcp_m1_matched_spatial_support_20261014/decision_receipt.json")
    macro_check=source("results/fcp_macro_direction_20261013/independent_verified.json")
    spatial_check=source("results/fcp_m1_matched_spatial_support_20261014/independent_verified.json")
    assert macro["n_macro_scope_cohorts"]==20
    assert macro["directional_estimates_and_sensitivities"]==1280
    assert macro["source_photo_classified"]==45953
    assert spatial["n_environmental_direction_rows"]==576
    assert macro_check["status"]=="PASS" and spatial_check["status"]=="PASS"
    assert macro_check["direct_source_genus_WLS_colour_coefficients_verified"]==352
    assert spatial_check["independent_direct_WLS_slope_reconstructions"]==288
    assert macro_check["source_SHA256"]["photo"]==PHOTO_SHA
    assert macro_check["source_SHA256"]["environment"]==ENV_SHA
    assert spatial_check["source_photo_sha256"]==PHOTO_SHA
    assert spatial_check["source_abiotic_sha256"]==ENV_SHA
    with (ROOT/"results/fcp_macro_direction_20261013/direction_ledger.csv").open(newline="",encoding="utf8") as f:
        rows=list(csv.DictReader(f))
    assert len(rows)==1280
    with (ROOT/"results/fcp_m1_matched_spatial_support_20261014/matched_spatial_support_differences.csv").open(newline="",encoding="utf8") as f:
        rows2=list(csv.DictReader(f))
    assert len(rows2)==576
    white=[r for r in macro["key_white_directional_observations"] if (r["mode"],r["quality"],r["region"],r["adjustment"],r["axis"],r["hue"])==("climate","all","world","genus_geo","wc_bio12","white")]
    assert len(white)==1
    slope=white[0]
    assert slope["n_photo_species"]==20694 and slope["n_genera"]==1934
    assert abs(slope["point_pp"]-2.254545235123587)<1e-12
    assert slope["family_95_interval_pp"][0]>0
    east=spatial["BIO12_white_by_region_and_weight"]["E90_E180|genus_equal"]
    assert east["n_matched_species"]==2886
    assert not east["difference_excludes_zero_after_24_maxT"]
    ci=east["paired_24_hue_axis_maxT_95_interval_pp"]
    assert ci[0]<0<ci[1]

    matched=source("results/source_receipts/matched_congeneric_vs_withinspecies_rainfall.json")["matched_source_photo_and_species"]
    assert (matched["n_common_original_photos"],matched["n_common_species"])==(21448,6143)
    assert matched["paired_genus_bootstrap95_within"][0]<0<matched["paired_genus_bootstrap95_within"][1]
    mnar=source("results/source_receipts/photo_colour_missingness_mnar.json")
    assert mnar["source_all_unclassifiable"]==54590
    soil=source("results/source_receipts/soil_addition_same_species.json")["soil_addition_on_identical_20786_species"]
    assert soil["CI95"][0]<0<soil["CI95"][1]
    photo_white=source("results/source_receipts/global_classified_white_fraction_and_compensation.json")["world_photo_white"]
    assert photo_white["original_white_photo_count"]==21574
    assert abs(photo_white["original_classifiable_photo_white_fraction"]-21574/45953)<1e-12
    assert not (ROOT/"results/polymorphism_h2_third_cohort_prospective_white_axis_20260917/result.json").exists()

    print("PASS mois: original blob SHA identity",len(SOURCE_SHA),"files;",
          len(rows),"macro estimates;",len(rows2),"matched spatial estimates.")
    print("This does NOT constitute a new original-photo model rerun.")

if __name__=="__main__":
    main()
