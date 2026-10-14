"""Audit independent original-ZIP reconstruction of all 576 macro spatial comparisons.

The complete source-table replay ran in an isolated working environment. Public
CI validates its aggregate receipt against the already-frozen source decision
and the exact original 576-result file; original source ZIPs are not in Git.
"""
from __future__ import annotations
import csv
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
INDEPENDENT=ROOT/"results/mois_original_matched_spatial_full576_replay_20261014.json"
SOURCE=ROOT/"results/fcp_m1_matched_spatial_support_20261014/decision_receipt.json"
SOURCE_CSV=ROOT/"results/fcp_m1_matched_spatial_support_20261014/matched_spatial_support_differences.csv"
REPLAYER=ROOT/"scripts/analysis/replay_mois_original_matched_spatial_576_20261014.py"


def load(path):
    return json.loads(path.read_text(encoding="utf8"))


def test_independent_real_source_spatial_replay_numerical_scope():
    x=load(INDEPENDENT)
    assert x["status"]=="PASS_ORIGINAL_576_POINTS_PERCENTILE_AND_MAXT_INTERVALS"
    assert x["original_source_ledger_rows"]==576
    assert x["bootstrap_resamples"]==999
    assert x["n_axes"]==6 and x["n_hues"]==4
    assert x["comparison_fields"]=={
        "points":576,"percentile_interval_endpoints":1152,"maxT_interval_endpoints":576}
    assert x["max_abs_difference_pp"]<1e-9
    assert x["photo_sha256"]=="2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0"
    assert x["environment_sha256"]=="6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca"
    assert len(x["main_results"])==4


def test_world_east_original_frozen_interval_and_taxon_equality():
    source=load(SOURCE)["BIO12_white_by_region_and_weight"]
    replay=load(INDEPENDENT)["main_results"]
    assert set(source)==set(replay)
    with SOURCE_CSV.open(newline="",encoding="utf8") as f:
        rows=list(csv.DictReader(f))
    assert len(rows)==576
    for key,original in source.items():
        x=replay[key]
        assert x["n_matched_species"]==original["n_matched_species"]
        assert x["n_matched_genera"]==original["n_matched_genera"]
        assert x["n_matched_photos"]==original["matched_photographic_individuals"]
        assert abs(x["stage_B_broad_white_BIO12_pp"]-original["same_species_with_broad_coordinate_pp"])<1e-12
        assert abs(x["stage_C_cell_white_BIO12_pp"]-original["same_species_within_genus_cell_pp"])<1e-12
        assert abs(x["B_minus_C_white_BIO12_pp"]-original["paired_adjustment_difference_pp"])<1e-12
        a=x["B_minus_C_white_BIO12_maxT95_pp"]
        b=original["paired_24_hue_axis_maxT_95_interval_pp"]
        assert all(abs(y-z)<1e-12 for y,z in zip(a,b))
        assert a[0]<0<a[1]
        assert original["difference_excludes_zero_after_24_maxT"] is False
        scope,weight=key.split("|")
        chosen=[q for q in rows if (q["region"],q["weight"],q["type"],
                q["comparison"],q["anchor"],q["hue"])==
                (scope,weight,"paired_difference",
                "spatial_adjustment_on_same_species","wc_bio12","white")]
        assert len(chosen)==1
        assert abs(float(chosen[0]["point_pp"])-x["B_minus_C_white_BIO12_pp"])<1e-12


def test_scientific_boundaries_and_private_recovery_still_hold():
    x=load(INDEPENDENT)
    assert x["interpretation"]["same_species_B_C_identical_photographic_hue_and_environment"]
    assert not x["interpretation"]["world_genus_equal_spatial_adjustment_difference_maxT95_excludes_zero"]
    assert not x["interpretation"]["eastern_genus_equal_spatial_adjustment_difference_maxT95_excludes_zero"]
    assert not x["interpretation"]["matched_geography_is_causal_intervention"]
    assert not x["interpretation"]["programme_wide_post_outcome_multiple_testing_controlled"]
    assert not x["interpretation"]["independent_biological_replication"]
    assert not x["interpretation"]["original_image_pixels_remeasured"]
    assert not x["interpretation"]["original_taxon_photo_ids_or_coordinates_exported"]
    assert x["interpretation"]["photo_rights_and_geoprivacy_review_required"]
    code=REPLAYER.read_text(encoding="utf8")
    assert "one_to_one" in code and "999" in code and "np.linalg.solve" in code
    assert "sha256" in code.lower() and "no_photo_or_taxon_locations_exported" in code
    assert "raw_inputs/" in (ROOT/".gitignore").read_text(encoding="utf8")
