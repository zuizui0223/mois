"""Report-only guard: does a single named genus drive the two main water slopes?

1943? No: original source correctly contains 1,934 qualifying nominal genera.
The restricted genus-moment NPZ is NOT a published or CI input. CI checks the
aggregate-only executed sensitivity receipt against original frozen coefficients.
"""
from __future__ import annotations
import csv
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
REPORT=ROOT/"results/mois_world_primary_moisture_single_genus_influence_20261014.json"
LEDGER=ROOT/"results/fcp_macro_direction_20261013/direction_ledger.csv"
SCRIPT=ROOT/"scripts/analysis/loo_mois_world_genus_moisture_20261014.py"


def report():
    return json.loads(REPORT.read_text(encoding="utf8"))


def test_primary_moisture_eight_original_slopes_and_all_exclusions_preserve_sign():
    z=report()
    assert z["status"]=="PASS_1934_INDEPENDENT_SINGLE_GENUS_DELETIONS_ON_SOURCE_LOCKED_TWO_MOISTURE_AXES"
    assert (z["source_original_measurement_denominator_species"],
            z["source_original_nominal_genera"])==(20694,1934)
    assert z["n_total_leave_one_genus_out_slope_values"]==15472
    assert z["source_private_genus_sufficient_moments_sha256"]=="77dff4acc5756610e463bec3e5877c711f1695d8984be0b675e86f7d8afb4aae"
    with LEDGER.open(newline="",encoding="utf8") as f:
        rows=list(csv.DictReader(f))
    assert len(rows)==1280
    outcomes=z["axis_hue_results"]
    assert set(outcomes)=={"wc_bio12","wc_vapr_annual_mean"}
    for axis,views in outcomes.items():
        assert set(views)=={"white","yellow_orange","red_pink","blue_purple"}
        for hue,v in views.items():
            matching=[r for r in rows if (
                r["mode"],r["quality"],r["region"],r["adjustment"],
                r["weight"],r["anchor"],r["hue"])==
                ("climate","all","world","genus_geo","genus_equal",axis,hue)]
            assert len(matching)==1
            ref=float(matching[0]["direction_pp_per_env_global_sd"])
            assert math.isclose(v["all_genera_original_point_pp"],ref,abs_tol=1e-11)
            assert v["leave_one_genus_out_min_pp"]<=ref<=v["leave_one_genus_out_max_pp"]
            assert v["point_sign_changes_over_1934_genus_exclusions"]==0
            assert v["leave_one_genus_out_min_pp"]*v["leave_one_genus_out_max_pp"]>0
            maxdelta=max(abs(v["leave_one_genus_out_min_pp"]-ref),
                         abs(v["leave_one_genus_out_max_pp"]-ref))
            assert math.isclose(maxdelta,v["largest_absolute_point_change_pp"],abs_tol=1e-10)
            assert 0<=v["p95_absolute_point_change_pp"]<=v["largest_absolute_point_change_pp"]
    assert outcomes["wc_bio12"]["white"]["leave_one_genus_out_min_pp"]>2.17
    assert outcomes["wc_vapr_annual_mean"]["white"]["leave_one_genus_out_min_pp"]>3.84


def test_exploratory_influence_diagnostic_never_claims_confidence_or_causality():
    z=report()
    assert z["new_statistical_significance_test"] is False
    assert z["confidence_intervals_reestimated"] is False
    assert z["all_eight_original_full_model_points_match_frozen_source_ledger"] is True
    assert z["all_eight_slope_signs_unchanged_in_all_single_genus_deletions"] is True
    assert z["licence_and_geoprivacy_review_required_before_releasing_source_moments"] is True
    assert z["private_moments_uploaded_to_public_repository"] is False
    assert z["no_source_photo_ids_species_genus_names_or_site_locations_exported"] is True
    for token in ("geography","strict-photo","phylogeny","MNAR"):
        assert token in z["interpretation"],token
    assert not (ROOT/"results/mois_anonymised_genus_sufficient_stats_INTERNAL_20261014.npz").exists()
    s=SCRIPT.read_text(encoding="utf8")
    assert "PRIVATE_MOMENT_SHA256" in s
    assert "np.linalg.solve" in s and "--private-moments" in s
