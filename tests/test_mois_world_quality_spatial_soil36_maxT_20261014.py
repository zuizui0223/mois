"""Verify all 76 independent original-source global quality/geography/soil bootstrap contrasts.

The original two source ZIPs were independently opened in the working source
recovery environment. Public CI checks saved source-scope results and sensitivity
limits only; no raw photographs/precise original locations are in this repo.
"""
from __future__ import annotations
import csv
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/"results/mois_original_world_strict_geocell_soil36_76maxT_source_audit_20261014.json"
SOURCE=ROOT/"results/fcp_macro_direction_20261013/direction_ledger.csv"
SCRIPT=ROOT/"scripts/analysis/replay_mois_world_quality_geocell_soil36_from_original_zips_20261014.py"


def test_original_world_three_additional_maxT_families_are_source_identical():
    data=json.loads(RECEIPT.read_text(encoding="utf8"))
    with SOURCE.open(newline="",encoding="utf8") as f:
        original=list(csv.DictReader(f))
    assert len(original)==1280
    assert data["n_reviewed_hue_environment_contrasts"]==76
    assert data["status"]=="PASS_ALL_76_FOUR_HUE_COEFFICIENTS_AND_BOOTSTRAP_CI_FAMILIES_MATCH_ORIGINAL"
    assert data["max_abs_point_and_CI_difference_pp"]<1e-9
    expected={
        "strict_photo_climate26":(24,5848,727,8538,727,2026101455),
        "within_genus_geocell_climate26":(24,9370,845,14268,1847,2026102253),
        "global_complete_soil36":(28,16801,1665,31925,1665,2026101397),
    }
    assert {s["scope"] for s in data["scopes"]}==set(expected)
    for scope in data["scopes"]:
        n,sp,g,photos,groups,seed=expected[scope["scope"]]
        assert (len(scope["rows"]),scope["n_species"],scope["n_genera"],scope["n_source_photo_rows"],
                scope["n_genus_subgroups"],scope["original_seed"])==(n,sp,g,photos,groups,seed)
        assert scope["n_original_bootstraps"]==999
        assert math.isfinite(scope["maxT_critical"]) and 2.5<scope["maxT_critical"]<3.5
        assert scope["independent_max_abs_point_and_CI_error_pp"]<1e-10
        for row in scope["rows"]:
            target=[v for v in original if all((
                v["mode"]==scope["source_mode"],v["quality"]==scope["photo_classification"],
                v["region"]=="world",v["adjustment"]==scope["spatial_adjustment"],
                v["weight"]=="genus_equal",v["anchor"]==row["axis"],
                v["hue"]==row["hue"]))]
            assert len(target)==1
            v=target[0]
            for x,y in [
                (row["coefficient_pp_per_source_sd"],v["direction_pp_per_env_global_sd"]),
                (row["bootstrap_95_pp"][0],v["cluster_ci_lo_pp"]),
                (row["bootstrap_95_pp"][1],v["cluster_ci_hi_pp"]),
                (row["within_scope_maxT95_pp"][0],v["maxT_ci_lo_pp"]),
                (row["within_scope_maxT95_pp"][1],v["maxT_ci_hi_pp"])]:
                assert abs(float(x)-float(y))<1e-12
            assert row["maxT_excludes_zero"]==(v["maxT_excludes_zero"]=="True")


def test_original_null_quality_and_locality_checks_not_upgraded_to_climatic_selection():
    data=json.loads(RECEIPT.read_text(encoding="utf8"))
    scopes={s["scope"]:s for s in data["scopes"]}
    for mode in ("strict_photo_climate26","within_genus_geocell_climate26"):
        assert scopes[mode]["source_family_maxT_positive_contrasts"]==[]
        for axis in ("wc_bio12","wc_vapr_annual_mean"):
            row=next(x for x in scopes[mode]["rows"] if x["axis"]==axis and x["hue"]=="white")
            assert row["within_scope_maxT95_pp"][0]<0<row["within_scope_maxT95_pp"][1]
    soil=scopes["global_complete_soil36"]
    assert soil["source_family_maxT_positive_contrasts"]==[
      "wc_bio12:white","wc_vapr_annual_mean:white","soil_available_water_proxy:white"]
    water=next(x for x in soil["rows"] if x["axis"]=="soil_available_water_proxy" and x["hue"]=="white")
    assert 0<water["within_scope_maxT95_pp"][0]<0.01
    assert not data["reports_selection_causality"]
    assert not data["program_wide_multi_hypothesis_selection_corrected"]
    assert data["not_an_independent_biological_sample"] is True
    assert not data["source_photo_image_pixels_remeasured"]
    assert data["original_source_id_photo_site_rights_not_released"]


def test_raw_source_replayer_has_explicit_origins_no_export_and_no_git_inputs():
    code=SCRIPT.read_text(encoding="utf8")
    for token in ("one_to_one","999","np.linalg.solve","sha256","source_rights_and_stable_reviewer_data_deposit_not_completed"):
        assert token in code,token
    assert "raw_inputs/" in (ROOT/".gitignore").read_text(encoding="utf8")
    assert (ROOT/"docs/MOIS_GEB_RESEARCH_ARTICLE_DRAFT_20261014.md").is_file()
    assert (ROOT/"docs/MOIS_GEB_SUPPORTING_INFORMATION_20261014.md").is_file()
