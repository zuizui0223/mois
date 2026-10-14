"""Defend the source-verified *within-scope* 24-hue maxT replay.

The original 999 bootstrap worlds and point estimates were already reconstructed
from checksum-frozen original measured source input tables outside CI. CI can
verify this numerical receipt and its scopes, but cannot rerun the unpublished
anonymous genus moment NPZ or original full photo ZIPs.
"""
from __future__ import annotations
import csv, json, math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
JSON_FILE=ROOT/"results/mois_world24_genus_moment_source_replay_audit_20261014.json"
LEDGER=ROOT/"results/fcp_macro_direction_20261013/direction_ledger.csv"
SCRIPT=ROOT/"scripts/analysis/replay_mois_world24_from_genus_moments_20261014.py"


def load():
    return json.loads(JSON_FILE.read_text(encoding="utf8"))


def test_source_24_world_maxT_ci_family_uses_original_results_exactly():
    r=load()
    assert r["n_tests"]==24
    assert r["bootstrap_draws"]==999 and r["bootstrap_seed"]==2026101430
    assert (r["source_genus_anonymous_group_count"],r["original_eligible_species"])==(1934,20694)
    assert math.isclose(r["maxT_critical_value"],3.0109108050453743,abs_tol=1e-11)
    assert r["maximum_absolute_point_and_ci_error_pp"]<1e-10
    assert r["all_24_points_and_percentile_ci_and_maxT_ci_match_original"]
    assert r["all_original_24_maxT_significance_flags_match"]
    with LEDGER.open(newline="",encoding="utf8") as f:
        rows=list(csv.DictReader(f))
    assert len(rows)==1280
    all24=r["source_rows"]
    assert len(all24)==24
    for q in all24:
        selected=[v for v in rows if all([
          v["mode"]=="climate",v["quality"]=="all",v["region"]=="world",
          v["adjustment"]=="genus_geo",v["weight"]=="genus_equal",
          v["anchor"]==q["axis"],v["hue"]==q["hue"]])]
        assert len(selected)==1,(q["axis"],q["hue"])
        orig=selected[0]
        assert int(orig["n_species"])==20694
        assert abs(q["point_pp"]-float(orig["direction_pp_per_env_global_sd"]))<1e-12
        for ours,names in ((q["nominal95_pp"],("cluster_ci_lo_pp","cluster_ci_hi_pp")),
                           (q["maxT95_pp"],("maxT_ci_lo_pp","maxT_ci_hi_pp"))):
            assert all(abs(a-float(orig[name]))<1e-12 for a,name in zip(ours,names))
        assert q["maxT_excludes_zero"]==(orig["maxT_excludes_zero"]=="True")
    assert sum(x["maxT_excludes_zero"] for x in all24)==2


def test_bootstrap_replay_is_not_studywide_selection_or_public_morph_archive():
    r=load()
    assert not r["original_source_programme_wide_selection_adjustment"]
    assert not r["independent_biological_sample_replication"]
    assert r["computed_from_original_measured_photo_environment_tables"]
    assert not r["original_image_pixels_used"]
    assert not r["anonymous_genus_sufficient_moments_available_in_public_repo"]
    assert not r["anonymous_moments_approved_for_public_redistribution"]
    assert not r["source_photo_original_IDs_or_coordinates_exported"]
    assert r["original_license_privacy_review_required"] and not r["permanent_archive_or_review_link_exists"]
    assert not (ROOT/"results/mois_anonymised_genus_sufficient_stats_INTERNAL_20261014.npz").exists()
    code=SCRIPT.read_text(encoding="utf8")
    assert "np.load" in code and "999" in code and "multinomial" in code
    assert "photo IDs" in code or "photo ID" in code
