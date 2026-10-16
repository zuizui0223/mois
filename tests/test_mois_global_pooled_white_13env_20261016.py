"""Source-safe evidence gate: pooled floral-white versus 13 environment predictors.

Original exact measured-photo/site ZIPs remain held separately under attribution
and location privacy controls. Public CI checks frozen source-recomputed
aggregate facts and synthetic linear algebra, NEVER pretending to replay ZIPs.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
RESULT=ROOT/"results/mois_global_overall_white_13env_20261016.json"
SOURCE=ROOT/"scripts/analysis/mois_global_pooled_white_13env_20261016.py"
PREVIOUS=ROOT/"results/mois_env_shapley_and_genus_residual_transfer_20261015.json"
REPORT=ROOT/"docs/MOIS_GLOBAL_OVERALL_WHITE_13ENV_AND_GENUS_SCALES_20261016.md"


def close(x,y,tol=1e-10):
    assert math.isclose(x,y,abs_tol=tol,rel_tol=0),(x,y)


def test_original_cohort_and_distinct_global_targets():
    a=json.loads(RESULT.read_text(encoding="utf8"))
    assert a["photo_source_archive_SHA256"]=="2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0"
    assert a["environment_source_archive_SHA256"]=="6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca"
    assert a["original_source_cohort"]=={
        "eligible_classifiable_photo_rows":26163,
        "eligible_species":16801,
        "eligible_genera":1665,
        "soil36_complete_photo_rows_before_genus_gate":31925,
    }
    b=a["predictor_blocks_13_distinct"]
    assert {x:len(k) for x,k in b.items()}=={
      "Geography":3,"Moisture":2,"Temperature and radiation":2,
      "Elevation":1,"Soil5":5
    }
    assert len({f for col in b.values() for f in col})==13
    g=a["pooled_overall_white_genus_equal"]
    s=a["pooled_overall_white_species_equal"]
    close(g["pooled_all13_R2"],.01999939742187049)
    close(s["pooled_all13_R2"],.012823012026220198)
    close(g["same_cohort_genus_identity_only_R2"],.44315326849102415)
    close(s["same_cohort_genus_identity_only_R2"],.34258997766250865)
    close(g["same_cohort_within_genus_all13_R2_of_remaining_within_white_variance"],.006202366064066545)
    close(s["same_cohort_within_genus_all13_R2_of_remaining_within_white_variance"],.003399461515192022)
    close(g["same_cohort_genus_identity_plus_all13_R2"],.4466070357614218)
    close(g["added_env_after_genus_R2_of_total_white_variance"],(0.4466070357614218-0.44315326849102415))
    close(g["same_cohort_genus_identity_plus_all13_R2"],
          g["same_cohort_genus_identity_only_R2"] +
          (1-g["same_cohort_genus_identity_only_R2"])*g["same_cohort_within_genus_all13_R2_of_remaining_within_white_variance"])
    for target in (g,s):
        close(sum(target["pooled_block_shapley_fraction_total_white_variance"].values()),
              target["pooled_all13_R2"])
        assert target["pooled_block_unique_partial_R2"]["Soil5"]>=0
    full=a["all_soil36_classifiable_species_without_congeneric_gate_species_equal"]
    assert (full["eligible_species"],full["observed_named_genera"],full["source_classifiable_photo_rows"])==(20786,4853,31925)
    close(full["all13_in_sample_pooled_white_R2"],.012806646762199388)
    assert abs(full["all13_in_sample_pooled_white_R2"]-s["pooled_all13_R2"])<.0001
    old=json.loads(PREVIOUS.read_text(encoding="utf8"))
    close(g["same_cohort_within_genus_all13_R2_of_remaining_within_white_variance"],
          old["shapley_explanatory_in_sample"]["in_sample_total_R2"])
    assert not any(t in RESULT.read_text(encoding="utf8") for t in ('inat_taxon_id','photo_id','decimalLatitude','gps_coordinates'))


def test_pooled_block_attribution_and_uncertainty():
    a=json.loads(RESULT.read_text(encoding="utf8"))
    g=a["pooled_overall_white_genus_equal"]
    shares=g["pooled_block_shapley_fraction_total_white_variance"]
    for name,expected in {"Geography":.006490063358,"Moisture":.005574291189,
                          "Soil5":.004021893324,"Elevation":.002372890403,
                          "Temperature and radiation":.001540259148}.items():
        close(shares[name],expected)
    boot=a["genus_equal_999_cluster_bootstrap"]
    assert (boot["draws"],boot["seed"])==(999,202610161210)
    assert boot["pooled_all13_R2_bootstrap_percentile95"][0]<g["pooled_all13_R2"]<boot["pooled_all13_R2_bootstrap_percentile95"][1]
    ci=boot["all13_coefficients"]
    assert len(ci)==13 and boot["simultaneous_13_conditional_coefficients_critical_maxT95"]>2
    close(ci["wc_bio12"]["point_pp_per_original_soil36_photo_site_SD"],3.7700812278013878)
    assert ci["wc_bio12"]["family_13_variable_maxT95"][0]>0
    assert ci["wc_vapr_annual_mean"]["family_13_variable_maxT95"][0]<0<ci["wc_vapr_annual_mean"]["family_13_variable_maxT95"][1]
    assert ci["wc_elevation_m"]["family_13_variable_maxT95"][0]<0<ci["wc_elevation_m"]["family_13_variable_maxT95"][1]
    assert ci["soil_available_water_proxy"]["family_13_variable_maxT95"][0]<0<ci["soil_available_water_proxy"]["family_13_variable_maxT95"][1]
    assert ci["lon_cos"]["family_13_variable_maxT95"][1]<0
    null=a["fixed_genus_size_species_white_permutation_null"]
    assert (null["draws"],null["seed"])==(1999,202610161214)
    close(null["genus_equal"]["mean"],.19063009721271282)
    close(null["species_equal"]["mean"],.09902315308989651)
    assert null["genus_equal"]["one_sided_p_plus_one"]==.0005
    assert null["species_equal"]["one_sided_p_plus_one"]==.0005
    assert null["genus_equal"]["central_95_envelope"][1]<g["same_cohort_genus_identity_only_R2"]


def test_source_replayer_and_ecological_scope_notes():
    script=SOURCE.read_text(encoding="utf8")
    text=REPORT.read_text(encoding="utf8")
    for word in ("--self-test","--photo-zip","--env-zip",
                 "source.original","block_shapley","genus_shuffle",
                 "within-genus", "Original source"):
        assert word in script or word in text,word
    assert "2.000%" in text and "1.282%" in text and "0.6202%" in text
    assert "20,786" in text and "44.315%" in text
    assert "not a genetic" in text.lower() or "not a global population" in text.lower() or "not a genetic" in script.lower()
