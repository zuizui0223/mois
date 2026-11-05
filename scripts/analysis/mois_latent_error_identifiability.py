#!/usr/bin/env python3
"""MOIS latent error audit before any Hui-Walter gold-standard-free model fit.

No labels fabricated, no images downloaded, no inference, no MNAR scenarios.
A 500px YOLO flower BOX/CIELAB value is not an independent binary test for
botanically visible white PETALS. This program does not estimate Se/Sp.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

SOURCE_SHA256 = "2d648c832929f8ba59301906d5ea8d100bdfec6449c6e2828f6b9a910ddd0a02"
SOURCE_N = 131072
SOURCE_STRICT9 = 48784
VALID_T1 = {"original_FCP_ROI9_white"}
VALID_T2 = {"independent_VLM_white", "independent_OWLv2_petal_white",
            "independent_nonshared_white_classifier"}
LABELS = {"white", "nonwhite", "unknown"}
PARAM_MODES = {"shared", "population_specific"}


def require(condition, message):
    if not bool(condition):
        raise ValueError("MOIS_HW_IDENTIFIABILITY: " + message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(4 << 20), b""):
            h.update(block)
    return h.hexdigest()


def identifiability(tests: int, populations: int, accuracy: str) -> dict:
    require(isinstance(tests, int) and tests >= 2, "At least two binary tests required")
    require(isinstance(populations, int) and populations >= 2,
            "At least two different-prevalence populations required for Hui-Walter design")
    require(accuracy in PARAM_MODES, "Unknown accuracy invariance assumption")
    observable = populations * (2**tests - 1)
    parameters = populations + 2*tests * (1 if accuracy == "shared" else populations)
    if observable < parameters:
        status = "UNDERIDENTIFIED"
    elif observable == parameters:
        status = "JUST_IDENTIFIED_COUNT_ONLY"
    else:
        status = "OVERIDENTIFIED_COUNT_ONLY"
    return {
        "tests": tests, "populations": populations,
        "accuracy_across_populations": accuracy,
        "observable_degrees": observable, "free_parameters": parameters,
        "degrees_after_parameter_count": observable-parameters,
        "count_status": status,
        "local_identifiability_proven": False,
        "conditional_independence_verified_by_counts": False,
        "flower_petal_truth_verified_by_counts": False,
        "prevalence_separation_guaranteed_by_strata_names": False,
        "latent_class_label_switching_resolved": False,
        "supports_environment_specific_error_claim": (
            status != "UNDERIDENTIFIED" and accuracy == "population_specific"),
        "not_a_sensitivity_or_specificity_estimate": True,
    }


def paired_count_audit(table: pd.DataFrame) -> dict:
    """Count same-photo pairs ONLY after both are declared binary white-petal tests."""
    fields = {"photo_id", "observation_id", "inat_taxon_id", "source_arm",
              "environment_stratum", "test1_label", "test2_label",
              "test1_method", "test2_method",
              "test1_model_sha256", "test2_model_sha256", "source_photo_sha256"}
    require(fields <= set(table.columns), "Missing joined paired-test identity or provenance")
    d=table.copy()
    for key in ("photo_id","observation_id","inat_taxon_id","source_arm",
                "environment_stratum","test1_method","test2_method",
                "test1_model_sha256","test2_model_sha256","source_photo_sha256"):
        d[key]=d[key].fillna("").astype(str).str.strip()
        require(d[key].ne("").all(), "Missing required source identity or test provenance: "+key)
    require(d.photo_id.is_unique, "duplicated photo ID; do not count as independent photographs")
    require(d.observation_id.is_unique,
            "Repeated observation UUID; use explicit observation-cluster likelihood instead")
    require(d.test1_method.isin(VALID_T1).all() and d.test2_method.isin(VALID_T2).all(),
            "Two tests must target same botanical white petals; YOLO BOX Lab is not a white test")
    require(d["test1_label"].isin(LABELS).all() and d["test2_label"].isin(LABELS).all(),
            "Invalid label; unknown/abstention must not be coerced to nonwhite")
    for method in ("test1_model_sha256","test2_model_sha256","source_photo_sha256"):
        require(d[method].str.fullmatch("[a-f0-9]{64}").all(),
                "Invalid method/image provenance SHA field: "+method)
    require((d.test1_model_sha256 != d.test2_model_sha256).all(),
            "Two tests may not silently reuse the same model weights")
    require(d.environment_stratum.nunique() >= 2,
            "Need two or more prespecified environment populations")
    require(d.source_arm.nunique() == 1,
            "Sampling source-reference and targeted source-arm cannot silently pool")
    # Same source photo record per line; separate model execution/image identity
    # evidence is an additional external provenance requirement.
    complete=d.test1_label.ne("unknown") & d.test2_label.ne("unknown")
    strata=[]
    for stratum,frame in d.groupby("environment_stratum",sort=True):
        p=frame.loc[complete.loc[frame.index]]
        counts={f"{x}{y}":int(((p.test1_label==x)&(p.test2_label==y)).sum())
                for x in ("white","nonwhite") for y in ("white","nonwhite")}
        strata.append({"environment_stratum":str(stratum),"original_photo_ids":len(frame),
                       "n_complete":len(p),"n_unknown_or_abstaining":len(frame)-len(p),
                       "cells":counts,
                       "observed_test1_white_frequency":(
                           float(p.test1_label.eq("white").mean()) if len(p) else None),
                       "observed_test2_white_frequency":(
                           float(p.test2_label.eq("white").mean()) if len(p) else None)})
    return {
        "schema":"mois_same_photo_binary_petal_two_test_quality_gate_v1",
        "original_photo_ids":len(d),
        "same_photo_paired_complete":int(complete.sum()),
        "unknown_abstentions":int((~complete).sum()),
        "source_arm":str(d.source_arm.iloc[0]),
        "strata":strata,
        "conditional_independence_verified":False,
        "strata_prevalence_distinct_proven":False,
        "observed_white_label_rates_can_prove_prevalence":False,
        "genuine_plant_open_flower_or_petal_verified":False,
        "independent_model_execution_and_photo_content_must_be_authenticated_separately":True,
        "latent_white_prevalence_estimated":False,
        "environment_specific_misclassification_estimated":False,
        "unknowns_imputed_as_nonwhite":False,
        "source_observation_clustering_ignored":False,
        "MEE_level_novelty_or_model_correctness_established":False,
    }


def audit_frozen_bulk128(parquet,output):
    """Use a real, SHA-locked numeric source to detect unavailable second test."""
    path=Path(parquet);out=Path(output)
    require(not out.exists(),"Refuse to overwrite the source verification receipt")
    require(sha(path)==SOURCE_SHA256,"Bulk128 actual measured joined source SHA mismatch")
    source=pd.read_parquet(path)
    require(len(source)==SOURCE_N,"Real bulk128 source count differs")
    print("MOIS_HW_SIGNED_SOURCE_COLUMN_NAMES",json.dumps(list(map(str,source.columns)),sort_keys=True),flush=True)
    require({"photo_id","observation_id","strict_original_roi9",
             "screen_status","floral_box_main_L_star","medium_box_proxy_not_strict9",
             "unverified_floral_box_not_petal_ROI"} <= set(source.columns),
            "Original paired same-source photo identity columns changed")
    for key in ("photo_id","observation_id"):
        vals=source[key].dropna().astype(str)
        require(len(vals)==SOURCE_N and vals.is_unique and vals.ne("").all(),
                "Exact "+key+" identity duplicates or missing values")
    strict=source.strict_original_roi9.astype(str).str.lower()
    require(strict.isin(("true","false")).all(),"Strict9 state must be binary success/unknown")
    n_success=int(strict.eq("true").sum())
    require(n_success==SOURCE_STRICT9,"Original strict9 success count drift")
    # Critical distinction: successful white palette measurement !=
    # second independently generated anatomical WHITE label.
    independent_columns={"independent_VLM_white","independent_OWLv2_petal_white",
                         "second_white_classifier_label"}
    require(not (independent_columns & set(source.columns)),
            "Unexpected second white classifier requires new signed method provenance")
    report={
        "schema":"mois_bulk128_original_two_methods_hui_walter_source_eligibility_v1",
        "input_sha256":SOURCE_SHA256,
        "original_distinct_photo_ids":SOURCE_N,
        "original_distinct_observation_ids":SOURCE_N,
        "original_strict9_success":n_success,
        "original_strict9_unknown":SOURCE_N-n_success,
        "medium500_flower_BOX_and_CIELAB_available":True,
        "medium500_is_independent_white_petal_binary_test":False,
        "full_original_ROI9_automated_colour_not_independently_botanically_verified":True,
        "second_independent_white_petal_classifier_outputs_in_frozen_source":False,
        "same_photo_paired_binary_white_tests":0,
        "Hui_Walter_fit_run":False,
        "Se_Sp_environment_specific_estimated":False,
        "corrected_white_moisture_association_estimated":False,
        "original_strict9_unknown_imputed_or_relabelled":False,
        "raw_original_photos_downloaded":0,
        "new_original_photo_ROI_inference":0,
        "next_source_requirement":"same original photo+observation ID, second independent model with white/nonwhite/abstain, SHA-authenticated and source-arm stratified",
        "source_columns":list(map(str,source.columns)),
        "identifiability_examples":[
            identifiability(2,2,"shared"),
            identifiability(2,2,"population_specific"),
            identifiability(3,2,"population_specific")
        ]
    }
    out.mkdir(parents=True)
    (out/"source_eligibility.json").write_text(
        json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("MOIS_HW_SOURCE_STOP_NO_INDEPENDENT_WHITE_CLASSIFIER",json.dumps({
        "photo_ids":SOURCE_N,"strict9_success":n_success,
        "eligible_second_classifier_white_binary":False,
        "Hui_Walter_fit":False},sort_keys=True),flush=True)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bulk-parquet",type=Path)
    p.add_argument("--audit-paired",type=Path)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    require(bool(a.bulk_parquet) != bool(a.audit_paired),
            "Choose exactly one numeric-only audit source")
    if a.bulk_parquet:
        audit_frozen_bulk128(a.bulk_parquet,a.output)
    else:
        path=Path(a.audit_paired)
        report=paired_count_audit(pd.read_csv(path,dtype=str))
        require(not a.output.exists(),"Refuse overwrite")
        a.output.mkdir(parents=True)
        report["input_sha256"]=sha(path)
        (a.output/"paired_two_white_tests_AUDIT_ONLY.json").write_text(
            json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
        print("MOIS_HW_PAIRED_COUNTS_ONLY_NO_LATENT_CLASS_FIT",json.dumps({
            "paired_complete":report["same_photo_paired_complete"],
            "unknown_abstentions":report["unknown_abstentions"]},sort_keys=True))


if __name__=="__main__":
    main()
