#!/usr/bin/env python3
"""Exact completed original-photo census: history + all 75 K1 + 3141 FCP U100.

Recheck published SHA-pinned real-photo Parquet, FCP first-photo candidate
probabilities, outcome-blind source dispositions, signed K1 76,753-photo roster,
and independently published all-75 K1 worldwide union.  Photo IDs, not taxa
or observations, are the deduplication grain.  Never count candidate-only
93.2m IDs, assume true open flowers, or turn technical unknowns into colours.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import tempfile

from scripts.acquisition.mois_completed_pixel_checkpoint import HUES, sha, require
from scripts.acquisition.mois_k1_compact_global_pixel_totals import (
    SCHEMA as K1_GLOBAL_SCHEMA,
)
from scripts.acquisition.mois_fcp_u100_collect_new_first_photo_ROIs import (
    SCHEMA as FCP_ROI_SCHEMA,
    PREFIX_PHOTOS_SHA,
    check_prefix_probabilities,
)
from scripts.acquisition.mois_fcp_u100_depth1_prior_k1_reuse_queue import (
    PREFIX_ARM,
)
from scripts.acquisition.mois_fcp_u100_depth1_preflight import (
    SCHEMA as FCP_PREFLIGHT_SCHEMA,
)

SCHEMA = "mois_global_historical_allK1_FCP3141_exact_original_photo_ID_union_v1"
K1_GLOBAL_SHA = "b75380101cbe7cc9935fbf33dd94a287452c4a08f90ee9865b9bce0c51c6134c"
FCP_ROI_RECEIPT_SHA = "5dbd4590ae96729cd9469fc8403a8ee4528d18e1fb5b928d5a108483ea29150b"
FCP_ROI_PARQUET_SHA = "ee79c60671dda65a3d2f4f9319af952dbe81ab805ce269fe32b3bb3448b0d6b8"
FCP_PREFLIGHT_SHA = "7920cdefb4dd1f5a35def3cdf3fb3c1d19e86e87d4f8db8dab08645e37542eef"
FCP_DISPOSITIONS_SHA = "dacf3a8df021bac2de9f6f02fa779b246e71e544be4e473296f90ad19f66e63f"
K1_ROSTER_SHA = "b620e56cc8386670d183f7a7a74f707cd22de5ccdbb55c1dace2c4a9617a0386"
UPSTREAM=93_209_351
HISTORICAL=929_138
K1=76_753
U100=3_141
ELIGIBLE=20_242
REUSE_OLD=15_483
REUSE_K1=1_618
SOURCES = (
    "prior_original_photo_ROI_reuse",
    "reserved_by_inflight_K1_original_photo",
    "new_original_photo_ROI_queued",
    "prior_independent_technical_photo_holdout",
    "unsupported_original_photo_format_technical_unknown",
)


def assert_sha(path, expected, label):
    require(sha(path)==expected, "Published original photo SHA mismatch: "+label)


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def union(k1_global, fcp_receipt, fcp_parquet, preflight_receipt,
          disposition_csv, prefix_csv, k1_source_roster, output, *,
          k1_global_sha=K1_GLOBAL_SHA,
          fcp_receipt_sha=FCP_ROI_RECEIPT_SHA,
          fcp_parquet_sha=FCP_ROI_PARQUET_SHA,
          preflight_sha=FCP_PREFLIGHT_SHA,
          disposition_sha=FCP_DISPOSITIONS_SHA,
          prefix_sha=PREFIX_PHOTOS_SHA,
          k1_roster_sha=K1_ROSTER_SHA,
          old_n=HISTORICAL,k1_n=K1,new_n=U100,
          eligible=ELIGIBLE,old_reuse=REUSE_OLD,k1_reuse=REUSE_K1,
          k1_hues_expected=444604, fcp_hues_expected=1201,
          k1_groups_expected=75):
    import pandas as pd
    output=Path(output)
    require(not output.exists(), "Refuse overwrite of immutable genuine global ROI union")
    for p,h,name in [
        (k1_global,k1_global_sha,"completed all K1 exact worldwide union"),
        (fcp_receipt,fcp_receipt_sha,"real FCP3141 completed photo receipt"),
        (fcp_parquet,fcp_parquet_sha,"FCP3141 original full-resolution ROI Parquet"),
        (preflight_receipt,preflight_sha,"FCP 20242 first-photo/K1 disjoint preflight"),
        (disposition_csv,disposition_sha,"FCP per-photo disposition roster"),
        (prefix_csv,prefix_sha,"FCP 20242-photo 1/N/M source roster"),
        (k1_source_roster,k1_roster_sha,"all frozen new K1 photo IDs"),
    ]: assert_sha(p,h,name)
    g=read_json(k1_global)
    f=read_json(fcp_receipt)
    p=read_json(preflight_receipt)
    require(g.get("schema")==K1_GLOBAL_SCHEMA
            and g.get("signed_historical_true_distinct_original_photo_IDs")==old_n
            and g.get("K1_original_new_photo_attempts")==k1_n
            and g.get("genuine_K1_completed_and_SHA_sealed_groups")==k1_groups_expected
            and g.get("source_K1_all_precolour_original_photo_IDs_disjoint_from_historical")==k1_n
            and g.get("zero_all_source_K1_vs_historical_photo_ID_overlap") is True
            and g.get("counts_reproduce_independent_full_929138_ID_union") is True
            and g.get("exact_worldwide_historical_plus_completed_K1_distinct_original_photo_IDs")==old_n+k1_n
            and g.get("exact_worldwide_historical_plus_completed_K1_four_hue_ROIs")==k1_hues_expected
            and g.get("additional_raw_image_downloads_or_ROI_inference_during_this_update")==0,
            "K1 all 75 actual photo groups/historical identity not independently signed")
    require(f.get("schema")==FCP_ROI_SCHEMA
            and f.get("original_93m_photo_ID_source_population")==UPSTREAM
            and f.get("eligible_global_first_FCP_nominal_taxon_IDs")==eligible
            and f.get("first_original_photo_candidates_retained")==eligible
            and f.get("previously_measured_original_image_ROIs_reused_elsewhere")==old_reuse
            and f.get("inflight_K1_original_photo_IDs_reserved_elsewhere")==k1_reuse
            and f.get("newly_true_image_measured_original_photo_attempts")==new_n
            and f.get("individual_photo_design_probabilities_1_over_N_over_M_preserved") is True
            and f.get("photo_colour_inferred_by_unchanged_frozen_FCP_ROI_v4") is True
            and f.get("technical_unknowns_kept_unclassified") is True
            and f.get("individual_botanical_native_wild_anthesis_verified") is False
            and f.get("globally_unbiased_white_flora_fraction_estimated") is False
            and f.get("actual_new_image_parquet_SHA256")==fcp_parquet_sha
            and f.get("new_images_downloaded_OR_ROI_inferred_by_this_data_only_collector")==0
            and f.get("previous_real_photo_pixel_Parquets_rewritten") is False,
            "FCP genuine new original-photo source receipt/probability/unknown status altered")
    require(p.get("schema")==FCP_PREFLIGHT_SCHEMA
            and p.get("eligible_nominal_FCP_U100_source_taxa")==eligible
            and p.get("already_measured_old_photo_IDs_to_reuse_NOT_RECOMPUTED")==old_reuse
            and p.get("already_running_K1_new_photo_IDs_to_reuse_NOT_RECOMPUTED")==k1_reuse
            and p.get("never_before_measured_FCP_U100_original_photo_ROI_target_IDs")==new_n
            and p.get("new_original_photos_overlap_with_K1")==0
            and p.get("source_original_observations_one_per_taxon")==eligible
            and p.get("all_photo_inclusion_probabilities_individually_source_checked") is True
            and p.get("new_original_image_bytes_downloaded_here")==0
            and p.get("flower_colour_models_run_here")==0
            and p.get("ready_for_inference_ONLY_AFTER_K1_75_WORKERS_AND_75_SIGNED_GROUP_RELEASES") is True,
            "Signed old/K1/new photo-ID source preflight receipt invalid")
    require(old_reuse+k1_reuse+new_n==eligible,
            "Global FCP first-photo outcome-blind selection denominator lost")
    prefix=check_prefix_probabilities(prefix_csv,prefix_sha=prefix_sha,
                                      expected_eligible=eligible)
    with gzip.open(k1_source_roster,"rt",encoding="utf-8",newline="") as infile:
        reader=csv.DictReader(infile)
        k1_ids=set()
        for row in reader:
            pid=row["photo_id"]
            require(pid.isdecimal() and pid not in k1_ids,"Original K1 precolour photo IDs repeated")
            k1_ids.add(pid)
        require(len(k1_ids)==k1_n,"Frozen actual K1 original IDs missing")
    chosen={}
    kinds={key:0 for key in SOURCES}
    with Path(disposition_csv).open(encoding="utf-8",newline="") as handle:
        reader=csv.DictReader(handle)
        require({"photo_id","observation_id","inat_taxon_id",
                 "depth1_joint_original_photo_inclusion_probability",
                 "reuse_or_new_photo_outcome","source_arm"}<=set(reader.fieldnames or []),
                "FCP per-photo old/K1/new evidence grain missing")
        taxon_ids=set()
        observation_ids=set()
        for row in reader:
            pid,obs,taxon=(row[x] for x in ("photo_id","observation_id","inat_taxon_id"))
            action=row["reuse_or_new_photo_outcome"]
            require(pid in prefix and pid not in chosen and action in kinds
                    and taxon not in taxon_ids and obs not in observation_ids
                    and row["source_arm"]==PREFIX_ARM,
                    "Repeated FCP taxon/observation/photo or unknown source action")
            old_obs,old_taxon,prob=prefix[pid]
            require(obs==old_obs and taxon==old_taxon
                    and math.isclose(float(row["depth1_joint_original_photo_inclusion_probability"]),prob,rel_tol=1e-12),
                    "Photo-specific FCP first-photo selection probability/UUID altered")
            require((pid in k1_ids)==(action=="reserved_by_inflight_K1_original_photo"),
                    "K1 already-selected original photo would be recoloured/double-counted")
            chosen[pid]=(obs,taxon,prob,action)
            taxon_ids.add(taxon)
            observation_ids.add(obs)
            kinds[action]+=1
    require(len(chosen)==eligible and kinds["prior_original_photo_ROI_reuse"]==old_reuse
            and kinds["reserved_by_inflight_K1_original_photo"]==k1_reuse
            and kinds["new_original_photo_ROI_queued"]==new_n
            and not kinds["prior_independent_technical_photo_holdout"]
            and not kinds["unsupported_original_photo_format_technical_unknown"],
            "FCP first-photo 15483 old/1618 K1/3141 genuinely new photo partition inconsistent")
    df=pd.read_parquet(fcp_parquet)
    mandatory={"photo_id","observation_id","inat_taxon_id","morph",
               "measurement_status","image_sha256","selection_probability","source_arm"}
    require(mandatory<=set(df.columns) and len(df)==new_n,
            "FCP real ROI terminal source columns or 3141 original photos missing")
    photos=df["photo_id"].fillna("").astype(str)
    obs=df["observation_id"].fillna("").astype(str)
    taxa=df["inat_taxon_id"].fillna("").astype(str)
    hues=df["morph"].fillna("").astype(str)
    states=df["measurement_status"].fillna("").astype(str)
    imagehash=df["image_sha256"].fillna("").astype(str)
    probabilities=pd.to_numeric(df["selection_probability"],errors="coerce")
    known=hues.isin(HUES)
    failures=states.eq("image_acquisition_failed")
    require(photos.str.fullmatch(r"[0-9]+").all() and photos.is_unique
            and obs.is_unique and taxa.is_unique
            and not set(photos)&k1_ids
            and set(photos)=={pid for pid,entry in chosen.items()
                              if entry[3]=="new_original_photo_ROI_queued"}
            and df["source_arm"].fillna("").eq(PREFIX_ARM).all()
            and hues.isin((*HUES,"mixed_uncertain")).all()
            and states.eq("classified_four_state_morph").eq(known).all()
            and ((imagehash.str.fullmatch("[0-9a-f]{64}"))|failures).all()
            and probabilities.notna().all() and probabilities.gt(0).all()
            and probabilities.le(1).all(),
            "FCP new measured photo IDs/terminal image SHA/unknown hue inconsistent")
    for pid,oo,tt,pp in zip(photos,obs,taxa,probabilities):
        original=chosen[pid]
        require(oo==original[0] and tt==original[1]
                and math.isclose(float(pp),original[2],rel_tol=1e-12),
                "FCP photo source original observation/taxon/selection probability changed")
    classified=int(known.sum())
    failed=int(failures.sum())
    uncertain=int((~known & ~failures).sum())
    require(classified+failed+uncertain==new_n
            and (classified,failed,uncertain)==(
                f["newly_true_image_FCP_four_hue_classified"],
                f["newly_true_image_acquisition_failed"],
                f["newly_true_image_unclassifiable_acquired"])
            and classified==fcp_hues_expected
            and all(int(hues.eq(h).sum())==f["new_original_photo_hue_counts"][h]
                    for h in (*HUES,"mixed_uncertain"))
            and int(hues.eq("mixed_uncertain").sum())==failed+uncertain,
            "FCP SHA-signed original technical nonresponse and automatic four-colour counts disagree")
    source_total=g["exact_worldwide_historical_plus_completed_K1_distinct_original_photo_IDs"]
    source_classified=g["exact_worldwide_historical_plus_completed_K1_four_hue_ROIs"]
    report={
        "schema":SCHEMA,
        "full_93m_candidate_photo_ids_metadata_only":UPSTREAM,
        "global_prior_original_photo_ID_attempts_from_all75_K1":source_total,
        "global_prior_original_four_hue_ROIs_from_all75_K1":source_classified,
        "additional_exact_FCP_U100_new_original_photo_IDs":new_n,
        "additional_original_photo_automated_four_hue_ROIs":classified,
        "additional_FCP_original_acquisition_failed":failed,
        "additional_FCP_acquired_but_four_hue_unclassifiable":uncertain,
        "exact_global_unique_original_photo_ID_attempts":source_total+new_n,
        "exact_global_automatic_four_hue_photo_ROIs":source_classified+classified,
        "eligible_U100_nominal_taxa_source_only_not_all_measured":eligible,
        "FCP_old_exact_original_photo_IDs_reused_without_new_inference":old_reuse,
        "FCP_K1_exact_photo_IDs_reused_without_new_inference":k1_reuse,
        "new_FCP_images_overlapping_historical_or_K1":0,
        "FCP_new_photo_original_observation_IDs_unique_within_arm":True,
        "K1_vs_FCP_photo_IDs_proved_disjoint_from_signed_full_source":True,
        "crossarm_independent_observation_IDs_not_imputed":True,
        "source_K1_complete_worldwide_union_SHA256":k1_global_sha,
        "source_FCP_measured_original_parquet_SHA256":fcp_parquet_sha,
        "source_FCP_actual_pixel_receipt_SHA256":fcp_receipt_sha,
        "source_FCP_20242_K1_disjoint_preflight_SHA256":preflight_sha,
        "source_FCP_disposition_roster_SHA256":disposition_sha,
        "source_K1_76753_original_photo_ID_roster_SHA256":k1_roster_sha,
        "source_original_first_photo_1_over_N_over_M_SHA256":prefix_sha,
        "previous_raw_image_downloads_or_image_inference_repeated_for_this_union":0,
        "true_wild_native_open_blooming_flower_colour_proportions_not_verified":True,
        "globally_unbiased_white_flora_rate_not_estimated":True,
    }
    output.parent.mkdir(parents=True,exist_ok=True)
    stage=Path(tempfile.mkdtemp(prefix="mois-real-K1-FCP-total-",dir=output.parent))
    try:
        (stage/"mois_K1_FCP_U100_source_audited_global_original_photo_union.json").write_text(
            json.dumps(report,sort_keys=True,indent=2)+"\n",encoding="utf-8")
        os.replace(stage,output)
    except BaseException:
        shutil.rmtree(stage,ignore_errors=True)
        raise
    print("MOIS_REAL_FCP_K1_WORLDWIDE_UNIQUE_ORIGINAL_PHOTO_PIXEL_ROIS",json.dumps({
        "genuine_photo_attempts":report["exact_global_unique_original_photo_ID_attempts"],
        "four_hue_ROIs":report["exact_global_automatic_four_hue_photo_ROIs"],
        "new_FCP_originals":new_n,"new_FCP_classifiable":classified,
        "reused_historical_and_K1_FCP_photo_IDs":old_reuse+k1_reuse,
        "repeated_photo_fetches":0,"new_photo_ID_overlaps":0,
    },sort_keys=True),flush=True)
    return report


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ("k1-global","fcp-receipt","fcp-parquet","preflight",
                 "dispositions","prefix","k1-roster","output"):
        ap.add_argument("--"+name,required=True,type=Path)
    a=ap.parse_args()
    union(a.k1_global,a.fcp_receipt,a.fcp_parquet,a.preflight,
          a.dispositions,a.prefix,a.k1_roster,a.output)


if __name__=="__main__":
    main()
