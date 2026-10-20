#!/usr/bin/env python3
"""Unify existing real FCP447500 photo pixels and first genuine 93m-source RGB wave.

No new photo download, morphology inference, geo, native-status or ecological
analysis. Retain first 350399 original MOIS RGB outcomes and then FCP's signed
447500 original pixels. Where the newly measured 93m p0002 observation sample
shares a photo ID with the FCP archive, preserve the existing FCP pixel result
and record disagreement, including source photo labels. All source attempts
and technical unknowns remain auditable even if an ID appears in both arms.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis.mois_actual_photo_cohort_season_audit_20261019 import (
    HUES, require, sha256,
)
from scripts.analysis.mois_append_full93m_first_original_RGB_photo_wave_20261019 import (
    COLUMNS, OUTPUT_PHOTOS as FIRST_PARQUET,
    OUTPUT_REPORT as FIRST_RECEIPT,
)
from scripts.acquisition.mois_reuse_complete_FCP447500_real_original_ROI_pixels import (
    OUT_ALL as FCP_PARQUET,
    OUT_AUDIT as FCP_RECEIPT,
)

FCP_UNION=594081
FIRST_UNION=368861
BASE=350399
FIRST_NEW=18462
TOTAL_FCP_CLASSIFIED=288747
FIRST_CLASSIFIED=171379
RESULT_PARQUET="mois_source_verified_FCP447500_and_93m_p0002_combined_unique_photo_ROIs.parquet"
RESULT_RECEIPT="mois_original_photo_RGB_FCP447500_plus_first93m_id_union_receipt.json"
RESULT_COUNTS="mois_original_photo_RGB_all_reused_source_arm_hue_counts.csv"


def verified_sources(fcp_parquet,fcp_receipt,first_parquet,first_receipt):
    fcp_report=json.loads(Path(fcp_receipt).read_text(encoding="utf-8"))
    first_report=json.loads(Path(first_receipt).read_text(encoding="utf-8"))
    require(
        fcp_report.get("schema")==
        "mois_import_source_signed_FCP447500_true_original_ROIs_on350399_v1"
        and fcp_report.get("union_total_genuinely_original_photo_IDs")==FCP_UNION
        and fcp_report.get("union_total_four_hue_classified")==TOTAL_FCP_CLASSIFIED
        and sha256(fcp_parquet)==fcp_report["output_SHA256"][FCP_PARQUET],
        "FCP447500 original measured photo ID union or immutable bytes changed"
    )
    require(
        first_report.get("schema")==
        "mois_genuine_full93m_first_probability_photo_RGB_append_to_350399_v1"
        and first_report.get("total_actual_original_photo_attempts")==FIRST_UNION
        and first_report.get("new_actual_original_photo_attempts")==FIRST_NEW
        and first_report.get("total_actual_four_hue_classified")==FIRST_CLASSIFIED
        and first_report.get("photo_ID_duplicate_count")==0
        and first_report.get("old_photo_IDs_altered") is False
        and sha256(first_parquet)==first_report["output_sha256"][FIRST_PARQUET],
        "First full93m actual ROI source/immutable original pixel SHA is not valid"
    )
    fcp=pd.read_parquet(fcp_parquet)
    first=pd.read_parquet(first_parquet)
    return fcp,first,fcp_report,first_report


def _check_pixels(frame,label):
    required={"photo_id","inat_taxon_id","species","observation_id","morph",
              "classified","white","arm","measurement_status","technical_failure"}
    require(required<=set(frame),
            f"{label}: original ROI missing source columns {required-set(frame)}")
    ids=frame.photo_id.astype(str)
    require(ids.str.fullmatch(r"\d+").all() and ids.is_unique,
            f"{label}: invalid or duplicate original photo ID")
    require(frame.classified.eq(frame.morph.isin(HUES)).all(),
            f"{label}: old original classifier four-hue state changed")
    require(frame.white.eq(frame.morph.eq("white")).all(),
            f"{label}: white photo ROI colour is not the original classified hue")
    require(frame.loc[~frame.classified,"white"].eq(False).all(),
            f"{label}: technically unknown photographs were called white")
    return ids


def append(fcp,first,*,expected_fcp=FCP_UNION,expected_first=FIRST_UNION,
           expected_base=BASE,expected_first_new=FIRST_NEW):
    fcp_ids=_check_pixels(fcp,"FCP source-authenticated original pixels")
    first_ids=_check_pixels(first,"first 93m source-authenticated original pixels")
    require(len(fcp)==expected_fcp and len(first)==expected_first
            and len(fcp)-expected_base==243682,
            "FCP and first original sources not at frozen 350399 base")
    require(
        fcp_ids.iloc[:expected_base].tolist()==first_ids.iloc[:expected_base].tolist(),
        "The 350399 old original photo IDs were not retained in both frozen source releases"
    )
    for col in ("inat_taxon_id","morph","classified","white"):
        q=fcp[col].iloc[:expected_base].fillna("").astype(str).reset_index(drop=True)
        r=first[col].iloc[:expected_base].fillna("").astype(str).reset_index(drop=True)
        require(q.equals(r),
                "Original 350399 source image pixels or accepted taxon changed: "+col)
    first_new=first.iloc[expected_base:].copy()
    require(len(first_new)==expected_first_new
            and first_new.arm.eq("global93m_probability_observation").all(),
            "New observation-probability original-photo append is not a genuine distinct source arm")
    require(not first_ids.iloc[expected_base:].isin(
                set(first_ids.iloc[:expected_base])).any(),
            "Previously measured original MOIS pixels re-occur in new probability source")
    old_ids=set(fcp_ids)
    first_new_ids=first_new.photo_id.astype(str)
    overlaps=first_new_ids.isin(old_ids)
    common=first_new.loc[overlaps,["photo_id","inat_taxon_id","morph",
                                   "classified","white"]].set_index("photo_id")
    compare=common.join(
        fcp.loc[fcp.photo_id.isin(common.index),
                ["photo_id","inat_taxon_id","morph","classified","white"]].
            set_index("photo_id"),lsuffix="_first93m",rsuffix="_FCP",
            how="inner",validate="one_to_one")
    require(len(compare)==int(overlaps.sum()),
            "Existing FCP colour and new original source photo ID overlap not one-to-one")
    taxon_disagreement=int((compare.inat_taxon_id_first93m.astype(str)!=
                            compare.inat_taxon_id_FCP.astype(str)).sum())
    both_classified=compare.classified_first93m & compare.classified_FCP
    hue_disagreement=int((both_classified &
                          compare.morph_first93m.ne(compare.morph_FCP)).sum())
    confidence_disagreement=int(compare.classified_first93m.ne(
        compare.classified_FCP).sum())

    genuinely_added=first_new.loc[~overlaps].copy()
    # Preserve the FCP union's existing full source provenance schema, with
    # new first93m fields mapped only where genuinely absent in FCP.
    for key in fcp.columns:
        if key not in genuinely_added:
            if key=="photo_colour_provenance":
                genuinely_added[key]="actual_original_93m_observation_p0002_ROI_v4"
            elif key=="source_original_fcp_morph":
                genuinely_added[key]=""
            elif key=="fcp_source_cohort":
                genuinely_added[key]=""
            elif key=="source_fcp_global_classifiable":
                genuinely_added[key]=False
            else:
                genuinely_added[key]=pd.NA
    new_columns=[z for z in genuinely_added.columns if z not in fcp]
    for key in new_columns:
        fcp[key]=pd.NA
    genuinely_added=genuinely_added[fcp.columns]
    output=pd.concat((fcp,genuinely_added),ignore_index=True,sort=False)
    _check_pixels(output,"combined FCP and independent 93m observations")
    require(len(output)==expected_fcp+len(genuinely_added)
            and output.photo_id.is_unique
            and output.iloc[:expected_fcp].reset_index(drop=True).equals(
                fcp.reset_index(drop=True)),
            "Existing 594081 genuinely measured photo original ROI pixels mutated")
    require(int(output.classified.sum())==
            int(fcp.classified.sum())+int(genuinely_added.classified.sum()),
            "Four-hue counts diverged from preserved original pixel source")
    previously_classifiable=set(fcp.loc[fcp.classified,"inat_taxon_id"].astype(str))
    newly_classifiable=set(genuinely_added.loc[
        genuinely_added.classified,"inat_taxon_id"].astype(str))
    report={
        "schema":"mois_original_photo_ROI_FCP447500_plus_first_full93m_p0002_ID_union_v1",
        "preexisting_FCP350399_source_union_original_photo_attempts":len(fcp),
        "first_full93m_p0002_actual_photo_attempts":len(first_new),
        "first_full93m_p0002_shared_with_FCP_completed_pixel_source":int(overlaps.sum()),
        "additional_original_photo_IDs_not_previously_in_FCP_or_MOIS":len(genuinely_added),
        "all_genuinely_attempted_unique_original_photo_IDs":len(output),
        "all_original_photo_four_hue_classified":int(output.classified.sum()),
        "additional_first93m_four_hue_classified":int(genuinely_added.classified.sum()),
        "first93m_photo_ID_overlap_accepted_taxon_disagreements":taxon_disagreement,
        "first93m_photo_ID_overlap_both_classified_hue_disagreements":hue_disagreement,
        "first93m_photo_ID_overlap_classifiability_disagreements":confidence_disagreement,
        "first93m_photo_IDs_all_accounted_for":True,
        "new_first93m_classified_nominal_taxa_not_in_FCP_or_old_MOIS":len(
            newly_classifiable-previously_classifiable),
        "all_old_FCP_and_MOIS_photo_pixels_preserved_immutable":True,
        "duplicate_source_photo_IDs_in_final_union":0,
        "original_color_model_unmodified":True,
        "original_licensed_photo_images_downloaded_for_this_union":False,
        "meta_candidate_photo_IDs_not_counted_as_measured":True,
        "new_native_status_or_environmental_ecology_analyses_run":False,
        "source_observation_probability_first_wave":0.0002,
        "photo_cohort_provenance":"FCP original pixels immutable, then first probability 93m original photo IDs absent in FCP",
        "unknown_image_ROIs_in_final_union":int((~output.classified).sum()),
    }
    return output,report


def run(fcp_parquet,fcp_receipt,first_parquet,first_receipt,out):
    out=Path(out)
    require(not out.exists() or not any(out.iterdir()),
            "Previously published genuine original photo ID union is immutable")
    fcp,first,fcp_report,first_report=verified_sources(
        fcp_parquet,fcp_receipt,first_parquet,first_receipt)
    union,report=append(fcp,first)
    out.mkdir(parents=True,exist_ok=True)
    union.to_parquet(out/RESULT_PARQUET,compression="zstd",index=False)
    coverage=[]
    for name,frame in union.groupby("arm",observed=True,sort=True):
        classified=frame.loc[frame.classified]
        coverage.append({
            "original_photo_source_arm":name,
            "genuine_source_photo_attempts":len(frame),
            "actual_original_four_hue_ROIs":len(classified),
            "unknown_or_technical_failed_photo_attempts":len(frame)-len(classified),
            "white":int(classified.morph.eq("white").sum()),
            "yellow_orange":int(classified.morph.eq("yellow_orange").sum()),
            "red_pink":int(classified.morph.eq("red_pink").sum()),
            "blue_purple":int(classified.morph.eq("blue_purple").sum()),
            "named_taxa_with_actual_four_hue":int(
                classified.inat_taxon_id.astype(str).loc[
                    classified.inat_taxon_id.astype(str).ne("")].nunique()),
            "unique_source_observations":int(frame.observation_id.astype(str).nunique())
        })
    pd.DataFrame(coverage).to_csv(out/RESULT_COUNTS,index=False)
    report["input_original_pixel_Sha256"]={
        FCP_PARQUET:sha256(fcp_parquet),
        FCP_RECEIPT:sha256(fcp_receipt),
        FIRST_PARQUET:sha256(first_parquet),
        FIRST_RECEIPT:sha256(first_receipt),
    }
    report["completed_output_sha256"]={
        RESULT_PARQUET:sha256(out/RESULT_PARQUET),
        RESULT_COUNTS:sha256(out/RESULT_COUNTS),
    }
    (out/RESULT_RECEIPT).write_text(
        json.dumps(report,indent=2,ensure_ascii=False,sort_keys=True,
                   allow_nan=False)+"\n",encoding="utf-8")
    print("MOIS_MERGED_FCP447500_REAL_RGB_AND_FIRST_FULL93M_PHOTOS",
          json.dumps({
              k:report[k] for k in (
                  "preexisting_FCP350399_source_union_original_photo_attempts",
                  "first_full93m_p0002_shared_with_FCP_completed_pixel_source",
                  "additional_original_photo_IDs_not_previously_in_FCP_or_MOIS",
                  "all_genuinely_attempted_unique_original_photo_IDs",
                  "all_original_photo_four_hue_classified",
                  "new_first93m_classified_nominal_taxa_not_in_FCP_or_old_MOIS",
                  "duplicate_source_photo_IDs_in_final_union")
          },sort_keys=True),flush=True)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("fcp-parquet","fcp-receipt","first-parquet","first-receipt","output"):
        p.add_argument("--"+name,type=Path,required=True)
    a=p.parse_args()
    run(a.fcp_parquet,a.fcp_receipt,a.first_parquet,a.first_receipt,a.output)


if __name__=="__main__":
    main()
