#!/usr/bin/env python3
"""Create truly photo-ID-unioned MOIS original-pixel acquisition lake from two real sources.

Old: SHA-authenticated 350,399 MOIS original ROIs + 447,500 FCP original
photo ROIs, with 203,818 overlapping photo IDs = 594,081 unique real pixels.
New: independently acquired full93m observation-probability p=.0002,
18,462 genuine original image attempts / 7,082 known four-hue ROIs.
The second cohort may overlap previously processed FCP photos even though it
does not overlap the first 350,399 MOIS photo IDs. Never double count or
silently change original colour outcomes. No ecological/native/pollination fit.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis.mois_actual_photo_cohort_season_audit_20261019 import (
    HUES, load_pixels, require, sha256,
)

FCP_TOTAL=594081
FCP_CLASSIFIED=288747
P0002_ATTEMPTED=18462
P0002_CLASSIFIED=7082
OLD_FCP_SCHEMA="mois_import_source_signed_FCP447500_true_original_ROIs_on350399_v1"
NEW_SCHEMA="mois_global93m_observation_baseline_real_RGB_photo_results_v1"
SOURCE_SCHEMA="mois_full_source_observation_sampling_cohort_v1"
OLD_FILE="mois_actual350399_plus_FCP447500_photo_ID_deduplicated_original_hues.parquet"
OUTPUT="mois_original_p0002_plus_FCP594081_unique_actual_photo_ROIs.parquet"
REPORT="mois_p0002_plus_FCP594081_actual_pixel_unioned_source_receipt.json"
ARM_TABLE="mois_p0002_FCP594081_original_photo_source_hue_coverage.csv"
UNION_COLUMNS=(
    "photo_id","inat_taxon_id","species","observation_id",
    "latitude","longitude","observed_on","morph","classified",
    "white","arm","measurement_status","technical_failure",
    "photo_colour_provenance","source_original_fcp_morph",
    "source_fcp_global_classifiable","fcp_source_cohort"
)


def validate_fcp_reuse(parquet, receipt):
    rec=json.loads(Path(receipt).read_text())
    require(rec.get("schema")==OLD_FCP_SCHEMA
            and rec.get("full_original_FCP_real_photo_attempts")==447500
            and rec.get("MOIS_previous_real_photo_attempts")==350399
            and rec.get("original_FCP_photo_IDs_already_in_MOIS")==203818
            and rec.get("new_genuinely_measured_FCP_photo_IDs")==243682
            and rec.get("union_total_genuinely_original_photo_IDs")==FCP_TOTAL
            and rec.get("union_total_four_hue_classified")==FCP_CLASSIFIED
            and rec.get("old_source_image_ROI_results_preserved") is True
            and rec.get("no_original_images_downloaded_or_redistributed") is True
            and rec.get("WCVP_native_or_environmental_analysis_run") is False,
            "FCP/MOIS prior authentic original-photo union receipt not verified")
    require(rec["output_SHA256"][Path(parquet).name]==sha256(parquet),
            "FCP/MOIS 594081 original photo pixel bytes do not match frozen Release")
    old=pd.read_parquet(parquet)
    require(len(old)==FCP_TOTAL and old.photo_id.is_unique
            and set(UNION_COLUMNS)<=set(old),
            "Previously measured 594081 photograph records incomplete")
    old["photo_id"]=old.photo_id.astype(str)
    old["inat_taxon_id"]=old.inat_taxon_id.fillna("").astype(str)
    old["observation_id"]=old.observation_id.fillna("").astype(str)
    require(old.photo_id.str.fullmatch(r"\d+").all()
            and old.classified.eq(old.morph.isin(HUES)).all()
            and old.white.eq(old.morph.eq("white")).all()
            and int(old.classified.sum())==FCP_CLASSIFIED,
            "Frozen FCP/MOIS original image colour labels have drifted")
    return old,rec


def validate_new_original(parquet,receipt,source):
    rec=json.loads(Path(receipt).read_text())
    manifest=json.loads(Path(source).read_text())
    require(rec.get("schema")==NEW_SCHEMA and
            manifest.get("schema")==SOURCE_SCHEMA,
            "Original p0002 real photo outcome/observation source schema changed")
    require(sha256(parquet)==rec["measured_photo_parquet_sha256"]
            and sha256(source)==rec["source_plan_sha256"],
            "Original full93m p0002 image or before-colour source SHA has changed")
    require(rec["attempted"]==P0002_ATTEMPTED
            and rec["classifiable"]==P0002_CLASSIFIED
            and rec["downloaded"]==18436
            and rec["acquisition_failed"]==26
            and rec["unclassifiable_acquired"]==11354
            and rec["full_original_upstream_photo_candidates"]==93209351
            and rec["observation_level_selection_probability"]==0.0002
            and rec["matching_prior_measured_photos_reused"]==71
            and rec["technical_unknowns_not_assumed_nonwhite"] is True
            and rec["ecological_analysis_run"] is False,
            "Expected genuinely measured 18462 p0002 original photo ROIs absent")
    require(manifest.get("selection_probability")==0.0002
            and manifest.get("selection_unit")=="observation_uuid"
            and manifest.get("sampled_photo_ids_before_prior_dedup")==18533
            and manifest.get("previously_measured_photo_ids_skipped")==71
            and manifest.get("new_photo_attempts_planned")==P0002_ATTEMPTED
            and manifest.get("selected_new_photos_unsupported_by_frozen_loader")==0
            and manifest.get("all_unselected_source_photos_eligible_at_higher_probability") is True,
            "Original source photo sampling was not the signed p0002 cohort")
    new=load_pixels(parquet,"global93m_probability_observation")
    require(len(new)==P0002_ATTEMPTED and new.photo_id.is_unique
            and new.photo_id.str.fullmatch(r"\d+").all()
            and int(new.classified.sum())==P0002_CLASSIFIED
            and new.source_arm.eq("global93m_probability_sample_observations").all()
            and new.observation_id.ne("").all(),
            "Photo IDs, observation IDs or source-arm labels missing in new real RGB")
    require(int(new.technical_failure.sum())==rec["acquisition_failed"],
            "Technical image acquisition failures were changed into colour labels")
    return new,rec


def combine(fcp, new):
    require(fcp.photo_id.is_unique and new.photo_id.is_unique,
            "Original image source photo ID duplications")
    old_ids=set(fcp.photo_id)
    overlapping=new[new.photo_id.isin(old_ids)].copy()
    extra=new.loc[~new.photo_id.isin(old_ids)].copy()
    old_lookup=fcp.set_index("photo_id")
    shared=overlapping.set_index("photo_id")
    incompatible_taxon=0
    incompatible_known_colour=0
    if len(shared):
        prior=old_lookup.loc[shared.index]
        incompatible_taxon=int(prior.inat_taxon_id.ne(shared.inat_taxon_id).sum())
        incompatible_known_colour=int(
            (prior.classified & shared.classified &
             prior.morph.ne(shared.morph)).sum())
    # Source priority: an original image already measured under earlier
    # FCP/MOIS gets its existing published colour and QC, never re-inferred
    # from the p0002 outcome, even if two independent runs labelled it differently.
    mapped=pd.DataFrame(index=extra.index)
    for k in UNION_COLUMNS:
        if k in extra:
            mapped[k]=extra[k].values
        else:
            mapped[k]=pd.NA
    mapped["photo_colour_provenance"]="MOIS_first93m_original_observation_probability_p0002_ROI_v4"
    mapped["source_original_fcp_morph"]=""
    mapped["source_fcp_global_classifiable"]=False
    mapped["fcp_source_cohort"]=""
    mapped["technical_failure"]=mapped.measurement_status.eq("image_acquisition_failed")
    mapped["classified"]=mapped.morph.isin(HUES)
    mapped["white"]=mapped.morph.eq("white")
    # Retain established schema/column order; do not overwrite the frozen
    # 594081-photo original FCP/MOIS table.
    compatible=mapped.reindex(columns=fcp.columns)
    for col in compatible:
        if col not in mapped.columns:
            compatible[col]=pd.NA
    combined=pd.concat([fcp,compatible],ignore_index=True,sort=False)
    require(len(combined)==len(fcp)+len(extra)
            and combined.photo_id.is_unique
            and combined.iloc[:len(fcp)][list(UNION_COLUMNS)].
                reset_index(drop=True).equals(
                    fcp[list(UNION_COLUMNS)].reset_index(drop=True)),
            "Original source 594081 photo hue metadata were overwritten")
    require(combined.classified.eq(combined.morph.isin(HUES)).all()
            and combined.white.eq(combined.morph.eq("white")).all(),
            "New source technically unknown hues were converted into colours")
    prev_taxa=set(fcp.loc[fcp.classified,"inat_taxon_id"])-{""}
    new_taxa=set(extra.loc[extra.classified,"inat_taxon_id"])-{""}
    old_obs=set(fcp.loc[fcp.observation_id.ne(""),"observation_id"])
    novel_obs=set(extra.loc[extra.observation_id.ne(""),"observation_id"])
    report={
        "schema":"mois_source_verified_FCP594081_plus_new93m_p0002_real_photo_ROI_union_v1",
        "preexisting_actual_photo_attempts_FCP_MOIS":len(fcp),
        "new_first_full93m_observation_probability_photo_attempts":len(new),
        "new_first93m_photo_IDs_already_in_FCP_MOIS":len(overlapping),
        "new_photo_IDs_not_previously_measured_in_FCP_MOIS":len(extra),
        "combined_unique_original_photo_attempts":len(combined),
        "previous_FCP_MOIS_four_hue_classified":int(fcp.classified.sum()),
        "new_unique_first93m_four_hue_classified":int(extra.classified.sum()),
        "combined_unique_original_four_hue_classified":int(combined.classified.sum()),
        "new_unique_unclassifiable_or_failed_ROIs":int((~extra.classified).sum()),
        "overlapping_original_photo_ID_taxon_conflicts_preserved_old":incompatible_taxon,
        "overlapping_original_photo_ID_both_classified_hue_conflicts_preserved_old":
            incompatible_known_colour,
        "added_nominal_taxa_with_genuine_four_hue_not_in_FCP_MOIS":len(new_taxa-prev_taxa),
        "new_unique_actual_image_observation_IDs":len(novel_obs),
        "new_unique_observation_IDs_not_in_old_FCP_MOIS":len(novel_obs-old_obs),
        "first93m_original_source_inclusion_probability":.0002,
        "all_original_photo_IDs_unambiguous":True,
        "old_genuine_FCP_MOIS_colour_rows_untouched":True,
        "source_metadata_without_real_RGB_never_inferred":True,
        "all_original_photo_classification_failures_remain_missing_not_nonwhite":True,
        "no_native_species_or_environmental_modelling_performed":True,
        "reused_FCP_is_same_platform_not_independent_photo_provider":True
    }
    return combined,report


def run(fcp_parquet,fcp_receipt,new_parquet,new_receipt,source_manifest,output):
    output=Path(output)
    require(not output.exists() or not any(output.iterdir()),
            "Cannot overwrite a previously signed original photo pixel union")
    old,_=validate_fcp_reuse(fcp_parquet,fcp_receipt)
    new,_=validate_new_original(new_parquet,new_receipt,source_manifest)
    joined,report=combine(old,new)
    output.mkdir(parents=True,exist_ok=True)
    joined.to_parquet(output/OUTPUT,index=False,compression="zstd")
    arms=[]
    for arm,g in joined.groupby("arm",observed=True,sort=True):
        known=g.loc[g.classified]
        arms.append({
            "source_arm":str(arm),"real_photo_attempts":len(g),
            "classifiable_four_hue_ROIs":len(known),
            "unknown_colour_ROIs":int((~g.classified).sum()),
            "acquisition_failed_when_retained":int(
                g.measurement_status.fillna("").eq("image_acquisition_failed").sum()),
            "nominal_taxa_with_classifiable_hue":int(known.inat_taxon_id.nunique()),
            **{h:int(known.morph.eq(h).sum()) for h in HUES}
        })
    pd.DataFrame(arms).to_csv(output/ARM_TABLE,index=False)
    report["original_input_sha256"]={
        "FCP_MOIS_594081_original_pixel_Parquet":sha256(fcp_parquet),
        "FCP_MOIS_594081_source_receipt":sha256(fcp_receipt),
        "full93m_p0002_real_pixel_Parquet":sha256(new_parquet),
        "full93m_p0002_real_photo_receipt":sha256(new_receipt),
        "full93m_before_hue_observation_source_manifest":sha256(source_manifest),
    }
    report["output_sha256"]={
        OUTPUT:sha256(output/OUTPUT),
        ARM_TABLE:sha256(output/ARM_TABLE),
    }
    (output/REPORT).write_text(json.dumps(
        report,ensure_ascii=False,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("MOIS_FCP594081_PLUS_FULL93M_FIRST18462_UNIQUE_ACTUAL_PHOTOS",
          json.dumps({k:report[k] for k in (
              "combined_unique_original_photo_attempts",
              "new_first93m_photo_IDs_already_in_FCP_MOIS",
              "new_photo_IDs_not_previously_measured_in_FCP_MOIS",
              "combined_unique_original_four_hue_classified",
              "added_nominal_taxa_with_genuine_four_hue_not_in_FCP_MOIS"
          )},sort_keys=True),flush=True)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("fcp-parquet","fcp-receipt","new-parquet","new-receipt",
                "source-manifest","output"):
        p.add_argument("--"+key,type=Path,required=True)
    a=p.parse_args()
    run(a.fcp_parquet,a.fcp_receipt,a.new_parquet,a.new_receipt,
        a.source_manifest,a.output)


if __name__=="__main__":
    main()
