#!/usr/bin/env python3
"""Append only genuine second p=.0004 original 93m-source RGB photos.

Existing 612483 source-photo pixels are immutable, technical unknowns never
become coloured flowers, and no native/abiotic/phenology model is fitted.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import pandas as pd
from scripts.analysis.mois_actual_photo_cohort_season_audit_20261019 import (
    HUES, require, sha256,
)

OLD_N=612483
OLD_KNOWN=295801
NEW_N=18119
SOURCE_FRAME=93209351
PRIOR_SOURCE_DRAW=36758
PRIOR_MEASURED=18639
P=.0004
BASE_FILE="mois_original_p0002_plus_FCP594081_unique_actual_photo_ROIs.parquet"
BASE_RECEIPT="mois_p0002_plus_FCP594081_actual_pixel_unioned_source_receipt.json"
OUT_PHOTOS="mois_FCP_and_first_second_93m_actual_original_photo_ROIs.parquet"
OUT_RECEIPT="mois_p0004_actual_photo_append_receipt.json"
OUT_QC="mois_p0004_actual_photo_source_arm_hue_quality.csv"


def combine(old,new,*,expected_old=OLD_N,expected_new=NEW_N):
    require(len(old)==expected_old and old.photo_id.is_unique,
            "Historical original image photo ID count/uniqueness changed")
    require(len(new)==expected_new and new.photo_id.is_unique,
            "Second original source photo IDs incomplete or duplicated")
    require(new.photo_id.astype(str).str.fullmatch(r"\d+").all(),
            "Non-numeric original photo IDs")
    require(not new.photo_id.isin(old.photo_id).any(),
            "Second p0004 photo duplicates a previously genuinely measured original ID")
    require(old.classified.eq(old.morph.isin(HUES)).all()
            and new.classified.eq(new.morph.isin(HUES)).all(),
            "Classifiability disagrees with original four-hue ROI label")
    require(old.white.eq(old.morph.eq("white")).all()
            and new.white.eq(new.morph.eq("white")).all(),
            "Unknown/source technically failed photographic ROIs are not white")
    frozen=old.copy(deep=True)
    old=old.copy(deep=True)  # Never add columns to the caller's source-signed historic frame.
    fresh=new.copy(deep=True)
    fresh["arm"]="global93m_second_p0004"
    fresh["photo_colour_provenance"]="genuine_93m_second_nested_p0004_observation_RGB"
    fresh["source_original_fcp_morph"]=""
    fresh["source_fcp_global_classifiable"]=False
    fresh["fcp_source_cohort"]=""
    fresh["technical_failure"]=fresh.measurement_status.eq("image_acquisition_failed")
    for col in old.columns:
        if col not in fresh:
            fresh[col]=pd.NA
    for col in fresh.columns:
        if col not in old:
            old[col]=pd.NA
    # The SHA-authenticated old photo table stores environmental coordinates
    # and climate fields as numeric values. The new per-photo CSV shard is read
    # with dtype=str, so these numeric fields arrive as numeric-looking text.
    # Make only those fields use the old numeric schema: Arrow must not silently
    # coerce 29.7127033333 into an incompatible mixed float/string object.
    # Neither old source row values nor source photo / observation IDs change.
    for col in frozen.columns:
        old_series = frozen[col]
        if (col not in fresh or
                not pd.api.types.is_numeric_dtype(old_series) or
                pd.api.types.is_bool_dtype(old_series)):
            continue
        raw = fresh[col]
        stripped = raw.astype("string").str.strip()
        given = stripped.notna() & stripped.ne("")
        parsed = pd.to_numeric(raw, errors="coerce")
        require(not bool((given & parsed.isna()).any()),
                "Non-numeric original source metadata in authenticated field: " + col)
        if pd.api.types.is_integer_dtype(old_series):
            nonmissing = parsed.dropna()
            require(bool(nonmissing.mod(1).eq(0).all()),
                    "Non-numeric fractional value in old integer metadata field: " + col)
            fresh[col] = parsed.astype("Int64")
        else:
            fresh[col] = parsed.astype("float64")
    fresh=fresh[old.columns]
    out=pd.concat([old,fresh],ignore_index=True,sort=False)
    common=["photo_id","inat_taxon_id","species","observation_id",
            "morph","classified","white","arm"]
    require(len(out)==expected_old+expected_new and out.photo_id.is_unique
            and out.iloc[:expected_old][common].reset_index(drop=True).equals(
                frozen[common].reset_index(drop=True)),
            "Existing FCP/MOIS first-93m source photo colour rows were altered")
    require(out.classified.eq(out.morph.isin(HUES)).all()
            and out.white.eq(out.morph.eq("white")).all(),
            "Final source includes an invented four-hue classification")
    prior_taxa=set(frozen.loc[frozen.classified,"inat_taxon_id"].astype(str))
    added_taxa=set(fresh.loc[fresh.classified,"inat_taxon_id"].astype(str))
    return out,{
        "old_original_photo_attempts":len(frozen),
        "second_nested_actual_image_attempts":len(fresh),
        "new_unique_original_photo_IDs":len(fresh),
        "all_unique_original_photo_attempts":len(out),
        "old_four_hue_classified":int(frozen.classified.sum()),
        "new_four_hue_classified":int(fresh.classified.sum()),
        "total_four_hue_classified":int(out.classified.sum()),
        "new_unknown_or_technical_failed":int((~fresh.classified).sum()),
        "new_acquisition_failed":int(fresh.technical_failure.sum()),
        "new_four_hue_nominal_taxa_not_in_old":len(added_taxa-prior_taxa),
        "old_original_image_petal_ROI_rows_immutable":True,
        "new_source_photo_IDs_disjoint_from_all_old_FCP_and_MOIS":True,
        "native_climate_pollination_models_not_run":True
    }


def run(base,base_receipt,second,second_receipt,second_manifest,output):
    output=Path(output)
    require(not output.exists() or not any(output.iterdir()),
            "Already signed original photo data cannot be overwritten")
    a=json.loads(Path(base_receipt).read_text())
    b=json.loads(Path(second_receipt).read_text())
    q=json.loads(Path(second_manifest).read_text())
    require(a.get("schema")==
            "mois_source_verified_FCP594081_plus_new93m_p0002_real_photo_ROI_union_v1"
            and a.get("combined_unique_original_photo_attempts")==OLD_N
            and a.get("combined_unique_original_four_hue_classified")==OLD_KNOWN
            and a.get("old_genuine_FCP_MOIS_colour_rows_untouched") is True
            and sha256(base)==a["output_sha256"][BASE_FILE],
            "612483 genuinely measured old original RGB photo Parquet not authenticated")
    require(q.get("schema")=="mois_full_source_observation_sampling_cohort_v1"
            and q.get("full_frame_photo_ids")==SOURCE_FRAME
            and q.get("selection_unit")=="observation_uuid"
            and q.get("selection_probability")==P
            and q.get("sampled_photo_ids_before_prior_dedup")==PRIOR_SOURCE_DRAW
            and q.get("previously_measured_photo_ids_skipped")==PRIOR_MEASURED
            and q.get("new_photo_attempts_planned")==NEW_N
            and q.get("selected_new_photos_unsupported_by_frozen_loader")==0
            and q.get("no_per_taxon_or_observation_photo_count_cap") is True,
            "Second p0004 pre-colour source observation sampling manifest modified")
    require(b.get("schema")=="mois_global93m_observation_baseline_real_RGB_photo_results_v1"
            and b.get("attempted")==NEW_N
            and b.get("matching_prior_measured_photos_reused")==PRIOR_MEASURED
            and b.get("observation_level_selection_probability")==P
            and b.get("technical_unknowns_not_assumed_nonwhite") is True
            and b.get("ecological_analysis_run") is False
            and b.get("downloaded")+b.get("acquisition_failed")==NEW_N
            and b.get("classifiable")+b.get("unclassifiable_acquired")+
                 b.get("acquisition_failed")==NEW_N
            and sha256(second)==b.get("measured_photo_parquet_sha256")
            and sha256(second_manifest)==b.get("source_plan_sha256"),
            "Original second 18119 source image terminal results incomplete/changed")
    old=pd.read_parquet(base)
    fresh=pd.read_parquet(second)
    fresh["classified"]=fresh.morph.isin(HUES)
    stated=fresh.classified_four_hue.astype(str).str.lower().eq("true")
    require(fresh.classified.eq(stated).all()
            and int(fresh.classified.sum())==b["classifiable"]
            and int(fresh.measurement_status.eq("image_acquisition_failed").sum())==
                b["acquisition_failed"]
            and fresh.source_arm.eq("global93m_probability_sample_observations").all()
            and (pd.to_numeric(fresh.selection_probability,errors="coerce")==P).all(),
            "Raw photo image model outputs differ from signed source or classifiability")
    fresh["white"]=fresh.morph.eq("white")
    out,stats=combine(old,fresh)
    require(stats["old_four_hue_classified"]==OLD_KNOWN,
            "Original FCP/first-image photo four-hue ROI count changed")
    output.mkdir(parents=True,exist_ok=True)
    out.to_parquet(output/OUT_PHOTOS,index=False,compression="zstd")
    rows=[]
    for arm,g in out.groupby("arm",observed=True):
        classified=g[g.classified]
        rows.append({"arm":arm,"attempted":len(g),"classifiable":len(classified),
            "technical_unknown":len(g)-len(classified),
            **{h:int(classified.morph.eq(h).sum()) for h in HUES}})
    pd.DataFrame(rows).to_csv(output/OUT_QC,index=False)
    stats.update({
        "schema":"mois_FCP_and_full93m_second_p0004_real_original_photo_append_v1",
        "full_original_photo_metadata_candidates":SOURCE_FRAME,
        "source_second_nested_selection_probability":P,
        "original_photo_source_sha256":{
            BASE_FILE:sha256(base),BASE_RECEIPT:sha256(base_receipt),
            "measured_photos.parquet":sha256(second),
            "receipt.json":sha256(second_receipt),
            "source.json":sha256(second_manifest)},
        "output_sha256":{
            OUT_PHOTOS:sha256(output/OUT_PHOTOS),
            OUT_QC:sha256(output/OUT_QC)}
    })
    (output/OUT_RECEIPT).write_text(json.dumps(
        stats,indent=2,sort_keys=True,ensure_ascii=False,allow_nan=False)+"\n")
    print("MOIS_SECOND_FULL93M_P0004_REAL_RGB_ADDED",
          json.dumps({k:stats[k] for k in (
              "new_unique_original_photo_IDs",
              "all_unique_original_photo_attempts",
              "new_four_hue_classified",
              "total_four_hue_classified",
              "new_four_hue_nominal_taxa_not_in_old",
          )},sort_keys=True),flush=True)
    return stats


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("base","base-receipt","second","second-receipt",
                "second-manifest","output"):
        p.add_argument("--"+key,type=Path,required=True)
    a=p.parse_args()
    run(a.base,a.base_receipt,a.second,a.second_receipt,
        a.second_manifest,a.output)


if __name__=="__main__":
    main()
