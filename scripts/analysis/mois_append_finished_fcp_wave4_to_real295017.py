#!/usr/bin/env python3
"""Admit only EXTRA completed original FCP photo pixels beyond MOIS real295017.

Reuses FCP's original exact-photo-ID / 256 terminal authenticity verifier as an
upstream gate. This module takes its verified completed-partition union and
the independently finished 295017-photo RGB source, requiring *all* existing
photos/colours/taxa/dates unchanged and accepting only unseen photo IDs with
genuine completed FCP original ROI-v4 terminal four-hue/unknown measurements.

Data-only transfer for MOIS worldwide interspecific/photographic colour study;
never import FCP W/D/50km polymorphism estimates or impute flowering anatomy.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import zipfile
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis.mois_seasonal_267k_climate_normals import (
    labelled_calendar, HUES, require,
)

OLD_ORIGINAL_PHOTOS = 295017
OLD_FOUR_HUES = 142538
OLD_DATED_ORIGINAL = 289170
OLD_DATED_FOUR_HUES = 140247
OLD_FCP_PARTITIONS = 180
OLD_RAW_SHA = "44e3cb6b4c6174dfb1035de208343c7c6a2e6c321e2557a509224c4d9f2ac7b3"
OLD_DATED_SHA = "b797d85464d2e3e39d38326196dd5ab3462d4e8e19b227dc62604b3bdf047a93"
OLD_RAW = "mois_real295017_original_photo_pixels_three_waves_FCP180.csv.gz"
OLD_DATED = "mois_real295017_original_photo_hues_plus_verified_observed_dates.csv.gz"
OLD_RECEIPT = "mois_real295017_signed_photo_union_date_and_colour_receipt.json"
FCP_SOURCE = "mois_fcp_wave4_PARTIAL_completed_pixel_original_photo_union.csv.gz"
FCP_RECEIPT = "mois_fcp_wave4_PARTIAL_real_completed_original_ROI_receipt.json"
FCP_RUN_ID = 38992543115
FCP_SOURCE_SNAPSHOT_HEAD = "ef448c9514e0c702e7887e00b288ee966728b0d5"
OUT_RAW = "mois_next_finished_fcp_original_photo_pixels_plus295017.csv.gz"
OUT_DATED = "mois_next_finished_fcp_original_photo_date_hues_plus295017.csv.gz"
OUT_SEASON = "mois_next_finished_fcp_original_date_hue_season_table.csv"
OUT_RECEIPT = "mois_next_finished_fcp_original_photo_new_pixel_date_hue_receipt.json"
DATE_TOKEN = "original_FCP_wave4_pixel_terminal_public_observed_on"


def sha(raw:bytes)->str:
    return hashlib.sha256(raw).hexdigest()


def frozen_base(archive:Path):
    with zipfile.ZipFile(archive) as z:
        raw=z.read(OLD_RAW)
        dated=z.read(OLD_DATED)
        old_receipt=json.loads(z.read(OLD_RECEIPT))
    require(sha(raw)==OLD_RAW_SHA and sha(dated)==OLD_DATED_SHA,
            "Original 295017 photo-ROI pixels or dates SHA changed")
    require(old_receipt["files_SHA256"][OLD_RAW]==OLD_RAW_SHA and
            old_receipt["files_SHA256"][OLD_DATED]==OLD_DATED_SHA,
            "Signed 295017 photo/data SHA differs")
    source=pd.read_csv(io.BytesIO(raw),compression="gzip",
                       dtype=str,keep_default_na=False)
    enriched=pd.read_csv(io.BytesIO(dated),compression="gzip",
                         dtype=str,keep_default_na=False)
    require(len(source)==len(enriched)==OLD_ORIGINAL_PHOTOS and
            source.photo_id.is_unique and enriched.photo_id.is_unique and
            source.photo_id.tolist()==enriched.photo_id.tolist(),
            "Signed 295017 original photo identities or order changed")
    require(tuple(enriched.columns[:len(source.columns)])==tuple(source.columns)
            and enriched[source.columns].equals(source),
            "Original 295017 photo source metadata/color was changed by date recovery")
    require(source.morph.isin(HUES).sum()==OLD_FOUR_HUES and
            enriched.qualified_date.str.len().eq(10).sum()==OLD_DATED_ORIGINAL
            and (enriched.morph.isin(HUES)&
                 enriched.qualified_date.str.len().eq(10)).sum()==OLD_DATED_FOUR_HUES,
            "Old authentic measured photo, date or hue source baseline drift")
    return source,enriched,old_receipt


def complete_FCP_source(path:Path,receipt_path:Path):
    raw=path.read_bytes()
    receipt=json.loads(receipt_path.read_text())
    require(receipt["schema"]==
              "mois_authentic_completed_original_FCP_wave4_ROI_PARTITIONS_source_v1"
            and receipt["FCP_live_original_actions_run"]==FCP_RUN_ID
            and receipt["FCP_original_run_head_sha"]==FCP_SOURCE_SNAPSHOT_HEAD
            and OLD_FCP_PARTITIONS<receipt["FCP_completed_ROI_partitions"]<=256
            and receipt["FCP_total_ROI_partitions"]==256
            and receipt["source_photo_level_csv_sha256"]==sha(raw),
            "No newer authentic original FCP worker partitions beyond the 180-photo source")
    require(receipt["missing_FCP_wave4_partitions_are_not_hues"] is True
            and receipt["FCP_genetic_or_within_species_D_W_50km_polymorphism_statistics_imported"] is False,
            "Original photo pixel workflow imported nonmeasured or genetic results")
    fcp=pd.read_csv(io.BytesIO(raw),compression="gzip",
                    dtype=str,keep_default_na=False)
    require(fcp.photo_id.is_unique and
            len(fcp)==receipt["combined_unique_photo_IDs"] and
            fcp.morph.isin(HUES).sum()==receipt["combined_four_hue_classified_photo_IDs"],
            "New source photo terminal mapping/count or original hue changed")
    return fcp,receipt,sha(raw)


def expand(base:pd.DataFrame,enriched:pd.DataFrame,fcp:pd.DataFrame,
           *, min_partitions:int,old_fcp_ids=None):
    """Return genuine photo pixel append and date-bound append; no overlap rewrites."""
    require(base.photo_id.is_unique and enriched.photo_id.is_unique and fcp.photo_id.is_unique,
            "Source original photo ID no longer unique")
    require(len(base)==len(enriched) and len(base)==OLD_ORIGINAL_PHOTOS,
            "Existing original photograph count changed")
    oldid=set(base.photo_id)
    if old_fcp_ids is None:
        old_fcp_ids=set(base.loc[
            base.source_original_measured.eq("fcp_wave4_partial_original_ROI_v4"),
            "photo_id"])
    require(old_fcp_ids.issubset(set(fcp.photo_id)),
            "The expanded original FCP source dropped previously completed authentic partitions")
    overlap=fcp[fcp.photo_id.isin(oldid)]
    new=fcp[~fcp.photo_id.isin(oldid)].copy()
    require(len(new)>0 and new.photo_id.is_unique,
            "No genuine new FCP pixel result beyond the full 295017 source")
    require(new.old_mois_or_new_fcp.eq("fcp_additional_ROI").all() and
            new.source_original_measured.eq("fcp_wave4_partial_original_ROI_v4").all(),
            "New source photographs were not independently completed FCP ROI-v4 pixels")
    unknown_reason_codes=(
        "not_evaluable_roi_or_flip_gate",
        "not_evaluable_ambiguous_palette_composition",
        "not_evaluable_no_biological_palette_mass",
    )
    require(new.morph.isin((*HUES,"mixed_uncertain")).all() and
            new.measurement_status.isin(("classified_four_state_morph",*unknown_reason_codes)).all() and
            new.morph.isin(HUES).eq(
                new.measurement_status.eq("classified_four_state_morph")).all() and
            new.morph.eq("mixed_uncertain").eq(
                new.measurement_status.isin(unknown_reason_codes)).all(),
            "New FCP photo pixel ROI outcome states are not source-coherent")

    compare=overlap.merge(base,on="photo_id",how="inner",validate="one_to_one",
                          suffixes=("_fcp","_old"))
    require(len(compare)==len(overlap),
            "FCP overlap failed to match existing original photo IDs")
    for key in ("inat_taxon_id","morph","observed_on"):
        require(compare[key+"_fcp"].fillna("").eq(
            compare[key+"_old"].fillna("")).all(),
            f"Original FCP terminal disagrees with historic MOIS {key} for a reused photo ID")
    differing_obs=compare.observation_id_fcp.ne(compare.observation_id_old)
    # A FCP UUID vs numeric observation identifier was already recorded for
    # one old photograph. Never count it as two source photo identities.
    mismatch_ids=int(differing_obs.sum())

    # Date must be an ORIGINAL current photo-source observation day. We do not
    # fetch a new date or infer a month from climate/region/season.
    dated=new.observed_on.str.fullmatch(r"\d{4}-\d{2}-\d{2}")
    invalid=new.loc[dated,"observed_on"].map(lambda x: _valid_day(x))
    require(invalid.all(), "New photographic observations have an invalid date")
    new["qualified_date"]=new.observed_on.where(dated,"")
    new["original_source_date_provenance"]=np.where(
        dated,DATE_TOKEN,"unknown_original_FCP_observed_on")
    require(tuple(base.columns)==tuple(fcp.columns) and
            tuple(enriched.columns)==tuple(new.columns),
            "Cannot append a source photo outcome with a changed photo/date schema")
    merged=pd.concat([base,new[list(base.columns)]],ignore_index=True)
    with_dates=pd.concat([enriched,new[list(enriched.columns)]],ignore_index=True)
    require(merged.photo_id.is_unique and with_dates.photo_id.is_unique and
            merged.iloc[:len(base)].reset_index(drop=True).equals(base.reset_index(drop=True)) and
            with_dates.iloc[:len(base)].reset_index(drop=True).equals(enriched.reset_index(drop=True)),
            "Original 295017 measured source photo rows were edited or reclassified")
    original_date=(enriched.qualified_date.str.len()==10)
    final_date=(with_dates.qualified_date.str.len()==10)
    summary={
        "prior_original_measured_photo_IDs":len(base),
        "prior_original_four_hue_photo_ROIs":int(base.morph.isin(HUES).sum()),
        "prior_source_dated_four_hue_photos":int(
            (enriched.morph.isin(HUES)&original_date).sum()),
        "additional_original_photo_pixels":len(new),
        "additional_original_four_hue_ROI_labels":int(new.morph.isin(HUES).sum()),
        "additional_new_photos_with_original_observed_date":int(dated.sum()),
        "additional_new_photo_hues_with_original_observed_date":
            int((new.morph.isin(HUES)&dated).sum()),
        "prior_measured_FCP180_original_photo_ids":len(old_fcp_ids),
        "FCP_finished_original_source_photo_ids_overlapping_295017":len(overlap),
        "FCP_original_photo_source_numeric_observation_ID_discrepancies":mismatch_ids,
        "expanded_distinct_measured_source_photo_IDs":len(merged),
        "expanded_four_hue_source_ROIs":int(merged.morph.isin(HUES).sum()),
        "expanded_dated_photo_attempts":int(final_date.sum()),
        "expanded_dated_original_four_hue_ROIs":
            int((with_dates.morph.isin(HUES)&final_date).sum()),
        "real_original_photo_pixels_and_source_hues_not_relabelled":True,
        "missing_hues_never_reclassified_as_nonwhite_or_flowering_absent":True,
        "old_99817_S3_API_historical_original_photo_dates_unmodified":True,
        "FCP_original_local_D_W_or_genetic_50km_estimands_NOT_imported":True,
    }
    require(summary["expanded_distinct_measured_source_photo_IDs"]==
            summary["prior_original_measured_photo_IDs"]+
            summary["additional_original_photo_pixels"] and
            summary["expanded_dated_original_four_hue_ROIs"]==
            summary["prior_source_dated_four_hue_photos"]+
            summary["additional_new_photo_hues_with_original_observed_date"],
            "Distinct real pixel original source/date and four-hue counts not additive")
    return merged,with_dates,summary


def _valid_day(x:str)->bool:
    try:
        return date.fromisoformat(x).isoformat()==x
    except (ValueError,TypeError):
        return False


def run(prior:Path,finished:Path,finished_receipt:Path,output:Path):
    require(not output.exists() or not any(output.iterdir()),
            "Never overwrite a prior actual photographed source result")
    raw,dated,previous_receipt=frozen_base(prior)
    fcp,src_proof,fcp_hash=complete_FCP_source(finished,finished_receipt)
    union,enriched,result=expand(
        raw,dated,fcp,
        min_partitions=src_proof["FCP_completed_ROI_partitions"])
    result.update({
        "schema":"mois_next_finished_FCP_wave4_ROI_plus_295017_source_photo_date_union_v1",
        "FCP_original_run":FCP_RUN_ID,
        "FCP_real_completed_photo_ROI_partitions":src_proof["FCP_completed_ROI_partitions"],
        "FCP_total_planned_partitions_not_assumed_measured":256,
        "FCP_completed_source_snapshot_full_SHA256":fcp_hash,
        "prior295017_source_original_photo_SHA256":OLD_RAW_SHA,
        "prior295017_public_date_source_SHA256":OLD_DATED_SHA,
        "photo_specific_biological_anthesis_not_verified":True,
        "real_observed_year_month_weather_not_recomputed":True,
        "overlapping_FCP_source_photo_pixels_agreed":True,
    })
    output.mkdir(parents=True,exist_ok=True)
    rawfile=output/OUT_RAW
    datefile=output/OUT_DATED
    for df,path in ((union,rawfile),(enriched,datefile)):
        df.to_csv(path,index=False,
                  compression={"method":"gzip","mtime":0},lineterminator="\n")
    result["expanded_real_original_photo_source_SHA256"]=sha(rawfile.read_bytes())
    result["expanded_real_date_hue_photo_source_SHA256"]=sha(datefile.read_bytes())

    panel=enriched.copy()
    panel["observed_on"]=panel.qualified_date
    cal=labelled_calendar(panel)
    require(int((cal.date.notna()&cal.classified).sum())==
            result["expanded_dated_original_four_hue_ROIs"],
            "Date and genuine original hue summation disagrees with seasonal parser")
    annual=[]
    for zone in ("north_extratropical","south_extratropical"):
        z=cal[cal.zone.eq(zone)&cal.date.notna()&cal.season.isin(("summer","winter"))]
        for season in ("summer","winter"):
            s=z[z.season.eq(season)]
            classified=s[s.classified]
            annual.append({
                "hemisphere":zone,"meteorological_season":season,
                "source_original_photo_attempts_with_date":len(s),
                "source_actual_classifiable_photo_hues_with_date":len(classified),
                "source_photo_white_hue_ROIs_with_date":int(classified.white.sum()),
                "photo_white_fraction_conditional_classifiable":
                    round(float(classified.white.mean()),8) if len(classified) else None,
            })
    seasonal=output/OUT_SEASON
    pd.DataFrame(annual).to_csv(seasonal,index=False,lineterminator="\n")
    result["actual_original_photo_season_table_SHA256"]=sha(seasonal.read_bytes())
    (output/OUT_RECEIPT).write_text(json.dumps(
        result,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("MOIS_NEW_FINISHED_FCP_PIXELS_PLUS_ORIGINAL295017",
          json.dumps({k:result[k] for k in (
            "FCP_real_completed_photo_ROI_partitions",
            "additional_original_photo_pixels","additional_original_four_hue_ROI_labels",
            "additional_new_photo_hues_with_original_observed_date",
            "expanded_distinct_measured_source_photo_IDs",
            "expanded_dated_original_four_hue_ROIs")},sort_keys=True),flush=True)
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--prior295017-authentic-archive",type=Path,required=True)
    p.add_argument("--finished-fcp-photo-source",type=Path,required=True)
    p.add_argument("--finished-fcp-photo-receipt",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    q=p.parse_args()
    run(q.prior295017_authentic_archive,q.finished_fcp_photo_source,
        q.finished_fcp_photo_receipt,q.output)
