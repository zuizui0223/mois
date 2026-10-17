#!/usr/bin/env python3
"""Restore actual iNaturalist observed dates for 1751 still-undated NUMERIC IDs.

Only original 295017 source photos, with source-confirmed photo_id and original
numeric observation_id, are queried via bounded 50-observation public API
requests. 4096 other undated photos use nonnumeric observation UUID strings:
they cannot be sent to the numeric ID API without an independently verified
original UUID→numeric observation mapping; they remain DATE UNKNOWN.

Existing historical 99817 S3+API date records never overwritten, original ROI
hues remain immutable. Public present observed_on is observation metadata,
not independently verified photograph capture/exif or open flower stage.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import time
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import requests

from scripts.acquisition.mois_historical_original_observed_on_inat_pilot_20261017 import (
    normalize_public_observation, public_observations_batch,
)
from scripts.analysis.mois_seasonal_267k_climate_normals import (
    HUES, labelled_calendar, require,
)

ORIGINAL_PHOTO_COHORT=295017
INITIAL_DATED_PHOTOS=289170
INITIAL_HUE_DATE=140247
TOTAL_MISSING_PHOTO_DATES=5847
NUMERIC_SOURCE_MISSING=1751
NUMERIC_HUE_MISSING=715
UUID_NONNUMERIC_MISSING=4096
PUBLIC_QUERY_N=50
SOURCE_ZIP_RECEIPT="mois_real295017_signed_photo_union_date_and_colour_receipt.json"
SOURCE_PHOTO="mois_real295017_original_photo_hues_plus_verified_observed_dates.csv.gz"
SOURCE_SHA="b797d85464d2e3e39d38326196dd5ab3462d4e8e19b227dc62604b3bdf047a93"
OUTPUT_SOURCE="mois_real295017_photo_source_plus_exact_remaining_numeric_API_dates.csv.gz"
OUTPUT_LEDGER="mois_original1751_numeric_missing_photo_date_public_iNat_result.csv"
OUTPUT_JSON="mois_original295017_remaining1751_numeric_photo_date_recovery_QC_receipt.json"
SLEEP=1.5


def digest(data:bytes)->str:
    return hashlib.sha256(data).hexdigest()


def load_original(path:Path)->pd.DataFrame:
    with zipfile.ZipFile(path) as z:
        raw=z.read(SOURCE_PHOTO)
        receipt=json.loads(z.read(SOURCE_ZIP_RECEIPT))
    require(digest(raw)==SOURCE_SHA and
            receipt["files_SHA256"][SOURCE_PHOTO]==SOURCE_SHA,
            "Exact real295017 original photo + existing historic qualified-date source drift")
    df=pd.read_csv(io.BytesIO(raw),compression="gzip",
                   dtype=str,keep_default_na=False)
    require(len(df)==ORIGINAL_PHOTO_COHORT and df.photo_id.is_unique and
            df.qualified_date.str.len().eq(10).sum()==INITIAL_DATED_PHOTOS and
            (df.morph.isin(HUES)&df.qualified_date.str.len().eq(10)).sum()==INITIAL_HUE_DATE,
            "Original photo/date/automatic hue baseline changed")
    return df


def choose_unknown_numerical(d:pd.DataFrame)->pd.DataFrame:
    missing=d[d.qualified_date.eq("")].copy()
    require(len(missing)==TOTAL_MISSING_PHOTO_DATES and
            missing.photo_id.is_unique and missing.observation_id.is_unique,
            "Actual missing original photo/day attempt baseline changed")
    numeric=missing[missing.observation_id.str.fullmatch(r"\d+")].copy()
    other=missing[~missing.observation_id.str.fullmatch(r"\d+")].copy()
    require(len(numeric)==NUMERIC_SOURCE_MISSING and len(other)==UUID_NONNUMERIC_MISSING
            and numeric.morph.isin(HUES).sum()==NUMERIC_HUE_MISSING and
            other.source_original_measured.eq("mois_wave1_new_original_ROI_v4").all(),
            "Actual numeric iNaturalist original photo/date source support drift")
    require(numeric.old_mois_or_new_fcp.value_counts().to_dict()=={
            "historical":726,"fcp_additional_ROI":1025},
            "Actual original source photo missing numeric date composition changed")
    require(numeric.photo_id.str.isdecimal().all() and
            numeric.inat_taxon_id.str.isdecimal().all() and
            numeric.observation_id.str.isdecimal().all(),
            "Source photo/taxon/observation identifiers aren't numerical canonical IDs")
    # Exact literal original photo+observation ID, never source climate/photo hue selection.
    return numeric.sort_values("photo_id",kind="stable").reset_index(drop=True)


def merge(d:pd.DataFrame,source:pd.DataFrame,
          retrieved:pd.DataFrame)->tuple[pd.DataFrame,dict]:
    require(len(source)==NUMERIC_SOURCE_MISSING and
            len(retrieved)==NUMERIC_SOURCE_MISSING and
            retrieved.original_photo_id.is_unique and
            retrieved.original_photo_id.tolist()==source.photo_id.tolist(),
            "Recovered iNat date photo provenance not source ID matched")
    require(retrieved.original_observation_id.tolist()==source.observation_id.tolist()
            and retrieved.source_taxon_id.tolist()==source.inat_taxon_id.tolist(),
            "Source original numeric observation/photo/taxon mapping changed")
    exact=retrieved.provenance_status.eq(
        "source_photo_and_observation_verified_with_current_public_date")
    accepted=retrieved.loc[exact,["original_photo_id","retrieved_observed_on"]].rename(
        columns={"original_photo_id":"photo_id",
                 "retrieved_observed_on":"extra_api_date"})
    merged=d.merge(accepted,on="photo_id",how="left",validate="one_to_one")
    require(len(merged)==ORIGINAL_PHOTO_COHORT and
            not ((merged.qualified_date.ne("")) &
                 merged.extra_api_date.notna()).any(),
            "Original 99817 historical S3/API and all other source dates must not be overwritten")
    enriched=merged[d.columns].copy()
    extra=merged.extra_api_date.notna()
    enriched.loc[extra,"qualified_date"]=merged.loc[extra,"extra_api_date"]
    enriched.loc[extra,"original_source_date_provenance"] = (
        "current_iNat_exact_numeric_observation_original_photo_date_verified")
    require(enriched.photo_id.tolist()==d.photo_id.tolist() and
            enriched[[c for c in d if c not in (
                "qualified_date","original_source_date_provenance")]].equals(
                    d[[c for c in d if c not in (
                        "qualified_date","original_source_date_provenance")]]),
            "Source original photograph pixels/automatic hues/taxa/location reclassified")
    before=(d.morph.isin(HUES)&d.qualified_date.str.len().eq(10)).sum()
    after=(enriched.morph.isin(HUES)&enriched.qualified_date.str.len().eq(10)).sum()
    extra_class=(source.morph.isin(HUES).to_numpy()&exact.to_numpy()).sum()
    require(int(after-before)==int(extra_class),
            "Date-enrichment classifiable original photograph count does not reconcile")
    return enriched,{
        "original_actual_photo_pixels":ORIGINAL_PHOTO_COHORT,
        "original_date_and_four_hue_photo_ids":int(before),
        "newly_verified_source_photo_numeric_original_obs_dates":int(exact.sum()),
        "newly_dated_classifiable_four_hue_original_source_photos":int(extra_class),
        "after_verified_source_photo_dates":int(enriched.qualified_date.str.len().eq(10).sum()),
        "after_date_and_actual_classifiable_four_hue_original_source_photos":int(after),
        "non_numeric_source_missing_date_photos_not_queried":UUID_NONNUMERIC_MISSING,
        "numeric_source_attempts_still_missing_day":int((~exact).sum()),
        "full_source_hue_classification_not_recomputed":True,
        "no_current_existing_verified_observed_on_day_changed":True,
        "photo_specific_open_flower_organ_not_verified":True,
    }


def run(archive:Path,output:Path):
    require(not output.exists() or not any(output.iterdir()),
            "Never overwrite a previous original photo-date source result")
    original=load_original(archive)
    selected=choose_unknown_numerical(original)
    rows=[];batch=[]
    with requests.Session() as session:
        for start in range(0,len(selected),PUBLIC_QUERY_N):
            if start:time.sleep(SLEEP)
            chunk=selected.iloc[start:start+PUBLIC_QUERY_N]
            ids=chunk.observation_id.tolist()
            found,receipt=public_observations_batch(ids,session=session)
            for source in chunk[["photo_id","observation_id","inat_taxon_id"]].to_dict(
                    orient="records"):
                rows.append(normalize_public_observation(
                    source,found.get(str(source["observation_id"]))))
            batch.append({"batch_number":len(batch)+1,**receipt})
            if len(batch)%6==0:
                print("MOIS_NUMERIC_MISSING_ORIGINAL_PHOTO_DATE_PROGRESS",
                      json.dumps({"api_batches":len(batch),"source_original_photos_done":len(rows)}),
                      flush=True)
    require(len(rows)==NUMERIC_SOURCE_MISSING and len(batch)==36,
            "Original source numeric photo/date API batches incomplete")
    records=pd.DataFrame(rows)
    enriched,counts=merge(original,selected,records)
    output.mkdir(parents=True,exist_ok=True)
    patch=output/OUTPUT_LEDGER
    records.to_csv(patch,index=False,lineterminator="\n")
    full=output/OUTPUT_SOURCE
    enriched.to_csv(full,index=False,lineterminator="\n",
                    compression={"method":"gzip","mtime":0})
    results={
       "schema":"mois_295017_source_missing_1751_original_numeric_observation_photo_ID_public_date_recovery_v1",
       "previous_real_295017_original_photo_plus_99817_historic_date_SHA256":SOURCE_SHA,
       "source_photo_date_target_blind_to_photo_hue_environment_and_geography":True,
       "actual_original_photo_numeric_observation_IDs_without_dated_photo":NUMERIC_SOURCE_MISSING,
       "genuine_missing_original_photo_hues_among_numeric_date_targets":NUMERIC_HUE_MISSING,
       "original_missing_photo_IDs_without_numeric_observation_not_queried":UUID_NONNUMERIC_MISSING,
       "actual_public_API_requests":len(batch),
       "max_original_observations_per_API_call":PUBLIC_QUERY_N,
       "public_API_query_provider_metadata_fetch_date_UTC":datetime.now(timezone.utc).date().isoformat(),
       "original_API_response_SHA256_and_receipts":batch,
       "source_status":{str(k):int(v) for k,v in records.provenance_status.value_counts().items()},
       "current_taxonomy_matches_original_source_for_date_eligible":int(
           records.loc[records.provenance_status.eq(
               "source_photo_and_observation_verified_with_current_public_date"),
               "taxon_id_now_matches_source"].astype(str).str.lower().eq("true").sum()),
       "original_photo_exact_date_row_CSV_SHA256":digest(patch.read_bytes()),
       "full_295017_source_plus_verified_day_csv_SHA256":digest(full.read_bytes()),
       "no_technical_unknown_classified_as_chromatic_or_flowering":True,
       "original_99817_historical_date_precedence_unchanged":True,
       "current_public_observation_date_not_independently_image_EXIF":True,
       "no_original_photo_pixel_reclassification_or_pigment_plasticity_inference":True,
       **counts,
    }
    (output/OUTPUT_JSON).write_text(json.dumps(results,ensure_ascii=False,
                                         indent=2,sort_keys=True)+"\n")
    print("MOIS_REAL_REMAINING_NUMERIC_ORIGINAL_PHOTO_DATES",json.dumps(
        {k:results[k] for k in (
          "newly_verified_source_photo_numeric_original_obs_dates",
          "newly_dated_classifiable_four_hue_original_source_photos",
          "after_verified_source_photo_dates",
          "after_date_and_actual_classifiable_four_hue_original_source_photos"
        )},sort_keys=True),flush=True)
    return results


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--signed-original295017-archive",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    run(a.signed_original295017_archive,a.output)
