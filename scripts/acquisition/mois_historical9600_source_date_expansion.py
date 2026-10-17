#!/usr/bin/env python3
"""Actual 9,600 ORIGINAL historical photo dates: bounded iNaturalist enrichment.

The original 267425 photo/hue table remains immutable. Select photo IDs blind
to colour/environment by exactly the prefix of the source-verified hash order.
Reuse the 200 exact-photo-verified outcomes from MOIS run 39007591041; query
only 9,400 more original observations (188 x 50 API records). Original current
observed_on is mutable public observation metadata, NOT camera capture time,
confirmed open flowering stage, or an observed-year climate record.
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
    read_signed, labelled_calendar, HUES, MOIS_ROWS, MOIS_HUES, PHOTO_DATA,
    PHOTO_RECEIPT, photo_season_table, require,
)

TARGET_N = 9600
PRIOR_N = 200
BATCH_N = 50
SLEEP_SECONDS = 1.5
HISTORICAL_ROWS = 100543
SOURCE_SHA = "fbb54b58dc91e40e2e166081212588acc18adb19af4bdbcd1c7fe1778bee28d9"
BLIND200_SOURCE_SHA = "2e2cb7098a8539407a54badd89dce2058da2b9d2c56c58ebd8f0f898496604e7"
PRIOR_CSV_SHA = "5a8a4b5e5af8f4a3c2187b1e38b7dbdbafb88de035e376c28c8a0da43f946dff"
PRIOR_CSV = "mois_historical_original_photo_observation_date_blind200_current_iNat.csv"
PRIOR_RECEIPT = "mois_historical_blind200_observation_date_QC_receipt.json"
OUTPUT_CSV = "mois_historical_9600_original_photo_dates_current_public.csv"
SEASON_CSV = "mois_267425_plus_historic9600_derived_season_photo_denominators.csv"
RESULT_JSON = "mois_historical_9600_original_photo_dates_and_photo_seasons_receipt.json"
SELECT_SALT = "mois_historical_dated_source_pilot_v1:"


def hexdigest(binary: bytes) -> str:
    return hashlib.sha256(binary).hexdigest()


def fixed_prefix(d: pd.DataFrame, n: int = TARGET_N) -> pd.DataFrame:
    require(0 < n <= TARGET_N, "Never exceed the explicitly bounded <10k API observations study")
    require(d.photo_id.is_unique and d.observation_id.is_unique and
            d.photo_id.notna().all() and d.observation_id.notna().all(),
            "Original source photo and observation identifier uniqueness was violated")
    keys = d.photo_id.map(lambda p: hexdigest((SELECT_SALT+str(p)).encode()))
    sampled = (d.assign(hash_selection_key=keys)
               .sort_values(["hash_selection_key", "photo_id"],kind="stable")
               .head(n).drop(columns="hash_selection_key").reset_index(drop=True))
    require(sampled.observation_id.is_unique and len(sampled)==n,
            "A source-blind observation was sampled twice")
    return sampled


def read_prior200(pilot_zip: Path, prefix200: pd.DataFrame) -> tuple[pd.DataFrame,dict]:
    with zipfile.ZipFile(pilot_zip) as z:
        receipt = json.loads(z.read(PRIOR_RECEIPT))
        csvdata = z.read(PRIOR_CSV)
    require(receipt["original_full_source_photo_SHA256"]==SOURCE_SHA and
            receipt["blind_sample_photo_count"]==PRIOR_N and
            receipt["blind_source_photo_ids_SHA256"]==BLIND200_SOURCE_SHA and
            receipt["date_recovery_photo_csv_sha256"]==PRIOR_CSV_SHA and
            hexdigest(csvdata)==PRIOR_CSV_SHA,
            "Actual previous blind 200 photo+observation+date source drifted")
    require(hexdigest(("\n".join(prefix200.photo_id)+"\n").encode())==BLIND200_SOURCE_SHA,
            "Original first 200 SHA-prefix photo IDs differ from the signed pilot")
    d=pd.read_csv(io.BytesIO(csvdata),dtype=str)
    require(len(d)==PRIOR_N and d.original_photo_id.tolist()==prefix200.photo_id.tolist() and
            d.original_observation_id.tolist()==prefix200.observation_id.tolist() and
            d.source_taxon_id.tolist()==prefix200.inat_taxon_id.tolist(),
            "Prior 200 historical records were reselected, shuffled, or replaced")
    require(d.provenance_status.eq(
        "source_photo_and_observation_verified_with_current_public_date").all(),
        "Prior verified dates were lost or retrospectively filled")
    return d,receipt


def normal_season_table(source:pd.DataFrame, selected:pd.DataFrame,
                        recovered:pd.DataFrame) -> tuple[pd.DataFrame,dict]:
    """Source photo composition, not confirmed anthesis, same-year weather or site effects."""
    require(source.photo_id.is_unique and selected.photo_id.is_unique,
            "Repeated source photograph not allowed")
    require(recovered.original_photo_id.is_unique and
            recovered.original_photo_id.tolist()==selected.photo_id.tolist(),
            "Date recovery does not preserve blind original photo-ID order")
    eligible = recovered.provenance_status.eq(
        "source_photo_and_observation_verified_with_current_public_date")
    admitted = recovered.loc[eligible,["original_photo_id","retrieved_observed_on"]].rename(
        columns={"original_photo_id":"photo_id",
                 "retrieved_observed_on":"verified_current_observation_date"})
    q = source.merge(admitted,on="photo_id",how="left",validate="one_to_one")
    require(len(q)==MOIS_ROWS and q.photo_id.is_unique,
            "Exact original photographic full-source union changed on date merge")
    had_date=q.observed_on.notna()
    require((~(had_date & q.verified_current_observation_date.notna())).all(),
            "New photo dates collide with original source dates")
    q["calendar_day_for_source_qualified_summary"]=q.observed_on.where(
        had_date,q.verified_current_observation_date)
    # Never silently overwrite or edit any original CSV field/measurement.
    table = q[["photo_id","inat_taxon_id","morph","observation_id",
               "latitude","longitude","old_mois_or_new_fcp",
               "calendar_day_for_source_qualified_summary"]].rename(
        columns={"calendar_day_for_source_qualified_summary":"observed_on"})
    labelled = labelled_calendar(table)
    old_count=int((source.observed_on.notna() & source.morph.isin(HUES)).sum())
    new_classified=int((selected.photo_id.isin(admitted.photo_id) &
                        selected.morph.isin(HUES)).sum())
    require(int((labelled.date.notna() & labelled.classified).sum())==
            old_count+new_classified,
            "A historic source photo or one of its original hue labels was recounted")
    overall=photo_season_table(labelled)
    overall["photo_source_group"]="full_combined_source_after_traced_date_recovery"
    original_only=photo_season_table(labelled[
        labelled.old_mois_or_new_fcp=="historical"])
    original_only["photo_source_group"]="historical_sample_with_verified_current_date"
    historical_not_selected=labelled[
        (labelled.old_mois_or_new_fcp=="historical")&labelled.date.isna()]
    require(len(historical_not_selected)==HISTORICAL_ROWS-int(eligible.sum()),
            "Historical source date remains unknown except exact matched photo IDs")
    grouped=pd.concat([overall,original_only],ignore_index=True)
    stats={
        "original_classified_with_original_observed_on":old_count,
        "new_classified_historical_source_dates_from_verified_original_photo_ID":new_classified,
        "expanded_classified_original_photo_day_records":old_count+new_classified,
        "historic_original_photo_dates_still_unverified":len(historical_not_selected),
        "n_original_source_full_photo_attempts":len(q),
        "n_historic_verified_new_date_ids":int(eligible.sum()),
    }
    return grouped,stats


def run(photo_zip:Path, prior200_zip:Path, output:Path,
        *, n:int=TARGET_N, postbatch_sleep:float=SLEEP_SECONDS) -> dict:
    require(n==TARGET_N,
            "This source-locked stage is exactly 9,600 old-photo IDs; use a separate source-cohort stage for other sizes")
    require(not output.exists() or not any(output.iterdir()),
            "Never rewrite a signed original photograph/date source")
    source,src_receipt,src_sha=read_signed(
        photo_zip,PHOTO_DATA,PHOTO_RECEIPT,"source_photo_level_csv_sha256",
        ["photo_id","observation_id","inat_taxon_id","old_mois_or_new_fcp",
         "observed_on","morph","latitude","longitude"])
    require(src_sha==SOURCE_SHA and len(source)==MOIS_ROWS and
            int(source.morph.isin(HUES).sum())==MOIS_HUES and
            src_receipt["combined_unique_photo_IDs"]==MOIS_ROWS,
            "Genuine original source photo IDs or measured photo hues changed")
    historical=source[source.old_mois_or_new_fcp=="historical"]
    require(len(historical)==HISTORICAL_ROWS and historical.observed_on.isna().all(),
            "Exact original source 100,543-photo historical frame was modified")
    selected=fixed_prefix(historical,n)
    prior200,prior_receipt=read_prior200(prior200_zip,selected.head(PRIOR_N))
    output.mkdir(parents=True,exist_ok=True)
    # Failures remain in action artifacts with the finite batches that succeeded.
    new_rows=[]
    evidence=[]
    with requests.Session() as client:
        rest=selected.iloc[PRIOR_N:]
        for offset in range(0,len(rest),BATCH_N):
            if offset:time.sleep(postbatch_sleep)
            chunk=rest.iloc[offset:offset+BATCH_N]
            ids=chunk.observation_id.tolist()
            retrieved,ledger=public_observations_batch(ids,session=client)
            for srcrow in chunk.to_dict(orient="records"):
                new_rows.append(normalize_public_observation(
                    srcrow,retrieved.get(str(srcrow["observation_id"]))))
            evidence.append({"chunk":1+offset//BATCH_N,
                             "first_photo_hash":hexdigest(str(chunk.photo_id.iloc[0]).encode()),
                             "original_observation_IDs_requested":len(ids),**ledger})
            # Preserve partial source QA if API fails, without treating it
            # as a completed 9600-photo cohort.
            if len(evidence)%12==0:
                (output/"partial_progress_nonfinal.json").write_text(
                    json.dumps({"processed_additional_original_photos":len(new_rows),
                        "authenticated_batches":evidence,
                        "scientific_status":"INCOMPLETE_DO_NOT_ADMIT_TO_RESULT"},indent=2)+"\n")
    require(len(new_rows)==TARGET_N-PRIOR_N and len(evidence)==188,
            "Not all 188 new-source iNaturalist observation batches returned")
    result_rows=pd.concat([prior200,pd.DataFrame(new_rows)],ignore_index=True)
    require(len(result_rows)==TARGET_N and result_rows.original_photo_id.is_unique and
            result_rows.original_photo_id.tolist()==selected.photo_id.tolist(),
            "Original 9600 photo order/completion changed")
    outcome_path=output/OUTPUT_CSV
    result_rows.to_csv(outcome_path,index=False,lineterminator="\n")
    season,stats=normal_season_table(source,selected,result_rows)
    season_path=output/SEASON_CSV
    season.to_csv(season_path,index=False,lineterminator="\n",float_format="%.8f")
    status_counts={str(k):int(v) for k,v in
                   result_rows.provenance_status.value_counts(dropna=False).items()}
    exact=int(result_rows.provenance_status.eq(
        "source_photo_and_observation_verified_with_current_public_date").sum())
    classifiable=int(selected.morph.isin(HUES).sum())
    header={
        "schema":"mois_9600_exact_original_historical_photo_public_dates_v1",
        "original_full_source_photo_sha256":src_sha,
        "prior200_source_photo_ID_selection_SHA256":BLIND200_SOURCE_SHA,
        "prior200_exact_photo_observation_date_csv_sha256":PRIOR_CSV_SHA,
        "sample_selected_entirely_from_original_photo_ID_hash":True,
        "first200_unchanged_original_signed_pilot_reused":True,
        "total_historical_source_original_photo_ids":HISTORICAL_ROWS,
        "total_previous_genuinely_measured_original_photo_ids":MOIS_ROWS,
        "total_previous_classified_photo_hues":MOIS_HUES,
        "blind_photo_source_selected_ids":TARGET_N,
        "original_blind_source_photo_ids_SHA256":hexdigest(
            ("\n".join(selected.photo_id)+"\n").encode()),
        "original_sample_classifiable_ROI_hues":classifiable,
        "additional_current_public_iNaturalist_observations_queried":TARGET_N-PRIOR_N,
        "new_public_iNaturalist_API_requests":len(evidence),
        "new_public_iNaturalist_API_response_SHA256_per_batch":evidence,
        "old200_current_public_dates_originally_retrieved_UTC":prior_receipt["current_public_API_fetch_date_UTC"],
        "new9400_current_public_dates_retrieved_UTC":datetime.now(timezone.utc).date().isoformat(),
        "exact_original_photo_observation_and_current_date_matches":exact,
        "n_current_original_taxonomy_matches":int(
            result_rows.taxon_id_now_matches_source.astype(str).str.lower().eq("true").sum()),
        "current_public_Flowers_observation_level_tags":int(
            result_rows.has_current_observation_level_flowering_annotation
            .astype(str).str.lower().eq("true").sum()),
        "technical_public_photo_date_identity_status":status_counts,
        "historical_source_colour_and_taxonomy_unchanged":True,
        "no_untested_dates_or_open_flower_states_imputed":True,
        "new_original_photo_date_CSV_sha256":hexdigest(outcome_path.read_bytes()),
        "original_photo_month_season_denominators_sha256":hexdigest(season_path.read_bytes()),
        **stats,
        "ecological_claims_not_identified":[
          "iNaturalist current public observed_on may have changed since original archival source",
          "Original photo date metadata is not independent photographic open-flower/petal verification",
          "Summer/winter photo white share varies with sampled species, geography and photographic effort",
          "This blind 9600-photo sample is not full 100543-photo date coverage",
          "No observed-year/month local climate anomalies, genetic polymorphism or pigment fitness fit",
          "No FCP D/W within-species or 50 km spatial colour estimands reused"
        ]
    }
    (output/RESULT_JSON).write_text(json.dumps(
        header,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("MOIS_REAL_9600_HISTORIC_PUBLIC_OBSERVED_ON",json.dumps({
        "observed_on_recovered":exact,
        "blind_photo_n":TARGET_N,
        "new_dated_automatic_hues":stats["new_classified_historical_source_dates_from_verified_original_photo_ID"],
        "full_combined_date_and_colour":stats["expanded_classified_original_photo_day_records"],
        "status":status_counts},sort_keys=True),flush=True)
    return header


if __name__ == "__main__":
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("--source-photo-archive",type=Path,required=True)
    a.add_argument("--prior200-archive",type=Path,required=True)
    a.add_argument("--output",type=Path,required=True)
    args=a.parse_args()
    run(args.source_photo_archive,args.prior200_archive,args.output)
