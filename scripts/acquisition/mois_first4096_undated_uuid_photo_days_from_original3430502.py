#!/usr/bin/env python3
"""Restore actual observation dates for the initial 4096 unique UUID-ID MOIS photos.

The initial real original first-wave 4096 MOIS images lost their observation
dates when their source numeric IDs were masked/replaced with nonnumeric UUID
values. The FCP authoritative Sept-2026 metadata ledger has 3,430,502 genuine
photo ID rows with original observed_on and source taxon. Stream it ONCE to
match exact original photo_id, source nominal taxon, and a valid real date.
Photo_id is the join, NOT an invented numeric observation ID inferred from UUID.

The 4096 authentic ROI outcomes remain frozen and not botanically reclassified
into open petals. No iNaturalist API call, photo re-download or weather model.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import zipfile
from datetime import date, datetime, timezone
from pathlib import Path

import pandas as pd
from scripts.analysis.mois_seasonal_267k_climate_normals import HUES, require

FCP_AUTHORITATIVE_FULL_SOURCE_SHA="58e037a210e23d3ff5d6f04d392665ada5ee9ac5149d0676672c239cf092df01"
FCP_AUTHORITATIVE_ROWS=3430502
FCP_SOURCE_SNAPSHOT="2026-09-27"
MOIS_295017_SOURCE_SHA="b797d85464d2e3e39d38326196dd5ab3462d4e8e19b227dc62604b3bdf047a93"
MOIS_SOURCE_RECEIPT="mois_real295017_signed_photo_union_date_and_colour_receipt.json"
MOIS_295017_DATE_SOURCE="mois_real295017_original_photo_hues_plus_verified_observed_dates.csv.gz"
TARGET_ORIGINAL_PHOTOS=4096
TARGET_DATE_KNOWN=0
TARGET_CLASSIFIED_HUES=1576
OUTPUT_CSV="mois_first4096_uuid_only_original_photo_source_date_exact3430502.csv"
OUTPUT_JSON="mois_first4096_uuid_observation_photo_day_original_source_QC_receipt.json"


def sha(raw:bytes)->str:
    return hashlib.sha256(raw).hexdigest()


def target_photo_ids(prior_zip:Path)->pd.DataFrame:
    with zipfile.ZipFile(prior_zip) as z:
        data=z.read(MOIS_295017_DATE_SOURCE)
        receipt=json.loads(z.read(MOIS_SOURCE_RECEIPT))
    require(sha(data)==MOIS_295017_SOURCE_SHA and
            receipt["files_SHA256"][MOIS_295017_DATE_SOURCE]==MOIS_295017_SOURCE_SHA,
            "Original 295017 source photographed pixel/date GZIP not frozen")
    d=pd.read_csv(io.BytesIO(data),compression="gzip",
                  dtype=str,keep_default_na=False)
    eligible=d[d.source_original_measured.eq(
        "mois_wave1_new_original_ROI_v4")].copy()
    require(len(d)==295017 and len(eligible)==TARGET_ORIGINAL_PHOTOS and
            eligible.photo_id.is_unique and eligible.qualified_date.eq("").all() and
            eligible.observation_id.str.isdecimal().sum()==0 and
            eligible.morph.isin(HUES).sum()==TARGET_CLASSIFIED_HUES,
            "Exact real first original4096 UUID-only undated four-hue photo cohort changed")
    require(eligible.photo_id.str.isdecimal().all() and
            eligible.inat_taxon_id.str.isdecimal().all(),
            "Original first4096 literal photographic identity/source taxon ID invalid")
    return eligible.sort_values("photo_id",kind="stable").reset_index(drop=True)


def valid_day(x:str)->str|None:
    try:
        day=date.fromisoformat(x)
        if day.year<1800 or day>datetime.now(timezone.utc).date():
            return None
        return day.isoformat()
    except (ValueError,TypeError):
        return None


def scan_full_source(source_csv:Path,target:pd.DataFrame):
    """Stream entire SHA-frozen original photo source and match ONLY target IDs."""
    ids=set(target.photo_id)
    expected=target.set_index("photo_id")
    outcomes={}
    source_ids_seen=set()
    shaobj=hashlib.sha256()
    n=0
    with source_csv.open("rb") as raw:
        class Hashing(io.RawIOBase):
            # Byte-transparent rolling SHA256 of the ORIGINAL data provider
            # metadata source; buffered I/O may read ahead, never modify bytes.
            def readable(self):return True
            def readinto(self, b):
                n_read=raw.readinto(b)
                if n_read:
                    shaobj.update(memoryview(b)[:n_read])
                return n_read
        # The native FCP exact complete source is an uncompressed CSV, not
        # photographic pixels. Digest the entire original downloaded file.
        with io.TextIOWrapper(io.BufferedReader(Hashing()),encoding="utf-8",newline="") as stream:
            reader=csv.DictReader(stream)
            required={"photo_id","inat_taxon_id","observation_id","observed_on",
                      "source_snapshot_date","latitude","longitude"}
            require(required.issubset(reader.fieldnames or []),
                    "FCP authoritative 3.43M original-photo metadata schema drift")
            for row in reader:
                n+=1
                pid=row["photo_id"]
                if pid not in ids:continue
                require(pid not in source_ids_seen,
                        "Original complete source repeats an initial4096 photo ID")
                source_ids_seen.add(pid)
                original=expected.loc[pid]
                taxon_matches=row["inat_taxon_id"]==original["inat_taxon_id"]
                day=valid_day(row["observed_on"])
                date_identity=row["source_snapshot_date"]==FCP_SOURCE_SNAPSHOT
                photosource_has_observation=bool(row["observation_id"])
                latloc=original["latitude"]
                lonloc=original["longitude"]
                try:
                    source_loc_diff=max(abs(float(latloc)-float(row["latitude"])),
                                        abs(float(lonloc)-float(row["longitude"])))
                except (ValueError,TypeError):
                    source_loc_diff=None
                status=("photo_ID_original_taxon_and_snapshot_date_verified"
                        if taxon_matches and day and date_identity and
                           photosource_has_observation
                        else "photo_ID_present_but_taxon_day_or_snapshot_mismatch")
                outcomes[pid]={
                    "original_photo_id":pid,
                    "original_first4096_observation_source_UUID":original["observation_id"],
                    "original_first4096_taxon_id":original["inat_taxon_id"],
                    "current_FCP_source_observation_ID":row["observation_id"],
                    "source_taxon_ID_matches_original_ROI_photo":taxon_matches,
                    "verified_original_full_source_observed_on":
                        day if status=="photo_ID_original_taxon_and_snapshot_date_verified"
                        else None,
                    "source_metadata_snapshot":row["source_snapshot_date"],
                    "photo_source_maximum_abs_geocoordinate_difference":source_loc_diff,
                    "original_source_photo_hue":original["morph"],
                    "date_identity_status":status,
                }
    digest=shaobj.hexdigest()
    require(n==FCP_AUTHORITATIVE_ROWS and digest==FCP_AUTHORITATIVE_FULL_SOURCE_SHA,
            "Full original FCP 3,430,502-photo metadata source SHA or row total differs")
    return outcomes,{"FCP_authoritative_3430502_metadata_SHA256":digest,
                     "FCP_source_rows_fully_scanned":n,
                     "matched_original_first4096_photo_IDs_in_source":len(source_ids_seen)}


def run(photo_source_archive:Path,full_source:Path,out:Path):
    require(not out.exists() or not any(out.iterdir()),
            "Refuse to overwrite completed source photo/date artifact")
    target=target_photo_ids(photo_source_archive)
    found,scan_receipt=scan_full_source(full_source,target)
    rows=[]
    for row in target.to_dict(orient="records"):
        pid=row["photo_id"]
        if pid in found:
            rows.append(found[pid])
        else:
            rows.append({
                "original_photo_id":pid,
                "original_first4096_observation_source_UUID":row["observation_id"],
                "original_first4096_taxon_id":row["inat_taxon_id"],
                "current_FCP_source_observation_ID":None,
                "source_taxon_ID_matches_original_ROI_photo":False,
                "verified_original_full_source_observed_on":None,
                "source_metadata_snapshot":None,
                "photo_source_maximum_abs_geocoordinate_difference":None,
                "original_source_photo_hue":row["morph"],
                "date_identity_status":"photo_ID_not_in_authoritative_FCP_original_source",
            })
    require(len(rows)==TARGET_ORIGINAL_PHOTOS and
            len({r["original_photo_id"] for r in rows})==TARGET_ORIGINAL_PHOTOS,
            "No original UUID-source-photo may be duplicated or imputed")
    d=pd.DataFrame(rows)
    accepted=d.date_identity_status.eq(
        "photo_ID_original_taxon_and_snapshot_date_verified")
    extra=int((accepted&d.original_source_photo_hue.isin(HUES)).sum())
    out.mkdir(parents=True,exist_ok=True)
    csvfile=out/OUTPUT_CSV
    d.to_csv(csvfile,index=False,lineterminator="\n",float_format="%.8f")
    receipt={
      "schema":"mois_first4096_measured_original_photo_UUID_source_date_restored_from_complete3430502_v1",
      "original_firstwave_measured_photo_IDs_without_date":TARGET_ORIGINAL_PHOTOS,
      "original_firstwave_classifiable_four_hue_photo_IDs_without_date":TARGET_CLASSIFIED_HUES,
      "original_full_source_inat_3430502_photo_IDs":FCP_AUTHORITATIVE_ROWS,
      "original_photo_IDs_found_in_full_source":scan_receipt[
           "matched_original_first4096_photo_IDs_in_source"],
      "original_photo_IDs_with_authenticated_full_source_calendar_day":int(accepted.sum()),
      "new_source_dated_original_four_hue_classifiable_photo_IDs":extra,
      "remaining_first4096_photo_IDs_without_source_valid_date":int((~accepted).sum()),
      "original_full_source_RGFCA_20260927_observation_UUID_to_numeric_ID_not_inferred":True,
      "original_four_hue_photo_ROIs_not_reclassified":True,
      "actual_observation_year_month_weather_exposure_NOT_derived":True,
      "photo_specific_open_flower_petal_verification_NOT_claimed":True,
      "original_FCP_source_monthly_metadata_snapshot":FCP_SOURCE_SNAPSHOT,
      "FCP_within_species_D_W_50km_not_imported":True,
      "source_all_4096_original_photo_date_result_SHA256":sha(csvfile.read_bytes()),
      **scan_receipt
    }
    (out/OUTPUT_JSON).write_text(
        json.dumps(receipt,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("MOIS_UUID4096_REAL_ORIGINAL_3430502_SOURCE_PHOTO_DATE_RECOVERY",
          json.dumps({k:receipt[k] for k in (
             "original_photo_IDs_found_in_full_source",
             "original_photo_IDs_with_authenticated_full_source_calendar_day",
             "new_source_dated_original_four_hue_classifiable_photo_IDs",
             "remaining_first4096_photo_IDs_without_source_valid_date"
          )},sort_keys=True),flush=True)
    return receipt


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original295017-source-photo-zip",type=Path,required=True)
    p.add_argument("--authoritative-original-full3430502-photo-metadata-csv",
                   type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    run(a.original295017_source_photo_zip,
        a.authoritative_original_full3430502_photo_metadata_csv,a.output)
