#!/usr/bin/env python3
"""MOIS: stream official monthly iNaturalist open metadata, restore historic dates.

Two original provider files are large (31.75 GiB total, 2026-09-27). They are
streamed ONCE, with pre-verified Content-Length+ETag and rolling compressed
SHA256, never written to the GitHub runner's disk. photos.csv.gz joins original
photo_id -> observation_uuid; observations.csv.gz joins UUID -> observed_on.
This public OPEN PHOTO dataset excludes some licenses; no unmatched photo is
assigned an invented date. Critically observation_uuid is NOT source numeric
observation_id. Prior 200 photo-specific API checks remain distinct corroboration.
No flower stage, pigmentation adaptation, or causal weather estimate is inferred.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import io
import json
import time
import zipfile
from datetime import date, datetime, timezone
from pathlib import Path

import pandas as pd
import requests

from scripts.analysis.mois_seasonal_267k_climate_normals import (
    HUES, MOIS_ROWS, MOIS_HUES, PHOTO_DATA, PHOTO_RECEIPT,
    read_signed, labelled_calendar, photo_season_table, require,
)
from scripts.analysis.mois_seasonal_spatial_identifiability_20261017 import eligible_site_photo_support

HISTORIC=100543
ORIGINAL_SHA="fbb54b58dc91e40e2e166081212588acc18adb19af4bdbcd1c7fe1778bee28d9"
PROBE="mois_official_inaturalist_monthly_open_s3_metadata_source_probe.json"
OLD_CSV="mois_historical_original_photo_observation_date_blind200_current_iNat.csv"
OLD_RECEIPT="mois_historical_blind200_observation_date_QC_receipt.json"
OLD_SHA="5a8a4b5e5af8f4a3c2187b1e38b7dbdbafb88de035e376c28c8a0da43f946dff"
FINAL_TABLE="mois_historical_100543_original_photo_id_s3_uuid_date_real_snapshot.csv.gz"
SEASON_TABLE="mois_full267425_after_official_S3_historical_date_recovery_season_photos.csv"
RECEIPT="mois_full_historical_photo_dates_official_monthly_s3_source_receipt.json"
PROGRESS_N=20_000_000


class DigestReader:
    """Only binary read; gzip.GzipFile consumes it as an input stream."""
    def __init__(self, raw):
        self.raw=raw
        self.sha256=hashlib.sha256()
        self.bytes=0
    def read(self,n=-1):
        v=self.raw.read(n)
        if v:
            self.bytes+=len(v)
            self.sha256.update(v)
        return v


def provider_snapshot_stream(client:requests.Session,entry:dict):
    """HTTP snapshot-source lease; returns response and hashing compressed reader."""
    url=entry["url"]
    probe_length=int(entry["content_length_bytes"])
    expected_etag=entry["etag"]
    head=client.head(url,timeout=30,allow_redirects=True)
    head.raise_for_status()
    require(int(head.headers.get("Content-Length","0"))==probe_length
            and head.headers.get("ETag")==expected_etag
            and head.headers.get("Last-Modified")==entry["last_modified"],
            "Official monthly open-data archive changed since the source probe")
    response=client.get(url,headers={"If-Match":expected_etag,
                                     "Accept-Encoding":"identity"},
                        stream=True,timeout=(40,180))
    response.raise_for_status()
    require(response.status_code==200 and
            response.headers.get("ETag")==expected_etag,
            "Source monthly object was switched between HEAD and GET")
    response.raw.decode_content=False
    return response,DigestReader(response.raw)


def scan_table(client:requests.Session,entry:dict,
               columns:tuple[str,...], callback, progress_name:str) -> dict:
    response,reader=provider_snapshot_stream(client,entry)
    total=0
    started=time.monotonic()
    try:
        with gzip.GzipFile(fileobj=reader,mode="rb") as gz:
            header=gz.readline().rstrip(b"\r\n").split(b"\t")
            require(len(set(header))==len(header) and
                    all(s.encode("ascii") in header for s in columns),
                    f"Unexpected official S3 {progress_name} metadata schema: {header[:20]}")
            indices={s:header.index(s.encode("ascii")) for s in columns}
            last=max(indices.values())
            for line in gz:
                total+=1
                fields=line.rstrip(b"\r\n").split(b"\t",last)
                require(len(fields)>last,
                        f"Missing required official metadata column on {progress_name} row {total}")
                callback(fields,indices)
                if total%PROGRESS_N==0:
                    print("MOIS_PROVIDER_OPEN_DATA_STREAM_PROGRESS",
                          json.dumps({"file":progress_name,"rows_seen":total,
                                      "compressed_gib_read":round(reader.bytes/1024**3,3),
                                      "elapsed_seconds":round(time.monotonic()-started,1)}),
                          flush=True)
    finally:
        response.close()
    require(reader.bytes==int(entry["content_length_bytes"]),
            "Official S3 compressed metadata length mismatches pinned source")
    return {"uncompressed_source_table_rows_seen":total,
            "compressed_bytes_read":reader.bytes,
            "compressed_original_SHA256":reader.sha256.hexdigest(),
            "etag":entry["etag"],"last_modified":entry["last_modified"],
            "source_verified_by_full_GZIP_CRC":True}


def valid_public_date(v:bytes,as_of:date)->str|None:
    try:
        s=v.decode("ascii")
        require(len(s)==10 and s[4]=="-" and s[7]=="-",
                "Unexpected provider observed_on date format")
        parsed=date.fromisoformat(s)
        return s if date(1800,1,1)<=parsed<=as_of else None
    except (UnicodeDecodeError,ValueError):
        return None


def scan_metadata(source:pd.DataFrame,probe:dict,client:requests.Session,
                  as_of:date) -> tuple[dict[str,dict],dict]:
    targets={str(p).encode("ascii") for p in source.photo_id}
    matched_by_photo={}
    uuid_to_photo={}
    duplicated_photo=0
    def one_photo(fields,at):
        nonlocal duplicated_photo
        pid=fields[at["photo_id"]]
        if pid not in targets:return
        uid=fields[at["observation_uuid"]]
        require(len(uid)>0, "Original photo maps to blank official observation UUID")
        if pid in matched_by_photo:
            duplicated_photo+=1
            require(matched_by_photo[pid]==uid,
                    "Original source photo maps to conflicting observation UUIDs")
        else:
            matched_by_photo[pid]=uid
            require(uid not in uuid_to_photo,
                    "Distinct historical original observations map to one provider observation UUID")
            uuid_to_photo[uid]=pid
    photo_receipt=scan_table(client,probe["source"]["photo_table"],
        ("photo_id","observation_uuid"),one_photo,"photos.csv.gz")
    print("MOIS_S3_PHOTO_UUID_JOIN",
          json.dumps({"original_historical_photos":len(targets),
              "source_original_photos_found_in_open_data":len(matched_by_photo),
              "duplicated_photo_rows":duplicated_photo}),flush=True)
    recovered={}
    def one_observation(fields,at):
        uid=fields[at["observation_uuid"]]
        pid=uuid_to_photo.get(uid)
        if pid is None:return
        require(pid not in recovered,
                "Official observation UUID repeated; cannot assign unique source date")
        today=valid_public_date(fields[at["observed_on"]],as_of)
        recovered[pid]={
            "current_public_snapshot_observation_uuid":uid.decode("ascii"),
            "current_public_snapshot_observed_on":today,
            "current_public_snapshot_taxon_id":fields[at["taxon_id"]].decode("ascii"),
            "current_public_quality_grade":fields[at["quality_grade"]].decode("ascii"),
        }
    obs_receipt=scan_table(client,probe["source"]["observation_table"],
        ("observation_uuid","observed_on","taxon_id","quality_grade"),
        one_observation,"observations.csv.gz")
    return {p.decode("ascii"):v for p,v in recovered.items()},{
        "photos":photo_receipt,"observations":obs_receipt,
        "photo_IDs_found_in_open_snapshot":len(matched_by_photo),
        "duplicated_photo_rows":duplicated_photo,
        "matched_observation_uuid_found_in_metadata":len(recovered),
    }


def read_prior(prior_zip:Path) -> pd.DataFrame:
    with zipfile.ZipFile(prior_zip) as z:
        csv_bytes=z.read(OLD_CSV)
        receipt=json.loads(z.read(OLD_RECEIPT))
    require(hashlib.sha256(csv_bytes).hexdigest()==OLD_SHA
            and receipt["date_recovery_photo_csv_sha256"]==OLD_SHA
            and receipt["n_original_photo_ID_matched_and_current_observation_date_present"]==200,
            "Source independently verified 200-photo date checks changed")
    return pd.read_csv(io.BytesIO(csv_bytes),dtype=str)


def assemble(historical:pd.DataFrame,all_photos:pd.DataFrame,
             snapshot:dict[str,dict], old200:pd.DataFrame, output:Path,
             probe:dict, full_stream:dict) -> dict:
    require(len(historical)==HISTORIC and historical.photo_id.is_unique
            and historical.observed_on.isna().all(),
            "Original 100543-photo historical unit changed")
    rows=[]
    for source in historical.to_dict(orient="records"):
        pid=str(source["photo_id"])
        meta=snapshot.get(pid)
        if meta is None:
            reason="no_complete_open_photo_to_observation_UUID_record"
            day=None
            taxon_now=None
            uuid=None
            quality=None
        else:
            day=meta["current_public_snapshot_observed_on"]
            reason=("matched_photo_UUID_observation_and_public_date"
                    if day is not None else "matched_original_photo_but_no_valid_observed_on")
            uuid=meta["current_public_snapshot_observation_uuid"]
            taxon_now=meta["current_public_snapshot_taxon_id"]
            quality=meta["current_public_quality_grade"]
        rows.append({
            "original_photo_id":pid,
            "original_numeric_observation_id":str(source["observation_id"]),
            "original_taxon_id":str(source["inat_taxon_id"]),
            "original_frozen_ROI_morph":source["morph"],
            "original_latitude":source["latitude"],
            "original_longitude":source["longitude"],
            "snapshot_observation_uuid":uuid,
            "snapshot_observed_on":day,
            "snapshot_taxon_id":taxon_now,
            "snapshot_quality_grade":quality,
            "now_taxon_matches_original":
                (None if taxon_now is None else taxon_now==str(source["inat_taxon_id"])),
            "source_photo_date_status":reason,
        })
    enriched=pd.DataFrame(rows)
    require(enriched.original_photo_id.is_unique and len(enriched)==HISTORIC,
            "Public snapshot exact photo ID join changed source denominator")
    recovered=enriched.loc[enriched.source_photo_date_status.eq(
        "matched_photo_UUID_observation_and_public_date"),[
            "original_photo_id","snapshot_observed_on"]].rename(
                columns={"original_photo_id":"photo_id",
                         "snapshot_observed_on":"snapshot_public_observation_date"})
    q=all_photos.merge(recovered,on="photo_id",how="left",validate="one_to_one")
    require(len(q)==MOIS_ROWS and
            not (q.observed_on.notna() & q.snapshot_public_observation_date.notna()).any(),
            "S3 date would override an original source observed_on")
    q["effective_observed_on"]=q.observed_on.where(
        q.observed_on.notna(),q.snapshot_public_observation_date)
    frame=q[["photo_id","observation_id","inat_taxon_id","morph",
             "latitude","longitude","old_mois_or_new_fcp","effective_observed_on"]].rename(
        columns={"effective_observed_on":"observed_on"})
    dated=labelled_calendar(frame)
    before=int((all_photos.observed_on.notna() & all_photos.morph.isin(HUES)).sum())
    new_hue_count=int(enriched.loc[
        enriched.snapshot_observed_on.notna(),"original_frozen_ROI_morph"].isin(HUES).sum())
    after=int((dated.date.notna() & dated.classified).sum())
    require(before==81867 and after==before+new_hue_count,
            "Full original source photo dated hue count does not reconcile")
    group=[]
    a=photo_season_table(dated)
    a["source_group"]="all_orig_267425_plus_public_S3_historical_date_only"
    group.append(a)
    historic_table=dated.loc[dated.old_mois_or_new_fcp.eq("historical")]
    a=photo_season_table(historic_table)
    a["source_group"]="historic_source_100543_open_S3_photo_dates"
    group.append(a)
    result_table=pd.concat(group,ignore_index=True)
    csv_path=output/SEASON_TABLE
    result_table.to_csv(csv_path,index=False,lineterminator="\n",float_format="%.8f")
    target_path=output/FINAL_TABLE
    enriched.to_csv(target_path,index=False,compression={"method":"gzip","mtime":0},
                    lineterminator="\n",float_format="%.8f")
    oldcheck=old200.merge(enriched,left_on="original_photo_id",
          right_on="original_photo_id",validate="one_to_one",how="left")
    assert len(oldcheck)==200
    previous_exact=oldcheck.retrieved_observed_on.notna()
    current_s3=oldcheck.snapshot_observed_on.notna()
    comparisons=int((previous_exact&current_s3).sum())
    disagreements=int(((oldcheck.retrieved_observed_on!=oldcheck.snapshot_observed_on)
         &(previous_exact&current_s3)).sum())
    source_s3_taxon_flags=enriched.now_taxon_matches_original.astype(str).str.lower()
    status={str(k):int(v) for k,v in enriched.source_photo_date_status.value_counts().items()}
    historical_source_hue=enriched.original_frozen_ROI_morph.isin(HUES)
    southern={}
    for label,sub in [("north",historic_table[historic_table.zone=="north_extratropical"]),
                      ("south",historic_table[historic_table.zone=="south_extratropical"])]:
        eligible=sub[sub.classified&sub.season.isin(["summer","winter"])]
        if len(eligible) and eligible.season.nunique()==2:
            counts=eligible.groupby("season").white.agg(["size","sum"])
            southern[label]={
                "summer_original_classified_photos":int(counts.loc["summer","size"]),
                "winter_original_classified_photos":int(counts.loc["winter","size"]),
                "raw_source_photo_white_summer_minus_winter_pp":round(100*(
                    counts.loc["summer","sum"]/counts.loc["summer","size"]-
                    counts.loc["winter","sum"]/counts.loc["winter","size"]),6),
            }
    paired={}
    for name,panel in (("historical",historic_table),("combined",dated)):
        q=panel[panel.classified&panel.date.notna()&
             panel.zone.isin(("north_extratropical","south_extratropical"))&
             panel.season.isin(("summer","winter"))&
             panel.latitude.between(-90,90)&panel.longitude.between(-180,180)]
        paired[name]={}
        for zone in ("north_extratropical","south_extratropical"):
            sub=q[q.zone.eq(zone)]
            paired[name][zone]=eligible_site_photo_support(
                sub,0.25,True,2)
    return {
        "full_historical_attempted_source_photo_count":HISTORIC,
        "full_historical_classifiable_original_hues":int(historical_source_hue.sum()),
        "photo_ids_with_verified_snapshot_observation_date":int(enriched.snapshot_observed_on.notna().sum()),
        "photo_ids_absent_or_date_invalid":int(enriched.snapshot_observed_on.isna().sum()),
        "white_original_photo_hues_with_verified_snapshot_date":int(
            (enriched.original_frozen_ROI_morph.eq("white") &
             enriched.snapshot_observed_on.notna()).sum()),
        "classifiable_original_hues_with_verified_snapshot_date":new_hue_count,
        "original_before_classifiable_dated":before,
        "new_after_classifiable_dated":after,
        "original_historical_source_taxonomy_matches_snapshot":int(
            source_s3_taxon_flags.eq("true").sum()),
        "snapshot_current_taxon_differs_from_original":int(source_s3_taxon_flags.eq("false").sum()),
        "original_200_current_API_photo_IDs_with_monthly_snapshot_dates":comparisons,
        "original_200_date_disagreements_between_current_API_and_20260927_snapshot":disagreements,
        "original_photo_current_open_metadata_identity_status":status,
        "observed_summer_winter_from_historical_original_photo_source":southern,
        "same_taxon_coarse_025_cell_same_season_year_summer_winter_min2_observation_ids":paired,
        "original_source_photos_and_four_hues_not_modified":True,
        "no_inferred_month_for_unmatched_or_technical_missing_photo":True,
        "recovered_historical_photo_table_SHA256":hashlib.sha256(target_path.read_bytes()).hexdigest(),
        "recovered_photo_season_source_table_SHA256":hashlib.sha256(csv_path.read_bytes()).hexdigest(),
    }


def run(photo_zip:Path,probe_path:Path,prior200_zip:Path,output:Path):
    require(not output.exists() or not any(output.iterdir()),
            "Refuse to overwrite a prior official source result tree")
    all_photos,prior,sha=read_signed(photo_zip,PHOTO_DATA,PHOTO_RECEIPT,
       "source_photo_level_csv_sha256",
       ["photo_id","observation_id","inat_taxon_id","morph","latitude","longitude",
        "observed_on","old_mois_or_new_fcp"])
    require(sha==ORIGINAL_SHA and len(all_photos)==MOIS_ROWS and
            all_photos.morph.isin(HUES).sum()==MOIS_HUES,
            "Genuine 267425 original photo union checksum or hue count changed")
    historical=all_photos[all_photos.old_mois_or_new_fcp.eq("historical")]
    require(len(historical)==HISTORIC and historical.observed_on.isna().all(),
            "Historical 100543 original photo date frame differs")
    probe=json.loads(probe_path.read_text())
    require(probe["schema"]=="mois_official_monthly_inat_open_data_metadata_source_probe_v1"
            and len(probe["source"])==2 and probe["estimated_input_transfer_GiB"]<50,
            "Source monthly provider ETags not locked or exceptionally large")
    prior200=read_prior(prior200_zip)
    output.mkdir(parents=True,exist_ok=True)
    (output/"full_source_scan_not_completed.json").write_text(json.dumps({
        "status":"INCOMPLETE_DO_NOT_REPORT_AS_AN_ECOLOGICAL_RESULT",
        "original_source_photo_CSV_sha256":sha,
        "provider_probe_SHA256":hashlib.sha256(probe_path.read_bytes()).hexdigest(),
        "official_monthly_photos_metadata_etag":probe["source"]["photo_table"]["etag"],
        "official_monthly_observations_metadata_etag":probe["source"]["observation_table"]["etag"],
        "requested_data_transfer_GiB":probe["estimated_input_transfer_GiB"],
        "no_photo_pixels_reclassified":True,
    },indent=2)+"\n")
    with requests.Session() as client:
        current,source_scan=scan_metadata(historical,probe,client,datetime.now(timezone.utc).date())
    analysis=assemble(historical,all_photos,current,prior200,output,probe,source_scan)
    result={
        "schema":"mois_historical_full_100543_exact_photo_to_S3_uuid_observed_date_v1",
        "source_photo_sha256":sha,
        "source_photo_rows":len(all_photos),
        "provider_probe_sha256":hashlib.sha256(probe_path.read_bytes()).hexdigest(),
        "provider_monthly_snapshot_source":probe["source"],
        "provider_GZIP_files_original_SHA256_and_full_row_counts":source_scan,
        "provider_snapshot_asof_current_UTC_date":datetime.now(timezone.utc).date().isoformat(),
        "original_numeric_observation_id_in_S3_metadata":False,
        "exact_original_photo_id_current_photo_to_observation_UUID_match":True,
        "retrospective_original_numeric_observation_id_verified_for_entire_100543":False,
        "past_200_original_observation_IDs_verified_independently_by_public_API":True,
        "historical_current_public_snapshot_observed_on_not_necessarily_archived_upload_date":True,
        "no_open_flower_petal_identity_verified":True,
        "no_observed_year_weather_or_causal_pigmentation_inference":True,
        "no_FCP_genetic_or_W_D_within_species_estimator":True,
        **analysis
    }
    (output/RECEIPT).write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    (output/"full_source_scan_not_completed.json").unlink()
    print("MOIS_FULL_HISTORICAL_OPEN_DATA_DATE_REJOIN",
          json.dumps({k:result[k] for k in (
              "photo_ids_with_verified_snapshot_observation_date",
              "classifiable_original_hues_with_verified_snapshot_date",
              "new_after_classifiable_dated",
              "original_200_date_disagreements_between_current_API_and_20260927_snapshot")},sort_keys=True),flush=True)
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-original-photo-zip",type=Path,required=True)
    p.add_argument("--provider-probe-receipt",type=Path,required=True)
    p.add_argument("--prior200-original-photo-date-zip",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    run(a.source_original_photo_zip,a.provider_probe_receipt,
        a.prior200_original_photo_date_zip,a.output)
