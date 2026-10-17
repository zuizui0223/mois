#!/usr/bin/env python3
"""MOIS: source-exact original photo day reconciliation, 302242 pixels/146052 hues.

Previously successful real-photo sources:
 A. 302242 measured original pixels, 146061 frozen ROI-v4 hue outcomes, incl.
    218 genuinely completed FCP wave4 partitions. 143770 initially date+hue.
 B. 1732 extra already verified public numeric observation/photo date records
    on the first 295017 original photo rows, including 706 classifiable ROIs.
 C. Original *full* 3430502-source FCP 2026-09-27 photo/observation UUID/date
    CSV; scan ALL source rows and exact-match the 4096 existing date-missing
    MOIS wave1 measured original photo IDs and source UUIDs/taxa/coordinates.
    There are 1576 already measured hues among these 4096 dated photos.

Never classify new photo pixels here, invent day, overwrite existing historical
S3/API dates, import FCP within-species W/D fitness/50km results or treat an
automatically classified ROI as a confirmed open flower.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
from datetime import date
from pathlib import Path
import zipfile

import pandas as pd

from scripts.analysis.mois_seasonal_267k_climate_normals import (
    HUES, labelled_calendar, require,
)
from scripts.analysis.mois_original_photo_community_season_geoscale import (
    matched_community, cell_cluster_interval,
)

BASE_N=302242
BASE_ROI=146061
BASE_DATED=296395
BASE_DATED_ROI=143770
FIRST_N=295017
NUMERIC_NEW=1732
NUMERIC_NEW_HUE=706
UUID_MISSING=4096
UUID_NEW_HUE=1576
FINAL_DATED=302223
FINAL_DATED_ROI=146052
TOTAL_ORIGINAL_SOURCE_ROWS=3430502
SOURCE_ARCHIVE_DAY="2026-09-27"
FINAL_SOURCE_MAX_DAY=date(2026,10,17)
BASE_CSV="mois_next_finished_fcp_original_photo_date_hues_plus295017.csv.gz"
BASE_RECEIPT="mois_next_finished_fcp_original_photo_new_pixel_date_hue_receipt.json"
BASE_SHA="b83871c6cd69cb624693364742c2a1d61003253d084af0cb6d80ab9f898c3959"
NUMERIC_CSV="mois_real295017_photo_source_plus_exact_remaining_numeric_API_dates.csv.gz"
NUMERIC_RECEIPT="mois_original295017_remaining1751_numeric_photo_date_recovery_QC_receipt.json"
NUMERIC_SHA="65c072a97f7a59b1933dbc0541041a29069d468bf6fc889933d5a922c743a11a"
SOURCE_FULL="rgfca_1000_exact_source_photo_IDs/rgfca_1000_complete_original_source_photo_IDs.csv"
SOURCE_RECEIPT="rgfca_1000_exact_source_photo_IDs/rgfca_1000_targeted_monthly_photo_ID_receipt.json"
SOURCE_SHA="58e037a210e23d3ff5d6f04d392665ada5ee9ac5149d0676672c239cf092df01"
OUT_CSV="mois_original302242_exact_observed_on_full_uuid_numeric_verified.csv.gz"
OUT_UUID="mois_original4096_photo_UUID_observation_day_20260927_exact.csv"
OUT_SEASON="mois_original146052_dated_hue_full_geoseason_original_photo_composition.csv"
OUT_RECEIPT="mois_original302242_146052_photo_hue_full_observed_on_provenance_receipt.json"
COLOUR=("white","yellow_orange","red_pink","blue_purple")


def sha(raw:bytes)->str:
    return hashlib.sha256(raw).hexdigest()


def load_already_signed(base_zip:Path, numeric_zip:Path)->tuple[pd.DataFrame,pd.DataFrame,dict]:
    with zipfile.ZipFile(base_zip) as z:
        raw=z.read(BASE_CSV)
        receipt=json.loads(z.read(BASE_RECEIPT))
    require(sha(raw)==BASE_SHA and
            receipt["expanded_real_date_hue_photo_source_SHA256"]==BASE_SHA and
            receipt["expanded_distinct_measured_source_photo_IDs"]==BASE_N and
            receipt["expanded_four_hue_source_ROIs"]==BASE_ROI and
            receipt["expanded_dated_original_four_hue_ROIs"]==BASE_DATED_ROI and
            receipt["FCP_real_completed_photo_ROI_partitions"]==218,
            "Actual 302242 completed-pixel photo source SHA/ROI/worker terminal changed")
    base=pd.read_csv(io.BytesIO(raw),compression="gzip",dtype=str,keep_default_na=False)
    require(len(base)==BASE_N and base.photo_id.is_unique and
            base.morph.isin(HUES).sum()==BASE_ROI and
            base.qualified_date.str.len().eq(10).sum()==BASE_DATED and
            (base.morph.isin(HUES)&base.qualified_date.str.len().eq(10)).sum()==BASE_DATED_ROI,
            "Baseline authenticated 302242 photo/date/hue counts drifted")

    with zipfile.ZipFile(numeric_zip) as z:
        raw=z.read(NUMERIC_CSV)
        numeric_receipt=json.loads(z.read(NUMERIC_RECEIPT))
    require(sha(raw)==NUMERIC_SHA and
            numeric_receipt["full_295017_source_plus_verified_day_csv_SHA256"]==NUMERIC_SHA and
            numeric_receipt["newly_verified_source_photo_numeric_original_obs_dates"]==NUMERIC_NEW and
            numeric_receipt["newly_dated_classifiable_four_hue_original_source_photos"]==NUMERIC_NEW_HUE,
            "Independent current original-photo numeric iNaturalist observation+day source changed")
    numeric=pd.read_csv(io.BytesIO(raw),compression="gzip",dtype=str,keep_default_na=False)
    require(len(numeric)==FIRST_N and numeric.photo_id.is_unique and
            base.iloc[:FIRST_N].photo_id.tolist()==numeric.photo_id.tolist(),
            "The original 295017 photo rows have been reordered/replaced after authentic image append")
    return base,numeric,{"base_302242_source_SHA256":BASE_SHA,
                          "numeric_1732_public_date_source_SHA256":NUMERIC_SHA}


def numeric_patch(base:pd.DataFrame, numeric:pd.DataFrame,*,first_n:int=FIRST_N):
    require(len(base)>=first_n and len(numeric)==first_n and
            base.iloc[:first_n].photo_id.tolist()==numeric.photo_id.tolist(),
            "Photo ID current measurement and numeric source date rows are not aligned")
    fixed_columns=[c for c in base.columns if c not in
                   ("qualified_date","original_source_date_provenance")]
    require(base.iloc[:first_n][fixed_columns].reset_index(drop=True).equals(
            numeric[fixed_columns].reset_index(drop=True)),
            "Current numeric date provider changed existing original pixels/taxon/hue/location")
    old=base.iloc[:first_n]
    changed=old.qualified_date.ne(numeric.qualified_date).to_numpy()
    idx=old.index[changed]
    require(len(idx)>0 and old.loc[idx,"qualified_date"].eq("").all() and
            numeric.loc[idx,"qualified_date"].str.fullmatch(r"\d{4}-\d{2}-\d{2}").all(),
            "Extra numeric observed-on would replace a previously recorded photo date")
    restored=int(numeric.loc[idx,"morph"].isin(HUES).sum())
    out=base.copy()
    out.loc[idx,"qualified_date"]=numeric.loc[idx,"qualified_date"].to_numpy()
    out.loc[idx,"original_source_date_provenance"]=numeric.loc[
        idx,"original_source_date_provenance"].to_numpy()
    return out,{"new_original_numeric_photo_observation_date_IDs":len(idx),
                "new_original_numeric_four_hue_date_photo_IDs":restored}


def missing_original_UUID_rows(d:pd.DataFrame)->pd.DataFrame:
    frame=d[d.qualified_date.eq("")&
            d.old_mois_or_new_fcp.eq("mois_new_original_ROI")].copy()
    require(len(frame)==UUID_MISSING and frame.photo_id.is_unique and
            frame.observation_id.is_unique and
            frame.observation_id.str.fullmatch(
                r"[a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{12}").all(),
            "Actual 4096 original missing-day photo source UUID population changed")
    require(int(frame.morph.isin(HUES).sum())==UUID_NEW_HUE,
            "Current source old 4096 genuinely measured ROI hue denominator changed")
    return frame


def valid_source_day(s:str)->bool:
    try:
        today=date.fromisoformat(s)
        return date(1800,1,1)<=today<=FINAL_SOURCE_MAX_DAY and today.isoformat()==s
    except (TypeError,ValueError):
        return False


def scan_all_full_FCP_source(source:Path, source_receipt:Path,
                             missing:pd.DataFrame) -> tuple[pd.DataFrame,dict]:
    receipt=json.loads(source_receipt.read_text())
    require(receipt["schema"]=="fcp_rgfca_targeted_original_20260927_photo_ID_census_v1"
            and receipt["selected_all_photo_IDs_csv_sha256"]==SOURCE_SHA
            and receipt["n_selected_original_photo_IDs"]==TOTAL_ORIGINAL_SOURCE_ROWS
            and receipt["source_snapshot_date"]==SOURCE_ARCHIVE_DAY,
            "Full original FCP source manifest changed")
    expected={r["photo_id"]:r for r in missing.to_dict(orient="records")}
    found={}
    digest=hashlib.sha256()
    row_count=0
    with source.open("rb") as raw:
        def lines():
            for line in raw:
                digest.update(line)
                yield line.decode("utf-8")
        r=csv.DictReader(lines())
        required={"photo_id","observation_id","inat_taxon_id",
                  "latitude","longitude","observed_on","source_snapshot_date"}
        require(required.issubset(set(r.fieldnames or [])),
                "Upstream official 3430502 photo ID metadata schema changed")
        for row in r:
            row_count+=1
            pid=row["photo_id"]
            if pid not in expected:
                continue
            require(pid not in found,"Duplicate target original photo ID in full 3430502-source census")
            old=expected[pid]
            require(row["observation_id"]==old["observation_id"]
                    and row["inat_taxon_id"]==old["inat_taxon_id"],
                    "Original photo ID matched wrong source observation UUID/taxon")
            require(abs(float(row["latitude"])-float(old["latitude"]))<0.00001 and
                    abs(float(row["longitude"])-float(old["longitude"]))<0.00001,
                    "Target originally selected photo locality changed")
            require(row["source_snapshot_date"]==SOURCE_ARCHIVE_DAY,
                    "Monthly original FCP snapshot day differs among original photographs")
            require(valid_source_day(row["observed_on"]),
                    "An original source 4096 UUID observation day missing or invalid")
            found[pid]={
                "original_photo_id":pid,
                "matched_source_original_observation_UUID":row["observation_id"],
                "matched_source_original_taxon_id":row["inat_taxon_id"],
                "verified_original_observed_on":row["observed_on"],
                "source_snapshot_date":row["source_snapshot_date"],
                "provenance":"exact_FCP_full_3430502_original_photo_observation_UUID_date",
            }
    require(row_count==TOTAL_ORIGINAL_SOURCE_ROWS and
            digest.hexdigest()==SOURCE_SHA and len(found)==len(expected),
            "Full original 3430502 source was incomplete, source changed or 4096 photos unmatched")
    matched=pd.DataFrame([found[pid] for pid in missing.photo_id])
    return matched,{
        "original_provider_source_photo_records_streamed":row_count,
        "original_3430502_source_full_CSV_SHA256":digest.hexdigest(),
        "original_photo_and_observation_UUID_and_taxon_and_locality_exact_matches":len(matched),
        "actual_original_photo_observation_dates_from_snapshot":len(matched),
        "original_source_snapshot":SOURCE_ARCHIVE_DAY,
    }


def apply_uuid_days(d:pd.DataFrame,uuid:pd.DataFrame)->tuple[pd.DataFrame,dict]:
    require(uuid.original_photo_id.is_unique and
            uuid.matched_source_original_observation_UUID.is_unique,
            "Previously undated original UUID photo/observation pair duplicated")
    require(uuid.verified_original_observed_on.map(valid_source_day).all(),
            "Original 4096 source UUID observed days changed")
    matched_ids=set(uuid.original_photo_id)
    rows=d[d.photo_id.isin(matched_ids)]
    require(len(rows)==len(uuid) and rows.qualified_date.eq("").all(),
            "Full FCP source date would change or overwrite an old confirmed day")
    cross=rows[["photo_id","observation_id","inat_taxon_id"]].merge(
        uuid,left_on="photo_id",right_on="original_photo_id",validate="one_to_one")
    require(cross.observation_id.eq(
              cross.matched_source_original_observation_UUID).all() and
            cross.inat_taxon_id.eq(cross.matched_source_original_taxon_id).all(),
            "Final date assembly breaks photo ID to original UUID/taxon identity")
    out=d.copy()
    dates=dict(zip(uuid.original_photo_id,uuid.verified_original_observed_on))
    mask=out.photo_id.isin(matched_ids)
    out.loc[mask,"qualified_date"]=out.loc[mask,"photo_id"].map(dates)
    out.loc[mask,"original_source_date_provenance"]=(
        "FCP_original3430502_20260927_photo_ID_and_observation_UUID_verified")
    return out,{
        "new_source_UUID_original_photo_day_IDs":len(rows),
        "new_UUID_source_classifiable_original_hues":int(rows.morph.isin(HUES).sum()),
    }


def run(base:Path,numeric:Path,full_csv:Path,full_receipt:Path,output:Path):
    require(not output.exists() or not any(output.iterdir()),
            "Never replace completed original pixels/photographic source-date receipts")
    original,ready,proof=load_already_signed(base,numeric)
    plus_numeric,ns=numeric_patch(original,ready)
    require(ns["new_original_numeric_photo_observation_date_IDs"]==NUMERIC_NEW and
            ns["new_original_numeric_four_hue_date_photo_IDs"]==NUMERIC_NEW_HUE,
            "Authenticated 1732 current original photo/observed date matches changed")
    missing=missing_original_UUID_rows(plus_numeric)
    source_uuid,scan=scan_all_full_FCP_source(full_csv,full_receipt,missing)
    complete,us=apply_uuid_days(plus_numeric,source_uuid)
    require(us["new_source_UUID_original_photo_day_IDs"]==UUID_MISSING and
            us["new_UUID_source_classifiable_original_hues"]==UUID_NEW_HUE and
            int(complete.photo_id.nunique())==len(complete)==BASE_N and
            int(complete.morph.isin(HUES).sum())==BASE_ROI and
            int(complete.qualified_date.str.len().eq(10).sum())==FINAL_DATED and
            int((complete.qualified_date.str.len().eq(10)&
                 complete.morph.isin(HUES)).sum())==FINAL_DATED_ROI,
            "No source original photo day/colour total may be counted without direct photo provenance")
    unchanged=[x for x in complete.columns if x not in
               ("qualified_date","original_source_date_provenance")]
    require(complete[unchanged].equals(original[unchanged]) and
            complete.photo_id.tolist()==original.photo_id.tolist(),
            "Original exact measured photograph/photo ROI taxonomy/pixel/location was modified")
    output.mkdir(parents=True,exist_ok=True)
    csvpath=output/OUT_CSV
    complete.to_csv(csvpath,index=False,compression={"method":"gzip","mtime":0},
                    lineterminator="\n")
    uuidpath=output/OUT_UUID
    source_uuid.to_csv(uuidpath,index=False,lineterminator="\n")
    panel=complete[["photo_id","observation_id","inat_taxon_id","morph","latitude",
                    "longitude","old_mois_or_new_fcp","qualified_date"]].rename(
                columns={"qualified_date":"observed_on"}).copy()
    panel["latitude"]=pd.to_numeric(panel.latitude,errors="coerce")
    panel["longitude"]=pd.to_numeric(panel.longitude,errors="coerce")
    panel=labelled_calendar(panel)
    require(int((panel.date.notna()&panel.classified).sum())==FINAL_DATED_ROI,
            "Original complete source season parser changes date+true measured hue count")
    table=[]
    local={}
    for zone in ("north_extratropical","south_extratropical"):
        q=panel[panel.zone.eq(zone)&panel.date.notna()]
        for season in ("winter","spring","summer","autumn"):
            t=q[q.season.eq(season)]
            c=t[t.classified]
            table.append({"hemisphere":zone,"season":season,
                          "actual_date_bearing_photo_attempts":len(t),
                          "actually_classifiable_original_photo_ROIs":len(c),
                          "observed_white_photo_ROIs":int(c.white.sum()),
                          "source_white_photo_fraction_given_classifiable":
                              round(float(c.white.mean()),8) if len(c) else None})
        matched=matched_community(panel,width=0.25,min_obs=5,min_taxa=2)
        one=matched[matched.zone.eq(zone)]
        local[zone]={
            "same_geocell_and_same_season_year_qualified_community_groups":int(len(one)),
            "source_actual_classifiable_photo_ROIs_in_matched_groups":
                int(one.summer_photos.sum()+one.winter_photos.sum()) if len(one) else 0,
            "summer_minus_winter_photo_white_pp_999_geocell_bootstrap":
                cell_cluster_interval(one,draws=999,seed=20261017) if len(one) else None,
        }
    season=output/OUT_SEASON
    pd.DataFrame(table).to_csv(season,index=False,lineterminator="\n")
    receipt={
        "schema":"mois_original302242_photo_146052_date_and_hue_all_exact_source_rejoin_v1",
        "original_302242_actual_measured_photo_pixels":BASE_N,
        "original_146061_automatically_classified_four_hue_ROIs":BASE_ROI,
        "previous_source_public_dated_original_photo_attempts":BASE_DATED,
        "previous_source_date_and_four_hue_ROIs":BASE_DATED_ROI,
        "newly_exact_original_numeric_observation_dates":NUMERIC_NEW,
        "newly_numeric_date_and_classified_source_ROIs":NUMERIC_NEW_HUE,
        "newly_FCP_full_source_photo_to_UUID_to_original_day":UUID_MISSING,
        "newly_UUID_source_date_and_original_four_hue_ROIs":UUID_NEW_HUE,
        "final_source_observation_dated_photo_attempts":FINAL_DATED,
        "final_dated_and_original_four_hue_classifiable_ROIs":FINAL_DATED_ROI,
        "remaining_original_photo_day_unknown":BASE_N-FINAL_DATED,
        "remaining_dated_colour_unknown":FINAL_DATED-FINAL_DATED_ROI,
        "remaining_classifiable_hue_photos_without_original_day":BASE_ROI-FINAL_DATED_ROI,
        "named_taxa_with_original_day_and_four_hue_ROI":
            int(panel.loc[panel.date.notna()&panel.classified,"inat_taxon_id"].nunique()),
        "genuine_pixel_baseline_and_numeric_date_signatures":proof,
        "authentic_original_full_FCP_photo_opportunity_source_scan":scan,
        "original_source_validated_photo_UUID_days_csv_SHA256":sha(uuidpath.read_bytes()),
        "final_source_data_csv_SHA256":sha(csvpath.read_bytes()),
        "source_season_group_table_SHA256":sha(season.read_bytes()),
        "source_0p25_degree_locality_seasonyear_matched_community":local,
        "all_302242_original_photo_taxon_ROI_pixel_coordinates_unchanged":True,
        "source_numeric_observation_date_and_99817_historic_S3_priority_preserved":True,
        "no_technical_unclassified_photo_assumed_chromatic_or_nonflowering":True,
        "no_photo_specific_open_flower_or_petal_qualification_done":True,
        "no_observed_year_temperature_solar_rain_or_VPD_response_fitted":True,
        "MOIS_global_taxon_colour_not_FCP_within_species_D_W_50km_result":True,
    }
    (output/OUT_RECEIPT).write_text(json.dumps(
        receipt,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("MOIS_AUTHENTIC_302242_DATE_REPAIR_ALL_OBSERVATION_ID_TYPES",
          json.dumps({k:receipt[k] for k in (
              "newly_exact_original_numeric_observation_dates",
              "newly_FCP_full_source_photo_to_UUID_to_original_day",
              "original_302242_actual_measured_photo_pixels",
              "final_source_observation_dated_photo_attempts",
              "final_dated_and_original_four_hue_classifiable_ROIs",
              "remaining_classifiable_hue_photos_without_original_day",
              "named_taxa_with_original_day_and_four_hue_ROI")},sort_keys=True),flush=True)
    return receipt


if __name__=="__main__":
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("--signed302242-pixel-and-date-archive",type=Path,required=True)
    a.add_argument("--numeric1732-photo-date-archive",type=Path,required=True)
    a.add_argument("--original-full3430502-photo-metadata-csv",type=Path,required=True)
    a.add_argument("--original-full3430502-source-receipt",type=Path,required=True)
    a.add_argument("--output",type=Path,required=True)
    z=a.parse_args()
    run(z.signed302242_pixel_and_date_archive,z.numeric1732_photo_date_archive,
        z.original_full3430502_photo_metadata_csv,z.original_full3430502_source_receipt,
        z.output)
