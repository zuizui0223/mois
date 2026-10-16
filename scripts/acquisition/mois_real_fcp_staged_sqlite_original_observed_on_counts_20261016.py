#!/usr/bin/env python3
"""Count ALL source ORIGINAL observed_on dates straight from immutable FCP SQLite.

This is a fast scientific checkpoint AFTER FCP's single durable four-table
provider original-source ingest has completed but BEFORE source 128 per-taxon
photo-ID manifests have been emitted. It does NOT download image pixels and
does NOT require any flower-colour/ROI or iNaturalist live API records.

Exactly FCP's globally deduplicated licensed, research-grade, original
2026-09-27 source SQLite: photos(photo_id,uuid,taxon,ambiguous),
observations(uuid,observed_on,lat,...). Count every ambiguous=0 photo
opportunity and classify its original observation date with Python's strict
Gregorian date parser and hemisphere, including missing/invalid dates.
The result gives the first true ALL-photo observed-date/season sampling
opportunity histogram as soon as the permanent staged source exists.

It is NOT proof of actual blooming flowers, and not yet available when only
FCP's independently verified PHOTO COUNT ledger is present. On production,
require the true complete FCP staged SQLite original SHA, exact source archive
SHA, taxon ledger SHA, and all 93,209,351 photo opportunities. Read-only
mode is immutable; never re-index or alter the original FCP source.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sqlite3
from collections import Counter,defaultdict
from pathlib import Path
from urllib.parse import quote

from scripts.acquisition.mois_full_source_dated_flower_candidates_20261016 import (
    source_observation_day,climatic_season,load_counts,sha,
    FRAME_SHA,COUNTS_SHA,FULL_TAXA,FULL_PHOTO_OPPORTUNITIES,
    FULL_WORK_PARTITIONS,SNAPSHOT
)

ARCHIVE_SHA="bfc1a97c436bf5989fccd43f12392953128489b76b92cea24da335ad73da4d8e"
STAGE_SCHEMA="fcp_staged_sqlite_photo_universe_v1"
STAGE_STATUS="SOURCE_SQLITE_STAGED_SHARDS_NOT_YET_EMITTED"
SCHEMA="mois_real_fcp_staged_complete_42111_photo_observed_on_original_climate_month_v1"
STATUS="REAL_FULL_42111_SOURCE_PHOTO_DATE_OPPORTUNITIES_NOT_REAL_FLOWERING_COLOURS"
DATE_STATES=(
    "ORIGINAL_SOURCE_DATE_KNOWN","SOURCE_DATE_NOT_VALID_ISO_DAY",
    "SOURCE_DATE_IMPOSSIBLE_CALENDAR_DAY","SOURCE_DATE_AFTER_FROZEN_SNAPSHOT"
)
PRODUCTION_MONTH_ZONE=("NORTH_EXTRATROPICAL","SOUTH_EXTRATROPICAL","TROPICAL")


def source_stage_gate(stage_root:Path,counts:Path,*,synthetic=False):
    receipt_path=stage_root/"fcp_staged_sqlite_source_receipt.json"
    if not receipt_path.is_file():
        raise ValueError("FCP four-table permanent staged source not yet ready")
    receipt=json.loads(receipt_path.read_text(encoding="utf8"))
    sql_file=stage_root/"_working_fcp_snapshot.sqlite"
    if (receipt.get("schema")!=STAGE_SCHEMA
        or receipt.get("status")!=STAGE_STATUS
        or receipt.get("sqlite_basename")!=sql_file.name
        or receipt.get("sqlite_bytes")!=sql_file.stat().st_size
        or receipt.get("sqlite_sha256")!=sha(sql_file)
        or receipt.get("historical_flower_colour_csv_used") is not False
        or receipt.get("source_image_pixels_opened") is not False
        or receipt.get("original_photo_ids_deduplicated_globally") is not True
        or receipt.get("full_source_census_complete") is not False
        or receipt.get("full_study_original_images_measured") is not False
        or receipt.get("source_snapshot_date")!=SNAPSHOT):
        raise ValueError("Original FCP full-provider staged metadata source not byte authenticated")
    if not synthetic:
        verified_archive=receipt.get("verified_original_dated_archive",{})
        if (receipt.get("source_species")!=FULL_TAXA
            or receipt.get("source_species_frame_sha256")!=FRAME_SHA
            or receipt.get("n_target_taxon_licensed_photos")!=FULL_PHOTO_OPPORTUNITIES
            or receipt.get("shards")!=FULL_WORK_PARTITIONS
            or verified_archive.get("source_archive_sha256")!=ARCHIVE_SHA
            or sha(counts)!=COUNTS_SHA):
            raise ValueError("Not genuine complete 42,111-original-taxa and 93,209,351 unique licensed source")
    baseline=load_counts(counts,synthetic=synthetic)
    if len(baseline)!=receipt["source_species"] or (
        sum(n for _,n in baseline.values())!=receipt["n_target_taxon_licensed_photos"]):
        raise ValueError("Source taxonomy/count ledger does not match original photo table")
    return receipt,sql_file,baseline


def audit_staged_original_dates(stage_root:Path,counts:Path,out:Path,*,
                                synthetic=False,min_public=10)->dict:
    if min_public<10 or out.exists() or out.with_name(out.name+".INCOMPLETE").exists():
        raise ValueError("Original full source date counts cannot be overwritten or publicly under-suppressed")
    source,sql_file,baseline=source_stage_gate(stage_root,counts,synthetic=synthetic)
    raw_count=valid=0
    date_states=Counter()
    domain_month=Counter()
    domain_local_season=Counter()
    annual_counts=Counter()
    taxon_date_presence=set()
    taxa_by_month=defaultdict(set)
    by_taxon=Counter()
    # URI read-only+immutable avoids SQLite journals and accidental writes to
    # frozen FCP source. Source source SHA was independently validated above.
    uri=sql_file.as_uri()+"?mode=ro&immutable=1"
    with sqlite3.connect(uri,uri=True) as db:
        db.execute("PRAGMA query_only=ON")
        actual_source=db.execute(
            "SELECT COUNT(*) FROM photos WHERE ambiguous=0").fetchone()[0]
        if actual_source!=source["n_target_taxon_licensed_photos"]:
            raise ValueError("Frozen original image ID source photo-count drift")
        # This indexed UUID join can stream at bounded memory rather than
        # loading all 93M PHOTO records into RAM. It may take hours and MUST
        # run on the operator's durable host containing the original SQLite.
        rows=db.execute("""SELECT p.taxon,o.observed_on,o.lat
            FROM photos AS p JOIN observations AS o ON p.uuid=o.uuid
            WHERE p.ambiguous=0""")
        for tid,observed,lat in rows:
            if type(tid) is not int or tid not in baseline:
                raise ValueError("Source photo taxon not in independent original species frame")
            raw_count+=1
            by_taxon[tid]+=1
            day,status=source_observation_day(observed)
            if status not in DATE_STATES:
                raise ValueError("Unknown original snapshot observation-date state")
            date_states[status]+=1
            if day is None:continue
            if lat is None or not math.isfinite(float(lat)) or not -90<=float(lat)<=90:
                raise ValueError("FCP previously source-screened original photo lacks valid latitude")
            zone,season=climatic_season(float(lat),day.month)
            valid+=1
            taxon_date_presence.add(tid)
            domain_month[(zone,day.month)]+=1
            domain_local_season[(zone,season)]+=1
            annual_counts[day.year]+=1
            taxa_by_month[(zone,day.month)].add(tid)
    if (raw_count!=actual_source or sum(date_states.values())!=raw_count
        or valid!=sum(domain_month.values()) or valid!=sum(domain_local_season.values())
        or valid!=sum(annual_counts.values())
        or set(by_taxon)!=set(t for t,(_,n) in baseline.items() if n>0)
        or any(by_taxon.get(t,0)!=expected_n for t,(_,expected_n) in baseline.items())):
        raise ValueError("All original full source photo/taxa/date group conservations failed")
    month_rows=[]
    for (zone,month),n in sorted(domain_month.items()):
        month_rows.append({
            "latitude_domain":zone,"original_calendar_month":month,
            **({"n_date_valid_original_photo_opportunities":n,
                "n_original_named_taxa_with_valid_date_photo":len(taxa_by_month[(zone,month)])}
               if n>=min_public else {"public_suppressed_below":min_public})
        })
    season_rows=[]
    for (zone,season),n in sorted(domain_local_season.items()):
        season_rows.append({
            "latitude_domain":zone,"source_local_season":season,
            **({"n_date_valid_original_photo_opportunities":n}
               if n>=min_public else {"public_suppressed_below":min_public})
        })
    result={
        "schema":SCHEMA,"status":STATUS,
        "source_original_fcp_staged_sqlite_receipt_sha256":
            sha(stage_root/"fcp_staged_sqlite_source_receipt.json"),
        "source_original_sqlite_sha256":source["sqlite_sha256"],
        "source_original_snapshot_taxa_frame_sha256":source["source_species_frame_sha256"],
        "source_complete_42111_photo_count_ledger_sha256":sha(counts),
        "source_original_snapshot_date":SNAPSHOT,
        "n_original_source_named_taxa":len(baseline),
        "n_original_photo_opportunities_with_globally_unique_ids":raw_count,
        "n_original_snapshot_observed_on_date_known":valid,
        "n_original_snapshot_observed_on_date_unknown_or_invalid":raw_count-valid,
        "n_original_named_taxa_with_valid_observed_on_photos":len(taxon_date_presence),
        "n_source_taxa_with_no_valid_original_observation_dates":
            len(baseline)-len(taxon_date_presence),
        "original_observed_on_status_counts":dict(sorted(date_states.items())),
        "original_snapshot_photo_observation_year_counts":dict(sorted(annual_counts.items())),
        "original_latitude_domain_observed_month_photo_opportunity_counts":month_rows,
        "original_latitude_domain_local_season_photo_opportunities":season_rows,
        "all_original_photo_ids_inclusion_cap":None,
        "source_original_date_not_created_or_uploaded_at":True,
        "all_source_photo_ids_permanently_manifested_in_this_report":False,
        "all_original_images_opened_or_flowering_roi_confirmed":False,
        "n_classifiable_flowering_photographs_estimated_here":0,
        "real_original_photo_white_fraction_estimated_here":False,
        "historic_100543_photo_colour_outcomes_reused":False,
        "snapshot_photo_count_denominator_not_wild_blooming_individuals":True,
        "synthetic_test_only":bool(synthetic),
        "note":"True original photo-ID opportunities with original observed_on, not confirmed flowering; only original ROI-v4 can establish classifiable white, and source all-photo ID manifests still require FCP 128-shard freeze."
    }
    out.parent.mkdir(parents=True,exist_ok=True)
    temp=out.with_name(out.name+".INCOMPLETE")
    temp.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf8")
    temp.replace(out)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fcp-complete-original-staged-sqlite-root",type=Path,required=True)
    p.add_argument("--independent-42111-original-photo-count-csv",type=Path,required=True)
    p.add_argument("--out-public-original-observed-on-photo-opportunities",type=Path,required=True)
    a=p.parse_args()
    r=audit_staged_original_dates(
        a.fcp_complete_original_staged_sqlite_root,
        a.independent_42111_original_photo_count_csv,
        a.out_public_original_observed_on_photo_opportunities)
    print("MOIS_ACTUAL_FULL_ORIGINAL_PHOTO_OBSERVED_ON_ONLY_NOT_FLOWER_PIXEL",
          json.dumps({k:r[k] for k in (
              "status","n_original_photo_opportunities_with_globally_unique_ids",
              "n_original_snapshot_observed_on_date_known",
              "n_original_snapshot_observed_on_date_unknown_or_invalid"
          )},sort_keys=True))


if __name__=="__main__":main()
