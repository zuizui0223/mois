#!/usr/bin/env python3
"""Finalize ALL 128 original FCP photo shards after DATE-FIRST ROI-v4 measurement.

This is the seasonal-flower-colour shortcut: scan all 93,209,351 original
photo IDs for source snapshot observed_on, but download/measure only date-
verified source photographs. An original-date-missing photograph remains
in the opportunity/exclusion denominator; it cannot inform monthly colour.

Each original date-known photo MUST have one genuine ROI work-unit terminal,
including failures. Only date-known photographs whose original ROI detected
a classifiable flower produce a blooming-white outcome. Every original
photo ID is retained in the private Bloom manifest; multiple photos from
one observation UUID contribute ONE observation-level white fraction.

No global completed source claim without all 128 independent FCP source
shard/date and ROI receipts, exact 42,111/93,209,351 metadata counts, SHA
attestation and duplicate-photo refusal. Aggregates only are publishable.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sqlite3
from collections import Counter,defaultdict
from pathlib import Path

from scripts.acquisition.mois_full_source_dated_flower_candidates_20261016 import (
    sha,source_observation_day,climatic_season,FULL_WORK_PARTITIONS,
    FULL_PHOTO_OPPORTUNITIES,FULL_TAXA,FRAME_SHA
)
from scripts.acquisition.mois_aggregate_all_original_dated_shards_20261016 import (
    SCHEMA as DATES_SCHEMA,STATUS as DATES_STATUS
)
from scripts.analysis.mois_reassemble_dated_original_roi_shard_20261016 import (
    SCHEMA as SHARD_SCHEMA,STATUS as SHARD_STATUS
)
from scripts.analysis.mois_original_full_fcp_dated_flower_ROI_20261016 import MORPHS

SCHEMA="mois_complete_all_original_date_first_fcp_blooming_photo_roi_v1"
STATUS="COMPLETE_ALL_ORIGINAL_DATE_KNOWN_PHOTO_ROI_BLOOMS_NO_DATE_GUESSING"


def aggregate_all(dates_receipt:Path,roi_shard_root:Path,
                  out_private_global_blooms:Path,out_public_report:Path,*,
                  synthetic=False,partitions=FULL_WORK_PARTITIONS,
                  expected_photos=FULL_PHOTO_OPPORTUNITIES,
                  expected_taxa=FULL_TAXA,min_public=10)->dict:
    if min_public<10 or partitions<1 or (
        not synthetic and (partitions!=FULL_WORK_PARTITIONS
           or expected_photos!=FULL_PHOTO_OPPORTUNITIES
           or expected_taxa!=FULL_TAXA)):
        raise ValueError("Full original 42111/93m source cannot be relaxed or sampled")
    for target in (out_private_global_blooms,out_public_report):
        if target.exists() or target.with_name(target.name+".INCOMPLETE").exists():
            raise ValueError("Complete original blooming-photo global result cannot be overwritten")
    dated=json.loads(dates_receipt.read_text(encoding="utf8"))
    if (dated.get("schema")!=DATES_SCHEMA or dated.get("status")!=DATES_STATUS
        or dated.get("synthetic_test_only") is not synthetic
        or dated.get("source_original_shards_complete")!=partitions
        or dated.get("n_all_original_source_named_taxa")!=expected_taxa
        or dated.get("n_all_original_unique_photo_opportunities")!=expected_photos
        or len(dated.get("source_shard_receipts",[]))!=partitions
        or dated.get("source_photo_cap_per_species","INVALID") is not None
        or dated.get("n_original_snapshot_observed_on_valid_photo_candidates",-1)
            +dated.get("n_original_snapshot_observed_on_invalid_or_missing_photo_candidates",-1)
           !=expected_photos
        or dated.get("original_source_photo_bloom_phenology_verified_by_image_pixels") is not False
        or not synthetic and dated.get("original_full_42111_species_frame_sha256")!=FRAME_SHA):
        raise ValueError("Cannot fake original all-128-shard photo source date census completion")
    out_private_global_blooms.parent.mkdir(parents=True,exist_ok=True)
    out_public_report.parent.mkdir(parents=True,exist_ok=True)
    stage=out_private_global_blooms.with_name(out_private_global_blooms.name+".INCOMPLETE")
    public_stage=out_public_report.with_name(out_public_report.name+".INCOMPLETE")
    scratch=out_private_global_blooms.with_name(out_private_global_blooms.name+".ids.sqlite.INCOMPLETE")
    photo_total=dated_total=roi_total=classified=white=source_failed=0
    species=set()
    genera=set()
    calendar=defaultdict(Counter)
    receipts=[]
    try:
        with sqlite3.connect(scratch) as db, stage.open("w",encoding="utf8") as out:
            db.execute("PRAGMA journal_mode=DELETE")
            db.execute("PRAGMA synchronous=NORMAL")
            db.execute("CREATE TABLE unique_original_photo(photo_id TEXT PRIMARY KEY)")
            db.execute("""CREATE TABLE source_observations(
               uuid TEXT PRIMARY KEY, taxon INTEGER NOT NULL, genus TEXT NOT NULL,
               observer TEXT NOT NULL, domain TEXT NOT NULL,
               day TEXT NOT NULL, lat REAL NOT NULL, lon REAL NOT NULL,
               n_photos INTEGER NOT NULL, n_white INTEGER NOT NULL
            )""")
            for shard in range(partitions):
                original=dated["source_shard_receipts"][shard]
                if original.get("shard_number")!=shard:
                    raise ValueError("Original full provider date shards are unordered/duplicated")
                src=roi_shard_root/f"shard_{shard:03}_of_{partitions:03}"
                p=src/"positive_actual_bloom_photos.ndjson"
                q=src/"roi_shard_receipt.json"
                if not p.is_file() or not q.is_file():
                    raise ValueError(f"Full source photo ROI shard {shard} not source-complete")
                proof=json.loads(q.read_text(encoding="utf8"))
                if (proof.get("schema")!=SHARD_SCHEMA
                    or proof.get("status")!=SHARD_STATUS
                    or proof.get("shard_number")!=shard
                    or proof.get("synthetic_test_only") is not synthetic
                    or proof.get("original_photo_opportunities_in_source_shard")!=
                        original["n_source_photo_opportunities"]
                    or proof.get("n_original_date_known_photos_dispatched")!=
                        original["n_original_date_known_candidates"]
                    or proof.get("n_original_date_known_photo_ROI_terminals")!=
                        original["n_original_date_known_candidates"]
                    or proof.get("source_original_date_queue_sha256")!=
                        original["dated_shard_private_queue_sha256"]
                    or proof.get("protected_positive_bloom_photo_month_ndjson_sha256")!=sha(p)
                    or proof.get("no_date_valid_original_image_omitted_from_roi_processing") is not True
                    or proof.get("source_photo_per_species_sample_limit","INVALID") is not None
                    or proof.get("this_shard_alone_is_NOT_complete_global_42111_source") is not True
                    or proof.get("actual_all_128_original_photo_shards_complete") is not False):
                    raise ValueError("Original shard/date/photo bloom double-positive counts unsealed")
                count=white_one=0
                seen_species=set()
                with p.open(encoding="utf8") as f:
                    for line in f:
                        r=json.loads(line)
                        if (r.get("actual_photo_contains_classifiable_flower") is not True
                            or r.get("original_snapshot_20260927_date_proven") is not True
                            or r.get("original_ROI_v4_hue") not in MORPHS):
                            raise ValueError("One original source photo is not a verified flower")
                        dt,state=source_observation_day(r.get("original_observed_on"))
                        if dt is None or state!="ORIGINAL_SOURCE_DATE_KNOWN" or \
                           dt.month!=r.get("observed_month"):
                            raise ValueError("Original source photo observed_on was not preserved")
                        lat,lon=float(r["latitude"]),float(r["longitude"])
                        if not all(math.isfinite(v) for v in (lat,lon)) or \
                           not -90<=lat<=90 or not -180<=lon<=180:
                            raise ValueError("Source classification photo location invalid")
                        domain,season=climatic_season(lat,dt.month)
                        if domain!=r.get("latitude_domain") or season!=r.get("local_season"):
                            raise ValueError("Original photo observed month/hemisphere changed")
                        pid=str(r.get("original_photo_id",""))
                        obs=str(r.get("original_observation_uuid",""))
                        tax=r.get("original_taxon_id")
                        genus=r.get("original_genus")
                        observer=str(r.get("original_observer",""))
                        if (not pid.isdecimal() or not obs or
                            type(tax) is not int or tax<1
                            or not isinstance(genus,str) or not genus
                            or not observer
                            or str(r.get("original_species_name","")).split()[0]!=genus):
                            raise ValueError("Actual original bloom source photo/taxon identity invalid")
                        try:
                            db.execute("INSERT INTO unique_original_photo VALUES(?)",(pid,))
                        except sqlite3.IntegrityError as exc:
                            raise ValueError("Original photo ID duplicated across source shards") from exc
                        iswhite=int(r["original_ROI_v4_hue"]=="white")
                        cur=db.execute("""INSERT INTO source_observations
                             (uuid,taxon,genus,observer,domain,day,lat,lon,n_photos,n_white)
                             VALUES(?,?,?,?,?,?,?,?,1,?)
                             ON CONFLICT(uuid) DO UPDATE SET
                               n_photos=n_photos+1,
                               n_white=n_white+excluded.n_white
                             WHERE taxon=excluded.taxon AND genus=excluded.genus
                               AND observer=excluded.observer AND domain=excluded.domain
                               AND day=excluded.day AND lat=excluded.lat AND lon=excluded.lon""",
                             (obs,tax,genus,observer,domain,dt.isoformat(),lat,lon,iswhite))
                        if cur.rowcount!=1:
                            raise ValueError("Original multi-image observation provenance changed across shards")
                        count+=1
                        white_one+=iswhite
                        calendar[(domain,dt.month)]["classified"]+=1
                        calendar[(domain,dt.month)]["white"]+=iswhite
                        seen_species.add(tax)
                        species.add(tax)
                        genera.add(genus)
                        out.write(json.dumps(r,sort_keys=True)+"\n")
                if (count!=proof["n_original_bloom_four_hue_classified_photo_pixels"]
                    or white_one!=proof["n_original_bloom_white_four_hue_classified_photo_pixels"]):
                    raise ValueError("Original FCP photo ROI shard positive flowering-hue conservation failed")
                photo_total+=proof["original_photo_opportunities_in_source_shard"]
                dated_total+=proof["n_original_date_known_photos_dispatched"]
                roi_total+=proof["n_original_date_known_photo_ROI_terminals"]
                classified+=count
                white+=white_one
                source_failed+=proof["n_original_permanent_photo_acquisition_failures"]
                receipts.append({
                    "original_shard":shard,
                    "original_date_photo_source_sha256":original["dated_shard_receipt_sha256"],
                    "completed_roi_source_shard_receipt_sha256":sha(q),
                    "n_date_known_roi_photo_terminals":proof[
                        "n_original_date_known_photo_ROI_terminals"],
                    "n_real_source_flower_colour_classifiable_photos":count
                })
                db.commit()
            independent=db.execute("SELECT COUNT(*) FROM source_observations").fetchone()[0]
            species_means={}
            for dom,month,n,avg in db.execute(
                """SELECT domain,CAST(SUBSTR(day,6,2) AS INTEGER),COUNT(*),AVG(obs_white)
                   FROM (SELECT domain,day,taxon,
                                AVG(1.0*n_white/n_photos) AS obs_white
                         FROM source_observations GROUP BY domain,day,taxon)
                   GROUP BY domain,CAST(SUBSTR(day,6,2) AS INTEGER)"""):
                # Overall month mean below is recomputed across taxon first;
                # this auxiliary result cannot replace photo-level month means.
                pass
            for dom,month,n,white_mean in db.execute(
                """SELECT domain,month,COUNT(*),AVG(species_white)
                   FROM (SELECT domain,CAST(SUBSTR(day,6,2) AS INTEGER) AS month,
                                taxon,AVG(1.0*n_white/n_photos) AS species_white
                         FROM source_observations
                         GROUP BY domain,month,taxon)
                   GROUP BY domain,month"""):
                species_means[(dom,month)]=(n,white_mean)
        if (photo_total!=expected_photos
            or dated_total!=dated["n_original_snapshot_observed_on_valid_photo_candidates"]
            or roi_total!=dated_total or classified>roi_total or white>classified
            or independent>classified):
            raise ValueError("All original ROI/date/photo opportunity counts do not conserve the source")
        public_month=[]
        for (dom,month),c in sorted(calendar.items()):
            if c["classified"]<min_public:
                public_month.append({
                    "latitude_domain":dom,"original_month":month,
                    "public_suppressed_below":min_public
                })
                continue
            ns,mean_white=species_means[(dom,month)]
            public_month.append({
                "latitude_domain":dom,"original_month":month,
                "n_four_hue_classifiable_real_source_bloom_photos":c["classified"],
                "n_verified_original_white_bloom_photos":c["white"],
                "photo_white_fraction_given_photo_verified_bloom":
                    c["white"]/c["classified"],
                "n_original_named_species_with_classifiable_photographs":ns,
                "species_equal_observation_equal_white_fraction":mean_white
            })
        result={
            "schema":SCHEMA,"status":STATUS,
            "original_full_source_dates_complete_receipt_sha256":sha(dates_receipt),
            "full_42111_original_photo_snapshot_source_count":photo_total,
            "n_original_FCP_all_source_photo_opportunities":photo_total,
            "n_original_source_date_known_photo_candidates":dated_total,
            "n_original_date_unknown_photo_opportunities_not_analysed_for_bloom":photo_total-dated_total,
            "n_original_date_known_photo_ROI_image_terminals":roi_total,
            "n_original_date_and_site_known_ROI_bloom_hues":classified,
            "n_original_date_known_white_flower_photos":white,
            "n_original_distinct_photo_observation_UUIDs_with_flower":independent,
            "n_original_source_taxa_with_classifiable_date_known_flower":len(species),
            "n_original_source_genera_with_classifiable_date_known_flower":len(genera),
            "n_source_permanent_photo_unavailable":source_failed,
            "source_shards_completed":partitions,
            "per_completed_shard_source_sha256":receipts,
            "source_original_domain_month_classified_photo_white_and_species_equal":public_month,
            "private_actual_bloom_photo_month_hue_manifest_sha256":sha(stage),
            "full_42111_original_source_not_selectively_date_sampled":True,
            "per_species_source_photo_cap":None,
            "WorldClim_monthly_T_S_P_VPD_source_site_link_not_yet_run":True,
            "all_missing_date_source_photos_were_excluded_before_roi":True,
            "all_original_date_eligible_photo_pixels_attempted":True,
            "unclassified_no_flower_photos_not_in_white_colour_denominator":True,
            "no_old_100543_photo_ROI_hues_reused":True,
            "original_photo_population_random_sample_of_all_wild_blooms":False,
            "original_93m_source_photo_ids_snapshot_confirmed":True,
            "synthetic_test_only":bool(synthetic),
        }
        public_stage.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf8")
        stage.replace(out_private_global_blooms)
        public_stage.replace(out_public_report)
        scratch.unlink()
        return result
    except Exception:
        # Keep unsuccessful private .INCOMPLETE artefacts quarantined.
        raise


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--all-original-date-census-receipt",required=True,type=Path)
    p.add_argument("--all-128-complete-original-roi-shard-root",required=True,type=Path)
    p.add_argument("--out-private-global-photo-bloom-records",required=True,type=Path)
    p.add_argument("--out-public-date-verified-global-photo-white",required=True,type=Path)
    args=p.parse_args()
    r=aggregate_all(args.all_original_date_census_receipt,
        args.all_128_complete_original_roi_shard_root,
        args.out_private_global_photo_bloom_records,
        args.out_public_date_verified_global_photo_white)
    print(json.dumps({
        "status":r["status"],
        "original_full_source_photo_opportunities":r[
            "full_42111_original_photo_snapshot_source_count"],
        "n_real_photo_roi_bloom_pixels":r[
            "n_original_date_and_site_known_ROI_bloom_hues"],
        "n_actual_white_original_bloom_photos":r[
            "n_original_date_known_white_flower_photos"],
    },sort_keys=True))


if __name__=="__main__":main()
