#!/usr/bin/env python3
"""Source-locked 128-shard all-dated photo candidate reconciliation.

Accept ONLY the entire original 2026-09-27 FCP 42,111-taxon census; exact
128 partitions, all 93,209,351 independent original photo opportunities.
Every shard's date-only NDJSON is source-SHA authenticated independently.

The output reports TRUE original date-known candidate photographs and monthly
original photo OBSERVATION effort. It reports ZERO actual flower-verified
photos because metadata do not contain per-photo confirmed blooming status.
It never opens pixels, consumes historical 100543 hues, or imposes a per-taxon
cap. Protected individual photo IDs/coordinates never enter public output.

A staged shard is not complete globe; missing even one shard hard fails.
"""
from __future__ import annotations
import argparse
import json
import hashlib
from collections import Counter
from pathlib import Path
from scripts.acquisition.mois_full_source_dated_flower_candidates_20261016 import (
    load_counts,sha,source_observation_day,climatic_season,SCHEMA as SRC_SCHEMA,
    STATUS as SRC_STATUS,FRAME_SHA,COUNTS_SHA,FULL_TAXA,
    FULL_PHOTO_OPPORTUNITIES,FULL_WORK_PARTITIONS,SNAPSHOT
)

SCHEMA="mois_all_original_42111_dated_photo_candidates_128_shards_v1"
STATUS="ALL_ORIGINAL_PHOTO_IDS_RECONCILED_DATES_ONLY_BLOOM_UNVERIFIED"


def aggregate(shards_root:Path, dated_root:Path, counts_csv:Path,
              output:Path, *, synthetic=False, partitions=FULL_WORK_PARTITIONS,
              expected_photo_total=FULL_PHOTO_OPPORTUNITIES,
              min_public_cell=10)->dict:
    if output.exists() or output.with_name(output.name+".INCOMPLETE").exists():
        raise ValueError("Never overwrite full source date validation result")
    if min_public_cell<10 or partitions<1 or (
        not synthetic and (partitions!=FULL_WORK_PARTITIONS or
                           expected_photo_total!=FULL_PHOTO_OPPORTUNITIES)):
        raise ValueError("No smaller biological source/cap or public count release")
    baseline=load_counts(counts_csv,synthetic=synthetic)
    if not synthetic and len(baseline)!=FULL_TAXA:
        raise ValueError("Original complete species opportunity count changed")
    photo_n=known_n=0
    months=Counter()
    seasons=Counter()
    date_states=Counter()
    taxa_all=set()
    taxa_dated=set()
    frame=None
    source_sqlite=None
    shard_receipts=[]
    protected_keys={
        "inat_taxon_id","species","observation_id","photo_id",
        "photo_url_large","photo_license","latitude","longitude",
        "observer_id","observed_on","observed_month","observed_year",
        "source_snapshot_date","latitude_domain","local_season",
        "original_photo_flowering_state","origin"
    }
    for part in range(partitions):
        shard=shards_root/f"shard_{part:03}_of_{partitions:03}"
        original_shard=shard/"shard_source_receipt.json"
        if not original_shard.is_file():
            raise ValueError(f"Cannot promote incomplete original photo-ID shard {part}")
        fcp=json.loads(original_shard.read_text(encoding="utf8"))
        receipt_path=dated_root/f"dated_candidates_shard_{part:03}_receipt.json"
        queue=dated_root/f"dated_candidates_shard_{part:03}.ndjson"
        if not receipt_path.is_file() or not queue.is_file():
            raise ValueError(f"Cannot promote unsealed original dated candidate shard {part}")
        receipt=json.loads(receipt_path.read_text(encoding="utf8"))
        if (receipt.get("schema")!=SRC_SCHEMA or receipt.get("status")!=SRC_STATUS
            or receipt.get("shard_number")!=part
            or receipt.get("n_fcp_shards_expected_for_global_completion")!=partitions
            or receipt.get("source_photo_cap_per_taxon","INVALID") is not None
            or receipt.get("metadata_date_known_is_not_blooming_photo_confirmation") is not True
            or receipt.get("source_original_counts_ledger_sha256")!=sha(counts_csv)
            or receipt.get("source_shard_receipt_sha256")!=sha(original_shard)
            or receipt.get("private_full_dated_source_photo_queue_sha256")!=sha(queue)
            or receipt.get("actual_flower_photographs_verified_from_pixels")!=0
            or receipt.get("photo_hue_classification_attempted") is not False
            or receipt.get("synthetic_test_only") is not synthetic
            or fcp.get("partition")!=part or fcp.get("n_shards")!=partitions
            or fcp.get("status")!="SOURCE_SHARD_COMPLETE_FROM_FROZEN_SQLITE_NO_IMAGES"):
            raise ValueError("One photo ID/month shard receipt was not independently source sealed")
        source_taxa=fcp.get("source_species_ids")
        if (not isinstance(source_taxa,list) or
            len(source_taxa)!=len(set(source_taxa)) or
            any(type(t) is not int or t not in baseline or t in taxa_all
                for t in source_taxa)):
            raise ValueError("Original source taxon/photo shards are not disjoint")
        taxa_all.update(source_taxa)
        n_source=sum(baseline[t][1] for t in source_taxa)
        if (n_source!=fcp.get("n_original_source_photo_ids") or
            n_source!=receipt.get("n_original_eligible_source_photographs_examined") or
            len(source_taxa)!=receipt.get("n_source_taxa_in_this_independent_shard")):
            raise ValueError("Per-shard photo counts differ from original independent census")
        ident=receipt.get("original_source_frame_sha256")
        if frame is None:frame=ident
        if frame!=ident or ident!=fcp.get("source_species_frame_sha256"):
            raise ValueError("Source frame changed within original photo shards")
        stage=fcp.get("stage_sqlite_sha256")
        if source_sqlite is None:source_sqlite=stage
        if stage!=source_sqlite:
            raise ValueError("Original provider SQLite photo census source changed between shards")
        photo_n+=n_source
        count_here=0
        taxa_here=set()
        with queue.open(encoding="utf8") as f:
            for line in f:
                r=json.loads(line)
                if set(r)!=protected_keys:
                    raise ValueError("Unverified original date candidate provenance field")
                t=r["inat_taxon_id"]
                if (type(t) is not int or t not in source_taxa
                    or r["species"]!=baseline[t][0]
                    or r["source_snapshot_date"]!=SNAPSHOT
                    or r["original_photo_flowering_state"]!="NOT_YET_PHOTO_VERIFIED"
                    or r["origin"]!="ORIGINAL_DATED_FCP_SHARD_NO_PHOTO_PIXELS"):
                    raise ValueError("Original dated photo candidate corrupted or wrong source taxon")
                d,state=source_observation_day(r["observed_on"])
                if d is None or state!="ORIGINAL_SOURCE_DATE_KNOWN" or (
                    r["observed_month"]!=d.month or r["observed_year"]!=d.year):
                    raise ValueError("Original source date not actual snapshot observation date")
                lat=float(r["latitude"])
                domain,season=climatic_season(lat,d.month)
                if domain!=r["latitude_domain"] or season!=r["local_season"]:
                    raise ValueError("Source original hemisphere season changed")
                if not str(r["photo_id"]).isdecimal() or not r["observation_id"]:
                    raise ValueError("Original photo identity not a real dated candidate")
                count_here+=1
                taxa_here.add(t)
                months[(domain,d.month)]+=1
                seasons[(domain,season)]+=1
        if (count_here!=receipt["n_original_photo_observation_dates_usable"]
            or len(taxa_here)!=receipt["n_taxa_with_usable_observed_date"]
            or n_source-count_here!=receipt["n_original_photo_observation_dates_excluded"]
            or sum(receipt["date_screen_state"].values())!=n_source
            or receipt["date_screen_state"].get("ORIGINAL_SOURCE_DATE_KNOWN",0)!=count_here):
            raise ValueError("Original dated queue and exclusion ledger counts inconsistent")
        known_n+=count_here
        taxa_dated.update(taxa_here)
        date_states.update(receipt["date_screen_state"])
        shard_receipts.append({
            "shard_number":part,
            "fcp_complete_shard_receipt_sha256":sha(original_shard),
            "dated_shard_receipt_sha256":sha(receipt_path),
            "dated_shard_private_queue_sha256":sha(queue),
            "n_source_photo_opportunities":n_source,
            "n_original_date_known_candidates":count_here,
        })
    if len(taxa_all)!=len(baseline) or taxa_all!=set(baseline):
        raise ValueError("All original 42,111 source taxa not present across sealed shards")
    if photo_n!=sum(n for _,n in baseline.values()) or photo_n!=expected_photo_total:
        raise ValueError("Actual source photo IDs do not reconcile to original photo opportunity census")
    if (not synthetic and frame!=FRAME_SHA) or (
        date_states["ORIGINAL_SOURCE_DATE_KNOWN"]!=known_n
        or known_n!=sum(months.values())
        or known_n!=sum(seasons.values())
        or sum(date_states.values())!=photo_n):
        raise ValueError("Original provider monthly/date conservation failed")
    result={
        "schema":SCHEMA,"status":STATUS,
        "original_full_42111_species_frame_sha256":frame,
        "original_full_photo_count_ledger_sha256":sha(counts_csv),
        "original_source_stage_sqlite_sha256":source_sqlite,
        "source_original_shards_complete":partitions,
        "n_all_original_source_named_taxa":len(taxa_all),
        "n_all_original_unique_photo_opportunities":photo_n,
        "n_original_snapshot_observed_on_valid_photo_candidates":known_n,
        "n_original_snapshot_observed_on_invalid_or_missing_photo_candidates":
            photo_n-known_n,
        "n_named_taxa_with_at_least_one_original_valid_dated_photo":len(taxa_dated),
        "n_named_taxa_without_original_valid_dated_photos":
            len(taxa_all)-len(taxa_dated),
        "date_quality_states":dict(sorted(date_states.items())),
        "latitude_zone_original_observed_month_counts":[
            {"latitude_domain":zone,"original_calendar_month":month,
             "n_original_dated_photo_candidates":n}
            if n>=min_public_cell else
            {"latitude_domain":zone,"original_calendar_month":month,
             "public_suppressed_below":min_public_cell}
            for (zone,month),n in sorted(months.items())],
        "latitude_zone_local_season_photo_opportunities":[
            {"latitude_domain":zone,"local_season":season,
             "n_original_dated_photo_candidates":n}
            if n>=min_public_cell else
            {"latitude_domain":zone,"local_season":season,
             "public_suppressed_below":min_public_cell}
            for (zone,season),n in sorted(seasons.items())],
        "source_shard_receipts":shard_receipts,
        "source_photo_cap_per_species":None,
        "source_observation_dates_from_immutable_provider_archive":True,
        "original_source_photo_bloom_phenology_verified_by_image_pixels":False,
        "n_actual_original_photo_ROI_flower_colour_measurements_this_stage":0,
        "unverified_source_photo_automatically_assigned_white_or_nonwhite":False,
        "old_original_100543_photo_cohort_reused_as_original_full_data":False,
        "photographed_blooming_plants_not_a_population_census":True,
        "synthetic_test_only":bool(synthetic)
    }
    output.parent.mkdir(parents=True,exist_ok=True)
    staged=output.with_name(output.name+".INCOMPLETE")
    staged.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf8")
    staged.replace(output)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fcp-sealed-shards-root",required=True,type=Path)
    p.add_argument("--mois-private-dated-shards-root",required=True,type=Path)
    p.add_argument("--original-42111-count-csv",required=True,type=Path)
    p.add_argument("--out-public-all-original-dated-photo-receipt",required=True,type=Path)
    a=p.parse_args()
    data=aggregate(a.fcp_sealed_shards_root,a.mois_private_dated_shards_root,
         a.original_42111_count_csv,a.out_public_all_original_dated_photo_receipt)
    print(json.dumps({k:data[k] for k in (
       "status","n_all_original_unique_photo_opportunities",
       "n_original_snapshot_observed_on_valid_photo_candidates",
       "n_actual_original_photo_ROI_flower_colour_measurements_this_stage"
    )},sort_keys=True))


if __name__=="__main__":main()
