#!/usr/bin/env python3
"""Stage ALL original 42,111-taxon FCP dated-photo candidates for bloom ROI.

This operates on one byte-sealed full-provider FCP source SHARD at a time.
No original 1-species-40-observation sampling cap, old 100543-photo cohort,
current inat API date, or downloaded plant image is permitted here.

iNaturalist Open Dataset's observations.tsv.gz includes observed_on and
observation_uuid; photos.tsv.gz has photo_id/observation_uuid but NO verified
per-image Flowers annotation. Thus dated source photographs here are ONLY
BLOOM-UNVERIFIED candidates. A separate original frozen image ROI outcome
must positively verify a flower and a hue before any photo enters a seasonal
white-flower numerator or denominator. Unclassifiable photos cannot be
silently called nonwhite or absent flowers.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path

FRAME_SHA="5d871bd1f1190c68fd513bb527e6c93de1606dc35ff8cfff5f6187850dc382cc"
COUNTS_SHA="ef26c16252dd785f2d08bde3b18938c7624d73a54386fa3d3974dceb09ed97a9"
SNAPSHOT="2026-09-27"
SNAPSHOT_DAY=date(2026,9,27)
FULL_TAXA=42111
FULL_PHOTO_OPPORTUNITIES=93209351
FULL_WORK_PARTITIONS=128
SOURCE_STATUS="SOURCE_SHARD_COMPLETE_FROM_FROZEN_SQLITE_NO_IMAGES"
SCHEMA="mois_original_42111_all_dated_photo_candidate_shard_v1"
STATUS="ORIGINAL_DATED_PHOTO_OPPORTUNITIES_BLOOM_UNVERIFIED_ROI_NOT_MEASURED"
FIELDS={
    "inat_taxon_id","species","observation_id","photo_id",
    "photo_url_large","photo_license","observer_id","observed_on",
    "latitude","longitude","source_snapshot_date","source_query"
}


def sha(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(4*1024*1024),b""):
            h.update(block)
    return h.hexdigest()


def source_observation_day(raw:str):
    """Only original recorded ISO date; no upload date or inferred month."""
    if not isinstance(raw,str) or not re.fullmatch(r"[0-9]{4}-[0-9]{2}-[0-9]{2}",raw):
        return None,"SOURCE_DATE_NOT_VALID_ISO_DAY"
    try:d=date.fromisoformat(raw)
    except ValueError:return None,"SOURCE_DATE_IMPOSSIBLE_CALENDAR_DAY"
    if d>SNAPSHOT_DAY:return None,"SOURCE_DATE_AFTER_FROZEN_SNAPSHOT"
    return d,"ORIGINAL_SOURCE_DATE_KNOWN"


def climatic_season(lat:float,month:int)->tuple[str,str]:
    if abs(lat)<23.5:
        return "TROPICAL",f"CALENDAR_Q{1+(month-1)//3}"
    season={12:"WINTER",1:"WINTER",2:"WINTER",3:"SPRING",4:"SPRING",
            5:"SPRING",6:"SUMMER",7:"SUMMER",8:"SUMMER",
            9:"AUTUMN",10:"AUTUMN",11:"AUTUMN"}[month]
    if lat<0:season={"SPRING":"AUTUMN","SUMMER":"WINTER",
                     "AUTUMN":"SPRING","WINTER":"SUMMER"}[season]
    return ("NORTH_EXTRATROPICAL" if lat>0 else
            "SOUTH_EXTRATROPICAL"),season


def load_counts(original_counts:Path,*,synthetic=False)->dict[int,tuple[str,int]]:
    if not synthetic and sha(original_counts)!=COUNTS_SHA:
        raise ValueError("Original 42111-taxon independent photo count ledger SHA mismatch")
    rows={}
    with original_counts.open(newline="",encoding="utf8") as f:
        reader=csv.DictReader(f)
        if reader.fieldnames!=[
            "inat_taxon_id","species","n_unique_source_photo_ids","snapshot_date"
        ]:
            raise ValueError("Not original taxon-specific FCP source photo count ledger")
        for row in reader:
            tid=row["inat_taxon_id"].strip()
            raw=row["n_unique_source_photo_ids"].strip()
            if (not tid.isdecimal() or not raw.isdecimal() or int(tid) in rows
                or not row["species"] or row["snapshot_date"]!=SNAPSHOT):
                raise ValueError("Duplicated, invalid or changed source-frame species/photo count")
            rows[int(tid)]=(row["species"],int(raw))
    if not synthetic and (len(rows)!=FULL_TAXA or
                          sum(n for _,n in rows.values())!=FULL_PHOTO_OPPORTUNITIES):
        raise ValueError("Source independent 42111/93209351 metadata count total drifted")
    return rows


def audit_shard(shard:Path,counts:Path,private_dated:Path,public_json:Path,*,
                synthetic=False,min_public_cell=10)->dict:
    """Read complete source taxon manifests (incl. zero-photo) in bounded chunks.

    One original FCP shard, NOT an independently eligible subset of all
    42,111 taxa. The 128 sealed original source shards still need global
    completion before claiming all 93M dated opportunities have been counted.
    """
    for out in (private_dated,public_json):
        if out.exists() or out.with_name(out.name+".INCOMPLETE").exists():
            raise ValueError("Cannot overwrite/interleave an original source-date shard")
    if min_public_cell<10:
        raise ValueError("Public small-source cell disclosure threshold may not be relaxed")
    baseline=load_counts(counts,synthetic=synthetic)
    receipt=json.loads((shard/"shard_source_receipt.json").read_text(encoding="utf8"))
    contract=json.loads((shard/"source_contract.json").read_text(encoding="utf8"))
    taxa=receipt.get("source_species_ids")
    if (receipt.get("schema")!="fcp_frozen_source_shard_emission_v1"
        or receipt.get("status")!=SOURCE_STATUS
        or not isinstance(taxa,list) or len(taxa)!=receipt.get("n_source_species")
        or any(not isinstance(v,int) or v not in baseline for v in taxa)
        or len(set(taxa))!=len(taxa)
        or contract.get("schema")!="fcp_uncapped_shard_source_contract_v1"
        or contract.get("photo_cap","INVALID") is not None
        or contract.get("outcome_blind") is not True
        or contract.get("created_d2")!=SNAPSHOT
        or contract.get("shard")!=receipt.get("partition")
        or contract.get("nshards")!=receipt.get("n_shards")
        or receipt.get("n_original_source_photo_ids")!=
            sum(baseline[tid][1] for tid in taxa)
        or receipt.get("historical_colour_outcomes_used") is not False
        or receipt.get("photo_pixels_opened") is not False
        or receipt.get("per_species_photo_cap","INVALID") is not None
        or receipt.get("source_species_frame_sha256")!=contract.get("source_species_frame_sha256")):
        raise ValueError("Not an intact, uncapped original dated metadata source shard")
    if not synthetic and (
        receipt.get("source_species_frame_sha256")!=FRAME_SHA
        or receipt.get("n_shards")!=FULL_WORK_PARTITIONS):
        raise ValueError("Source taxonomic/photo 42111-shard contract drifted")
    count_all=dated=0
    reason=Counter()
    domain_month=Counter()
    domain_season=Counter()
    by_tax=Counter()
    seen_photos=set()
    private_dated.parent.mkdir(parents=True,exist_ok=True)
    public_json.parent.mkdir(parents=True,exist_ok=True)
    staged=private_dated.with_name(private_dated.name+".INCOMPLETE")
    staged_public=public_json.with_name(public_json.name+".INCOMPLETE")
    try:
        with staged.open("w",encoding="utf8") as output:
            for tid in taxa:
                folder=shard/"species"/str(tid)
                check=json.loads((folder/"checkpoint.json").read_text(encoding="utf8"))
                file=folder/"photo_manifest.ndjson"
                if (check.get("taxon_id")!=tid or check.get("status")!="complete"
                    or check.get("source_photo_rows")!=baseline[tid][1]
                    or check.get("photo_cap","INVALID") is not None
                    or check.get("colour_pixels_opened") is not False):
                    raise ValueError("Original complete source taxon checkpoint changed")
                source_name=baseline[tid][0]
                digest=hashlib.sha256()
                n=0
                with file.open("rb") as inp:
                    for line in inp:
                        digest.update(line)
                        n+=1;count_all+=1
                        record=json.loads(line)
                        if (set(record)!=FIELDS | {"attribution","positional_accuracy_m"}
                            or record["inat_taxon_id"]!=tid
                            or record["species"]!=source_name
                            or record["source_snapshot_date"]!=SNAPSHOT
                            or not str(record["photo_id"]).isdecimal()
                            or not str(record["observation_id"]).strip()
                            or not str(record["photo_url_large"]).startswith(
                               "https://inaturalist-open-data.s3.amazonaws.com/photos/")
                            or not record["photo_license"]):
                            raise ValueError("Original photo-ID/location/date source row incompatible with FCP schema")
                        pid=str(record["photo_id"])
                        if pid in seen_photos:
                            raise ValueError("Original provider photo ID duplicated within shard")
                        seen_photos.add(pid)
                        day,state=source_observation_day(record["observed_on"])
                        reason[state]+=1
                        if day is None:continue
                        lat=float(record["latitude"])
                        lon=float(record["longitude"])
                        if not (math.isfinite(lat) and math.isfinite(lon)
                                and -90<=lat<=90 and -180<=lon<=180):
                            raise ValueError("Date-valid original source photo has invalid WGS84 site")
                        domain,season=climatic_season(lat,day.month)
                        dated+=1
                        by_tax[tid]+=1
                        domain_month[(domain,day.month)]+=1
                        domain_season[(domain,season)]+=1
                        # PRIVATE, immutable original snapshot provenance:
                        # never imply original photo contains a verified flower
                        row={
                            "inat_taxon_id":tid,
                            "species":source_name,
                            "observation_id":record["observation_id"],
                            "photo_id":pid,
                            "photo_url_large":record["photo_url_large"],
                            "photo_license":record["photo_license"],
                            "latitude":lat,"longitude":lon,
                            "observer_id":record["observer_id"],
                            "observed_on":day.isoformat(),
                            "observed_month":day.month,
                            "observed_year":day.year,
                            "source_snapshot_date":SNAPSHOT,
                            "latitude_domain":domain,
                            "local_season":season,
                            "original_photo_flowering_state":"NOT_YET_PHOTO_VERIFIED",
                            "origin":"ORIGINAL_DATED_FCP_SHARD_NO_PHOTO_PIXELS",
                        }
                        output.write(json.dumps(row,ensure_ascii=False,sort_keys=True)+"\n")
                if (n!=baseline[tid][1] or check.get("metadata_sha256")!=digest.hexdigest()):
                    raise ValueError("Original source taxon photo manifest SHA/count mismatch")
        if count_all!=receipt["n_original_source_photo_ids"] or sum(reason.values())!=count_all:
            raise ValueError("Original source shard enumerated-photo and calendar exclusions do not partition source")
        expected={
            "schema":SCHEMA,"status":STATUS,
            "original_source_frame_sha256":receipt["source_species_frame_sha256"],
            "source_shard_receipt_sha256":sha(shard/"shard_source_receipt.json"),
            "source_original_counts_ledger_sha256":sha(counts),
            "shard_number":receipt["partition"],
            "n_fcp_shards_expected_for_global_completion":receipt["n_shards"],
            "n_original_eligible_source_photographs_examined":count_all,
            "n_original_photo_observation_dates_usable":dated,
            "n_original_photo_observation_dates_excluded":count_all-dated,
            "date_screen_state":dict(sorted(reason.items())),
            "n_taxa_with_usable_observed_date":len(by_tax),
            "n_taxa_with_at_least_one_original_photo":sum(
                baseline[t][1]>0 for t in taxa),
            "n_source_taxa_in_this_independent_shard":len(taxa),
            "source_photo_cap_per_taxon":None,
            "metadata_rejected_by_observer_added_Flowers_annotation":0,
            "source_annotation_term_Flowers_available_in_original_four_tables":False,
            "actual_flower_photographs_verified_from_pixels":0,
            "photo_hue_classification_attempted":False,
            "is_an_all_fcp_128_shard_complete_census":False,
            "source_observation_month_domain_counts":[
                {"latitude_domain":d,"month":m,"n_original_date_known_candidate_photos":n}
                if n>=min_public_cell else
                {"latitude_domain":d,"month":m,"public_suppressed_below":min_public_cell}
                for (d,m),n in sorted(domain_month.items())
            ],
            "source_local_season_domain_counts":[
                {"latitude_domain":d,"season":s,"n_original_date_known_candidate_photos":n}
                if n>=min_public_cell else
                {"latitude_domain":d,"season":s,"public_suppressed_below":min_public_cell}
                for (d,s),n in sorted(domain_season.items())
            ],
            "private_full_dated_source_photo_queue_sha256":sha(staged),
            "snapshot_observed_on_only_not_current_API_date":True,
            "metadata_date_known_is_not_blooming_photo_confirmation":True,
            "site_photo_coordinates_and_ids_private_only":True,
            "synthetic_test_only":bool(synthetic),
        }
        staged_public.write_text(json.dumps(expected,indent=2,sort_keys=True)+"\n",encoding="utf8")
        os.replace(staged,private_dated)
        os.replace(staged_public,public_json)
        return expected
    except Exception:
        for f in (staged,staged_public):
            if f.exists():f.unlink()
        raise


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sealed-original-fcp-shard",type=Path,required=True)
    p.add_argument("--independent-42111-photo-counts",type=Path,required=True)
    p.add_argument("--out-private-dated-photo-queue",type=Path,required=True)
    p.add_argument("--out-public-source-receipt",type=Path,required=True)
    args=p.parse_args()
    result=audit_shard(args.sealed_original_fcp_shard,
        args.independent_42111_photo_counts,args.out_private_dated_photo_queue,
        args.out_public_source_receipt)
    print(json.dumps({
        "status":result["status"],
        "n_original_photo_opportunities_examined":result[
            "n_original_eligible_source_photographs_examined"],
        "n_dated_candidates_not_yet_flowering_verified":result[
            "n_original_photo_observation_dates_usable"],
        "n_photo_flower_labels_verified":0,
    },sort_keys=True))


if __name__=="__main__":main()
