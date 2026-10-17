#!/usr/bin/env python3
"""Select NEW real RGB photo pixels beyond all 309439 measured originals.

The 2026-09-27 FCP source has 3430502 genuine numeric iNaturalist photo IDs,
but photo IDs are NOT colour classifications or proof of open flowering.
Select 4096 with a deterministic hash of photo_id (balanced by 1000 original
taxon IDs, no repeated observer in THIS new cohort) without reading hue,
month/weather, biomass, colour, outcome or flower annotations. Exclude ALL
309439 measured original MOIS photo IDs before source selection.

A separate original ROI-v4 SAM/YOLO worker, not this source-only stage, must
actually retrieve pixels before any source photo is counted classified.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import heapq
import json
import math
import re
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path

from scripts.analysis.mois_full3430502_source_third4096 import (
    FULL_SHA, SOURCE_N, FIELDS, LINK_RE, save, digest,
)
from scripts.analysis.mois_reassemble_original_new4096_measurements import (
    HUES,
)

MEASURED_N=309439
MEASURED_HUES=149662
MEASURED_DATED_HUES=149653
MEASURED_SHA="546266e2bbaa340a2523864c28c5a877842eb4b0904a7e0dff90870c6519a47e"
SELECTION_N=4096
SOURCE_TAXA=1000
N_SHARDS=16
PHOTO_SHARD_SIZE=256
CANDIDATES_PER_TAXON=48
SALT="MOIS_NEXT4096_AFTER_309439_ORIGINAL_3430502_COMPLETE_PIXEL_WAVE4_SOURCE_V1"
SOURCE_OUT="mois_full3430502_next4096_original_photo_IDs_AFTER309439_WAVE4.csv"
FREEZE_OUT="mois_first4096_original_source_only_freeze.json"


def require(condition,why):
    if not condition:
        raise ValueError(why)


def measure_id(pid:str)->str:
    return "MOIS-"+hashlib.sha256((SALT+"|"+pid).encode()).hexdigest()[:24].upper()


def validate_site_day(row:dict)->bool:
    try:
        y=float(row["latitude"]);x=float(row["longitude"])
        d=date.fromisoformat(row["observed_on"])
        return (math.isfinite(y) and math.isfinite(x) and -90<=y<=90 and
                -180<=x<=180 and date(1800,1,1)<=d<=date(2026,10,17))
    except (TypeError,KeyError,ValueError):
        return False


def source_candidate(row:dict, prior_ids:set[str],counts:Counter):
    """Photo-id and provider fields only; no phenotype, sunlight, colour or month rank."""
    counts["source_original_photo_IDs_seen"]+=1
    pid=row.get("photo_id","")
    if not pid.isdecimal() or not row.get("inat_taxon_id","").isdecimal():
        raise ValueError("Original 2026-09-27 iNaturalist source photo/taxon IDs not numeric")
    if pid in prior_ids:
        counts["original_measured_photo_IDs_encountered"]=counts.get(
            "original_measured_photo_IDs_encountered",0)+1
        return None
    counts["original_eligible_unmeasured_photo_IDs_seen"]+=1
    if (row.get("source_snapshot_date")!="2026-09-27" or
        not all(row.get(k,"") for k in (
            "species","observer_id","observation_id","photo_license")) or
        not validate_site_day(row)):
        counts["invalid_source_metadata_no_pixel_inference"]+=1
        return None
    url=row.get("photo_url_large","")
    m=LINK_RE.fullmatch(url)
    if m is None or m.group(1)!=pid:
        counts["invalid_original_licensed_photo_URL_no_pixel_inference"]+=1
        return None
    h=int(hashlib.sha256((SALT+"|"+pid).encode()).hexdigest()[:16],16)
    entry={k:row.get(k,"") for k in (
        "photo_id","inat_taxon_id","species","observer_id","observation_id",
        "latitude","longitude","observed_on","photo_license",
        "photo_url_large","positional_accuracy_m","source_snapshot_date"
    )}
    entry["selection_hash"]=f"{h:016x}"
    entry["screen_fold"]="not_previously_measured_any_MOIS_309439_original_ID"
    return -h,pid,entry


def choose(source_rows, prior_ids:set[str], *, n:int=SELECTION_N,
           expected_taxa:int=SOURCE_TAXA):
    require(0<n<=SELECTION_N and expected_taxa>=1,"Invalid original photo sampling size")
    pools=defaultdict(list)
    counts=Counter()
    for row in source_rows:
        candidate=source_candidate(row,prior_ids,counts)
        if candidate is None:continue
        _,pid,record=candidate
        tax=record["inat_taxon_id"]
        heap=pools[tax]
        if len(heap)<CANDIDATES_PER_TAXON:
            heapq.heappush(heap,candidate)
        elif candidate>heap[0]:
            heapq.heapreplace(heap,candidate)
    require(len(pools)==expected_taxa,
            "Original source no longer has eligible photo IDs in all target taxa")
    sortedpools={tax:[item[-1] for item in sorted(
        heap,key=lambda v:(-v[0],int(v[1])))] for tax,heap in pools.items()}
    order=sorted(sortedpools,
        key=lambda tax:hashlib.sha256((SALT+"|TAXON|"+tax).encode()).hexdigest())
    pointers={tax:0 for tax in order}
    ids=set();obs=set();people=set();chosen=[]
    # Up to 48 source attempts/taxon. No future photo colours or climate labels
    # influence rank. Taxon-level balancing is a source design, not ecosystem
    # frequency/sampling-effort weighting in a universal flora.
    for turn in range(CANDIDATES_PER_TAXON):
        for tax in order:
            if len(chosen)>=n:break
            pool=sortedpools[tax]
            cursor=pointers[tax]
            while cursor<len(pool) and (
                pool[cursor]["observer_id"] in people or
                pool[cursor]["observation_id"] in obs or
                pool[cursor]["photo_id"] in ids):
                cursor+=1
            pointers[tax]=cursor+1
            if cursor>=len(pool):continue
            r=dict(pool[cursor])
            assert r["photo_id"] not in prior_ids
            chosen.append(r)
            people.add(r["observer_id"])
            obs.add(r["observation_id"])
            ids.add(r["photo_id"])
        if len(chosen)>=n:break
    require(len(chosen)==n and len(ids)==len(obs)==len(people)==n,
            "No fully source-valid unique photo/observation/observer 4096-ROI cohort")
    require(len({r["inat_taxon_id"] for r in chosen})==expected_taxa,
            "Original photo wave4 omits a source named taxon")
    for ix,row in enumerate(chosen):
        row["measurement_id"]=measure_id(row["photo_id"])
        row["image_filename"]=row["measurement_id"]+".jpg"
        row["wave_sequence"]=str(ix)
        row["shard"]=str(ix%N_SHARDS)
    return chosen,dict(counts)


def run(full:Path,old:Path,outdir:Path):
    require(not outdir.exists() or not any(outdir.iterdir()),
            "Never rewrite source-photo blind selection receipt")
    require(digest(full)==FULL_SHA,
            "Provider original 3430502 original-photo source SHA changed")
    require(digest(old)==MEASURED_SHA,
            "MOIS already genuinely measured 309439 original photo SHA changed")
    prior_ids=set()
    old_hues=0;old_dated_hues=0
    with gzip.open(old,"rt",encoding="utf-8",newline="") as f:
        for row in csv.DictReader(f):
            pid=row["photo_id"]
            require(pid and pid not in prior_ids,"Original measured photo ID duplicate")
            prior_ids.add(pid)
            hue=row["morph"] in HUES
            old_hues+=hue
            old_dated_hues+=hue and bool(row["qualified_date"])
    require(len(prior_ids)==MEASURED_N and old_hues==MEASURED_HUES and
            old_dated_hues==MEASURED_DATED_HUES,
            "Genuine measured original photo/date/hue source denominator changed")
    with full.open(newline="",encoding="utf-8") as stream:
        reader=csv.DictReader(stream)
        require({"photo_id","inat_taxon_id","observer_id","observation_id",
                 "observed_on","photo_url_large","source_snapshot_date"}.
                issubset(reader.fieldnames or ()),
                "Original complete iNaturalist photo source changed columns")
        chosen,counts=choose(reader,prior_ids)
    require(counts["source_original_photo_IDs_seen"]==SOURCE_N,
            "Full 3.43-million original photo-ID source not scanned completely")
    outdir.mkdir(parents=True,exist_ok=True)
    sourcefile=outdir/SOURCE_OUT
    save(sourcefile,chosen)
    shards={}
    for shard in range(N_SHARDS):
        f=outdir/f"mois_source_only_shard_{shard:02d}_256.csv"
        rows=[r for r in chosen if int(r["shard"])==shard]
        require(len(rows)==PHOTO_SHARD_SIZE,
                "Original source 4096 photographic IDs not partitioned exactly")
        save(f,rows)
        shards[str(shard)]=digest(f)
    receipt={
        "schema":"mois_original_full3430502_next4096_after309439_wave4_preoutcome_v1",
        "status":"SOURCE_ONLY_4096_RESERVED_NO_NEW_PIXEL_HUES",
        "source_photo_ids_original":SOURCE_N,
        "prior_original_measured_photo_ids":MEASURED_N,
        "prior_original_classified_four_hues":MEASURED_HUES,
        "prior_original_dated_classified_hues":MEASURED_DATED_HUES,
        "source_full_original_3430502_SHA256":FULL_SHA,
        "source_measured_original_309439_SHA256":MEASURED_SHA,
        "number_source_photo_ID_rows_scanned":counts["source_original_photo_IDs_seen"],
        "number_previously_measured_photo_IDs_encountered":
            counts["original_measured_photo_IDs_encountered"],
        "number_unmeasured_original_photo_ID_opportunities":
            counts["original_eligible_unmeasured_photo_IDs_seen"],
        "number_source_photo_metadata_unusable":counts.get("invalid_source_metadata_no_pixel_inference",0),
        "number_source_photo_url_license_unusable":
            counts.get("invalid_original_licensed_photo_URL_no_pixel_inference",0),
        "first_wave_original_source_photo_ids":SELECTION_N,
        "first_wave_taxa":len({r["inat_taxon_id"] for r in chosen}),
        "first_wave_distinct_source_observers":len({r["observer_id"] for r in chosen}),
        "first_wave_distinct_original_observation_IDs":len(
            {r["observation_id"] for r in chosen}),
        "first_wave_original_shards":N_SHARDS,
        "photos_per_original_shard":PHOTO_SHARD_SIZE,
        "first_wave_original_4096_csv_sha256":digest(sourcefile),
        "per_shard_sha256":shards,
        "no_photo_colour_or_fitness_or_bloom_label_used":True,
        "no_phenology_or_climate_or_geography_in_rank":True,
        "source_observer_disjoint_only_within_new_wave4":True,
        "new_photo_pixels_actually_measured":False,
        "FCP_genetic_withinspecies_D_W_local_50km_estimands_imported":False,
    }
    (outdir/FREEZE_OUT).write_text(json.dumps(
        receipt,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("MOIS_NEW_4096_BLIND_ORIGINAL_RGB_PHOTO_SOURCE_RESERVED",
          json.dumps({k:receipt[k] for k in (
              "prior_original_measured_photo_ids",
              "number_source_photo_ID_rows_scanned",
              "number_unmeasured_original_photo_ID_opportunities",
              "first_wave_original_source_photo_ids",
              "first_wave_taxa",
              "first_wave_distinct_source_observers")},sort_keys=True),flush=True)
    return receipt


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-3430502-original-photo-csv",type=Path,required=True)
    p.add_argument("--already-measured-309439-original-photo-csv-gz",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    run(a.source_3430502_original_photo_csv,
        a.already_measured_309439_original_photo_csv_gz,a.output)
