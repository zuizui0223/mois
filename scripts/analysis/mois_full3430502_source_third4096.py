#!/usr/bin/env python3
"""MOIS original 3,430,502-photo census: measure beyond the first 106,376 pool.

Outcome-blind source accounting only. Identifies all genuinely unmeasured,
photo-ID-disjoint original photos while allowing the FCP *photographs* to
be reused between studies. Assigns a reproducible third 4,096-photo
cross-species work wave OUTSIDE the original first 106,376 selected photo IDs.
No hue outcome, iNaturalist flowering annotation, or rainfall value is used.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import heapq
import json
import re
from collections import Counter,defaultdict
from pathlib import Path

FULL_SHA="58e037a210e23d3ff5d6f04d392665ada5ee9ac5149d0676672c239cf092df01"
FIRST_SELECTED_SHA="171d3fffdc518a6b025c430c6b809a6c1fe03654f4da913e47d77ea5bf40795b"
SOURCE_N=3430502
ALREADY_MEASURED_N=248312
SCREEN_SELECTED_N=106376
ALREADY_MEASURED_OVERLAP_N=9688
UNMEASURED_FULL_N=3420814
UNMEASURED_OUTSIDE_FIRST_SAMPLE_N=3315910
WAVE_N=4096
N_SHARDS=16
PHOTOS_PER_SHARD=256
TOPK_PER_SPECIES=32
SALT="MOIS_SOURCE_3430502_FULL_1000_TAXA_WAVE3_HASH_PHOTOS_20261017"
LINK_RE=re.compile(r"^https://inaturalist-open-data\.s3\.amazonaws\.com/photos/([0-9]+)/large\.(?:jpg|jpeg|png)$",re.I)
FIELDS=(
    "measurement_id","image_filename","photo_license","photo_url_large",
    "photo_id","inat_taxon_id","species","observation_id","observer_id",
    "latitude","longitude","observed_on","positional_accuracy_m",
    "screen_fold","selection_hash","source_snapshot_date","wave_sequence","shard"
)


def digest(p:Path)->str:
    h=hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""):
            h.update(chunk)
    return h.hexdigest()


def save(path:Path,rows:list[dict]):
    with path.open("x",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=FIELDS,lineterminator="\n")
        w.writeheader()
        w.writerows({key:r.get(key,"") for key in FIELDS} for r in rows)


def measure_id(pid):
    return "MOIS-"+hashlib.sha256((SALT+"|"+pid).encode()).hexdigest()[:24].upper()


def plan(full:Path,selected:Path,prior:Path,prior_receipt:Path,wave2:Path,wave2_receipt:Path,outdir:Path)->dict:
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError("No overwriting a colour-blind true original source wave")
    if digest(full)!=FULL_SHA or digest(selected)!=FIRST_SELECTED_SHA:
        raise ValueError("Untrusted original 20260927 full iNaturalist metadata")
    rec=json.loads(prior_receipt.read_text())
    if (rec["final_photo_ids"]!=ALREADY_MEASURED_N or
        digest(prior)!=rec["result_source_photo_csv_sha256"]):
        raise ValueError("Original 248312 photograph measured source changed")
    seen_old=set()
    with gzip.open(prior,"rt",newline="",encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["photo_id"] in seen_old:
                raise ValueError("Repeated original classified photo ID")
            seen_old.add(r["photo_id"])
    if len(seen_old)!=ALREADY_MEASURED_N:
        raise ValueError("Original measured ID count changed")

    screened=set()
    with selected.open(newline="",encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["photo_id"] in screened:
                raise ValueError("Repeated previously screened original source ID")
            screened.add(r["photo_id"])
    if len(screened)!=SCREEN_SELECTED_N:
        raise ValueError("Previously selected 106376 source photo frame changed")

    source2=json.loads(wave2_receipt.read_text(encoding="utf-8"))
    if (source2.get("schema")!="mois_full3430502_original_source_outside_first106376_wave2_v1"
        or source2.get("first_wave_original_source_photo_ids")!=4096
        or source2.get("first_wave_original_4096_csv_sha256")!=digest(wave2)
        or source2.get("no_photo_colour_or_fitness_or_bloom_label_used") is not True):
        raise ValueError("Previous original second wave not SHA authenticated")
    wave2_ids=set()
    with wave2.open(newline="",encoding="utf-8") as inp:
        for row in csv.DictReader(inp):
            pid=row["photo_id"]
            if pid in wave2_ids or pid in seen_old or pid in screened:
                raise ValueError("Second wave duplicates old or photo-screen ID")
            wave2_ids.add(pid)
    if len(wave2_ids)!=4096:
        raise ValueError("Previous second 4096-photo source is incomplete")
    heaps=defaultdict(list)
    taxa=Counter()
    counts=Counter()
    encountered_screen=set()
    with full.open(newline="",encoding="utf-8") as f:
        rdr=csv.DictReader(f)
        required={"inat_taxon_id","species","photo_id","observation_id",
                  "observer_id","observed_on","latitude","longitude",
                  "photo_url_large","photo_license","source_snapshot_date"}
        if not required.issubset(rdr.fieldnames or ()):
            raise ValueError("Signed full original photo source schema changed")
        for r in rdr:
            counts["all"]+=1
            pid=r["photo_id"]
            if not pid.isdecimal() or not r["inat_taxon_id"].isdecimal():
                raise ValueError("Non-numeric original source photo/taxon ID")
            if pid in screened:
                encountered_screen.add(pid)
                counts["in_first_106376"]+=1
            if pid in seen_old:
                counts["already_measured"]+=1
                continue
            counts["not_previously_measured"]+=1
            if pid in screened:continue
            counts["not_measured_outside_first_photo_screen"]+=1
            if pid in wave2_ids:
                counts["previously_reserved_wave2"]+=1
                continue
            taxa[r["inat_taxon_id"]]+=1
            if (r["source_snapshot_date"]!="2026-09-27"
                or not r["observer_id"] or not r["observation_id"]
                or not r["species"] or not r["latitude"] or not r["longitude"]
                or not r["observed_on"]):
                counts["source_identity_invalid_for_image_sampling"]+=1
                continue
            u=LINK_RE.fullmatch(r["photo_url_large"])
            if u is None or u.group(1)!=pid or not r["photo_license"]:
                counts["photo_source_license_or_url_unusable"]+=1
                continue
            rank=int(hashlib.sha256((SALT+"|"+pid).encode()).hexdigest()[:16],16)
            heap=heaps[r["inat_taxon_id"]]
            item=(-rank,pid,{
                "photo_id":pid,"inat_taxon_id":r["inat_taxon_id"],
                "species":r["species"],"observation_id":r["observation_id"],
                "observer_id":r["observer_id"],
                "latitude":r["latitude"],"longitude":r["longitude"],
                "observed_on":r["observed_on"],
                "positional_accuracy_m":r.get("positional_accuracy_m",""),
                "photo_license":r["photo_license"],
                "photo_url_large":r["photo_url_large"],
                "source_snapshot_date":r["source_snapshot_date"],
                "screen_fold":"not_in_first_106376_source_wave",
                "selection_hash":f"{rank:016x}",
            })
            if len(heap)<TOPK_PER_SPECIES:
                heapq.heappush(heap,item)
            elif item>heap[0]:
                heapq.heapreplace(heap,item)
    if (counts["all"]!=SOURCE_N or counts["already_measured"]!=ALREADY_MEASURED_OVERLAP_N
        or counts["not_previously_measured"]!=UNMEASURED_FULL_N
        or counts["not_measured_outside_first_photo_screen"]!=UNMEASURED_OUTSIDE_FIRST_SAMPLE_N
        or len(encountered_screen)!=SCREEN_SELECTED_N or len(taxa)!=1000
        or counts["previously_reserved_wave2"]!=4096):
        raise ValueError("Expanded 3.43M original photo source/candidate census did not reproduce")
    pools={tax:[row for rank,pid,row in sorted(
        heap,key=lambda t:(-t[0],int(t[1])))] for tax,heap in heaps.items()}
    if len(pools)!=1000:
        raise ValueError("Not all 1000 source taxa have original image opportunities")
    cursors={tax:0 for tax in pools}
    chosen=[];global_observers=set();species=sorted(pools,key=int)
    for iteration in range(12):
        for tax in species:
            if len(chosen)>=WAVE_N:break
            source=pools[tax]
            idx=cursors[tax]
            while idx<len(source) and source[idx]["observer_id"] in global_observers:
                idx+=1
            cursors[tax]=idx+1
            if idx>=len(source):continue
            row=source[idx].copy()
            if row["photo_id"] in seen_old or row["photo_id"] in screened or row["photo_id"] in wave2_ids:
                raise ValueError("Third wave reused old, source-screen or second-wave photo")
            global_observers.add(row["observer_id"])
            chosen.append(row)
        if len(chosen)>=WAVE_N:break
    if (len(chosen)!=WAVE_N or len(global_observers)!=WAVE_N
        or len({r["inat_taxon_id"] for r in chosen})!=1000):
        raise ValueError("Observer-disjoint four-photo-per-taxon source wave3 unavailable")
    original=[]
    for idx,r in enumerate(chosen):
        row=dict(r)
        row["wave_sequence"]=idx
        row["shard"]=idx%N_SHARDS
        row["measurement_id"]=measure_id(r["photo_id"])
        row["image_filename"]=row["measurement_id"]+".jpg"
        original.append(row)
    outdir.mkdir(parents=True,exist_ok=True)
    allfile=outdir/"mois_original_full_3430502_outside_first_screen_WAVE3_4096.csv"
    save(allfile,original)
    shards={}
    for k in range(N_SHARDS):
        path=outdir/f"mois_source_only_shard_{k:02d}_256.csv"
        x=[r for r in original if r["shard"]==k]
        if len(x)!=PHOTOS_PER_SHARD:
            raise ValueError("Invalid source photo wave2 executor partition")
        save(path,x)
        shards[str(k)]=digest(path)
    top=outdir/"mois_original_full_3430502_by_taxon_unmeasured_capacities.csv"
    with top.open("x",newline="",encoding="utf-8") as f:
        w=csv.writer(f,lineterminator="\n")
        w.writerow(("inat_taxon_id","unmeasured_photo_ids_outside_initial_106376"))
        w.writerows((tax,taxa[tax]) for tax in sorted(taxa,key=int))
    receipt={
      "schema":"mois_full3430502_original_source_outside_first106376_wave3_v1",
      "status":"SOURCE_ONLY_4096_RESERVED_NO_NEW_PIXEL_HUES",
      "authoritative_original_source_snapshot":"2026-09-27",
      "original_complete_source_photo_ids":counts["all"],
      "previously_measured_photo_id_overlap":counts["already_measured"],
      "new_original_unmeasured_source_photo_ids":counts["not_previously_measured"],
      "previous_original_106376_photo_screen_ids":len(encountered_screen),
      "new_unmeasured_source_photo_ids_outside_first_screen":counts["not_measured_outside_first_photo_screen"],
      "previous_wave2_photo_ids_excluded":counts["previously_reserved_wave2"],
      "unmeasured_outside_first_screen_after_wave2_exclusion":counts["not_measured_outside_first_photo_screen"]-counts["previously_reserved_wave2"],
      "previous_wave2_original_source_sha256":digest(wave2),
      "original_taxa_with_new_source_photos":len(taxa),
      "full_3p43M_original_source_file_sha256":FULL_SHA,
      "first_106376_source_file_sha256":FIRST_SELECTED_SHA,
      "prior_248312_measured_photo_ids_sha256":digest(prior),
      "full_unmeasured_source_taxon_counts_sha256":digest(top),
      "first_wave_original_source_photo_ids":len(chosen),
      "first_wave_taxa":len(species),
      "first_wave_distinct_source_observers":len(global_observers),
      "first_wave_original_shards":N_SHARDS,
      "photos_per_original_shard":PHOTOS_PER_SHARD,
      "first_wave_original_4096_csv_sha256":digest(allfile),
      "per_shard_sha256":shards,
      "source_invalid_id_date_site":counts["source_identity_invalid_for_image_sampling"],
      "source_unavailable_original_photo_license_or_URL":counts["photo_source_license_or_url_unusable"],
      "no_photo_colour_or_fitness_or_bloom_label_used":True,
      "FCP_photo_source_overlap_allowed_but_FCP_polymorphism_models_not_imported":True,
      "all_3420814_unmeasured_photos_ALREADY_classified":False
    }
    (outdir/"mois_first4096_original_source_only_freeze.json").write_text(
        json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print("MOIS_FULL_ORIGINAL_3430502_SOURCE_3420814_UNMEASURED_SECOND4096_FROZEN",
          json.dumps({"from_source":counts["all"],"new_unmeasured":counts["not_previously_measured"],
                      "outside_prior_screen":counts["not_measured_outside_first_photo_screen"],
                      "reserved_wave3":len(chosen)},sort_keys=True),flush=True)
    return receipt


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original-full3430502",required=True,type=Path)
    p.add_argument("--original-first106376",required=True,type=Path)
    p.add_argument("--already-measured-248312",required=True,type=Path)
    p.add_argument("--already-measured-receipt",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    plan(a.original_full3430502,a.original_first106376,
         a.already_measured_248312,a.already_measured_receipt,a.output)
