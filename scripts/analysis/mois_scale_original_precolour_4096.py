#!/usr/bin/env python3
"""Prepare a scalable, photo-ID-deduplicated colour-blind MOIS original ROI queue.

Reuse the EXACT original frozen 106,376-photo pre-colour 1,000-taxon
opportunity sample. Never read/choose hue outcomes from previous MOIS/FCP
sources. IDs are used solely to avoid needless re-downloading/reclassification.
The original available 3,430,502 photo IDs form a larger future pool; the
106,376 sample is a checkpoint, not a cap on source expansion.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

SOURCE_SHA="171d3fffdc518a6b025c430c6b809a6c1fe03654f4da913e47d77ea5bf40795b"
ORIGINAL_SELECTED=106376
ORIGINAL_TAXA=1000
BASELINE_PHOTOS=248312
ELIGIBLE_NEW=104904
WAVE_PHOTOS=4096
SHARDS=16
PHOTOS_PER_SHARD=256
SELECT_FIELDS=(
 "measurement_id","image_filename","photo_license","photo_url_large",
 "photo_id","inat_taxon_id","species","observation_id","observer_id",
 "latitude","longitude","screen_fold","selection_hash","source_snapshot_date",
 "wave_sequence","shard")
PHOTO_RE=re.compile(r"^https://inaturalist-open-data\.s3\.amazonaws\.com/photos/([0-9]+)/large\.(?:jpg|jpeg|png)$",re.I)
SEED="MOIS-expand-original-106376-first4096-source-ID-only-20261017-v1"


def hash_path(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""):
            h.update(chunk)
    return h.hexdigest()


def measurement(pid:str)->str:
    return "MOIS-"+hashlib.sha256((SEED+"|"+pid).encode()).hexdigest()[:24].upper()


def write_rows(path:Path,headers,rows):
    with path.open("x",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=headers,lineterminator="\n")
        writer.writeheader()
        writer.writerows({key:r.get(key,"") for key in headers} for r in rows)


def plan(source:Path,prior_gzip:Path,prior_receipt:Path,outdir:Path,
         *,selected=WAVE_PHOTOS)->dict:
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError("Cannot overwrite a precolour original photograph queue")
    if selected!=WAVE_PHOTOS:
        raise ValueError("First execution size is frozen at 4096, not species-colour selected")
    if hash_path(source)!=SOURCE_SHA:
        raise ValueError("Selected source ID photograph roster no longer matches frozen FCP SHA")
    doc=json.loads(prior_receipt.read_text(encoding="utf-8"))
    if (doc["final_photo_ids"]!=BASELINE_PHOTOS or
        hash_path(prior_gzip)!=doc["result_source_photo_csv_sha256"]):
        raise ValueError("Original MOIS measured data identity cannot be reconstructed")
    used=set()
    with gzip.open(prior_gzip,"rt",newline="",encoding="utf-8") as f:
        rdr=csv.DictReader(f)
        for r in rdr:
            pid=r["photo_id"]
            if pid in used:raise ValueError("Original MOIS photo union has duplicates")
            used.add(pid)
    if len(used)!=BASELINE_PHOTOS:
        raise ValueError("Original photo ID denominator drift")

    groups=defaultdict(list)
    available=[]
    old_overlap=0
    all_source=set()
    with source.open(newline="",encoding="utf-8") as f:
        rdr=csv.DictReader(f)
        need={"inat_taxon_id","species","photo_id","observation_id","observer_id",
              "latitude","longitude","photo_url_large","photo_license",
              "source_snapshot_date","screen_fold","selection_hash"}
        if not need.issubset(rdr.fieldnames or ()):
            raise ValueError("Original 20260927 source metadata columns unavailable")
        for r in rdr:
            pid=r["photo_id"]
            if pid in all_source:raise ValueError("Source photo-ID repeated")
            all_source.add(pid)
            if pid in used:
                old_overlap+=1
                continue
            if (not pid.isdecimal() or not r["inat_taxon_id"].isdecimal()
                or not r["observer_id"] or not r["observation_id"]
                or r["source_snapshot_date"]!="2026-09-27"):
                raise ValueError("Original photo metadata incomplete")
            url=PHOTO_RE.fullmatch(r["photo_url_large"])
            if url is None or url.group(1)!=pid or not r["photo_license"]:
                raise ValueError("Original iNaturalist source photo URL/license invalid")
            # Source opportunities, not flowering detections.
            available.append(r)
            groups[r["inat_taxon_id"]].append(r)
    if len(all_source)!=ORIGINAL_SELECTED or len(groups)!=ORIGINAL_TAXA:
        raise ValueError("Frozen source species and photo inclusion changed")
    if len(available)!=ELIGIBLE_NEW or len(used&all_source)!=old_overlap:
        raise ValueError("Previous classified photo exclusion does not reconcile")
    for rows in groups.values():
        rows.sort(key=lambda r:(r["selection_hash"],int(r["photo_id"])))
    species=sorted(groups,key=int)
    cursors={tax:0 for tax in species}
    global_observers=set()
    wave=[]
    for turn in range(12):
        if len(wave)>=selected:break
        for tax in species:
            pool=groups[tax]
            idx=cursors[tax]
            while idx<len(pool) and pool[idx]["observer_id"] in global_observers:
                idx+=1
            if idx>=len(pool):
                cursors[tax]=idx
                continue
            r=pool[idx]
            cursors[tax]=idx+1
            global_observers.add(r["observer_id"])
            wave.append(r)
            if len(wave)>=selected:break
    if len(wave)!=selected or len({r["photo_id"] for r in wave})!=selected:
        raise ValueError("Cannot form required colour-blind 4096 original source photos")
    if len(global_observers)!=selected or len({r["inat_taxon_id"] for r in wave})!=ORIGINAL_TAXA:
        raise ValueError("First wave must have 1000 taxa and 4096 independent photographers")
    outputs=[]
    for n,r in enumerate(wave):
        row=dict(r)
        row["wave_sequence"]=n
        row["shard"]=n%SHARDS
        row["measurement_id"]=measurement(r["photo_id"])
        row["image_filename"]=row["measurement_id"]+".jpg"
        outputs.append(row)
    outdir.mkdir(parents=True,exist_ok=True)
    selected_file=outdir/"mois_source_only_4096_original_photo_selection.csv"
    write_rows(selected_file,SELECT_FIELDS,outputs)
    shard_sha={}
    for k in range(SHARDS):
        subset=[r for r in outputs if r["shard"]==k]
        if len(subset)!=PHOTOS_PER_SHARD:
            raise ValueError("Unequal fixed independent original photo compute shard")
        target=outdir/f"mois_source_only_shard_{k:02d}_256.csv"
        write_rows(target,SELECT_FIELDS,subset)
        shard_sha[str(k)]=hash_path(target)
    # Preserve whole remaining 106,376-photo outcome-blind resource (all
    # licences, original IDs, coordinates). This is a durable work backlog.
    full=outdir/"mois_unmeasured_104904_original_source_queue.csv.gz"
    with full.open("wb") as raw:
        import io
        with gzip.GzipFile(filename="",fileobj=raw,mode="wb",mtime=0) as gf:
            with io.TextIOWrapper(gf,encoding="utf-8",newline="") as f:
                w=csv.DictWriter(f,fieldnames=(
                    "inat_taxon_id","species","photo_id","observation_id",
                    "observer_id","latitude","longitude","photo_url_large",
                    "photo_license","source_snapshot_date","screen_fold",
                    "selection_hash"),lineterminator="\n")
                w.writeheader()
                w.writerows({key:r[key] for key in w.fieldnames} for r in available)
    receipt={
       "schema":"mois_original_photo_precolour_1000taxa_unmeasured_queue_first4096_v1",
       "status":"SOURCE_ONLY_4096_RESERVED_NO_NEW_PIXEL_HUES",
       "authoritative_original_source_snapshot":"2026-09-27",
       "original_FCP_precolour_sample_photo_ids":ORIGINAL_SELECTED,
       "original_photo_source_species":ORIGINAL_TAXA,
       "already_measured_MOIS_original_photo_IDs_excluded_from_remeasurement":old_overlap,
       "actual_unmeasured_1000taxon_preselected_photo_ids":len(available),
       "first_wave_original_source_photo_ids":len(wave),
       "first_wave_distinct_source_observers":len(global_observers),
       "first_wave_taxa":len({r["inat_taxon_id"] for r in wave}),
       "first_wave_original_shards":SHARDS,
       "photos_per_original_shard":PHOTOS_PER_SHARD,
       "per_shard_sha256":shard_sha,
       "original_source_selection_sha256":SOURCE_SHA,
       "existing_measured_photo_union_sha256":hash_path(prior_gzip),
       "full_unmeasured_source_queue_sha256":hash_path(full),
       "first_wave_original_4096_csv_sha256":hash_path(selected_file),
       "no_photo_colour_or_fitness_or_bloom_label_used":True,
       "historical_MOIS_original_photo_colour_rows_not_imported_into_fresh_evaluator":True,
       "future_source_expansion_potential_original_1000_species_IDs":3430502,
       "source_3430502_or_full_93million_photo_ids_measured_yet":False,
       "subsequent_batch_no_biological_photo_cap":True,
    }
    (outdir/"mois_first4096_original_source_only_freeze.json").write_text(
        json.dumps(receipt,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("MOIS_NEW_SOURCE_PHOTOS_4096_FROZEN",json.dumps({
        k:receipt[k] for k in (
        "already_measured_MOIS_original_photo_IDs_excluded_from_remeasurement",
        "actual_unmeasured_1000taxon_preselected_photo_ids",
        "first_wave_original_source_photo_ids",
        "first_wave_taxa")},sort_keys=True))
    return receipt


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--frozen-original-precolour-106376",required=True,type=Path)
    p.add_argument("--already-measured-248312",required=True,type=Path)
    p.add_argument("--already-measured-receipt",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    plan(a.frozen_original_precolour_106376,a.already_measured_248312,
         a.already_measured_receipt,a.output)
