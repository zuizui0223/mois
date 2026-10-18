#!/usr/bin/env python3
"""MOIS fifth original ROI-v4 4096-photo source, strictly outcome-blind.

Choose photograph IDs from the exact 2026-09-27 original 3430502 FCP
photo-ID source, excluding already genuinely measured MOIS 313535 photo IDs
and the separate targeted4096 source reservation. 4096 reserved IDs are NOT
yet measured colour outcomes. Fixed 1000-taxon, globally observer-distinct
frozen source with deterministic hash before any pixel decoding.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import heapq
import json
import math
from collections import Counter,defaultdict
from datetime import date
from pathlib import Path

from scripts.analysis.mois_full3430502_source_third4096 import (
    FULL_SHA,SOURCE_N,FIELDS,LINK_RE,save,
)
from scripts.analysis.mois_reassemble_original_new4096_measurements import HUES

N_PRIOR=313535
N_CLASSIFIED=151239
N_DATED=151230
N_NEW=4096
N_SHARDS=16
SHARD_SIZE=256
TAXA=1000
MAX_PER_TAXON=64
SALT="MOIS_GLOBAL_FULL3430502_BLIND_FIFTH4096_AFTER_REAL313535_AND_TARGETED4096_20261018_V1"
SOURCE="mois_new4096_original_photo_source_after313535_outcome_blind.csv"
FREEZE="mois_new4096_photo_source_after313535_frozen_receipt.json"
PREVIOUS="mois_full313535_original_pixel_ROIs_and_observation_dates.csv.gz"
PREV_RECEIPT="mois_fourth_4096_actual_original_ROI_pixels_and_source_dates_decision.json"
TARGET="mois_targeted4096_photo_original_census_blind_selection.csv"
TARGET_RECEIPT="mois_first4096_original_source_only_freeze.json"
KNOWN_ORIGINAL_309439_SHA="546266e2bbaa340a2523864c28c5a877842eb4b0904a7e0dff90870c6519a47e"


def need(ok,why):
    if not ok:raise ValueError(why)


def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):h.update(chunk)
    return h.hexdigest()


def parse_old(old,oldreceipt):
    r=json.loads(oldreceipt.read_text())
    need(r["schema"]=="mois_actual_fourth4096_genuine_original_RGB_ROI_309439plus4096_v1"
         and r["initial_measured_original_photo_GZIP_SHA256"]==KNOWN_ORIGINAL_309439_SHA
         and r["expanded_genuine_original_photo_ids"]==N_PRIOR
         and r["expanded_original_four_hue_ROI_photo_ids"]==N_CLASSIFIED
         and r["expanded_actual_date_and_four_hue_ROI_photo_ids"]==N_DATED
         and r["source_original_photo_pixels_actually_processed"] is True
         and r["current_fourth4096_source_SHA256"]==digest(old),
         "Prior 313535 true pixel-record photo/date/colour signatures not verified")
    ids=set();known=dated=0
    with gzip.open(old,"rt",newline="",encoding="utf-8") as inp:
        reader=csv.DictReader(inp)
        need({"photo_id","inat_taxon_id","morph","qualified_date"}.issubset(
            reader.fieldnames or ()), "Original measured photo schema changed")
        for row in reader:
            pid=row["photo_id"]
            need(pid.isdecimal() and pid not in ids,
                 "Original previously processed photo ID duplicated")
            ids.add(pid)
            labelled=row["morph"] in HUES
            known+=labelled
            dated+=labelled and len(row["qualified_date"])==10
    need((len(ids),known,dated)==(N_PRIOR,N_CLASSIFIED,N_DATED),
         "Original previously image-measured colour and date totals changed")
    return ids,r


def reserved_target(path,receipt):
    r=json.loads(receipt.read_text())
    need(r["schema"]==
         "mois_3430502_targeted_opposite_season_after309439_and_reserved4096_v1"
         and r["source_photo_ids_original"]==SOURCE_N
         and r["first_wave_original_source_photo_ids"]==4096
         and r["first_wave_original_shards"]==16
         and r["no_photo_colour_or_fitness_or_bloom_label_used"] is True
         and r["new_photo_pixels_actually_measured"] is False
         and r["first_wave_original_4096_csv_sha256"]==digest(path),
         "Targeted source-only 4096 photo reservation not authenticated")
    excluded=set()
    with path.open(newline="",encoding="utf-8") as f:
        for x in csv.DictReader(f):
            need(x["photo_id"].isdecimal() and x["photo_id"] not in excluded,
                 "Targeted original photo reservation photo IDs duplicated")
            excluded.add(x["photo_id"])
    need(len(excluded)==4096,"Targeted 4096 reserved photos disappeared")
    return excluded,r


def valid(row):
    pid=row.get("photo_id","")
    if not (pid.isdecimal() and row.get("inat_taxon_id","").isdecimal()
            and row.get("source_snapshot_date")=="2026-09-27"
            and all(row.get(f,"") for f in (
                "species","observer_id","observation_id","photo_license",
                "photo_url_large","latitude","longitude","observed_on"))):
        return False
    match=LINK_RE.fullmatch(row["photo_url_large"])
    if match is None or match.group(1)!=pid:
        return False
    try:
        lat=float(row["latitude"]);lon=float(row["longitude"])
        day=date.fromisoformat(row["observed_on"])
        return (math.isfinite(lat) and math.isfinite(lon)
                and -90<=lat<=90 and -180<=lon<=180
                and date(1800,1,1)<=day<=date(2026,10,17))
    except (ValueError,TypeError):
        return False


def choose(rows,prior,reserved):
    heaps=defaultdict(list);counts=Counter()
    for r in rows:
        counts["provider_original_photo_ID_rows"]+=1
        pid=r.get("photo_id","")
        need(pid.isdecimal(),"Original 3.43M provider source ID malformed")
        if pid in prior:
            counts["photo_already_real_RGB_measured"]+=1
            continue
        if pid in reserved:
            counts["photo_reserved_for_season_target"]+=1
            continue
        if not valid(r):
            counts["source_metadata_unusable"]+=1
            continue
        counts["candidate_photo_IDs_with_valid_source"]+=1
        digest_key=hashlib.sha256((SALT+"|"+pid).encode()).hexdigest()
        rank=int(digest_key[:16],16)
        entry={k:r.get(k,"") for k in FIELDS}
        entry.update({
            "screen_fold":"outside_all_original_313535_measured_and_targeted4096",
            "selection_hash":f"{rank:016x}",
        })
        tax=r["inat_taxon_id"]
        item=(-rank,pid,entry)
        heap=heaps[tax]
        if len(heap)<MAX_PER_TAXON:
            heapq.heappush(heap,item)
        elif item>heap[0]:
            heapq.heapreplace(heap,item)
    need(counts["provider_original_photo_ID_rows"]==SOURCE_N,
         "Full source photograph provider original 3430502 IDs incomplete")
    need(len(heaps)==TAXA,"Taxon source changed: no candidate in all 1000 species")
    order=sorted(heaps,key=lambda k:hashlib.sha256(
        (SALT+"|TAXON|"+k).encode()).hexdigest())
    pools={t:[item[-1] for item in sorted(v,key=lambda w:(-w[0],int(w[1])))]
           for t,v in heaps.items()}
    cursor={t:0 for t in order}
    selected=[];ids=set();observers=set();obs=set()
    for repeat in range(MAX_PER_TAXON):
        for tax in order:
            if len(selected)>=N_NEW:break
            q=pools[tax];i=cursor[tax]
            while i<len(q) and (q[i]["photo_id"] in ids or
                                q[i]["observation_id"] in obs or
                                q[i]["observer_id"] in observers):
                i+=1
            cursor[tax]=i+1
            if i>=len(q):continue
            row=dict(q[i])
            selected.append(row)
            ids.add(row["photo_id"])
            obs.add(row["observation_id"])
            observers.add(row["observer_id"])
        if len(selected)>=N_NEW:break
    need(len(selected)==len(ids)==len(obs)==len(observers)==N_NEW
         and len({r["inat_taxon_id"] for r in selected})==TAXA,
         "Failed to find full 4096 genuinely new observer-distinct photo IDs")
    need(not (ids&prior) and not (ids&reserved),
         "Fifth original wave overlaps prior measured or reserved season source")
    for i,r in enumerate(selected):
        mid="MOIS-"+hashlib.sha256((SALT+"|"+r["photo_id"]).encode()
                                  ).hexdigest()[:24].upper()
        r.update({"measurement_id":mid,"image_filename":mid+".jpg",
                  "wave_sequence":str(i),"shard":str(i%N_SHARDS)})
    return selected,dict(counts)


def run(full,previous,previousreceipt,seasonplan,seasonreceipt,out):
    need(not out.exists() or not any(out.iterdir()),
         "Never overwrite frozen original RGB-photo source reservations")
    need(digest(full)==FULL_SHA,
         "Original 3,430,502 source photo IDs SHA drift")
    old,rprev=parse_old(previous,previousreceipt)
    excluded,rseason=reserved_target(seasonplan,seasonreceipt)
    need(not (old & excluded),
         "Season-targeted reserved 4096 photo IDs overlap prior 313535 measured")
    with full.open(newline="",encoding="utf-8") as f:
        selected,counts=choose(csv.DictReader(f),old,excluded)
    out.mkdir(parents=True,exist_ok=True)
    source=out/SOURCE
    save(source,selected)
    shards={}
    for k in range(N_SHARDS):
        photo=[r for r in selected if int(r["shard"])==k]
        need(len(photo)==SHARD_SIZE,"Per-worker photo source must be 256")
        shard=out/f"mois_source_only_shard_{k:02d}_256.csv"
        save(shard,photo)
        shards[str(k)]=digest(shard)
    receipt={
        "schema":"mois_next4096_after_true313535_photo_ID_blind_source_v1",
        "status":"SOURCE_ONLY_4096_RESERVED_NO_NEW_PIXEL_HUES",
        "source_photo_ids_original":SOURCE_N,
        "original_source_3430502_SHA256":FULL_SHA,
        "prior_full313535_source_SHA256":digest(previous),
        "prior_real_original_photo_IDs":N_PRIOR,
        "prior_four_hue_classifiable":N_CLASSIFIED,
        "prior_original_dated_four_hue":N_DATED,
        "separate_targeted4096_reserved_source_SHA256":digest(seasonplan),
        "targeted4096_excluded_as_reserved_NOT_independent_world_sample":True,
        "source_photo_IDs_excluded_old_and_targeted":True,
        "number_source_photo_ID_rows_scanned":counts[
            "provider_original_photo_ID_rows"],
        "eligible_photo_IDs_source_valid_before_actual_RGB":counts[
            "candidate_photo_IDs_with_valid_source"],
        "new_source_photo_IDs":N_NEW,
        "first_wave_original_source_photo_ids":N_NEW,
        "first_wave_original_shards":N_SHARDS,
        "photos_per_original_shard":SHARD_SIZE,
        "first_wave_taxa":len({r["inat_taxon_id"] for r in selected}),
        "first_wave_distinct_source_observers":len({
            r["observer_id"] for r in selected}),
        "first_wave_distinct_original_observation_IDs":len({
            r["observation_id"] for r in selected}),
        "first_wave_original_4096_csv_sha256":digest(source),
        "per_shard_sha256":shards,
        "new_photo_pixels_actually_measured":False,
        "no_photo_colour_or_fitness_or_bloom_label_used":True,
        "rank_blind_to_colour_climate_phenology_and_roi_outcome":True,
        "different_source_recruitment_mechanisms_not_identical_flora_prevalence":True,
        "no_FCP_intraspecific_colour_polymorphism_estimators":True,
    }
    (out/FREEZE).write_text(json.dumps(
        receipt,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("MOIS_FIFTH4096_ORIGINAL_IMAGE_PIXEL_SOURCE_ONLY",
          json.dumps({k:receipt[k] for k in (
              "prior_real_original_photo_IDs",
              "new_source_photo_IDs",
              "first_wave_taxa",
              "first_wave_distinct_source_observers",
              "source_photo_IDs_excluded_old_and_targeted")},sort_keys=True),
          flush=True)
    return receipt

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--full3430502-original-photo-source",type=Path,required=True)
    p.add_argument("--actually-measured313535-source",type=Path,required=True)
    p.add_argument("--measured313535-original-source-receipt",type=Path,required=True)
    p.add_argument("--season4096-reservation",type=Path,required=True)
    p.add_argument("--season4096-reservation-receipt",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    run(a.full3430502_original_photo_source,a.actually_measured313535_source,
        a.measured313535_original_source_receipt,
        a.season4096_reservation,a.season4096_reservation_receipt,a.output)
