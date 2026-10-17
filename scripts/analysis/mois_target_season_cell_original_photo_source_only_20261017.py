#!/usr/bin/env python3
"""Reserve 4096 new original photo IDs to test photo-season coverage, pixels later.

The 2026-09-27 source contains 3,430,502 actual PHOTO-IDENTIFIER candidates,
not real classified flowers. Exclude all 309,439 previously pixel-processed IDs
and all 4096 currently reserved in the independently running fourth wave.
Target opposite photographed seasons in an exact named taxon × 0.25-degree
cell × meteorological season-year, using ONLY prior photo-classifiability
presence, not white/color labels, and unseen candidates' taxon/site/day.

Reserve 16x256 genuine blind worker shards; only shard0 is processed in this
bounded stage. All other 3840 source IDs stay SOURCE ONLY, never counted
measured, no source photographic propensity/prevalence inference.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import heapq
import io
import json
import math
import zipfile
from collections import Counter,defaultdict
from datetime import date
from pathlib import Path

from scripts.analysis.mois_full3430502_source_third4096 import FULL_SHA,SOURCE_N,LINK_RE,FIELDS
from scripts.analysis.mois_reassemble_original_new4096_measurements import HUES

PHOTO_BASE_SHA="546266e2bbaa340a2523864c28c5a877842eb4b0904a7e0dff90870c6519a47e"
FULL_BASE_PHOTOS=309439
FULL_BASE_CLASSIFIED=149662
FULL_BASE_DATE_HUES=149653
PRIOR_WAVE4_RESERVED_SHA="65fc139909f94448e042b6539c6749e70d1ab21bb8b13f50383bb7417180337a"
RESERVE_N=4096
SELECT_ONLY_PIXEL_SHARD=0
N_SHARDS=16
SHARD_SIZE=256
SALT="MOIS_TARGETED_PREVIOUSLY_CLASSIFIABLE_PHOTO_OPPOSITE_SEASON_CELL_YEAR_SOURCE_20261017_V1"
PHOTO_SRC="mois_targeted4096_photo_original_census_blind_selection.csv"
FREEZE_OUT="mois_first4096_original_source_only_freeze.json"
COLS=("measurement_id","image_filename","photo_license","photo_url_large",
      "photo_id","inat_taxon_id","species","observation_id","observer_id",
      "latitude","longitude","observed_on","positional_accuracy_m",
      "screen_fold","selection_hash","source_snapshot_date","wave_sequence",
      "shard","target_priority","prior_classifiable_summer_obs","prior_classifiable_winter_obs")

def require(ok,why):
    if not ok:raise ValueError(why)

def sha(path:Path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):h.update(b)
    return h.hexdigest()

def save_csv(path:Path,records):
    with path.open("x",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=COLS,lineterminator="\n")
        w.writeheader()
        w.writerows(records)

def locus(day:str,lat,lon,tax:str):
    try:
        a=float(lat);b=float(lon)
        d=date.fromisoformat(day)
        if not (-90<=a<=90 and -180<=b<=180) or not (
            date(1800,1,1)<=d<=date(2026,10,17)):return None
        if abs(a)<=23.5:return None
        winter=(d.month in (12,1,2) if a>0 else d.month in (6,7,8))
        summer=(d.month in (6,7,8) if a>0 else d.month in (12,1,2))
        if not winter and not summer:return None
        season="winter" if winter else "summer"
        season_year=d.year+int(d.month==12)
        return (str(tax),math.floor((a+90)/.25),
                math.floor((b+180)/.25),season_year),season
    except (ValueError,TypeError):
        return None

def known_support(old:Path):
    require(sha(old)==PHOTO_BASE_SHA,"Prior actually measured 309439 photo/date SHA drift")
    already=set()
    observed=defaultdict(lambda:{"summer":set(),"winter":set()})
    oldphoto=0; oldcolour=0; oldday=0
    with gzip.open(old,"rt",newline="",encoding="utf-8") as f:
        for row in csv.DictReader(f):
            pid=row["photo_id"]
            require(pid not in already and pid.isdecimal(),"Measured original photo ID duplicated")
            already.add(pid);oldphoto+=1
            iscolour=row["morph"] in HUES
            oldcolour+=iscolour
            oldday+=iscolour and len(row["qualified_date"])==10
            if not iscolour or len(row["qualified_date"])!=10:continue
            p=locus(row["qualified_date"],row["latitude"],
                    row["longitude"],row["inat_taxon_id"])
            if p is not None:
                group,season=p
                observed[group][season].add(row["observation_id"])
    require((oldphoto,oldcolour,oldday)==
            (FULL_BASE_PHOTOS,FULL_BASE_CLASSIFIED,FULL_BASE_DATE_HUES),
            "Original prior photo/ROI/date baseline changed")
    counts={key:{s:len(groups[s]) for s in ("summer","winter")}
            for key,groups in observed.items()}
    pair1=sum(q["summer"]>=1 and q["winter"]>=1 for q in counts.values())
    pair2=sum(q["summer"]>=2 and q["winter"]>=2 for q in counts.values())
    return already,counts,{"source_taxon_cell_season_year_with_at_least_1_photo_both_summer_winter":pair1,
                           "source_taxon_cell_season_year_with_2_observation_IDs_both_seasons":pair2}

def check_reserved_fourth(plan:Path):
    with zipfile.ZipFile(plan) as z:
        names=z.namelist()
        source=next((n for n in names if n.endswith(
            "mois_full3430502_next4096_original_photo_IDs_AFTER309439_WAVE4.csv")),None)
        frozen=next((n for n in names if n.endswith(FREEZE_OUT)),None)
        require(bool(source) and bool(frozen),
                "Fourth original 4096 photo-ID reservation source unavailable")
        src=z.read(source)
        rec=json.loads(z.read(frozen))
    require(rec["first_wave_original_4096_csv_sha256"]==PRIOR_WAVE4_RESERVED_SHA
            and hashlib.sha256(src).hexdigest()==PRIOR_WAVE4_RESERVED_SHA
            and rec["source_measured_original_309439_SHA256"]==PHOTO_BASE_SHA
            and rec["first_wave_original_source_photo_ids"]==4096
            and rec["new_photo_pixels_actually_measured"] is False,
            "Never overlap fourth parallel real ROI worker source reserve")
    fourth=list(csv.DictReader(io.StringIO(src.decode())))
    ids={r["photo_id"] for r in fourth}
    require(len(ids)==4096,"Fourth frozen photo-ID source duplicated")
    return ids

def identity(photo:str):
    return hashlib.sha256((SALT+"|"+photo).encode()).hexdigest()

def push_smallest(heap:list,rank:str,row:dict,limit:int):
    key=(-int(rank[:16],16),row["photo_id"],row)
    if len(heap)<limit:heapq.heappush(heap,key)
    elif key>heap[0]:heapq.heapreplace(heap,key)

def scan(full:Path,used:set[str],old_pairs:dict):
    require(sha(full)==FULL_SHA,"Original 3430502 source-photo census SHA changed")
    top=defaultdict(list);general=defaultdict(list)
    stats=Counter()
    with full.open(newline="",encoding="utf-8") as f:
        reader=csv.DictReader(f)
        require(set(("photo_id","inat_taxon_id","observer_id","observation_id",
                     "observed_on","source_snapshot_date","photo_url_large")).issubset(
                    reader.fieldnames or ()),
                "Original licensed source metadata/photo URL contract drift")
        for r in reader:
            stats["source_photo_ID_candidates_scanned"]+=1
            pid=r["photo_id"]
            if not pid.isdecimal() or not r["inat_taxon_id"].isdecimal():
                raise ValueError("Invalid source numeric original photo/taxon")
            if pid in used:
                stats["prior_measured_or_reserved_photo_IDs_excluded"]+=1
                continue
            try:
                dy=date.fromisoformat(r["observed_on"])
                a=float(r["latitude"]);b=float(r["longitude"])
                valid=(-90<=a<=90 and -180<=b<=180 and
                    date(1800,1,1)<=dy<=date(2026,10,17))
            except (ValueError,TypeError,KeyError):
                valid=False
            match=LINK_RE.fullmatch(r.get("photo_url_large",""))
            if not (valid and r.get("source_snapshot_date")=="2026-09-27" and
                    r.get("observer_id") and r.get("observation_id") and
                    r.get("photo_license") and r.get("species") and
                    match and match.group(1)==pid):
                stats["source_photo_metadata_invalid_or_unlicensed"]+=1
                continue
            loc=locus(r["observed_on"],a,b,r["inat_taxon_id"])
            c=old_pairs.get(loc[0]) if loc else None
            season=loc[1] if loc else None
            if c:
                opposite="winter" if season=="summer" else "summer"
                if c[opposite]>=1 and c[season]<2:
                    priority=0 if c[season]>=1 else 1
                elif c[season]>=1:
                    priority=2
                else:priority=3
            else:priority=3
            rec={k:r.get(k,"") for k in COLS}
            rec.update({
                "target_priority":str(priority),
                "prior_classifiable_summer_obs":str(c["summer"] if c else 0),
                "prior_classifiable_winter_obs":str(c["winter"] if c else 0),
                "screen_fold":"unmeasured_and_not_reserved_original_photo",
                "selection_hash":identity(pid)[:16]})
            if priority<3:
                key=(priority,loc[0],loc[1])
                push_smallest(top[key],identity(pid),rec,4)
                stats["strict_source_candidate_potential_repeat_photo_IDs"]+=1
            else:
                push_smallest(general[rec["inat_taxon_id"]],
                              identity(pid),rec,24)
    require(stats["source_photo_ID_candidates_scanned"]==SOURCE_N,
            "Actual entire original 3430502 photo source not scanned")
    return top,general,dict(stats)

def pick(top:dict,general:dict,reserved:set[str],n:int=RESERVE_N):
    """Balance target site/taxon/year groups before weaker/unpaired fallbacks."""
    chosen=[];ids=set();obs=set();users=set()
    stage=Counter()
    def admit(x):
        if x["photo_id"] in ids or x["observer_id"] in users or x["observation_id"] in obs:
            return False
        ids.add(x["photo_id"]);obs.add(x["observation_id"]);users.add(x["observer_id"])
        chosen.append(dict(x))
        stage[str(x["target_priority"])]+=1
        return True

    selected_groups=sorted(top, key=lambda key:(
        key[0],hashlib.sha256((SALT+"|GROUP|"+repr(key)).encode()).hexdigest()))
    sorted_candidates={key:[x[2] for x in sorted(top[key],reverse=True)]
                       for key in selected_groups}
    # One source candidate per existing target group before seeking another
    # from that group, so a species photographed by one user doesn't dominate.
    for turn in range(4):
        for key in selected_groups:
            if len(chosen)>=n:break
            pool=sorted_candidates[key]
            if turn<len(pool):admit(pool[turn])
        if len(chosen)>=n:break
    by_tax=sorted(general,key=lambda k:hashlib.sha256(
        (SALT+"|GENUS|"+k).encode()).hexdigest())
    rest={tax:[x[2] for x in sorted(general[tax],reverse=True)] for tax in by_tax}
    for turn in range(24):
        for tax in by_tax:
            if len(chosen)>=n:break
            pool=rest[tax]
            if turn<len(pool):admit(pool[turn])
        if len(chosen)>=n:break
    require(len(chosen)==n and len(ids)==len(obs)==len(users)==n,
            "No genuine 4096 independent-observer source photos available after prior exclusions")
    require(set(ids).isdisjoint(reserved),"Targeted photo IDs overlap fourth worker reservation")
    # Priority 0/1 is overrepresented by design; DO NOT infer population
    # summer white prevalence or general adaptive climate slopes from this.
    for ix,row in enumerate(chosen):
        row.update({
            "measurement_id":"MOIS-"+identity(row["photo_id"])[:24].upper(),
            "image_filename":"MOIS-"+identity(row["photo_id"])[:24].upper()+".jpg",
            "wave_sequence":str(ix),"shard":str(ix//SHARD_SIZE)})
    return chosen,dict(stage)

def run(full:Path,base:Path,fourth:Path,output:Path):
    require(not output.exists() or not any(output.iterdir()),
            "Do not rewrite an earlier source-locked original-photo study")
    old,known,baseline=known_support(base)
    reserved=check_reserved_fourth(fourth)
    require(not (old & reserved),"Frozen fourth photo source overlaps previously measured")
    candidates,general,scanreceipt=scan(full,old|reserved,known)
    chosen,stage=pick(candidates,general,reserved)
    output.mkdir(parents=True,exist_ok=True)
    result=output/PHOTO_SRC
    save_csv(result,chosen)
    shards={}
    for k in range(N_SHARDS):
        arr=chosen[k*SHARD_SIZE:(k+1)*SHARD_SIZE]
        require(len(arr)==SHARD_SIZE,"Incomplete target source original photo shard")
        path=output/f"mois_source_only_shard_{k:02d}_256.csv"
        save_csv(path,arr)
        shards[str(k)]=sha(path)
    receipt={
        "schema":"mois_3430502_targeted_opposite_season_after309439_and_reserved4096_v1",
        "status":"SOURCE_ONLY_4096_RESERVED_NO_NEW_PIXEL_HUES",
        "source_photo_ids_original":SOURCE_N,
        "prior_original_measured_photo_ids":FULL_BASE_PHOTOS,
        "prior_original_classified_four_hues":FULL_BASE_CLASSIFIED,
        "prior_original_dated_classified_hues":FULL_BASE_DATE_HUES,
        "prior_full_original_source_SHA256":FULL_SHA,
        "prior_measured_original_309439_source_SHA256":PHOTO_BASE_SHA,
        "prior_fourth_photo_id_reserved_SHA256":PRIOR_WAVE4_RESERVED_SHA,
        "first_wave_original_source_photo_ids":RESERVE_N,
        "first_wave_original_shards":N_SHARDS,
        "photos_per_original_shard":SHARD_SIZE,
        "first_wave_distinct_source_observers":len({r["observer_id"] for r in chosen}),
        "first_wave_distinct_original_observation_IDs":len({r["observation_id"] for r in chosen}),
        "first_wave_taxa":len({r["inat_taxon_id"] for r in chosen}),
        "first_wave_original_4096_csv_sha256":sha(result),
        "per_shard_sha256":shards,
        "first_finished_photo_pixel_shard_to_process":SELECT_ONLY_PIXEL_SHARD,
        "source_selection_existing_classifiable_photograph_site_year_support":baseline,
        "candidate_site_season_categories":stage,
        "metadata_source_scan":scanreceipt,
        "new_photo_pixels_actually_measured":False,
        "no_new_source_photo_colour_or_bloom_label_used":True,
        "no_photo_colour_or_fitness_or_bloom_label_used":True,
        "prior_outcome_classifiability_presence_used_to_target_repeated_sites":True,
        "current_recruited_source_is_not_representative_global_photo_prevalence":True,
        "only_exact_source_photo_id_metadata_used_for_new_photo_rank":True,
        "source_photo_pixels_stage_no_pixels_yet":True,
        "FCP_within_species_D_W_or_50km_estimands_imported":False,
    }
    (output/FREEZE_OUT).write_text(json.dumps(
        receipt,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("MOIS_TARGET_NEW_ORIGINAL_OPPOSITE_SEASON_PHOTO_ID_SOURCE_ONLY",
          json.dumps({k:receipt[k] for k in (
              "first_wave_original_source_photo_ids","first_wave_taxa",
              "candidate_site_season_categories",
              "source_selection_existing_classifiable_photograph_site_year_support")},
              sort_keys=True),flush=True)
    return receipt


if __name__=="__main__":
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("--full3430502-source",type=Path,required=True)
    a.add_argument("--signed309439-original",type=Path,required=True)
    a.add_argument("--fourth4096-reserved-source-zip",type=Path,required=True)
    a.add_argument("--output",type=Path,required=True)
    q=a.parse_args()
    run(q.full3430502_source,q.signed309439_original,q.fourth4096_reserved_source_zip,q.output)
