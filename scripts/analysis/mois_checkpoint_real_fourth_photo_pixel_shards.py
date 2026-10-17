#!/usr/bin/env python3
"""Authenti­c finished original photo ROI shards of fourth4096: checkpoint ONLY.

This preserves all 309439 original photographic photo IDs/hues/dates, verifies
the independent source-locked 4096 photograph sample and each completed actual
256-photo model worker terminal, then admits ONLY finished image RGB results.
The unmeasured remaining source IDs DO NOT count as classified photo data.

A checkpoint is NOT a second independent cohort: the final fourth 4096 study
will *replace* this partial-source science denominator after all 16 finish.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import io
import json
import zipfile
from pathlib import Path

import pandas as pd
import numpy as np
from scripts.analysis.mois_reassemble_original_new4096_measurements import (
    HUES, BASE_COLS, METHOD, MEASURED_SCHEMA
)
from scripts.analysis.mois_seasonal_267k_climate_normals import (
    labelled_calendar,require
)

PRIOR_PHOTOS=309439
PRIOR_CLASSIFIED=149662
PRIOR_DATED_CLASSIFIED=149653
PRIOR_SOURCE="mois_original309439_fullFCP256_actual_photo_hues_and_source_dates.csv.gz"
PRIOR_SHA="546266e2bbaa340a2523864c28c5a877842eb4b0904a7e0dff90870c6519a47e"
PREPIXEL_SELECTION="mois_full3430502_next4096_original_photo_IDs_AFTER309439_WAVE4.csv"
SELECTION_SHA="65fc139909f94448e042b6539c6749e70d1ab21bb8b13f50383bb7417180337a"
FREEZE="mois_first4096_original_source_only_freeze.json"
RESULT_PHOTOS="mois_fourth_partial_finished_real_RGB_photo_hues_plus309439.csv.gz"
RESULT_RECEIPT="mois_fourth_partial_finished_pixel_source_checkpoint_receipt.json"
RESULT_SEASON="mois_fourth_partial_finished_pixel_source_hemisphere_season.csv"

def digest(raw:bytes)->str:
    return hashlib.sha256(raw).hexdigest()

def read_original(original_zip:Path):
    with zipfile.ZipFile(original_zip) as z:
        raw=z.read(PRIOR_SOURCE)
    require(digest(raw)==PRIOR_SHA,
            "Original 309439 actually measured source photo GZIP altered")
    old=pd.read_csv(io.BytesIO(raw),compression="gzip",
                    dtype=str,keep_default_na=False)
    require(len(old)==PRIOR_PHOTOS and old.photo_id.is_unique and
            tuple(old.columns)==(*BASE_COLS,"qualified_date",
                                 "original_source_date_provenance") and
            int(old.morph.isin(HUES).sum())==PRIOR_CLASSIFIED and
            int((old.morph.isin(HUES)&old.qualified_date.str.len().eq(10)).sum())==
                 PRIOR_DATED_CLASSIFIED,
            "Original real309439 source photo/date/hues differ")
    return old

def read_selection(plan_zip:Path):
    with zipfile.ZipFile(plan_zip) as z:
        snap=json.loads(z.read(FREEZE))
        data=z.read(PREPIXEL_SELECTION)
        shards={int(k):z.read(f"mois_source_only_shard_{int(k):02d}_256.csv")
                for k in range(16)}
    require(snap["schema"]=="mois_original_full3430502_next4096_after309439_wave4_preoutcome_v1" and
            snap["source_measured_original_309439_SHA256"]==PRIOR_SHA and
            snap["first_wave_original_source_photo_ids"]==4096 and
            snap["first_wave_original_4096_csv_sha256"]==SELECTION_SHA and
            snap["no_photo_colour_or_fitness_or_bloom_label_used"] and
            snap["new_photo_pixels_actually_measured"] is False and
            digest(data)==SELECTION_SHA and
            len(shards)==16 and
            all(digest(buf)==snap["per_shard_sha256"][str(k)]
                for k,buf in shards.items()),
            "Pre-pixel source photograph selection+old-source exclusion drift")
    source=list(csv.DictReader(io.StringIO(data.decode())))
    require(len(source)==4096 and len({r["photo_id"] for r in source})==4096
            and len({r["observer_id"] for r in source})==4096 and
            len({r["observation_id"] for r in source})==4096,
            "Fourth-wave original source IDs/observer/source observation no longer unique")
    by_worker={}
    for k,buf in shards.items():
        group=list(csv.DictReader(io.StringIO(buf.decode())))
        require(len(group)==256 and all(int(r["shard"])==k for r in group),
                "Fourth-wave 256-photo complete original source shard missing")
        by_worker[k]={r["measurement_id"]:r for r in group}
    return snap,by_worker

def complete_worker(source:dict,root:Path,k:int,freeze:dict):
    csvfiles=list(root.rglob(f"mois_original_256_fresh_ROI_shard_{k:02d}.csv"))
    receipts=list(root.rglob(f"mois_original_256_ROI_shard_{k:02d}_receipt.json"))
    require(len(csvfiles)==len(receipts)==1,
            "A completed RGB worker has no authentic original photo terminal/receipt")
    receipt=json.loads(receipts[0].read_text())
    require(receipt["schema"]==MEASURED_SCHEMA and receipt["shard"]==k and
            receipt["original_source_ids"]==receipt["terminal_rows"]==256 and
            receipt["source_sha256"]==freeze["per_shard_sha256"][str(k)] and
            receipt["source_frozen_4096_selection_sha256"]==SELECTION_SHA and
            receipt["source_model_attestation"]==METHOD and
            receipt["fcp_within_species_polymorphism_estimand_computed"] is False and
            receipt["pixels_ephemeral_no_original_image_bytes_persisted"] is True and
            receipt["original_new_photo_hue_terminal_SHA256"]==
                digest(csvfiles[0].read_bytes()),
            "Original unchanged model/photo pixels/shard SHA or source signature failed")
    terminal=list(csv.DictReader(csvfiles[0].open(newline="",encoding="utf-8-sig")))
    expected=source[k]
    seen=set();rows=[]
    for row in terminal:
        mid=row["measurement_id"]
        require(mid in expected and mid not in seen,
                "Unselected or repeated image original photo measurement")
        seen.add(mid)
        old=expected[mid]
        for field in ("photo_id","inat_taxon_id","species",
                      "observation_id","observer_id"):
            require(row[field]==old[field],
                    "Actual original measured photo ID, taxon, observation/observer changed")
        yes=row["morph"] in HUES
        require(row["morph"] in (*HUES,"mixed_uncertain") and
                (row["measurement_status"]=="classified_four_state_morph")==yes
                and row["classified_four_hue"]==("true" if yes else "false"),
                "Original visible hue-v4 ROI colour/unknown outcome tampered")
        new={x:"" for x in BASE_COLS}
        for field in BASE_COLS:
            if field in row:new[field]=row[field]
        for field in ("latitude","longitude","photo_license","observed_on"):
            new[field]=old[field]
        new["source_original_measured"]="mois_wave4_new_original_ROI_v4"
        new["old_mois_or_new_fcp"]="mois_new_original_ROI"
        new["qualified_date"]=old["observed_on"]
        new["original_source_date_provenance"]=(
            "fcp_full3430502_original_photo_observed_on_with_MOIS_wave4_unchanged_ROI_v4")
        rows.append((mid,new))
    require(len(rows)==256 and len(seen)==256 and
            sum(r[1]["morph"] in HUES for r in rows)==receipt["four_hue_classified"] and
            sum(r[1]["morph"] not in HUES for r in rows)==receipt["unclassified"],
            "Whole 256-photo shard not source-terminal complete")
    by_mid={k:v for k,v in rows}
    return [by_mid[mid] for mid in expected],receipt

def gather_workers(plan:dict,root:Path,freeze:dict):
    allrec=list(root.rglob("mois_original_256_ROI_shard_*_receipt.json"))
    require(2<=len(allrec)<=16,"Need at least 2 actually finished photo RGB shards, not source IDs")
    indices=[]
    for f in allrec:
        import re
        m=re.fullmatch(r"mois_original_256_ROI_shard_(\d\d)_receipt.json",f.name)
        require(m is not None,"Unrecognized finished 256-photo source ROI receipt")
        k=int(m[1])
        require(0<=k<16 and k not in indices,"Duplicate original source RGB worker")
        indices.append(k)
    all_photos=[];receipts=[]
    for k in sorted(indices):
        fresh,proof=complete_worker(plan,root,k,freeze)
        all_photos.extend(fresh);receipts.append({
            "shard":k,"original_real_photo_attempts":len(fresh),
            "classified_four_hue_actual_ROIs":proof["four_hue_classified"],
            "terminal_original_pixel_model_SHA":proof["original_new_photo_hue_terminal_SHA256"]})
    return all_photos,receipts

def run(original_zip:Path,plan_zip:Path,workers:Path,output:Path):
    require(not output.exists() or not any(output.iterdir()),
            "Never overwrite an actually measured checkpoint source")
    old=read_original(original_zip)
    freeze,source=read_selection(plan_zip)
    new,proof=gather_workers(source,workers,freeze)
    additions=pd.DataFrame(new,columns=old.columns).fillna("")
    require(len(new)==256*len(proof) and additions.photo_id.is_unique and
            not additions.photo_id.isin(old.photo_id).any(),
            "Source-only or already measured photo IDs cannot be new pixel data")
    require(additions.qualified_date.str.fullmatch(r"\d{4}-\d{2}-\d{2}").all(),
            "Original source observed date was not preserved in actual new pixels")
    result=pd.concat([old,additions],ignore_index=True)
    require(result.photo_id.is_unique and len(result)==PRIOR_PHOTOS+len(additions)
            and result.iloc[:PRIOR_PHOTOS].reset_index(drop=True).equals(old),
            "Existing 309439 genuine photo/classifier/source day records overwritten")
    nclass=int(additions.morph.isin(HUES).sum())
    nhues=int(result.morph.isin(HUES).sum())
    dated=int((result.morph.isin(HUES)&result.qualified_date.str.len().eq(10)).sum())
    require(nclass>0 and nhues==PRIOR_CLASSIFIED+nclass and
            dated==PRIOR_DATED_CLASSIFIED+nclass,
            "Original new photo RGB images were not actually classified+dated")
    output.mkdir(parents=True,exist_ok=True)
    actual=output/RESULT_PHOTOS
    result.to_csv(actual,index=False,lineterminator="\n",
                  compression={"method":"gzip","mtime":0})
    panel=result[["photo_id","observation_id","inat_taxon_id","morph",
         "latitude","longitude","qualified_date"]].rename(
         columns={"qualified_date":"observed_on"}).copy()
    for name in ("latitude","longitude"):
        panel[name]=pd.to_numeric(panel[name],errors="coerce")
    q=labelled_calendar(panel)
    seas=[]
    for zone in ("north_extratropical","south_extratropical"):
        for month in ("winter","spring","summer","autumn"):
            r=q[q.zone.eq(zone)&q.season.eq(month)&q.classified&q.date.notna()]
            seas.append({"zone":zone,"season":month,"classifiable_original_photo_ROIs":len(r),
              "white_original_photo_ROIs":int(r.white.sum()),
              "photographed_white_share":round(float(r.white.mean()),8) if len(r) else None})
    path=output/RESULT_SEASON
    pd.DataFrame(seas).to_csv(path,index=False,lineterminator="\n")
    receipt={
       "schema":"mois_fourth_original_actual_finished_partial_photo_pixel_checkpoint_v1",
       "source_original_real_309439_photo_SHA256":PRIOR_SHA,
       "source_original_prepixel_4096_photo_ID_SHA256":SELECTION_SHA,
       "completed_authentic_ROI_shards":len(proof),
       "actual_finished_shard_source_receipts":proof,
       "original_4096_photo_IDs_planned_not_measured":4096,
       "original_actual_photo_IDs_with_real_pixel_measurement_added":len(additions),
       "original_four_hue_ROI_source_photo_labels_added":nclass,
       "original_technical_colour_unknown_pixels_added":len(additions)-nclass,
       "original_photo_hue_plus_observation_day_added":nclass,
       "measured_original_photo_IDs_before":PRIOR_PHOTOS,
       "measured_original_photo_IDs_after_checkpoint":len(result),
       "original_classifiable_four_hue_ROIs_after_checkpoint":nhues,
       "original_date_plus_four_hue_ROIs_after_checkpoint":dated,
       "unmeasured_fourth_cohort_original_photo_IDs_remaining":4096-len(additions),
       "new_original_photo_source_full_GZIP_SHA256":digest(actual.read_bytes()),
       "new_original_photo_season_source_CSV_SHA256":digest(path.read_bytes()),
       "previous_photo_ID_hue_taxon_and_source_day_UNCHANGED":True,
       "NOT_AN_INDEPENDENT_REPLICATION_AND_WILL_BE_SUPERSEDED_BY_FULL_FOURTH4096":True,
       "original_ROI_label_DOES_NOT_mean_photo_specific_open_petals":True,
       "no_global_unbiased_moisture_genus_maxT_fitted_from_partial_source":True,
       "FCP_intraspecific_D_W_genetic_50km_not_imported":True
    }
    (output/RESULT_RECEIPT).write_text(json.dumps(
        receipt,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("MOIS_AUTHENTIC_NEW_FOURTH_PARTIAL_ROI_PIXELS",
          json.dumps({k:receipt[k] for k in (
             "completed_authentic_ROI_shards",
             "original_actual_photo_IDs_with_real_pixel_measurement_added",
             "original_four_hue_ROI_source_photo_labels_added",
             "measured_original_photo_IDs_after_checkpoint",
             "original_date_plus_four_hue_ROIs_after_checkpoint")},sort_keys=True),flush=True)
    return receipt

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original309439-photo-archive",type=Path,required=True)
    p.add_argument("--blind-original4096-source-archive",type=Path,required=True)
    p.add_argument("--complete-photo-RGB-shard-artifacts",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    run(args.original309439_photo_archive,args.blind_original4096_source_archive,
        args.complete_photo_RGB_shard_artifacts,args.output)
