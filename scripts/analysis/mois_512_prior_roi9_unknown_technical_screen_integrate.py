#!/usr/bin/env python3
"""Integrate 170 frozen medium results + 342 real bounded detector photos.

EVERY one of 512 distinct old original source photo/taxon/observation IDs
is a previous STRICT ROI9 UNKNOWN. This diagnoses detector coverage only;
even a high-confidence flower box is NOT validated open-petal colour or
proof a different full-resolution segmentation would succeed.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib,json
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.acquisition.mois_existing_roi_unknown_342_medium_flower_diagnostic import (
 SCHEMA as SHARD_SCHEMA,MODEL_SHA
)
from scripts.analysis.mois_512_prior_unknown_technical_screen_plan import (
 SCHEMA as SOURCE_SCHEMA,ROI_REVIEW_SHA,BULK_PARQUET_SHA
)

SCHEMA="mois_genuine_old512_originalROI9_unknown_two_medium_screen_methods_full_id_diagnostic_v1"
N_TOTAL=512
N_CACHE=170
N_NEW=342
N_WORKERS=8

def require(ok,msg):
    if not bool(ok):raise ValueError("MOIS_512_OLD_UNKNOWN_BOX_AUDIT: "+msg)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(4<<20),b""):h.update(block)
    return h.hexdigest()

def cached_as_rows(df):
    d=df.copy()
    d["technical_diagnostic_class"]=d.screen_status.map({
        "candidate":"existing_cached_flower_box_ROI9_failed",
        "negative":"existing_cached_detector_negative_not_botanical_absence",
        "error":"existing_cached_error",
        "uncertain":"existing_cached_uncertain"
    })
    require(d.technical_diagnostic_class.notna().all(),
            "The previously cached medium500 screen has unknown detection status")
    d["500px_pixel_decoded"]=d.medium_actual_RGB_decoded.astype(str).eq("True")
    d["sample_class"]="previously_500px_screened_NO_EXTRA_RGB_DOWNLOAD"
    d["largest_box_area_fraction"]=np.nan
    d["flower_box_count_002"]=np.nan
    d["diagnostic_source_model_SHA256"]="BULK128_PRIOR_SOURCE_MODEL_CONTRACT"
    return d

def combine(source,worker_root,out,*,expected_total=N_TOTAL,expected_cache=N_CACHE,
            expected_new=N_NEW,expected_workers=N_WORKERS):
    source,worker_root,out=map(Path,(source,worker_root,out))
    require(not out.exists(),"Refuse to overwrite existing review results")
    meta=json.loads((source/"manifest.json").read_text())
    require(meta.get("schema")==SOURCE_SCHEMA
            and meta.get("source_review_photo_rows")==expected_total
            and meta.get("existing_500px_photo_screen_results_reused")==expected_cache
            and meta.get("existing_unknown_photos_needing_small_bounded_500px_recheck")==expected_new
            and meta.get("new_photo_IDs_selected_for_world_census")==0
            and meta.get("input_review_csv_SHA256")==ROI_REVIEW_SHA
            and meta.get("existing_bulk128_RGB_pixel_parquet_SHA256")==BULK_PARQUET_SHA,
            "Frozen old-photo ROI unknown technical source plan drifted")
    p=source/"already_screened_500px_170_no_refetch.csv"
    q=source/"needs_bounded_medium500_screen_old_photoIDs.csv"
    require(sha(p)==meta["cached_medium_result_SHA256"]
            and sha(q)==meta["new_diagnostic_photo_queue_SHA256"],
            "Previously processed 170 photo status or 342 photo source ID SHA changed")
    cache=pd.read_csv(p,dtype={
        "photo_id":str,"observation_id":str,"inat_taxon_id":str},low_memory=False)
    qdf=pd.read_csv(q,dtype={
        "photo_id":str,"observation_id":str,"inat_taxon_id":str},low_memory=False)
    require(len(cache)==expected_cache and len(qdf)==expected_new,
            "Previously screened vs originally unknown technical review counts changed")
    parts=[];receipts=[]
    for i in range(expected_workers):
        folder=worker_root/f"mois-oldROI9-medium-diagnostic-{i}"
        f=folder/f"shard_{i:02d}_original_medium_diagnostic.csv"
        receipt=folder/f"shard_{i:02d}_source_receipt.json"
        require(f.is_file() and receipt.is_file(),
                "An actual old-photo medium diagnostic worker was not successfully downloaded")
        r=json.loads(receipt.read_text())
        require(r.get("schema")==SHARD_SCHEMA
                and r.get("worker")==i
                and r.get("model_SHA256")==MODEL_SHA
                and r.get("original_source_manifest_SHA256")==sha(source/"manifest.json")
                and r.get("individual_original_source_worker_SHA256")==
                    meta["worker_source_SHA256"][str(i)]
                and r.get("measurement_CSV_SHA256")==sha(f)
                and r.get("no_raw_image_bytes_archived")
                and r.get("prior_ROI9_unknowns_never_imputed_as_white_or_pigmented"),
                "Source SHA/YOLO detector photo content/source identity changed")
        d=pd.read_csv(f,dtype={
            "photo_id":str,"observation_id":str,"inat_taxon_id":str},low_memory=False)
        require(len(d)==r["original_old_photo_ids_attempted_for_first_medium_diagnostic"]
                and len(d)==meta["eight_worker_photo_counts"][i]
                and d.photo_id.is_unique
                and d.screen_status.isin(("candidate","negative","error")).all(),
                "Old previous image detector worker changed or lost photo IDs")
        d["sample_class"]="bounded_old_original_photo_FIRST_500px_diagnostic"
        d["diagnostic_source_model_SHA256"]=MODEL_SHA
        parts.append(d);receipts.append({"worker":i,"sha":sha(receipt),"decoded":r["decoded_medium500"]})
    new=pd.concat(parts,ignore_index=True)
    for df in (cache,new,qdf):
        for key in ("photo_id","observation_id","inat_taxon_id"):
            df[key]=df[key].fillna("").astype(str)
    proof=new.merge(qdf[["photo_id","observation_id","inat_taxon_id"]],
                    on=["photo_id","observation_id","inat_taxon_id"],
                    how="inner",validate="one_to_one")
    require(len(proof)==expected_new and set(new.photo_id)==set(qdf.photo_id),
            "Old technical unknown photo original observation/taxon mismatch")
    cache=cached_as_rows(cache)
    require(not set(cache.photo_id)&set(new.photo_id),
            "An already screened 500px photo was redownloaded by new diagnostic worker")
    ccols=["photo_id","observation_id","inat_taxon_id","source_design",
            "latitude_zone","local_season","postoutcome_source_gap",
            "previous_original_photo_attempts","source_photo_opportunities",
            "reason","screen_status","technical_diagnostic_class",
            "500px_pixel_decoded","largest_box_area_fraction","flower_box_count_002",
            "sample_class","diagnostic_source_model_SHA256"]
    d=pd.concat((cache.reindex(columns=ccols),new.reindex(columns=ccols)),
                ignore_index=True)
    require(len(d)==expected_total and d.photo_id.is_unique and
            d.observation_id.is_unique and d.inat_taxon_id.is_unique,
            "512 original-photo detection measurements unexpectedly replicated")
    require(d.screen_status.isin(("candidate","negative","uncertain","error")).all()
            and d["500px_pixel_decoded"].astype(str).isin(("True","False")).all(),
            "Technical detector state or authentic source image decode missing")
    d["detector_candidate_002"]=d.screen_status.eq("candidate")
    d["detector_negative"]=d.screen_status.eq("negative")
    d["old_ROI9_missing"]=True
    d["strict_original9_pigment_fraction"]=np.nan
    summary=(d.groupby(["sample_class","source_design","latitude_zone",
                        "local_season","postoutcome_source_gap"],
                observed=True,dropna=False)
             .agg(old_photo_ids=("photo_id","size"),
                  original_nominal_taxa=("inat_taxon_id","nunique"),
                  real_medium_rgb_decoded=("500px_pixel_decoded","sum"),
                  flower_box_candidate=("detector_candidate_002","sum"),
                  no_flower_box_detected=("detector_negative","sum"),
                  original_roi9_still_unknown=("old_ROI9_missing","sum"))
             .reset_index())
    def stable_rank(pid):
        return hashlib.sha256(
            ("MOIS_512_POST_BOX_BLIND_BOTANICAL_REVIEW|"+pid).encode()).hexdigest()
    # Old images WITH detected flower boxes should be evaluated for possible
    # improved segmentation and actual petal mask; old detector negatives
    # independently inspected for false-negative screening.
    manual=d.copy()
    manual["rank"]=manual.photo_id.map(stable_rank)
    manual["review_priority"]=manual.technical_diagnostic_class.map(
        lambda x:"inspect_strictROI9_failure_with_detected_large_flower_box"
          if x=="detectable_large_flower_box_ROI9_failed" else
          ("inspect_detector_negative_no_botanical_absence_claim"
           if "negative" in x else
           "inspect_flower_stage_and_original_photo_ROI9_gate"))
    manual=manual.sort_values(
        ["review_priority","postoutcome_source_gap","source_design",
         "latitude_zone","rank"],kind="stable")
    # These are source-photo cases awaiting botanical inspection; do not
    # relabel them, and do not call detector-positive a true petal mask.
    manual=manual[[
        "photo_id","observation_id","inat_taxon_id","sample_class","source_design",
        "latitude_zone","local_season","postoutcome_source_gap",
        "screen_status","technical_diagnostic_class","largest_box_area_fraction",
        "review_priority"
    ]]
    out.mkdir(parents=True)
    table=out/"all512_oldROI9_unknown_photoID_detector_status.csv"
    d.to_csv(table,index=False)
    summary_file=out/"review512_originalROI9_unknown_detector_strata.csv"
    summary.to_csv(summary_file,index=False)
    review_file=out/"review512_detector_to_botanist_next_action.csv"
    manual.to_csv(review_file,index=False)
    counts=Counter(d.screen_status)
    cases=Counter(d.technical_diagnostic_class)
    byarm={k:{
        "photos":len(g),"candidate":int(g.detector_candidate_002.sum()),
        "negative":int(g.detector_negative.sum()),
        "unknown_ROI9":len(g)
    } for k,g in d.groupby("source_design",observed=True)}
    result={
        "schema":SCHEMA,
        "old_original_photo_IDs_with_all_strictROI9_unknown":len(d),
        "distinct_source_taxa_reviewed_technically_not_botanically":int(d.inat_taxon_id.nunique()),
        "old_photo_ID_500px_detector_outputs_reused_without_fetch":expected_cache,
        "old_photo_IDs_newly_diagnosed_500px_but_already_original_fullROI_failed":expected_new,
        "real_medium500_RGB_decoded_combined":int(d["500px_pixel_decoded"].sum()),
        "detector_statuses_512":dict(counts),
        "technical_oldROI_failure_box_classes":dict(cases),
        "by_source_design":byarm,
        "genuine_new_strict9_ROI_labels_generated":0,
        "original_strict9_photo_outcomes_changed":0,
        "new_world_photo_IDs_or_plant_taxa_acquired":0,
        "old_photo_RGB_persisted_or_redistributed":False,
        "flower_detector_candidate_is_not_true_corroborated_botanical_petal":True,
        "ROI9_failure_causes_not_yet_manually_adjudicated":True,
        "source_review_manifest_SHA256":sha(source/"manifest.json"),
        "reused170_photo_status_SHA256":sha(p),
        "source342_unmeasured_medium_photo_queue_SHA256":sha(q),
        "workers_signed_SHA_receipts":receipts,
        "joint_diagnostic_source_photo_SHA256":sha(table),
        "next_action_original_photo_review_SHA256":sha(review_file),
        "region_season_technical_screen_summary_SHA256":sha(summary_file)
    }
    path=out/"mois_512_oldROI9_unknown_detector_technical_diagnostic_receipt.json"
    path.write_text(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+"\n")
    print("MOIS_512_OLD_ROI9_UNKNOWN_REAL_IMAGE_DETECTOR_DIAGNOSIS",json.dumps({
        "512_photo_count":len(d),
        "medium_RGB_decoded":result["real_medium500_RGB_decoded_combined"],
        "screen_candidates":counts["candidate"],"screen_negatives":counts["negative"],
        "screen_errors":counts["error"]+counts["uncertain"],
        "previous_screen_original_photo_IDs_reused":expected_cache,
        "old_photo_IDs_first_medium_diagnostic":expected_new
    },sort_keys=True),flush=True)
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("source","workers","output"):
        p.add_argument("--"+k,required=True,type=Path)
    x=p.parse_args()
    combine(x.source,x.workers,x.output)

if __name__=="__main__":
    main()
