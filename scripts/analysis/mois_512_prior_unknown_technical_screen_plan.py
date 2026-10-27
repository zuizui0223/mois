#!/usr/bin/env python3
"""Reuse bulk128 medium500 photo results before technical ROI-failure diagnosis.

The 512 source-photo / 512 source-taxon blind REVIEW roster contains only
original-image strict9 unknowns. Photo-ID overlap with the earlier bulk128
131072-image 500px detector already provides 170 real detector outcomes; do
NOT fetch these again. Freeze the remaining 342 (or data-verified count) old
photo IDs for bounded 500px diagnostic-only screening. This is NOT a new
photo ID acquisition and never upgrades ROI unknown to flower pigment.
"""
from __future__ import annotations

import argparse,csv,hashlib,json,math
from collections import Counter
from pathlib import Path
import pandas as pd

SCHEMA="mois_512_prior_strictROI9_unknown_existing_medium500_reuse_and_bounded_screen_v1"
ROI_REVIEW_SHA="03a36332e0475af3f5393c67fa0bf7fd5269c146455a4b174e9b36b4eb1d46c7"
BULK_PARQUET_SHA="2d648c832929f8ba59301906d5ea8d100bdfec6449c6e2828f6b9a910ddd0a02"
EXPECTED_REVIEW=512
EXPECTED_ALREADY_MEDIUM=170
WORKER_SHARDS=8
GAP=("P1_source_rare_no_colour_1to9",
     "P2_source_has_photos_no_colour_10plus")
DESIGN=("historical_K1_FCP_source_heterogeneous",
        "bias_refill_2048","p2_first_4096","bulk128_131072")

def require(ok,message):
    if not bool(ok):raise ValueError("MOIS_512_REVIEW_REUSE: "+message)

def sha(path):
    d=hashlib.sha256()
    with Path(path).open("rb") as f:
        for part in iter(lambda:f.read(4<<20),b""):
            d.update(part)
    return d.hexdigest()

def fix_id(df):
    x=df.copy()
    for k in ("photo_id","observation_id","inat_taxon_id"):
        x[k]=x[k].fillna("").astype(str).str.strip()
    require(x.photo_id.str.fullmatch(r"\d+").all()
            and x.photo_id.is_unique
            and x.observation_id.ne("").all()
            and x.observation_id.is_unique
            and x.inat_taxon_id.ne("").all()
            and x.inat_taxon_id.is_unique,
            "Source-review photo, observation, taxon identities must be independent")
    return x

def make_plan(review,bulk,*,workers=WORKER_SHARDS,strict=True):
    """Exact source ID joins; status of the 170 premeasured images is reusable."""
    review=fix_id(review)
    bulk=bulk.copy()
    for k in ("photo_id","observation_id","inat_taxon_id"):
        bulk[k]=bulk[k].fillna("").astype(str)
    require(bulk.photo_id.is_unique and
            {"strict_original_roi9","screen_status","medium_actual_RGB_decoded",
             "floral_box_main_L_star"}<=set(bulk.columns),
            "Original bulk128 medium photo statuses and exact source ROI9 missing")
    expected={"photo_id","observation_id","inat_taxon_id",
        "source_design","latitude_zone","local_season",
        "postoutcome_source_gap","previous_original_photo_attempts",
        "source_photo_opportunities","reason",
        "botanical_organ_or_flower_stage_verified"}
    require(expected<=set(review.columns),
            "Canonical 512 source review lacks authentic selection/data fields")
    require(review.source_design.isin(DESIGN).all()
            and review.postoutcome_source_gap.isin(GAP).all()
            and review.botanical_organ_or_flower_stage_verified.astype(str).eq("False").all(),
            "Diagnostic original ROI unknown was already botanically verified or priority drifted")
    ids=set(review.photo_id)
    p=bulk.loc[bulk.photo_id.isin(ids)].copy()
    require(p.observation_id.is_unique and p.inat_taxon_id.is_unique
            and p.screen_status.isin(("candidate","negative","uncertain","error")).all(),
            "Cached medium500 original photos are not unique or detector status invalid")
    joined=review.merge(p[["photo_id","observation_id","inat_taxon_id",
                             "strict_original_roi9","screen_status",
                             "medium_actual_RGB_decoded","floral_box_main_L_star"]],
                       on=["photo_id","observation_id","inat_taxon_id"],
                       how="left",validate="one_to_one",indicator=True)
    overlap=joined["_merge"].eq("both")
    require(int(overlap.sum())==len(p),
            "Same photo ID mismatched original observation/taxon in old cached detector")
    require(not joined.loc[overlap,"strict_original_roi9"].astype(bool).any(),
            "A review-selected supposed strict9 failure was actually coloured in frozen bulk")
    if strict:
        require(len(review)==EXPECTED_REVIEW and len(p)==EXPECTED_ALREADY_MEDIUM,
                "Expected independent 512 photo review or 170 bulk cached screens changed")
    cached=joined.loc[overlap].copy()
    todo=joined.loc[~overlap].copy()
    require(len(todo)+len(cached)==len(review)
            and todo.photo_id.is_unique and
            not set(todo.photo_id)&set(cached.photo_id),
            "Prior photo screening unexpectedly duplicated or omitted")
    # Pin 8 work groups: geographically stratified deterministic hashed
    # balanced assignment, so every job has 42-43 old unknown photos.
    todo["screen_rank"]=todo.photo_id.map(lambda v:hashlib.sha256(
        ("MOIS_512_20261027_DIAG_V1|"+v).encode()).hexdigest())
    todo=todo.sort_values(["source_design","latitude_zone","local_season",
                           "postoutcome_source_gap","screen_rank"],kind="stable")
    todo["worker"]=0
    workload=[0]*workers
    indices=list(todo.index)
    for i,row_index in enumerate(indices):
        k=min(range(workers),key=lambda j:(workload[j],j))
        todo.loc[row_index,"worker"]=k
        workload[k]+=1
    require(max(workload)-min(workload)<=1,"Review source 8-worker counts not balanced")
    return cached,todo,{"source_review_photo_rows":len(review),
        "existing_500px_photo_screen_results_reused":len(cached),
        "existing_unknown_photos_needing_small_bounded_500px_recheck":len(todo),
        "eight_worker_photo_counts":workload,
        "unverified_photo_petals_not_a_true_pigment_label":True,
        "prior_9palette_colour_outcomes_all_remain_UNKNOWN":True,
        "new_photo_IDs_selected_for_world_census":0}

def run(review_csv,bulk_parquet,output,*,verify_sha=True):
    review_csv,bulk_parquet,output=map(Path,(review_csv,bulk_parquet,output))
    require(not output.exists(),"Refuse overwriting source-reconciled review")
    if verify_sha:
        require(sha(review_csv)==ROI_REVIEW_SHA
                and sha(bulk_parquet)==BULK_PARQUET_SHA,
                "Source signed review or previously measured medium500 pixel bytes changed")
    r=pd.read_csv(review_csv,dtype={"photo_id":str,
             "observation_id":str,"inat_taxon_id":str},low_memory=False)
    bulk=pd.read_parquet(bulk_parquet,columns=[
        "photo_id","observation_id","inat_taxon_id",
        "strict_original_roi9","screen_status","medium_actual_RGB_decoded",
        "floral_box_main_L_star"
    ])
    cached,todo,counts=make_plan(r,bulk)
    output.mkdir(parents=True)
    (output/"workers").mkdir()
    recols=["photo_id","observation_id","inat_taxon_id","source_design",
        "latitude_zone","local_season","postoutcome_source_gap",
        "previous_original_photo_attempts","source_photo_opportunities",
        "reason","screen_status","medium_actual_RGB_decoded",
        "floral_box_main_L_star"]
    cached[recols].to_csv(output/"already_screened_500px_170_no_refetch.csv",index=False)
    queuecols=["photo_id","observation_id","inat_taxon_id","source_design",
        "latitude_zone","local_season","postoutcome_source_gap",
        "previous_original_photo_attempts","source_photo_opportunities",
        "reason"]
    todo[queuecols].to_csv(output/"needs_bounded_medium500_screen_old_photoIDs.csv",index=False)
    worker_hashes={}
    for worker in range(WORKER_SHARDS):
        path=output/"workers"/f"shard_{worker:02d}.csv"
        todo.loc[todo.worker.eq(worker),queuecols].to_csv(path,index=False)
        worker_hashes[str(worker)]=sha(path)
    m={"schema":SCHEMA,**counts,
        "design":{"roi9_unknown_photo_IDs_only":True,"candidate_conclusion":
          "Technical detector-box diagnostics, NOT botanist-confirmed corolla"},
        "input_review_csv_SHA256":sha(review_csv),
        "existing_bulk128_RGB_pixel_parquet_SHA256":sha(bulk_parquet),
        "cached_medium_result_SHA256":sha(output/"already_screened_500px_170_no_refetch.csv"),
        "new_diagnostic_photo_queue_SHA256":sha(output/"needs_bounded_medium500_screen_old_photoIDs.csv"),
        "worker_source_SHA256":worker_hashes,
        "no_photo_bytes_downloaded_at_queue_creation":True}
    (output/"manifest.json").write_text(json.dumps(m,sort_keys=True,indent=2)+"\n")
    print("MOIS_512_EXISTING_ROI9_UNKNOWN_MEDIUM_SCREEN_REUSE",json.dumps(
        counts,sort_keys=True),flush=True)
    return m

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("review-csv","bulk-parquet","output"):
        p.add_argument("--"+key,required=True,type=Path)
    x=p.parse_args()
    run(x.review_csv,x.bulk_parquet,x.output)

if __name__=="__main__":
    main()
