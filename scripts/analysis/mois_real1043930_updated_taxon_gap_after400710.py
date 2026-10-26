#!/usr/bin/env python3
"""Recompute MOIS 42,111-taxon data-collection priorities after 10,169 authentic
previously lost original nine-fraction ROI vectors were SHA-restored.

Use the OLD source photo/taxon opportunity and already-attempted photo register
without rescanning the 93.2m image-ID index. Verify exact 1,043,930 source
photo IDs, previously frozen 390,541 ROI9 counts per taxon, and true new
400,710 ROI9 from old signed wave2/cached/tail original archived photo data.
Do NOT treat categories or old unknowns as new measured pixels.
"""
from __future__ import annotations

import argparse,hashlib,json
from pathlib import Path
import pandas as pd
import numpy as np

from scripts.analysis.mois_real1043930_taxon_gap_priority import (
    priority,profile,sha as file_sha
)

SCHEMA="mois_original1043930_recovered400710_source_taxon_priorities_v2"
NEW_SHA="df0f562b58bdc7ba06e901530eb809f513a763e2b0524d3b195e99ee56e930ef"
PREVIOUS_STRICT=390541
STRICT=400710
N=1043930
TAXA=42111

def require(cond,msg):
    if not cond:raise ValueError("MOIS_UPDATED_TAXON_PRIORITIES: "+msg)

def audit(trait,receipt,old_register,old_report,output):
    trait,receipt,old_register,old_report,output=map(
        Path,(trait,receipt,old_register,old_report,output))
    require(not output.exists(),"New original-photo priority result must be immutable")
    r=json.loads(receipt.read_text(encoding="utf-8"))
    old=json.loads(old_report.read_text(encoding="utf-8"))
    require(r.get("schema")=="mois_real1043930_oldwaves_10169_original_roi9_restoration_v1"
            and r.get("new_source_photo9_lake_SHA256")==NEW_SHA==file_sha(trait)
            and r.get("full_original_photo_attempts")==N
            and r.get("new_original_strict9_photo_count")==STRICT
            and r.get("old_9component_strict_photo_count")==PREVIOUS_STRICT
            and r.get("source_photo9_restored_original_wave2_cached8k_tail16k")==10169
            and r.get("new_original_image_downloads")==0,
            "New 400710-photo source not restored solely from original stored numeric ROI")
    require(old.get("schema")=="mois_original1043930_taxon_photo_gap_and_source_priority_v1"
            and old.get("source_frame_taxa")==TAXA
            and old.get("strict9_photo_source")==PREVIOUS_STRICT
            and old.get("original_source_photo_ID_candidates")==93209351,
            "Previous 42k frame composition not verified")
    before=pd.read_csv(old_register,compression="gzip",
                       dtype={"inat_taxon_id":str},low_memory=False)
    require(len(before)==TAXA and before.inat_taxon_id.is_unique
            and int(before.n_strict_roi9.sum())==PREVIOUS_STRICT
            and int(before.n_image_photo_attempts.sum())==
                N-old["source_attempted_photos_but_outside_42111_taxon_frame"],
            "Previous original frame and exact original 1m photo IDs no longer match")
    photo=pd.read_parquet(trait,columns=[
        "photo_id","inat_taxon_id","source_arm_original",
        "original_four_hue","continuous9_strict","continuous9_exploratory_ambiguous"])
    photo.photo_id=photo.photo_id.astype(str)
    photo.inat_taxon_id=photo.inat_taxon_id.fillna("").astype(str)
    require(len(photo)==N and photo.photo_id.is_unique
            and int(photo.continuous9_strict.sum())==STRICT
            and int(photo.continuous9_exploratory_ambiguous.sum())==36224
            and int(photo.original_four_hue.isin(
                ("white","yellow_orange","red_pink","blue_purple")).sum())==459235,
            "New raw flower ROI9 photo source grain, hue or precision changed")
    oldcount=before.set_index("inat_taxon_id")["n_strict_roi9"]
    newcounts=photo.loc[photo.continuous9_strict].groupby(
         "inat_taxon_id",observed=True).size()
    current=before.set_index("inat_taxon_id").copy()
    current["previous_original9"]=current.n_strict_roi9
    current["n_strict_roi9"]=newcounts.reindex(
          current.index,fill_value=0).astype(int)
    current["roi9_from_frozen_oldwave_newly_restored"]=(
        current.n_strict_roi9-current.previous_original9)
    require(current.roi9_from_frozen_oldwave_newly_restored.ge(0).all(),
            "One 390541 old original nine vector was lost after restore")
    # Structural source photo attempt counts MUST remain exactly unchanged.
    imagecounts=photo.groupby("inat_taxon_id",observed=True).size()
    require(np.array_equal(
          current.n_image_photo_attempts.to_numpy(int),
          imagecounts.reindex(current.index,fill_value=0).to_numpy(int)),
          "Same exact photo-ID original attempt census changed during 10169 ROI9 recovery")
    outside=photo.loc[~photo.inat_taxon_id.isin(current.index)]
    require(len(outside)==old[
        "source_attempted_photos_but_outside_42111_taxon_frame"],
        "Older non-frame source photo ID population changed")
    added_outside=int(outside.continuous9_strict.sum())
    added_inside=int(current.roi9_from_frozen_oldwave_newly_restored.sum())
    require(added_inside+added_outside==STRICT-PREVIOUS_STRICT,
            "Source photo original new 10169 strict ROI9 must belong to exactly one source taxon")
    current["priority"]=[
      priority(n,k) for n,k in zip(
            current.n_unique_source_photo_ids,current.n_strict_roi9)]
    current["number_of_strict9_photo_ROIs_to_reach_5"]=(
        5-current.n_strict_roi9).clip(lower=0)
    stats=profile(current.reset_index())
    require(sum(v["source_frame_taxa"] for v in stats)==TAXA
            and sum(v["previous_strict_roi9_photo_count"] for v in stats)==STRICT-added_outside,
            "No source frame taxon may be silently removed by recovery")
    missing=current.n_strict_roi9.eq(0)
    opportunities=current.n_unique_source_photo_ids.gt(0)
    rare=current.n_unique_source_photo_ids.between(1,99)
    v1_missing=int(old["source_taxa_without_any_valid9"])
    result={
      "schema":SCHEMA,
      "full_source_93m_photo_IDs":93209351,
      "source_42111_taxa":TAXA,
      "real_original_photo_attempts":N,
      "prior_original_strict9_photo_count":PREVIOUS_STRICT,
      "strict9_after_no_image_downloads":STRICT,
      "strict9_added_from_old_worker_source":STRICT-PREVIOUS_STRICT,
      "strict9_added_inside_42111_original_taxa":added_inside,
      "strict9_added_OUTSIDE_42111_source_taxon_frame":added_outside,
      "source_taxa_with_any_strict9_NOW":int((~missing).sum()),
      "source_taxa_without_strict9_NOW":int(missing.sum()),
      "source_taxa_newly_covered_since_prior390541":v1_missing-int(missing.sum()),
      "source_taxa_with_original_photo_opportunities_but_no9_NOW":int(
          (missing&opportunities).sum()),
      "eligible_source_photo_opportunities_in_taxa_still_without9":int(
          current.loc[missing&opportunities,"n_unique_source_photo_ids"].sum()),
      "source_1to99_photo_taxa_without9_NOW":int((missing&rare).sum()),
      "source_1to99_photo_taxa_with9_NOW":int((~missing&rare).sum()),
      "valid9_1to4_photo_taxa_NOW":int(current.n_strict_roi9.between(1,4).sum()),
      "valid9_5to9_photo_taxa_NOW":int(current.n_strict_roi9.between(5,9).sum()),
      "valid9_10plus_photo_taxa_NOW":int(current.n_strict_roi9.ge(10).sum()),
      "taxa_with_photo_attempts_ge5_but_zero9":int((
         current.n_image_photo_attempts.ge(5)&missing).sum()),
      "taxa_with_photo_attempts_eq1_and_zero9":int((
         current.n_image_photo_attempts.eq(1)&missing).sum()),
      "priority_distribution_updated":stats,
      "prior_v1_priority_report_sha256":file_sha(old_report),
      "new_400710_original_photo9_lake_sha256":file_sha(trait),
      "no_new_photo_image_downloads":True,
      "no_legacy_four_hue_mutations":True,
      "all_original_93m_indexed_photos_remain_eventual_targets":True,
      "individual_new_photo_IDs_already_selected":False,
      "next_step":"Original93m source observation-deduplicated photo-ID priority queue for taxa still missing strict9",
    }
    output.mkdir(parents=True)
    register=output/"mois_42111_taxon_priority_AFTER_400710_original9.csv.gz"
    current.reset_index().sort_values(
        ["priority","n_unique_source_photo_ids","inat_taxon_id"]).to_csv(
        register,index=False,compression="gzip")
    pd.DataFrame(stats).to_csv(output/"mois_taxon_priority_AFTER_400710_summary.csv",index=False)
    (output/"mois_real1043930_updated_400710_taxon_gap_decision.json").write_text(
        json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+"\n")
    print("MOIS_REAL1043930_400710_SOURCE_TAXON_PRIORITY_AFTER_NO_IMAGE_RECOVERY",
          json.dumps({k:v for k,v in result.items() if k!="priority_distribution_updated"},
                     sort_keys=True),flush=True)
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("photo","receipt","old-register","old-report","output"):
        p.add_argument("--"+name,required=True,type=Path)
    q=p.parse_args()
    audit(q.photo,q.receipt,q.old_register,q.old_report,q.output)

if __name__=="__main__":
    main()
