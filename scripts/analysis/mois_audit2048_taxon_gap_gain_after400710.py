#!/usr/bin/env python3
"""Data-only MOIS taxon-coverage gain after exactly 2048 REAL new original ROIs.

Compare the SHA-frozen post-restoration 42,111-taxon 400,710 strict9 register
with the independently source-verified new 256 reference + 1792 targeted ROIs.
Keep all ROI unknowns and 64 previous-but-separate YOLO-box Lab outcomes out
of the nine-palette denominator. Never infer a global flora colour frequency.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

FRAME=42_111
BASE_PHOTOS=1_043_930
BASE_STRICT9=400_710
NEW_PHOTOS=2_048
NEW_STRICT9=682
PILOT_PHOTOS=64

PRIORITY_SHA="2dcc365c3b771a83477f1885703a2f573568bb3a7c996dcc11a194182068e8cc"
PRIORITY_RECEIPT_SHA="fe1d27302dffb6fe4fe06b38368fc5e90ee6012f4f6492098738ed6735713fe5"
NEW_PARQUET_SHA="dea4699566430aa8fe33802e20d997646357ef46fd3429454959ba54fbd535a9"
NEW_RECEIPT_SHA="aa5ee19ac329d1d87a96b4fa5c7aec67296aa0da8244b60d09ae20c10522d925"
PILOT_SHA="c8b9be91bb27edbe66b52725e34fba0eab8dd5ec5617bca1b99bdaf893ab36ed"
HUES=("white","yellow_orange","red_pink","blue_purple")
SCHEMA="mois_original2048_trueROI9_bias_first_taxon_coverage_gain_v1"

def require(condition, message):
    if not bool(condition):
        raise ValueError("MOIS_REAL2048_TAXON_GAP: "+message)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(4<<20),b""):
            h.update(chunk)
    return h.hexdigest()

def priority(source_photo_opportunities, strict9):
    n=int(source_photo_opportunities)
    v=int(strict9)
    require(n>=0 and v>=0,"Negative original source photo or ROI count")
    if n==0:return "P0_no_original_source_photo"
    if v==0 and n<=9:return "P1_source_rare_no_colour_1to9"
    if v==0:return "P2_source_has_photos_no_colour_10plus"
    if v<=4:return "P3_thin_colour_1to4"
    if v<=9:return "P4_thin_colour_5to9"
    return "P5_existing_10plus_strict9"

def apply_gain(frame, measured):
    """Keep the original 42k denominator; never cast outsider IDs to P0."""
    required_frame={"inat_taxon_id","n_strict_roi9","n_unique_source_photo_ids",
                    "n_image_photo_attempts","priority"}
    required_measured={"photo_id","observation_id","inat_taxon_id",
                       "continuous9_strict","sampling_design_cohort",
                       "latitude_region","recorded_local_season","original_four_hue",
                       "image_acquisition_failed"}
    require(required_frame<=set(frame.columns),"Original full-source taxon register columns absent")
    require(required_measured<=set(measured.columns),"Original measured pixel sample columns absent")
    before=frame.copy()
    before["inat_taxon_id"]=before.inat_taxon_id.astype(str)
    observed=measured.copy()
    for col in ("photo_id","observation_id","inat_taxon_id"):
        observed[col]=observed[col].astype(str)
    require(before.inat_taxon_id.is_unique and observed.photo_id.is_unique
            and observed.observation_id.is_unique,
            "Real photo, original observation or source taxon registry duplicated")
    require(observed.sampling_design_cohort.isin(("reference","targeted")).all(),
            "Source selection cohort missing, unknown or changed")
    require(observed.continuous9_strict.dtype == bool or
            str(observed.continuous9_strict.dtype)=="boolean",
            "Numeric ROI9 validity must be a true boolean, not a coerced string")
    require((~observed.continuous9_strict|
             observed.original_four_hue.isin(HUES)).all(),
            "Technically unknown four-hue images must never gain a strict9 label")
    require(not observed.loc[observed.image_acquisition_failed,"continuous9_strict"].any(),
            "Unacquired original photo cannot become a strict nine-palette colour")
    grouped=observed.groupby("inat_taxon_id",observed=True).agg(
        new_photo_attempts=("photo_id","size"),
        new_strict9=("continuous9_strict","sum"))
    known=before.set_index("inat_taxon_id",verify_integrity=True)
    outsiders=grouped.loc[~grouped.index.isin(known.index)]
    added=grouped.reindex(known.index,fill_value=0)
    known["new_photo_attempts"]=added.new_photo_attempts.astype(int)
    known["new_strict9"]=added.new_strict9.astype(int)
    require(known.loc[known.new_photo_attempts.gt(0),
                      "n_unique_source_photo_ids"].gt(0).all(),
            "New original photo assigned to taxon with zero upstream 93m photos")
    known["n_strict9_after"]=known.n_strict_roi9+known.new_strict9
    known["n_attempted_after"]=known.n_image_photo_attempts+known.new_photo_attempts
    known["priority_after"]=[
        priority(n,v) for n,v in zip(known.n_unique_source_photo_ids,
                                   known.n_strict9_after)]
    require(known.n_strict9_after.ge(known.n_strict_roi9).all(),
            "Old source nine-palette data lost while adding new photos")
    known["newly_colour_supported"]=known.n_strict_roi9.eq(0)&known.n_strict9_after.gt(0)
    known["newly_reach_5"]=known.n_strict_roi9.lt(5)&known.n_strict9_after.ge(5)
    known["newly_reach_10"]=known.n_strict_roi9.lt(10)&known.n_strict9_after.ge(10)
    return known.reset_index(),outsiders,observed

def audit(old_priority,old_receipt,new_parquet,new_receipt,pilot64,output,
          *,verify_shas=True,expected_frame=FRAME,expected_new=NEW_PHOTOS):
    files=list(map(Path,(old_priority,old_receipt,new_parquet,new_receipt,pilot64)))
    old_priority,old_receipt,new_parquet,new_receipt,pilot64=files
    output=Path(output)
    require(not output.exists(),"Refuse to overwrite published or staged source audit")
    if verify_shas:
        expected=[PRIORITY_SHA,PRIORITY_RECEIPT_SHA,NEW_PARQUET_SHA,
                  NEW_RECEIPT_SHA,PILOT_SHA]
        for p,digest in zip(files,expected):
            require(sha(p)==digest,"Original SHA-pinned input changed: "+p.name)
    old_summary=json.loads(old_receipt.read_text(encoding="utf-8"))
    new_summary=json.loads(new_receipt.read_text(encoding="utf-8"))
    require(old_summary.get("schema")==
            "mois_original1043930_recovered400710_source_taxon_priorities_v2"
            and old_summary.get("source_42111_taxa")==expected_frame
            and old_summary.get("strict9_after_no_image_downloads")==BASE_STRICT9,
            "Not the final 400710-original nine-palette source taxonomy baseline")
    require(new_summary.get("schema")==
            "mois_20261026_actual2048_novel_geo_season_balance_original_FCP_ROI9_v1"
            and new_summary.get("new_independent_original_photo_ids_attempted")==expected_new
            and new_summary.get("new_strict_original9_vectors")==NEW_STRICT9
            and new_summary.get("new_separate_reference_original_photos")==256
            and new_summary.get("new_separate_targeted_original_photos")==1792
            and new_summary.get("two_reference_and_targeted_design_arms_not_naively_poolable"),
            "New real measured cohort is not exact frozen reference+gap-targeted source")
    frame=pd.read_csv(old_priority,compression="gzip",dtype={"inat_taxon_id":str},
                      low_memory=False)
    measured=pd.read_parquet(new_parquet)
    require(len(frame)==expected_frame and len(measured)==expected_new,
            "Source frame or new original photograph census incomplete")
    require(int(frame.n_strict_roi9.sum())==
            BASE_STRICT9-int(old_summary["strict9_added_OUTSIDE_42111_source_taxon_frame"]),
            "Old original ROI9 source-frame counts differ from signed legacy recovery")
    require(int(frame.n_strict_roi9.gt(0).sum())==
            old_summary["source_taxa_with_any_strict9_NOW"],
            "Old original source taxonomic coverage signature changed")
    require(measured.continuous9_strict.sum()==NEW_STRICT9
            and measured.image_acquisition_failed.sum()==8,
            "New 2048 original measured ROI9 technical outcomes changed")
    known,outsiders,new=apply_gain(frame,measured)
    pilot=pd.read_csv(pilot64,dtype={"photo_id":str,"observation_id":str},
                      low_memory=False)
    require(len(pilot)==PILOT_PHOTOS and pilot.photo_id.is_unique
            and pilot.observation_id.is_unique,
            "Separately measured 64 original Lab pilot provenance not unique")
    require(not set(pilot.photo_id).intersection(new.photo_id)
            and not set(pilot.observation_id).intersection(new.observation_id),
            "New 2048 original photo or observation IDs overlap 64 existing Lab pilot")
    joined=new.merge(frame[["inat_taxon_id","priority","n_strict_roi9"]],
                     on="inat_taxon_id",how="left",validate="many_to_one")
    joined["priority"]=joined.priority.fillna("OUTSIDE_ORIGINAL_42111_SOURCE_TAXA")
    grouped=[]
    for (p,role),g in joined.groupby(["priority","sampling_design_cohort"],
                                    observed=True,dropna=False,sort=True):
        strict=g[g.continuous9_strict]
        group_before=g[g.n_strict_roi9.eq(0)]
        grouped.append({
            "priority_before":str(p),"cohort":str(role),
            "actual_original_photo_attempts":len(g),
            "actual_source_observations":int(g.observation_id.nunique()),
            "new_strict9":len(strict),
            "new_strict9_success_rate":len(strict)/len(g),
            "unique_source_taxa_attempted":int(g.inat_taxon_id.nunique()),
            "unique_source_taxa_newly_coloured":int(
                strict.loc[strict.n_strict_roi9.eq(0),"inat_taxon_id"].nunique()),
            "attempts_in_taxa_with_prior_zero":int(len(group_before)),
            "failed_image_acquisitions":int(g.image_acquisition_failed.sum()),
        })
    original_prior=Counter(frame.priority)
    current_prior=Counter(known.priority_after)
    priorities=sorted(set(original_prior)|set(current_prior))
    all_rows=[]
    for name in priorities:
        all_rows.append({"priority":name,
            "source_taxa_before":original_prior[name],
            "source_taxa_after":current_prior[name],
            "net_taxa":current_prior[name]-original_prior[name]})
    novel=known.loc[known.newly_colour_supported,
                    ["inat_taxon_id","priority","priority_after","n_strict_roi9",
                     "new_photo_attempts","new_strict9","n_strict9_after"]]
    photo_summary=(
        joined.groupby(["sampling_design_cohort","latitude_region",
                        "recorded_local_season"],observed=True,dropna=False)
        .agg(new_original_photos=("photo_id","size"),
             new_strict9=("continuous9_strict","sum"),
             image_acquisition_failed=("image_acquisition_failed","sum"),
             distinct_nominal_source_taxa=("inat_taxon_id","nunique"))
        .reset_index()
    )
    result={
       "schema":SCHEMA,
       "original_upstream_photo_ID_opportunities":93_209_351,
       "original_42111_source_taxa":expected_frame,
       "before_original_photo_attempts":BASE_PHOTOS,
       "before_strict9_worldwide_original_photo_outcomes":BASE_STRICT9,
       "new_source_disjoint_real_photo_attempts":expected_new,
       "new_actual_strict9_photo_ROIs":NEW_STRICT9,
       "after_original_photo_attempts_excluding_separate64Lab":BASE_PHOTOS+expected_new,
       "after_source_original_strict9_outcomes_excluding_separate64Lab":BASE_STRICT9+NEW_STRICT9,
       "source_frame_taxa_with_strict9_before":int(known.n_strict_roi9.gt(0).sum()),
       "source_frame_taxa_with_strict9_after":int(known.n_strict9_after.gt(0).sum()),
       "newly_strict9_covered_source_frame_taxa":int(known.newly_colour_supported.sum()),
       "new_source_taxa_reaching_5_strict9":int(known.newly_reach_5.sum()),
       "new_source_taxa_reaching_10_strict9":int(known.newly_reach_10.sum()),
       "new_photo_attempts_OUTSIDE_original42111_frame":int(outsiders.new_photo_attempts.sum()),
       "new_strict9_OUTSIDE_original42111_frame":int(outsiders.new_strict9.sum()),
       "new_original_photo_ID_overlap_with_prior_64_real_Lab_pilot":0,
       "new_original_observation_ID_overlap_with_prior_64_real_Lab_pilot":0,
       "reference_original_attempts":256,
       "reference_strict9":int(joined.loc[
           joined.sampling_design_cohort.eq("reference"),"continuous9_strict"].sum()),
       "targeted_original_attempts":1792,
       "targeted_strict9":int(joined.loc[
           joined.sampling_design_cohort.eq("targeted"),"continuous9_strict"].sum()),
       "individual_wild_native_petal_organ_status_independently_verified":False,
       "targeted_selected_using_existing_ROI9_gaps_not_worldwide_probability_sample":True,
       "region_and_season_baseline_joint_support_before_does_not_exist_in_source_receipts":True,
       "new_pixels_re_downloaded_or_re_inferred_by_this_analysis":0,
       "source_SHA256":{x.name:sha(x) for x in files},
    }
    require(result["reference_strict9"]==88 and result["targeted_strict9"]==594,
            "Actual per-design-cohort colour terminal counts mutated")
    require((result["new_strict9_OUTSIDE_original42111_frame"]+
             int(known.new_strict9.sum()))==NEW_STRICT9,
            "Original-frame and outside-frame measured colours do not add up")
    require(result["source_frame_taxa_with_strict9_after"]>=
            result["source_frame_taxa_with_strict9_before"],"Net taxon coverage decreased")
    output.mkdir(parents=True)
    pd.DataFrame(grouped).to_csv(
        output/"real2048_newROI9_success_by_prior_taxon_gap_and_cohort.csv",index=False)
    pd.DataFrame(all_rows).to_csv(
        output/"real2048_taxon_priority_BEFORE_AFTER.csv",index=False)
    novel.to_csv(output/"real2048_newly_covered_93m_source_taxa.csv",index=False)
    photo_summary.to_csv(output/"real2048_region_local_season_genuine9_increment.csv",
                         index=False)
    (output/"real2048_taxon_coverage_gap_delta_receipt.json").write_text(
        json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+"\n",
        encoding="utf-8")
    print("MOIS_REAL2048_REDUCED_GLOBAL_SOURCE_TAXON9_MISSINGNESS",
          json.dumps({k:v for k,v in result.items() if k not in (
              "source_SHA256","targeted_selected_using_existing_ROI9_gaps_not_worldwide_probability_sample",
              "individual_wild_native_petal_organ_status_independently_verified")},
              sort_keys=True),flush=True)
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for x in ("old-priority","old-receipt","new-parquet","new-receipt",
              "pilot64","output"):
        p.add_argument("--"+x,type=Path,required=True)
    v=p.parse_args()
    audit(v.old_priority,v.old_receipt,v.new_parquet,v.new_receipt,
          v.pilot64,v.output)

if __name__=="__main__":
    main()
