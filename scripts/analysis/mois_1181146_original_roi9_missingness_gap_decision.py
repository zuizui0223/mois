#!/usr/bin/env python3
"""Original photo ROI9 missingness by region, local season and sampling design.

Uses the SHA-frozen 1,181,146 one-photo-per-row four-arm cross-source union,
not raw images or the separate 500px box colour proxy. A fixed, deterministic
photo-and-taxon-distinct review sample diagnoses which prior original images
failed. Post-outcome source-taxon gap priority is only an operational triage
variable: never a predictor for unbiased global biological inference.
"""
from __future__ import annotations
import argparse,hashlib,json
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
import pandas as pd

SCHEMA="mois_real1181146_fourarm_strictroi9_geoseason_missingness_triage_v1"
PHOTO_SHA="3fbd0f0ec7eb347b3bbc7ab509ec9312ec6d797ed932dcf7bb425959898c3835"
TAXON_SHA="ea602238f3cd152fdd2e36e9ce4bcbc77b795a35dce802e9305f6f6cddd77fb8"
ARMS=("historical_K1_FCP_source_heterogeneous","bias_refill_2048","p2_first_4096","bulk128_131072")
ZERO_PRIORITIES=("P0_no_original_source_photo",
    "P1_source_rare_no_colour_1to9","P2_source_has_photos_no_colour_10plus")
HUES=("white","yellow_orange","red_pink","blue_purple")
EXPECTED=(1181146,451643,42111,30048,12063)
SEED="MOIS_TRUE1181146_SOURCE_ROI9_UNKNOWN_REVIEW_20261027"

def require(ok,msg):
    if not bool(ok):raise ValueError("MOIS_1181146_ROI9_TRIAGE: "+msg)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(4<<20),b""):h.update(b)
    return h.hexdigest()

def latitude_band(lat):
    x=pd.to_numeric(lat,errors="coerce").to_numpy(float)
    a=np.full(x.shape,"unknown",dtype=object)
    good=np.isfinite(x)&(x>=-90)&(x<=90)
    for condition,label in (
        (x< -66.5,"S_antarctic"),((x>=-66.5)&(x< -55),"S_subpolar"),
        ((x>=-55)&(x< -23.5),"S_temperate"),
        ((x>=-23.5)&(x<0),"S_tropical"),
        ((x>=0)&(x<23.5),"N_tropical"),
        ((x>=23.5)&(x<55),"N_temperate"),
        ((x>=55)&(x<66.5),"N_subpolar"),
        (x>=66.5,"N_arctic")
    ):a[good&condition]=label
    return pd.Series(a,index=lat.index)

def recorded_local_season(dates,lat):
    parsed=pd.to_datetime(dates.fillna("").astype(str),format="%Y-%m-%d",errors="coerce")
    month=parsed.dt.month
    x=pd.to_numeric(lat,errors="coerce").to_numpy(float)
    n=np.isfinite(x)&(x>=0)&(x<=90)
    s=np.isfinite(x)&(x<0)&(x>=-90)
    a=np.full(len(month),"unknown",dtype=object)
    for codes,label in (((12,1,2),"winter"),((3,4,5),"spring"),
                        ((6,7,8),"summer"),((9,10,11),"autumn")):
        good=month.isin(codes).to_numpy()
        a[n&good]=label
        reverse={"winter":"summer","spring":"autumn",
                 "summer":"winter","autumn":"spring"}[label]
        a[s&good]=reverse
    return pd.Series(a,index=dates.index)

def source_action(row):
    n=int(row.n_unique_source_photo_ids)
    attempts=int(row.after_photo_attempts)
    good=int(row.after_strict9)
    priority=str(row.priority_after)
    require(n>=0 and attempts>=0 and good>=0,"Negative original photo opportunity or ROI9")
    if good>0:return "SUPPORTED_NOT_ZERO_ROI9"
    require(priority in ZERO_PRIORITIES,"Unexpected no-colour priority")
    if n==0:return "P0_NO_SOURCE_OPPORTUNITY"
    if attempts>=5:return "REVIEW_EXISTING_REPEATED_ROI9_FAILURES_BEFORE_FETCH"
    if attempts>=1:return "REVIEW_FEW_OLD_UNKNOWN_PHOTOS_AND_SOURCE_IF_NEEDED"
    return "PREFER_FIRST_SOURCE_PHOTO_IF_AVAILABLE"

def review_queue(frame,quota=512):
    """Blind reproducible across source-design/region/season, one photo/taxon."""
    need=frame.loc[
        ~frame.strict9 &
        frame.priority_after.isin(ZERO_PRIORITIES[1:]) &
        frame.inat_taxon_id.ne("") &
        frame.observation_id.ne(""),
        ["photo_id","observation_id","inat_taxon_id","source_design",
         "latitude_zone","local_season","priority_after",
         "after_photo_attempts","n_unique_source_photo_ids",
         "original_roi_four_hue"]
    ].copy()
    if need.empty:return pd.DataFrame()
    need["stratum"]=(need.source_design.astype(str)+"|"+
                    need.latitude_zone.astype(str)+"|"+
                    need.local_season.astype(str)+"|"+
                    need.priority_after.astype(str))
    need["rank"]=need.photo_id.map(
        lambda x:hashlib.sha256((SEED+"|"+x).encode()).hexdigest())
    need=need.sort_values(["stratum","rank"],kind="stable")
    need=need.drop_duplicates(["stratum","inat_taxon_id"])
    groups={s:list(part.itertuples(index=False))
            for s,part in need.groupby("stratum",sort=True)}
    positions={k:0 for k in groups}
    selected=[];used=set()
    while len(selected)<quota:
        progress=False
        for key in sorted(groups):
            entries=groups[key];i=positions[key]
            while i<len(entries) and entries[i].inat_taxon_id in used:i+=1
            positions[key]=i
            if i>=len(entries):continue
            r=entries[i]
            used.add(r.inat_taxon_id)
            selected.append({
                "photo_id":r.photo_id,"observation_id":r.observation_id,
                "inat_taxon_id":r.inat_taxon_id,"source_design":r.source_design,
                "latitude_zone":r.latitude_zone,"local_season":r.local_season,
                "postoutcome_source_gap":r.priority_after,
                "previous_original_photo_attempts":int(r.after_photo_attempts),
                "source_photo_opportunities":int(r.n_unique_source_photo_ids),
                "reason":"FOUR_HUE_KNOWN_BUT_NO_STRICT9" if
                    r.original_roi_four_hue in HUES else
                    "FOUR_HUE_UNKNOWN_AND_NO_STRICT9",
                "botanical_organ_or_flower_stage_verified":False
            })
            positions[key]=i+1;progress=True
            if len(selected)==quota:break
        if not progress:break
    output=pd.DataFrame(selected)
    require(output.photo_id.is_unique and output.inat_taxon_id.is_unique,
            "Review photos/taxa must be independent across strata")
    return output

def audit(source_photo,source_taxa,out,*,verify_sha=True,verify_census=True,
          review_size=512,restored_geo_receipt=None):
    source_photo,source_taxa,out=map(Path,(source_photo,source_taxa,out))
    require(not out.exists() and 16<=review_size<=2048,
            "Immutable output exists or review quota unsupported")
    if verify_sha:
        require(sha(source_taxa)==TAXON_SHA,
                "Authoritative 42,111-source taxon gap register SHA changed")
        if restored_geo_receipt is None:
            require(sha(source_photo)==PHOTO_SHA,
                    "Frozen original-image 9-palette source SHA changed")
        else:
            restored=json.loads(Path(restored_geo_receipt).read_text())
            from scripts.analysis.mois_restore_bulk128_photo_geography_from_original_source import (
                SCHEMA as RESTORED_SCHEMA, WORLD_SHA as OLD_WORLD_SHA,
                SOURCE_MANIFEST_SHA, SOURCE_TAR_SHA)
            require(restored.get("schema")==RESTORED_SCHEMA
                    and restored.get("input_world_original_parquet_SHA256")==OLD_WORLD_SHA
                    and restored.get("signed_tier10_source_worker_manifest_SHA256")==SOURCE_MANIFEST_SHA
                    and restored.get("signed_tier10_source_worker_archive_SHA256")==SOURCE_TAR_SHA
                    and restored.get("source_georestored_original_parquet_SHA256")==sha(source_photo)
                    and restored.get("same_original_image_ID_and_original_9palette_unmodified")
                    and restored.get("source_original_bulk_photo_geopoints_available",0)>0,
                    "Restored bulk128 source photo geography lacks cryptographic original provenance")
    cols=["photo_id","observation_id","inat_taxon_id","source_design",
          "source_subarm","observed_on","latitude","longitude",
          "original_roi_four_hue","strict9","strict9_white"]
    d=pd.read_parquet(source_photo,columns=cols)
    q=pd.read_csv(source_taxa,compression="gzip",
                  dtype={"inat_taxon_id":str},low_memory=False)
    required={"inat_taxon_id","n_unique_source_photo_ids",
              "after_photo_attempts","after_strict9","priority_after"}
    require(required<=set(q.columns),"Original source-taxon gap fields absent")
    if verify_census:
        require((len(d),int(d.strict9.sum()),len(q),
                 int(q.after_strict9.gt(0).sum()),
                 int(q.after_strict9.eq(0).sum()))==EXPECTED,
                "One of original photo-ID / real numeric ROI9 / 42111-taxon totals drifted")
    for k in ("photo_id","observation_id","inat_taxon_id"):
        d[k]=d[k].fillna("").astype(str).str.strip()
    q.inat_taxon_id=q.inat_taxon_id.fillna("").astype(str).str.strip()
    require(d.photo_id.is_unique,
            "Duplicate photo IDs in authentic original ROI9 ledger")
    require(d.photo_id.str.fullmatch(r"\d+").all(),
            "Exact original photo IDs were malformed or nonnumeric")
    require(q.inat_taxon_id.is_unique,
            "Nominal source taxon IDs duplicated in authoritative source register")
    require(d.source_design.isin(ARMS).all(),
            "Original source acquisition-arm labels unexpected: "+
            repr(d.loc[~d.source_design.isin(ARMS),"source_design"].value_counts(dropna=False).head(8).to_dict()))
    require(str(d.strict9.dtype) in ("bool","boolean"),
            "Original valid-nine-palette column is not boolean: "+str(d.strict9.dtype))
    require(d.strict9.notna().all(),
            "Original strict9 photo success column contains undefined booleans")
    ok=d.strict9.astype(bool)
    whites=pd.to_numeric(d.strict9_white,errors="coerce")
    require((~ok|d.original_roi_four_hue.isin(HUES)).all() and
            (~ok|whites.between(0,1,inclusive="both")).all() and
            whites[~ok].isna().all(),
            "Original strict9 blank/unknown cannot be relabeled as a measured white flower")
    d["latitude_zone"]=latitude_band(d.latitude)
    d["local_season"]=recorded_local_season(d.observed_on,d.latitude)
    d["four_hue_known"]=d.original_roi_four_hue.isin(HUES)
    require(int(d.four_hue_known.sum())>=int(ok.sum()),
            "Four-hue colour cannot be sparser than source strict9")
    j=d.merge(q[["inat_taxon_id","priority_after","after_strict9",
                 "after_photo_attempts","n_unique_source_photo_ids"]],
              on="inat_taxon_id",how="left",validate="many_to_one",sort=False)
    j["priority_after"]=j.priority_after.fillna("OUTSIDE_42111_SOURCE_FRAME")
    reports={}
    for label,keys in (
       ("source_design",["source_design"]),
       ("source_subarm",["source_design","source_subarm"]),
       ("region_season",["source_design","latitude_zone","local_season"]),
       ("postoutcome_taxon_priority",["source_design","priority_after"])
    ):
        t=(j.groupby(keys,observed=True,dropna=False)
           .agg(original_photo_attempts=("photo_id","size"),
                original_observation_support=("observation_id","nunique"),
                distinct_source_taxa=("inat_taxon_id","nunique"),
                four_hue_known=("four_hue_known","sum"),
                original9_true=("strict9","sum")).reset_index())
        t["original9_unknown"]=t.original_photo_attempts-t.original9_true
        t["strict9_measurement_rate"]=t.original9_true/t.original_photo_attempts
        reports[label]=t
    q=q.copy()
    q["action"]=q.apply(source_action,axis=1)
    q["remaining_source_photo_opportunity_upper_bound"]=(
        q.n_unique_source_photo_ids-q.after_photo_attempts).clip(lower=0)
    reviews=review_queue(j,quota=review_size)
    if verify_census:
        p=Counter(q.loc[q.after_strict9.eq(0),"priority_after"])
        require(p=={ZERO_PRIORITIES[0]:16,ZERO_PRIORITIES[1]:2862,
                    ZERO_PRIORITIES[2]:9185},"Source-zero ROI9 taxon priority changed")
        require(len(reviews)==review_size,
                "Existing original photo bank could not supply reproducible QC audit")
    reports["manual_roi_failure_review_photoIDs"]=reviews
    reports["zero_colour_taxon_next_action"]=q.loc[q.after_strict9.eq(0),[
        "inat_taxon_id","priority_after","n_unique_source_photo_ids",
        "after_photo_attempts","after_strict9",
        "remaining_source_photo_opportunity_upper_bound","action"
    ]]
    out.mkdir(parents=True)
    paths={}
    for name,table in reports.items():
        dest=out/(name+".csv")
        table.to_csv(dest,index=False)
        paths[dest.name]=sha(dest)
    next_actions=Counter(q.loc[q.after_strict9.eq(0),"action"])
    receipt={
        "schema":SCHEMA,
        "original_93m_photo_opportunities":93209351,
        "exact_source_original_photo_IDs":len(d),
        "real_strict9_original_photos":int(ok.sum()),
        "original9_unknown_photo_attempts":int((~ok).sum()),
        "source_nominal_taxa":len(q),
        "source_nominal_taxa_strict9_supported":int(q.after_strict9.gt(0).sum()),
        "source_nominal_taxa_zero_strict9":int(q.after_strict9.eq(0).sum()),
        "source_P0_no_photo_taxa":int(q.priority_after.eq(ZERO_PRIORITIES[0]).sum()),
        "source_P1_rare_photos_zero9_taxa":int(q.priority_after.eq(ZERO_PRIORITIES[1]).sum()),
        "source_P2_photo_supported_zero9_taxa":int(q.priority_after.eq(ZERO_PRIORITIES[2]).sum()),
        "original_independent_observations_with_known_ID":int(d.loc[
            d.observation_id.ne(""),"observation_id"].nunique()),
        "original_photos_with_unknown_latitude_zone":int(d.latitude_zone.eq("unknown").sum()),
        "original_photos_with_unknown_local_season":int(d.local_season.eq("unknown").sum()),
        "zero9_source_taxa_action_counts":dict(next_actions),
        "manual_unknown_original_photo_review_targets":len(reviews),
        "one_diagnostic_review_photo_per_taxon":True,
        "no_new_S3_image_downloads":True,
        "no_original_ROI9_inferences":True,
        "technical_unknowns_preserved":True,
        "taxon_gap_priority_is_POST_outcome_not_an_ecological_predictor":True,
        "different_source_reference_and_targeted_arms_NOT_unbiased_world_prevalence":True,
        "manual_botanical_petals_and_native_status_NOT_yet_verified":True,
        "source_input_SHAs":{source_photo.name:sha(source_photo),
                            source_taxa.name:sha(source_taxa)},
        "output_SHAs":paths
    }
    dest=out/"mois_real1181146_original_roi9_bias_action_receipt.json"
    dest.write_text(json.dumps(receipt,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("MOIS_REAL1181146_SOURCE_ROI9_BIAS_AND_P2_MISSINGNESS",
          json.dumps({k:receipt[k] for k in (
              "real_strict9_original_photos","original9_unknown_photo_attempts",
              "source_P2_photo_supported_zero9_taxa",
              "zero9_source_taxa_action_counts",
              "manual_unknown_original_photo_review_targets"
          )},sort_keys=True),flush=True)
    return receipt

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--photo-parquet",required=True,type=Path)
    p.add_argument("--taxon-register",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--review-target",type=int,default=512)
    p.add_argument("--restored-geo-receipt",type=Path)
    a=p.parse_args()
    audit(a.photo_parquet,a.taxon_register,a.output,
          review_size=a.review_target,restored_geo_receipt=a.restored_geo_receipt)

if __name__=="__main__":
    main()
