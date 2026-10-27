#!/usr/bin/env python3
"""Photo-ID and taxon-ID exact cross-arm union of real MOIS original ROI9 data.

Inputs (all SHA-pinned existing original numerical photo outcomes, NO new S3):
  * 1,043,930 original photo/observation IDs from six historical/K1/FCP arms,
    together with their 400,710 restored genuine numerical nine fractions;
  * 2,048 independent geo-season original photos (682 real nine fractions);
  * 4,096 taxon-once P2 original photos (1,467 real nine fractions);
  * 128-shard 131,072 actual independent original photo outcomes
    (48,784 true original ROI9, beside 500px box-Lab proxy only);
  * the immutable 42,111-photo-opportunity source taxa register.

The original photo is the master grain. Repeated original observation UUIDs
*within* legacy source may represent several real photos of one observation;
count these separately as photos, but flag them for effective observation
support. Across the three new designs, source observation UUIDs must be
independent, and no original photo ID may be reacquired.

This is a frozen source-ledger union, not the live mutable main branch. Other
concurrent FCP/93m runs require their own photo-ID/observation-ID join before
claiming a newer global total. None of the nine numeric composition outcomes
is a botanically verified petal spectrum or wild-native flora frequency.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

SCHEMA="mois_20261027_crossarm_original1181146_source_photo_strict9_taxon_union_v1"
PALETTE=("white","yellow","orange","red","pink","magenta","purple","blue","bronze")
FINAL_COLOURS=tuple("strict9_"+k for k in PALETTE)
HUES=("white","yellow_orange","red_pink","blue_purple")
SOURCE="original_1043930"
ARMS=("original_1043930","bias_refill_2048","p2_first_4096","bulk128_131072")
COUNTS={"original_1043930":1043930,"bias_refill_2048":2048,
        "p2_first_4096":4096,"bulk128_131072":131072}
VALIDS={"original_1043930":400710,"bias_refill_2048":682,
        "p2_first_4096":1467,"bulk128_131072":48784}
SHAS={
  "baseline_nine":"df0f562b58bdc7ba06e901530eb809f513a763e2b0524d3b195e99ee56e930ef",
  "baseline_obs":"9a4f3dae6fbcee967b71cb29576e71d52c48e1f1f790b7d3ea4750c4bcb0e6d4",
  "geo2048":"dea4699566430aa8fe33802e20d997646357ef46fd3429454959ba54fbd535a9",
  "p2":"bb2af1a15c7f2cb0091b19bcb58c6d921594f0cf4311b852d7b5f06bb752cbab",
  "bulk128":"2d648c832929f8ba59301906d5ea8d100bdfec6449c6e2828f6b9a910ddd0a02",
  "priority":"2dcc365c3b771a83477f1885703a2f573568bb3a7c996dcc11a194182068e8cc",
}
STRICT_SOURCE_TAXA=21235
SOURCE_TAXA=42111
MASTER_PHOTOS=93209351

def require(ok,reason):
    if not bool(ok):
        raise ValueError("MOIS_CROSSARM_1181146_REALROI9: "+reason)

def sha(path):
    digest=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(4<<20),b""):
            digest.update(chunk)
    return digest.hexdigest()

def text_id(column,*,photo=False):
    v=column.fillna("").astype(str).str.strip()
    v=v.mask(v.str.lower().isin(("none","nan","null")),"")
    if photo:
        require(v.str.fullmatch(r"\d+").all(),"Photo ID is missing or non-numeric")
    return v

def bool_field(col,name):
    require(str(col.dtype) in ("bool","boolean")
            and not col.isna().any(),
            name+" must be a genuine measured boolean, not strings or unknowns")
    return col.astype(bool)

def source_base(lake,original_obs):
    """Full baseline atlas is one photo per row; source observations are a join key."""
    atlas=pd.read_parquet(lake,columns=[
        "photo_id","inat_taxon_id","continuous9_strict",
        "original_four_hue",*("continuous9_"+p for p in PALETTE)])
    obs=pd.read_parquet(original_obs,columns=[
        "photo_id","observation_id","inat_taxon_id","date","latitude",
        "longitude","arm"])
    require(len(atlas)==len(obs)==COUNTS[SOURCE],
            "1043930 original image atlas or official source UUID union is truncated")
    for x in (atlas,obs):
        x["photo_id"]=text_id(x.photo_id,photo=True)
        x["inat_taxon_id"]=text_id(x.inat_taxon_id)
        require(x.photo_id.is_unique,"Repeated historical original photo IDs")
    d=atlas.merge(obs,on="photo_id",how="inner",validate="one_to_one",
                  suffixes=("_lake","_obs"),sort=False)
    require(len(d)==COUNTS[SOURCE],"Baseline strict original9 IDs not present in official original observation registry")
    original=d.inat_taxon_id_lake
    fromobs=d.inat_taxon_id_obs
    require((original.eq(fromobs)|original.eq("")|fromobs.eq("")).all(),
            "Baseline photo has conflicting source plant taxon IDs")
    d["inat_taxon_id"]=original.where(original.ne(""),fromobs)
    d["observed_on"]=d.date.fillna("").astype(str)
    d["source_subarm"]=d.arm.fillna("").astype(str)
    d["source_design"]="historical_K1_FCP_source_heterogeneous"
    d["original_roi_four_hue"]=d.original_four_hue.fillna("").astype(str)
    d["strict9"]=bool_field(d.continuous9_strict,"Baseline ROI9")
    return d.rename(columns={**{
        "continuous9_"+p:"strict9_"+p for p in PALETTE
    }})[["photo_id","observation_id","inat_taxon_id","source_subarm",
        "source_design","observed_on","latitude","longitude",
        "original_roi_four_hue","strict9",*FINAL_COLOURS]]

def supplement(frame,arm):
    mapping={
        "bias_refill_2048":("continuous9_strict","original_four_hue","continuous9_"),
        "p2_first_4096":("strict_original9","source_original_four_hue","strict9_"),
        "bulk128_131072":("strict_original_roi9","original_four_hue","roi9_"),
    }
    require(arm in mapping,"Supplemental source arm must be pre-registered")
    field,hue,prefix=mapping[arm]
    needed={"photo_id","observation_id","inat_taxon_id",field,hue,
            *(prefix+p for p in PALETTE)}
    require(needed<=set(frame.columns),
            arm+": original nine-fraction numerical photo evidence incomplete")
    d=frame.copy()
    for key in ("photo_id","observation_id","inat_taxon_id"):
        d[key]=text_id(d[key],photo=(key=="photo_id"))
    require(d.photo_id.is_unique and d.observation_id.is_unique and
            d.observation_id.ne("").all(),
            arm+": independent observation and original photo IDs not unique")
    d["strict9"]=bool_field(d[field],arm+": original9")
    d["original_roi_four_hue"]=d[hue].fillna("").astype(str)
    d["source_subarm"]=(
        d["design_cohort"].fillna("").astype(str)
        if "design_cohort" in d else (
        d["sampling_design_cohort"].fillna("").astype(str)
        if "sampling_design_cohort" in d else arm))
    d["source_design"]=arm
    d["observed_on"]=(
        d["observed_on"].fillna("").astype(str)
        if "observed_on" in d else "")
    for coord in ("latitude","longitude"):
        if coord not in d:
            d[coord]=np.nan
    require(d.source_subarm.ne("").all(),arm+": lost original source allocation")
    d=d.rename(columns={prefix+p:"strict9_"+p for p in PALETTE})
    return d[["photo_id","observation_id","inat_taxon_id","source_subarm",
        "source_design","observed_on","latitude","longitude",
        "original_roi_four_hue","strict9",*FINAL_COLOURS]]

def validated_arm(frame,name,*,n=None,strict_n=None):
    d=frame.copy()
    require(d.photo_id.is_unique and d.photo_id.str.fullmatch(r"\d+").all()
            and (n is None or len(d)==n),
            name+": original photo grain, numeric photo ID or source count failed")
    if name!=SOURCE:
        require(d.observation_id.is_unique and d.observation_id.ne("").all(),
                name+": duplicate or missing original independent observation UUID")
    require(d.strict9.dtype==bool
            and (strict_n is None or int(d.strict9.sum())==strict_n),
            name+": strict numerical ROI9 source-specific count changed")
    hues=d.original_roi_four_hue.isin(HUES)
    require((~d.strict9|hues).all(),"Strict9 measured in a four-hue-unknown source image")
    f=d[list(FINAL_COLOURS)].apply(pd.to_numeric,errors="coerce").to_numpy(np.float32)
    finite=np.isfinite(f).all(axis=1)
    real=d.strict9.to_numpy(bool)
    valid=finite&(f>=0).all(axis=1)&(f<=1+1e-5).all(axis=1)&(
        np.abs(np.where(finite[:,None],f,0).sum(axis=1)-1)<=.001)
    require((~real|valid).all(),name+": real nine-palette numeric vector invalid")
    # Original ambiguous-palette vectors may be retained by source for
    # explorations, but they are NOT eligible as strict flower phenotypes.
    f[~real,:]=np.nan
    for i,col in enumerate(FINAL_COLOURS):
        d[col]=f[:,i]
    return d

def merge_arms(baseline,geo,p2,bulk,source_priority,*,verify_source_counts=True):
    frames={
        SOURCE:validated_arm(baseline,SOURCE,n=COUNTS[SOURCE] if verify_source_counts else None,
                             strict_n=VALIDS[SOURCE] if verify_source_counts else None),
        "bias_refill_2048":validated_arm(geo,"bias_refill_2048",
                         n=COUNTS["bias_refill_2048"] if verify_source_counts else None,
                         strict_n=VALIDS["bias_refill_2048"] if verify_source_counts else None),
        "p2_first_4096":validated_arm(p2,"p2_first_4096",
                         n=COUNTS["p2_first_4096"] if verify_source_counts else None,
                         strict_n=VALIDS["p2_first_4096"] if verify_source_counts else None),
        "bulk128_131072":validated_arm(bulk,"bulk128_131072",
                         n=COUNTS["bulk128_131072"] if verify_source_counts else None,
                         strict_n=VALIDS["bulk128_131072"] if verify_source_counts else None),
    }
    # Preserve any genuine historical observation that contains multiple
    # photos in its own arm, but never pretend these are independent visits.
    seen={}
    observed_among_new=set()
    old_obs=set(text_id(frames[SOURCE].observation_id))-{""}
    result=[]
    obs_overlap_to_old=Counter()
    for arm in ARMS:
        current=frames[arm]
        for pid,ob,taxon in zip(current.photo_id,current.observation_id,current.inat_taxon_id):
            require(pid not in seen,
                    "Exact original photo ID appears in two image inference arms: "+pid)
            seen[pid]=(ob,taxon,arm)
        if arm!=SOURCE:
            these=set(current.observation_id)
            require(not these.intersection(observed_among_new),
                    arm+": previously attempted another new source observation UUID")
            obs_overlap_to_old[arm]=len(these.intersection(old_obs))
            observed_among_new.update(these)
        result.append(current)
    combined=pd.concat(result,ignore_index=True)
    require(combined.photo_id.is_unique and len(combined)==sum(len(v) for v in frames.values()),
            "Cross-arm append introduced duplicate original photo rows")
    priority=source_priority.copy()
    priority.inat_taxon_id=text_id(priority.inat_taxon_id)
    require(priority.inat_taxon_id.is_unique and len(priority)==(
        SOURCE_TAXA if verify_source_counts else len(priority)),
        "42,111 original nominal source taxon registration not complete")
    reference=set(priority.inat_taxon_id)
    grp=combined.groupby("inat_taxon_id",observed=True).agg(
        after_photo_attempts=("photo_id","size"),
        after_strict9=("strict9","sum"))
    source_before=frames[SOURCE].groupby("inat_taxon_id",observed=True).agg(
        before_photo_attempts=("photo_id","size"),
        before_strict9=("strict9","sum"))
    q=priority.set_index("inat_taxon_id").copy()
    q["before_photo_attempts"]=source_before.before_photo_attempts.reindex(q.index,fill_value=0).astype(int)
    q["before_strict9"]=source_before.before_strict9.reindex(q.index,fill_value=0).astype(int)
    q["after_photo_attempts"]=grp.after_photo_attempts.reindex(q.index,fill_value=0).astype(int)
    q["after_strict9"]=grp.after_strict9.reindex(q.index,fill_value=0).astype(int)
    require((q.after_strict9>=q.before_strict9).all()
            and (q.after_photo_attempts>=q.before_photo_attempts).all(),
            "A preexisting source taxon original photo outcome was lost")
    if verify_source_counts:
        require(q.before_strict9.gt(0).sum()==STRICT_SOURCE_TAXA
                and (q.before_strict9.to_numpy(int)==
                     q.n_strict_roi9.to_numpy(int)).all(),
                "Baseline original9 photo taxon counts differ from original SHA source")
    q["new_strict9"]=q.after_strict9-q.before_strict9
    q["newly_colour_supported"]=q.before_strict9.eq(0)&q.after_strict9.gt(0)
    q["newly_reach_5"]=q.before_strict9.lt(5)&q.after_strict9.ge(5)
    q["newly_reach_10"]=q.before_strict9.lt(10)&q.after_strict9.ge(10)
    def updated_priority(n,v):
        n=int(n);v=int(v)
        if n==0:return "P0_no_original_source_photo"
        if v==0 and n<10:return "P1_source_rare_no_colour_1to9"
        if v==0:return "P2_source_has_photos_no_colour_10plus"
        if v<5:return "P3_thin_colour_1to4"
        if v<10:return "P4_thin_colour_5to9"
        return "P5_existing_10plus_strict9"
    q["priority_after"]=[
        updated_priority(n,v) for n,v in
        zip(q.n_unique_source_photo_ids,q.after_strict9)]
    by_arm=[]
    for arm in ARMS:
        d=frames[arm]
        by_arm.append({"source_design":arm,
            "true_original_photo_attempts":len(d),
            "photo_ID_unique":int(d.photo_id.nunique()),
            "original_observation_ID_present":int(d.observation_id.ne("").sum()),
            "unique_observation_UUIDs_known":int(d.loc[
                d.observation_id.ne(""),"observation_id"].nunique()),
            "strict_original_ROI9":int(d.strict9.sum()),
            "strict9_technical_unknown":int((~d.strict9).sum()),
            "new_observation_UUIDs_shared_with_historical_baseline":obs_overlap_to_old[arm],
        })
    outside=combined.loc[~combined.inat_taxon_id.isin(reference)]
    outside_photo=source_before.loc[~source_before.index.isin(reference)]
    legacy_strict_not_in_register=int(outside_photo.before_strict9.sum())
    after_strict_not_in_register=int(outside.strict9.sum())
    require(after_strict_not_in_register>=legacy_strict_not_in_register,
            "Orphan historical 9-palette photo count changed")
    counts=Counter()
    for arm in ARMS:
        counts[arm]=len(frames[arm])
    result={
        "schema":SCHEMA,
        "original_global_upstream_photo_IDs":MASTER_PHOTOS,
        "nominal_source_taxa_in_original_frame":int(len(q)),
        "baseline_frozen_master_original_photo_attempts":int(len(frames[SOURCE])),
        "strict_original_roi9_before":int(frames[SOURCE].strict9.sum()),
        "actual_crossarm_unique_original_photo_IDs":int(len(combined)),
        "actual_crossarm_strict_original_roi9":int(combined.strict9.sum()),
        "photo_IDs_duplicated_across_source_arms":0,
        "baseline_source_photo_observations_known":int(frames[SOURCE].observation_id.ne("").sum()),
        "baseline_legacy_photo_rows_reusing_same_observation":int(
            frames[SOURCE].observation_id.ne("").sum() -
            frames[SOURCE].loc[frames[SOURCE].observation_id.ne(""),"observation_id"].nunique()),
        "new_source_independent_observations_conflicting_across_new_arms":0,
        "supplementary_original_observation_overlap_to_historical_baseline":
            dict(obs_overlap_to_old),
        "original_source_nominal_taxa_with_strict9_before":int(q.before_strict9.gt(0).sum()),
        "original_source_nominal_taxa_with_strict9_after":int(q.after_strict9.gt(0).sum()),
        "net_newly_measured_nominal_source_taxa":int(q.newly_colour_supported.sum()),
        "nominal_source_taxa_without_any_strict9_after":int(q.after_strict9.eq(0).sum()),
        "nominal_source_taxa_newly_reaching_five":int(q.newly_reach_5.sum()),
        "nominal_source_taxa_newly_reaching_ten":int(q.newly_reach_10.sum()),
        "baseline_source_strict9_photo_rows_outside_42111_nominal_source_taxa":legacy_strict_not_in_register,
        "after_source_strict9_photo_rows_outside_42111_nominal_source_taxa":after_strict_not_in_register,
        "prior_source_frame_photo_opportunities_unchanged":True,
        "128shard_medium500_lab_not_admitted_as_strict9":True,
        "separate_64_real_photo_box_lab_pilot_not_admitted_as_strict9":True,
        "no_original_image_downloading_or_colour_inference_for_this_union":True,
        "not_the_live_current_main_waves_without_new_crossarm_digest_recheck":True,
        "technical_unknowns_kept_as_missing":True,
        "mixed_reference_and_targeted_cohorts_not_used_as_global_population_frequency":True,
        "biological_open_flower_petal_wild_native_status_independently_verified":False,
        "source_designs":by_arm,
    }
    if verify_source_counts:
        require(result["actual_crossarm_unique_original_photo_IDs"]==1181146
                and result["actual_crossarm_strict_original_roi9"]==451643
                and result["original_source_nominal_taxa_with_strict9_before"]==21235,
                "Old+new source numerical photo and taxon totals do not reconcile")
    return combined,q.reset_index(),result

def run(baseline_nine,baseline_obs,new2048,p2,bulk128,priority,output,
        *,verify_sha=True):
    paths={
        "baseline_nine":Path(baseline_nine),"baseline_obs":Path(baseline_obs),
        "geo2048":Path(new2048),"p2":Path(p2),
        "bulk128":Path(bulk128),"priority":Path(priority)}
    output=Path(output)
    require(not output.exists(),"Do not overwrite a frozen cross-arm original photo union")
    if verify_sha:
        for key,p in paths.items():
            require(sha(p)==SHAS[key],
                    "Untrusted original photo numerical SHA: "+key)
    old=source_base(paths["baseline_nine"],paths["baseline_obs"])
    recent={
        "bias_refill_2048":pd.read_parquet(paths["geo2048"]),
        "p2_first_4096":pd.read_parquet(paths["p2"]),
        "bulk128_131072":pd.read_parquet(paths["bulk128"]),
    }
    catalogue=pd.read_csv(paths["priority"],compression="gzip",
                          dtype={"inat_taxon_id":str},low_memory=False)
    combined,frame,report=merge_arms(old,*(supplement(recent[arm],arm)
        for arm in ARMS[1:]),catalogue,verify_source_counts=True)
    output.mkdir(parents=True)
    files={}
    def keep_csv(df,name,compressed=False):
        name=name+".csv"+(".gz" if compressed else "")
        dest=output/name
        df.to_csv(dest,index=False,compression="gzip" if compressed else None)
        files[dest.name]=sha(dest)
    photo_file=output/"mois_1181146_unique_original_photo_continuous9_by_arm.parquet"
    combined.to_parquet(photo_file,index=False,compression="zstd")
    files[photo_file.name]=sha(photo_file)
    keep_csv(frame,"mois_42111_original_taxon_9palette_before_after",compressed=True)
    keep_csv(frame.loc[frame.newly_colour_supported,[
        "inat_taxon_id","n_unique_source_photo_ids","before_strict9","after_strict9",
        "new_strict9","priority_after"]],
        "mois_newly_measured_source_taxa_after_131072_bulk")
    keep_csv(pd.DataFrame(report["source_designs"]),"mois_crossarm_original_roi9_four_data_arms")
    status=frame.groupby(["priority","priority_after"],observed=True).agg(
        source_taxa=("inat_taxon_id","size"),
        previously_no9_source_taxa_now_coloured=("newly_colour_supported","sum"),
        authentic_new_strict9_photo_ROIs=("new_strict9","sum")
    ).reset_index()
    keep_csv(status,"mois_42111_source_taxon_gap_priority_transitions")
    report["source_sha256"]={k:sha(v) for k,v in paths.items()}
    report["signed_data_outputs_SHA256"]=files
    (output/"mois_crossarm_1181146_strict9_taxon_global_receipt.json").write_text(
        json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("MOIS_AUTHENTIC_1181146_CROSSARM_PHOTO9_AND_TAXON_COVERAGE",
          json.dumps({k:report[k] for k in (
              "actual_crossarm_unique_original_photo_IDs",
              "actual_crossarm_strict_original_roi9",
              "original_source_nominal_taxa_with_strict9_after",
              "net_newly_measured_nominal_source_taxa",
              "nominal_source_taxa_without_any_strict9_after",
              "baseline_legacy_photo_rows_reusing_same_observation"
          )},sort_keys=True),flush=True)
    return report

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for k in ("baseline-nine","baseline-obs","new2048","p2",
              "bulk128","priority","output"):
        parser.add_argument("--"+k,type=Path,required=True)
    x=parser.parse_args()
    run(x.baseline_nine,x.baseline_obs,x.new2048,x.p2,x.bulk128,
        x.priority,x.output)

if __name__=="__main__":
    main()
