#!/usr/bin/env python3
"""Decide whether MOIS should acquire more 93m-source original photographs.

Read authenticated 1,043,930 original image-pixel photo ID outcomes,
the FULL frozen 42,111-taxon/93,209,351 photo-opportunity register,
and a finished independent full-ID source composition audit.

This creates a 42,111-row source taxon prioritization and source-stratified
readiness profiles. It NEVER fabricates species prevalence, source photo IDs
not present in the frozen register, automatic ROI recovery, or colour features.
No image downloading; no global 93m scan is initiated. Recommended stage-1
candidate 'slots' are a quota to be matched to the exact archived PHOTO IDs
and distinct OBSERVATIONS in a separate source-led acquisition workflow.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd

SCHEMA="mois_original1043930_taxon_photo_gap_and_source_priority_v1"
FRAME=42111
FRAME_WITH_PHOTOS=42095
PHOTO_UNIVERSE=93209351
MEASURED=1043930
STRICT9=390541
FROZEN_COLOURS_SHA="b20bbefe15020c21c3b5b0ef7edf3f8e08148a4ee7dda8bf772a89bf6b126f447"
FROZEN_COUNTS_SHA="ef26c16252dd785f2d08bde3b18938c7624d73a54386fa3d3974dceb09ed97a9"

def require(cond,msg):
    if not cond:
        raise ValueError("MOIS_TAXON_GAP: "+msg)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(4<<20),b""):h.update(chunk)
    return h.hexdigest()

def priority(nphotos,valid):
    """Source-eligible taxa only, independent of the colour outcome itself."""
    nphotos=int(nphotos);valid=int(valid)
    if nphotos==0:return "P0_no_original_source_photo"
    if valid==0 and nphotos<=9:return "P1_source_rare_no_colour_1to9"
    if valid==0:return "P2_source_has_photos_no_colour_10plus"
    if valid<5:return "P3_thin_colour_1to4"
    if valid<10:return "P4_thin_colour_5to9"
    return "P5_existing_10plus_strict9"

def source_bin(nphotos):
    n=int(nphotos)
    if n==0:return "00_no_source"
    if n<=9:return "01_1to9"
    if n<=99:return "02_10to99"
    if n<=999:return "03_100to999"
    if n<=9999:return "04_1000to9999"
    return "05_10000plus"

def conservative_photo_slots(nsource,nattempted,max_images=3):
    """Photo OPPORTUNITY cap, not confirmed distinct/new image selection."""
    return min(int(nsource),max_images) if int(nsource)>0 else 0

def profile(d):
    result=[]
    for cat,g in d.groupby("priority",observed=True):
        result.append({
            "priority":str(cat),
            "source_frame_taxa":int(len(g)),
            "source_original_photo_opportunities":int(g.n_unique_source_photo_ids.sum()),
            "previously_attempted_source_photo_count":int(g.n_image_photo_attempts.sum()),
            "previous_strict_roi9_photo_count":int(g.n_strict_roi9.sum()),
            "taxa_with_one_or_fewer_attempted_images":int(g.n_image_photo_attempts.le(1).sum()),
            "taxa_with_5plus_attempted_images_but_zero_strict9":int(
                (g.n_image_photo_attempts.ge(5)&g.n_strict_roi9.eq(0)).sum()),
            "source_candidate_slots_first_one_per_taxon":int(
                np.minimum(g.n_unique_source_photo_ids.to_numpy(dtype=int),1).sum()),
            "source_candidate_slots_first_three_per_taxon":int(
                np.minimum(g.n_unique_source_photo_ids.to_numpy(dtype=int),3).sum()),
            "not_actual_new_photo_IDs_unless_source_exactly_deduplicated":True,
        })
    return result

def run(photo,receipt,upstream,structure,output):
    photo,receipt,upstream,structure,output=map(
        Path,(photo,receipt,upstream,structure,output))
    require(not output.exists(),"Source-photo prioritization report must be immutable")
    old=json.loads(receipt.read_text(encoding="utf-8"))
    audited=json.loads(structure.read_text(encoding="utf-8"))
    require(old.get("schema")=="mois_real1043930_original9_with_FCP447500_independent_source_reuse_v1"
            and old.get("latest_real_photo_ID_attempts")==MEASURED
            and old.get("strict_original9_vectors_now")==STRICT9
            and old.get("output_photo_lake_SHA256")==FROZEN_COLOURS_SHA==sha(photo),
            "Source original 1m original photo-ID colour data changed")
    require(audited.get("schema")==
            "mois_real1043930_original_photo_population_structure_and_missingness_v1"
            and audited.get("source_metadata_distinct_taxa")==FRAME
            and audited.get("source_metadata_total_photoIDs")==PHOTO_UNIVERSE
            and audited.get("original_image_photo_attempts")==MEASURED
            and audited.get("strict9_original_flower_ROI_photos")==STRICT9
            and audited.get("source_SHA256",{}).get("original_nine_source")==
             FROZEN_COLOURS_SHA,
            "Original source composition audit not finished or not same photo pixels")
    require(sha(upstream)==FROZEN_COUNTS_SHA,
            "Full species frame 2026-09-27 exact 93m-photo count source differs")
    source=pd.read_csv(upstream,dtype={"inat_taxon_id":str},low_memory=False)
    require(len(source)==FRAME and source.inat_taxon_id.is_unique
            and int(source.n_unique_source_photo_ids.sum())==PHOTO_UNIVERSE,
            "Full unique source taxon or candidate photo counts changed")
    measured=pd.read_parquet(photo,columns=[
        "photo_id","inat_taxon_id","source_arm_original",
        "continuous9_strict","continuous9_exploratory_ambiguous",
        "original_four_hue"])
    measured.photo_id=measured.photo_id.astype(str)
    measured.inat_taxon_id=measured.inat_taxon_id.fillna("").astype(str)
    require(len(measured)==MEASURED and measured.photo_id.is_unique
            and int(measured.continuous9_strict.sum())==STRICT9,
            "Measured 1043930 original photo-ID cohort unavailable")
    check=pd.DataFrame({"inat_taxon_id":measured.inat_taxon_id,
           "n_image_photo_attempts":1,
           "n_strict_roi9":measured.continuous9_strict.astype(int),
           "n_exploratory_ambiguous9":
               measured.continuous9_exploratory_ambiguous.astype(int),
           "n_four_hue":measured.original_four_hue.isin(
               ("white","yellow_orange","red_pink","blue_purple")).astype(int)})
    sums=check.groupby("inat_taxon_id",observed=True).sum()
    frame=source.set_index("inat_taxon_id").join(sums,how="left")
    for k in sums.columns:
        frame[k]=frame[k].fillna(0).astype(int)
    frame["source_photo_opportunity_bin"]=[
        source_bin(n) for n in frame.n_unique_source_photo_ids]
    frame["priority"]=[
        priority(n,v) for n,v in zip(
            frame.n_unique_source_photo_ids,frame.n_strict_roi9)]
    frame["photo_census_attempted_in_existing1043930"]=frame.n_image_photo_attempts.gt(0)
    frame["minimum_positive_source_photos_if_missing"]=frame.n_unique_source_photo_ids.gt(0)
    frame["raw_source_candidate_slots_at_most3_per_taxon"]=[
        conservative_photo_slots(n,int(a)) for n,a in zip(
            frame.n_unique_source_photo_ids,frame.n_image_photo_attempts)]
    frame["additional_photo_request_not_an_actual_photoID_queue"]=True
    frame["known_source_observation_cluster_eligibility_from_taxon_totals"]=False
    frame["number_of_strict9_photo_ROIs_to_reach_5"]=(
        5-frame.n_strict_roi9).clip(lower=0)
    outside=measured.loc[~measured.inat_taxon_id.isin(frame.index)].copy()
    outside["known_four_hue"]=outside.original_four_hue.isin(
        ("white","yellow_orange","red_pink","blue_purple"))
    out_by_arm=outside.groupby("source_arm_original",dropna=False,observed=True).agg(
        n_images=("photo_id","size"),
        strict9=("continuous9_strict","sum"),
        old_four=("known_four_hue","sum"),
    ).reset_index()
    require(int(len(outside))==
            audited["source_photo_attempts_OUTSIDE_42111_upstream_taxon_frame"],
            "Measured photo IDs outside original upstream taxon frame disagree")
    require(int(outside.continuous9_strict.sum())==
            audited["source_strict9_OUTSIDE_upstream_taxon_frame"]==0,
            "No true measured strict colours allowed outside this particular source frame")
    stats=profile(frame.reset_index())
    zero=frame.n_strict_roi9.eq(0)
    available=frame.n_unique_source_photo_ids.gt(0)
    low=frame.n_strict_roi9.between(1,4)
    source_rare=frame.n_unique_source_photo_ids.between(1,99)
    report={
       "schema":SCHEMA,
       "source_frame_taxa":FRAME,
       "original_source_photo_ID_candidates":PHOTO_UNIVERSE,
       "source_taxa_with_eligible_photos":int(available.sum()),
       "source_taxa_with_zero_photo_opportunities":int((~available).sum()),
       "strict9_photo_source":STRICT9,
       "taxa_with_at_least_one_valid9":int((~zero).sum()),
       "source_taxa_without_any_valid9":int(zero.sum()),
       "source_taxa_with_photos_but_no_valid9":int((zero&available).sum()),
       "source_photo_opportunities_in_taxa_with_photos_but_no_valid9":int(
           frame.loc[zero&available,"n_unique_source_photo_ids"].sum()),
       "share_of_full_source_photo_opportunities_in_taxa_without_valid9_percent":
           100*float(frame.loc[zero&available,"n_unique_source_photo_ids"].sum())/PHOTO_UNIVERSE,
       "source_rare_1to99_taxa_with_photos":int(source_rare.sum()),
       "source_rare_1to99_photo_opportunities":int(frame.loc[
           source_rare,"n_unique_source_photo_ids"].sum()),
       "source_rare_1to99_taxa_with_valid9":int((source_rare&~zero).sum()),
       "source_rare_1to99_taxa_without_valid9":int((source_rare&zero).sum()),
       "source_rare_1to99_taxa_missing9_percentage":100*float(
           (source_rare&zero).sum())/source_rare.sum(),
       "source_taxa_with_1to4_valid9":int(low.sum()),
       "source_taxa_with_5to9_valid9":int(frame.n_strict_roi9.between(5,9).sum()),
       "source_taxa_with_10plus_valid9":int(frame.n_strict_roi9.ge(10).sum()),
       "source_taxa_ever_attempted_but_zero_valid9":int(
          ((frame.n_image_photo_attempts>0)&zero).sum()),
       "zero_colour_taxa_with_5plus_images_already_attempted":int(
          ((frame.n_image_photo_attempts>=5)&zero).sum()),
       "zero_colour_taxa_with_10plus_images_already_attempted":int(
          ((frame.n_image_photo_attempts>=10)&zero).sum()),
       "zero_colour_taxa_with_only_1_attempted_photo":int(
          ((frame.n_image_photo_attempts==1)&zero).sum()),
       "source_attempted_photos_but_outside_42111_taxon_frame":int(len(outside)),
       "source_outside_frame_strict9_photo_count":int(
          outside.continuous9_strict.sum()),
       "source_outside_frame_photo_counts_by_original_arm":
          out_by_arm.to_dict(orient="records"),
       "collection_priority_distribution":stats,
       "individual_source_photo_ID_queues_created":False,
       "all_original_93209351_images_still_eligible_for_future_pass":True,
       "next_recommended_operation":(
         "first source-photo-ID-exact, observation-disjoint priority queues "
         "for ~20.9k no-nine-colour taxa, then light ROI failure quality/rescue "
         "and bounded photos; avoid indiscriminately expanding common taxa"),
       "source_and_selection_restrictions":(
         "Source eligible counts per taxon != probability sample of flora; "
         "multi-photo observations are correlated; legacy arms have unknown "
         "global inclusion probabilities; green/petal CIELAB not inferable from 9 bins"),
       "source_attestation_SHA256":{
          "real_1043930_source_9palette":sha(photo),
          "full_93209351_taxon_photo_opportunity_counts":sha(upstream),
          "source_frame_audit_receipt":sha(structure)},
       "new_photo_downloads":0,
       "new_image_ROI_inferences":0,
    }
    require(int((~zero).sum())==audited["upstream_frame_taxa_with_strict9_colour"]
            and int(zero.sum())==audited[
                "upstream_frame_taxa_NO_valid_original_nine_colour"]
            and int(frame.loc[zero&available,"n_unique_source_photo_ids"].sum())==
                 audited["upstream_total_photoIDs_in_taxa_with_no_strict9_colour"],
            "Original data structure audit and source taxon priority disagree")
    require(sum(v["source_frame_taxa"] for v in stats)==FRAME,
            "Original full-frame taxa missing from priority output")
    output.mkdir(parents=True)
    register=output/"mois_42111_source_taxon_photo_ROI_priority_register.csv.gz"
    frame.reset_index().sort_values(
        ["priority","n_unique_source_photo_ids","inat_taxon_id"],
        ascending=[True,True,True]).to_csv(
             register,index=False,compression="gzip")
    pd.DataFrame(stats).to_csv(
        output/"mois_42111_source_taxon_priority_group_summary.csv",index=False)
    pd.DataFrame(out_by_arm).to_csv(
        output/"mois_photo_attempts_outside_42111_source_taxa_by_arm.csv",index=False)
    (output/"mois_real1043930_taxon_gap_and_next_acquisition_decision.json").write_text(
        json.dumps(report,indent=2,sort_keys=True,ensure_ascii=False,allow_nan=False)+"\n")
    print("MOIS_REAL_1043930_SOURCE_TAXON_GAP_NEXT_PHOTO_DECISION",
          json.dumps({k:v for k,v in report.items() if k not in (
             "collection_priority_distribution",
             "source_attestation_SHA256")},
             sort_keys=True,allow_nan=False),flush=True)
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("photo","receipt","upstream","structure","output"):
        p.add_argument("--"+name,type=Path,required=True)
    a=p.parse_args()
    run(a.photo,a.receipt,a.upstream,a.structure,a.output)

if __name__=="__main__":
    main()
