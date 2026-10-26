#!/usr/bin/env python3
"""Collect genuine up-to-4096 P2-only original photo ROI-v4 outputs data-only.

Never load old full-million pixels or re-download original images. New-source
PHOTO IDs and observation UUIDs were frozen, source SHA-verified and limited
to a different previously zero-strict9 taxon per photo. Reject unsealed source,
undone worker units, duplicate original photos and any flower-colour inference
from technical unknowns. *Not* an unbiased worldwide floristic proportion.
"""
from __future__ import annotations

import argparse,hashlib,json,math
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.acquisition.mois_completed_pixel_checkpoint import collect_wave
from scripts.analysis.mois_reassemble_original_new4096_measurements import METHOD
from scripts.acquisition.mois_p2_first_4096_source import SCHEMA as SOURCE_SCHEMA,ARM

SCHEMA="mois_P2_real_original_photo_one_taxon_once_source_ROI9_v1"
BASE_PHOTOS=1_045_978
BASE_STRICT9=401_392
BASE_TAXA=21_320
UPSTREAM=42_111
HUES=("white","yellow_orange","red_pink","blue_purple")
PALETTE=("white","yellow","orange","red","pink","magenta",
         "purple","blue","bronze")
FRACTIONS=tuple("flower_fraction_"+x for x in PALETTE)

def require(cond,message):
    if not bool(cond):raise ValueError("MOIS_P2_REAL_ROI: "+message)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):h.update(b)
    return h.hexdigest()

def attach_provenance(rows,source):
    df=pd.DataFrame(rows)
    pfile=Path(source)/"selected_P2_photo_source_provenance.csv"
    m=json.loads((Path(source)/"manifest.json").read_text())
    require(sha(pfile)==m["photo_source_provenance_SHA256"],
            "Precoloured P2 original source photo-ID/probability manifest changed")
    p=pd.read_csv(pfile,dtype={"photo_id":str,"inat_taxon_id":str,
             "original_observation_id":str},low_memory=False)
    require(len(p)==len(df)==m["new_photo_attempts_planned"]
            and p.photo_id.is_unique and p.inat_taxon_id.is_unique
            and p.original_observation_id.is_unique
            and set(p.source_gap_priority)=={"P2"}
            and set(p.prior_strict9.astype(int))=={0},
            "Source prior zero9 per exact original taxon lost")
    df.photo_id=df.photo_id.astype(str)
    df.inat_taxon_id=df.inat_taxon_id.astype(str)
    df.observation_id=df.observation_id.astype(str)
    df=df.merge(
        p.drop(columns=["inat_taxon_id"]).rename(columns={
            "original_observation_id":"original_observation_id_source"}),
        on="photo_id",how="left",validate="one_to_one",sort=False)
    require(df.photo_id.is_unique and df.observation_id.is_unique
            and df.inat_taxon_id.is_unique
            and df.observation_id.eq(df.original_observation_id_source).all()
            and df.source_arm.eq(ARM).all()
            and df.source_gap_priority.eq("P2").all()
            and df.prior_strict9.astype(int).eq(0).all(),
            "Original P2 selected taxon/photo/observation IDs differ from ROI worker")
    return df

def classify_nine(df):
    require(set(FRACTIONS)<=set(df.columns),
            "Frozen original FCP YOLO11/EfficientSAM numeric nine-palette columns absent")
    fractions=df[list(FRACTIONS)].apply(pd.to_numeric,errors="coerce").to_numpy(float)
    finite=np.isfinite(fractions).all(axis=1)
    inrange=(np.where(np.isfinite(fractions),fractions,-1)>=0).all(axis=1)
    masses=np.where(np.isfinite(fractions),fractions,0).sum(axis=1)
    capped=(np.where(np.isfinite(fractions),fractions,0)<=1+1e-6).all(axis=1)
    raw_valid=finite&inrange&capped&(masses>0)
    known=df.morph.isin(HUES).to_numpy(bool)
    strict=known&raw_valid
    failed=df.measurement_status.eq("image_acquisition_failed").to_numpy(bool)
    require(not (strict&failed).any()
            and int(known.sum())<=len(df)
            and df.morph.isin((*HUES,"mixed_uncertain")).all(),
            "Technical ROI unknown or failed acquisition misclassified")
    result=df.copy()
    result["source_original_four_hue"]=result.morph.astype(str)
    result["strict_original9"]=strict
    result["technical_four_hue_unknown"]=~known
    result["photo_acquisition_failed"]=failed
    result["strict9_white_fraction"]=np.nan
    result["strict9_palette_entropy"]=np.nan
    renorm=np.full(fractions.shape,np.nan)
    renorm[raw_valid,:]=fractions[raw_valid,:]/masses[raw_valid,None]
    result.loc[strict,"strict9_white_fraction"]=renorm[strict,0]
    entropy=-np.sum(np.where(np.isfinite(renorm)&(renorm>0),
        renorm*np.log(np.maximum(renorm,1e-12)),0),axis=1)/math.log(9)
    result.loc[strict,"strict9_palette_entropy"]=entropy[strict]
    for i,name in enumerate(PALETTE):
        result["strict9_"+name]=np.where(strict,renorm[:,i],np.nan)
    require(result.loc[~result.strict_original9,
                       ["strict9_"+x for x in PALETTE]].isna().all().all(),
            "Original ROI9 unknown gained invented pigment-colour pixels")
    return result

def summarize(df):
    out=[]
    for (zone,season),g in df.groupby(["latitude_region","recorded_local_season"],
                                      observed=True,sort=True):
        good=g[g.strict_original9]
        out.append({
            "latitude_region":zone,"recorded_local_season":season,
            "new_original_photos":len(g),
            "new_original_independent_observations":int(g.observation_id.nunique()),
            "new_distinct_original_P2_taxa":int(g.inat_taxon_id.nunique()),
            "strict_original9_colour_count":len(good),
            "P2_taxa_newly_strict9_coloured":int(good.inat_taxon_id.nunique()),
            "strict9_success_rate":len(good)/len(g),
            "old_P2_taxa_with_zero_strict9_photo_labelled_as_such":True,
            "new_four_hue_classifiable":int(g.morph.isin(HUES).sum()),
            "new_acquisition_failed":int(g.photo_acquisition_failed.sum()),
            "new_four_hue_unknown":int(g.technical_four_hue_unknown.sum()),
            "strict9_measured_white_fraction_mean":(
                float(good.strict9_white_fraction.mean()) if len(good) else None),
        })
    return out

def collect(source,workers,out):
    source,workers,out=map(Path,(source,workers,out))
    require(not out.exists(),"Immutable P2 original ROI9 dataset must not be overwritten")
    m=json.loads((source/"manifest.json").read_text())
    require(m.get("schema")=="mois_full_source_observation_sampling_cohort_v1"
            and m.get("sampling_schema")==SOURCE_SCHEMA
            and m.get("full_frame_photo_ids")==93_209_351
            and m.get("P2_only_no_existing_strict9_source_taxa")
            and m.get("source_P5_high_colour_taxa_targeted")==0
            and m.get("one_photo_per_source_taxon_for_priority_P2")
            and m.get("all_source_photo_ID_and_observation_IDs_new_and_unique")
            and m.get("remaining_original_P2_source_taxa_before_selection")==17_727
            and 512<=m.get("new_photo_attempts_planned",0)<=4096,
            "Wrong unclassified P2 true original photo source")
    photos,sealed=collect_wave(source,workers,METHOD)
    n=m["new_photo_attempts_planned"]
    require(sealed["wave_complete"] and sealed["pending_photo_attempts"]==0
            and sealed["attempted"]==n==len(photos)
            and sealed["planned_photo_attempts"]==n,
            "Not all original P2 photo model worker terminals completed")
    d=classify_nine(attach_provenance(photos,source))
    nine=int(d.strict_original9.sum())
    known=int(d.morph.isin(HUES).sum())
    failed=int(d.photo_acquisition_failed.sum())
    require(nine<=known and known+int(d.technical_four_hue_unknown.sum())==n
            and failed<=n-known and failed+sealed["downloaded"]==n
            and sealed["classifiable"]==known
            and sealed["acquisition_failed"]==failed
            and sealed["unclassifiable_acquired"]==n-known-failed,
            "Actual classified/technical P2 original image counts do not reconcile")
    # Exactly one photo / one observation / one prior-uncoloured original taxon.
    # Each valid colour thus adds a NEW observed source taxon, not another P5 image.
    newtaxa=int(d.loc[d.strict_original9,"inat_taxon_id"].nunique())
    require(newtaxa==nine and int(d.prior_strict9.sum())==0,
            "A P2 colour added only depth rather than a genuinely uncoloured source taxon")
    out.mkdir(parents=True)
    file=out/"mois_realP2_original_photo_ROI9_and_unknowns.parquet"
    d.to_parquet(file,index=False,compression="zstd")
    csv_geo=out/"mois_realP2_new_strict9_by_region_local_season.csv"
    pd.DataFrame(summarize(d)).to_csv(csv_geo,index=False)
    effort_table=d.groupby("prior_attempt_effort_bin",observed=True).agg(
        photo_attempts=("photo_id","size"),
        strict9_new_taxa=("strict_original9","sum"),
        acquired_failures=("photo_acquisition_failed","sum"),
        photographic_four_hue_known=("technical_four_hue_unknown",lambda x:int((~x).sum()))
    ).reset_index()
    effort_file=out/"mois_realP2_prior_photo_failure_effort_and_new_colour.csv"
    effort_table.to_csv(effort_file,index=False)
    details={
        "schema":SCHEMA,
        "original_photo_attempts_before_new_P2":BASE_PHOTOS,
        "original_strict9_before_new_P2":BASE_STRICT9,
        "original_frame_taxa_with_strict9_before":BASE_TAXA,
        "full_original_frame_taxa":UPSTREAM,
        "source_P2_taxa_without_strict9_before":UPSTREAM-BASE_TAXA,
        "new_P2_photo_IDs_and_unique_original_observations_attempted":n,
        "new_P2_nominal_original_taxa_attempted":n,
        "new_P2_four_hue_known":known,
        "new_P2_strict9_original_ROIs":nine,
        "new_P2_taxa_first_nine_palette_colours":newtaxa,
        "new_P2_original_images_downloaded":int(sealed["downloaded"]),
        "new_P2_failed_downloads":failed,
        "new_P2_classification_or_ROI_unknown":n-known,
        "combined_attempted_original_photos_excluding_different64Lab":BASE_PHOTOS+n,
        "combined_strict9_original_photo_ROIs":BASE_STRICT9+nine,
        "combined_42111_source_taxa_with_strict9":BASE_TAXA+newtaxa,
        "remaining_42111_source_taxa_without_strict9":UPSTREAM-BASE_TAXA-newtaxa,
        "targeted_P5_old_supported_source_taxa":0,
        "additional_photo_image_inference_in_collector":0,
        "all_new_source_original_photo_or_observation_IDs_distinct_from_previous_1045978_and64":True,
        "true_botanically_verified_anthesis_wild_or_native_status":False,
        "global_photo_or_worldwide_flora_sampling_probability_identified":False,
        "original_photo_inclusion_probability_is_conditional_on_P2_prior_gap":True,
        "original_source_manifest_sha256":sha(source/"manifest.json"),
        "original_source_prior_photo_ID_roster_sha256":sha(
            source/"selected_P2_photo_source_provenance.csv"),
        "actual_new_photo_ROI9_parquet_sha256":sha(file),
        "region_season_9palette_counts_sha256":sha(csv_geo),
        "effort9_counts_sha256":sha(effort_file),
        "original_genuine_measured_source_terminal_SHA_audit":sealed,
    }
    (out/"mois_realP2_source_colours_and_taxon_coverage_receipt.json").write_text(
        json.dumps(details,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("MOIS_REAL_ORIGINAL_P2_FIRST9_TAXA_GAIN",json.dumps({
        k:details[k] for k in (
            "new_P2_photo_IDs_and_unique_original_observations_attempted",
            "new_P2_four_hue_known","new_P2_strict9_original_ROIs",
            "new_P2_taxa_first_nine_palette_colours",
            "combined_strict9_original_photo_ROIs",
            "combined_42111_source_taxa_with_strict9",
            "remaining_42111_source_taxa_without_strict9",
            "targeted_P5_old_supported_source_taxa",
        )},sort_keys=True),flush=True)
    return details

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for x in ("source","workers","output"):
        ap.add_argument("--"+x,type=Path,required=True)
    x=ap.parse_args()
    collect(x.source,x.workers,x.output)

if __name__=="__main__":
    main()
