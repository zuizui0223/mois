#!/usr/bin/env python3
"""Incremental species-first 9-palette Atlas update after signed third image wave.

This avoids rebuilding all 1.18 million earlier nine-colour ROI vectors.
Wave3 was selected exclusively from prior no-strict9 taxa; each independent
observation can add at most one newly colour-supported source taxon. Preserve
uncertainty, never equate sampled photographic hue with species fixation.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
from scripts.analysis.mois_species_first_palette9 import (
    COLUMNS,PALETTE,class_labels,sha256
)
from scripts.analysis.mois_species_first_incremental_first256 import require
from scripts.analysis.mois_species_first_real1280_cumulative import choose_third_wave

SCHEMA="mois_species_first_wave3_incremental_compact_catalogue_v1"
PRIOR_ATLAS_SHA="e1438342218275078620cdd503d3a4b5e774eaa7c6cce0c160e8969eab875c52"
PRIOR_SAVED_PHOTO_ATTEMPTS=1182426
PRIOR_STRICT_ROI9=452097
HUES={"white","yellow_orange","red_pink","blue_purple"}


def observed_data(source_dir,actual_dir):
    source_dir,actual_dir=map(Path,(source_dir,actual_dir))
    m=json.loads((source_dir/"manifest.json").read_text())
    require(m.get("schema")=="mois_full_source_observation_sampling_cohort_v1"
            and m.get("source_arm")=="MOIS_SPECIES_FIRST_THIRD_NEW_ZERO9_1024"
            and m.get("new_photo_attempts_planned")==1024
            and m.get("work_units")==4
            and m.get("source_previous_real1280_original_attempts")==PRIOR_SAVED_PHOTO_ATTEMPTS,
            "Third-wave source manifest does not match signed second-wave cumulative")
    source=[];pixel=[];signatures={}
    for n in range(4):
        p=source_dir/"shards"/f"shard-{n:03d}.csv"
        terminal=actual_dir/f"shard_{n:03d}"/"photos.csv"
        receipt=actual_dir/f"shard_{n:03d}"/"receipt.json"
        require(sha256(p)==m["photo_shard_sha256"][str(n)],"Signed new source 256-photo shard changed")
        r=json.loads(receipt.read_text())
        require(r.get("schema")=="mois_full_original93m_real_RGB_photo_unit_v1"
                and r.get("attempted")==256 and r.get("shard")==n
                and r.get("photos_sha256")==sha256(terminal)
                and r.get("source_sha256")==sha256(p)
                and r.get("photo_colours_measured") is True,
                "Genuine original ROI-v4 raw photos/receipt SHA changed")
        a=pd.read_csv(p,dtype=str,keep_default_na=False)
        b=pd.read_csv(terminal,dtype=str,keep_default_na=False)
        require(len(a)==len(b)==256 and a.photo_id.is_unique and b.photo_id.is_unique
                and set(a.photo_id)==set(b.photo_id),
                "Third-wave source/photo terminal incomplete or duplicated")
        for col in ("inat_taxon_id","observation_id"):
            aa=a.set_index("photo_id")[col].sort_index()
            bb=b.set_index("photo_id")[col].sort_index()
            require(aa.equals(bb),"Third-wave original observation/taxon mismatch")
        source.append(a);pixel.append(b);signatures[str(n)]=sha256(receipt)
    observed=pd.concat(pixel,ignore_index=True)
    chosen=pd.concat(source,ignore_index=True)
    require(len(observed)==1024 and observed.photo_id.is_unique and
            observed.observation_id.is_unique and observed.inat_taxon_id.is_unique and
            chosen.photo_id.is_unique and chosen.inat_taxon_id.is_unique,
            "Wave3 photos/observation/species not globally independent")
    cols=["flower_fraction_"+x for x in PALETTE]
    fractions=observed[cols].apply(pd.to_numeric,errors="coerce").to_numpy(float)
    measured=(observed.morph.isin(HUES).to_numpy(bool)&
              np.isfinite(fractions).all(axis=1)&(fractions>=0).all(axis=1)&
              (fractions<=1.00001).all(axis=1)&
              (np.abs(np.nan_to_num(fractions).sum(axis=1)-1)<=0.001))
    labels=np.full(len(observed),"",dtype=object)
    labels[measured]=class_labels(fractions[measured])
    frame=observed[["photo_id","observation_id","inat_taxon_id"]].copy()
    frame["strict9"]=measured
    frame["dominant_label9"]=labels
    for i,k in enumerate(COLUMNS):
        frame[k]=np.where(measured,fractions[:,i],np.nan)
    return frame,chosen,signatures


def add_to_catalogue(atlas,third,*,expected_photos=1024):
    prior=atlas.copy()
    require(len(prior)==42111 and prior.inat_taxon_id.is_unique,
            "Full source taxon ledger absent")
    new=third.copy()
    for key in ("inat_taxon_id","photo_id","observation_id"):
        new[key]=new[key].fillna("").astype(str)
    require(len(new)==expected_photos and new.photo_id.is_unique and new.inat_taxon_id.is_unique,
            "Third-wave authentic photo selection changed")
    prior=prior.set_index("inat_taxon_id")
    require(new.inat_taxon_id.isin(prior.index).all(),
            "Third-wave source taxon outside prior 42111 census")
    statuses=prior.loc[new.inat_taxon_id]
    require(statuses.after_strict9.eq(0).all() and
            statuses.hue_status.eq("no_strict9").all(),
            "Source taxa already had a strict9 colour before this wave")
    prior.loc[new.inat_taxon_id,"after_photo_attempts"] = (
        prior.loc[new.inat_taxon_id,"after_photo_attempts"].astype(int).to_numpy()+1)
    positive=new.loc[new.strict9].copy()
    if not positive.empty:
        prior.loc[positive.inat_taxon_id,"after_strict9"]=1
        prior.loc[positive.inat_taxon_id,"independent_roi9_observations"]=1
        positive["hue_status"]=np.where(positive.dominant_label9.ne(""),
                                      "provisional_single","mixed_image_unresolved")
        prior.loc[positive.inat_taxon_id,"hue_status"]=positive.hue_status.to_numpy()
        prior.loc[positive.inat_taxon_id,"dominant_hue9"]=positive.dominant_label9.to_numpy()
        prior.loc[positive.inat_taxon_id,"observed_hue9_set"]=positive.dominant_label9.to_numpy()
        prior.loc[positive.inat_taxon_id,"distinct_palette_votes"]=(
            positive.dominant_label9.ne("").astype(int).to_numpy())
        prior.loc[positive.inat_taxon_id,"ambiguous_observation_rois"]=(
            positive.dominant_label9.eq("").astype(int).to_numpy())
        for c in PALETTE:
            prior.loc[positive.inat_taxon_id,"observed_mean_"+c]=(
                positive["strict9_"+c].to_numpy())
        prior.loc[positive.inat_taxon_id,"next_measure_priority"]=(
            np.where(positive.dominant_label9.ne(""),2,1))
    prior=prior.reset_index()
    prior["routine_hue_assigned"]=prior.dominant_hue9.fillna("").astype(str).ne("")
    require(len(prior)==42111 and prior.inat_taxon_id.is_unique,
            "Source taxon conservation lost in compact append")
    return prior


def run(prior_atlas,source_dir,actual_dir,global_options,
        prior_photos,first_delta,second_delta,output,*,quota=1024):
    output=Path(output)
    require(not output.exists(),"Refuse overwrite compact genuine photo atlas update")
    require(sha256(prior_atlas)==PRIOR_ATLAS_SHA,
            "Source signature for 42,111-taxon previous hue atlas drifted")
    prior=pd.read_csv(prior_atlas,compression="infer",dtype={"inat_taxon_id":str})
    frame,source,proof=observed_data(source_dir,actual_dir)
    updated=add_to_catalogue(prior,frame)
    # Verify no previously attempted image or observation was accidentally fetched.
    prev=pd.read_parquet(prior_photos,columns=["photo_id","observation_id","inat_taxon_id"])
    d1=pd.read_parquet(first_delta,columns=["photo_id","observation_id","inat_taxon_id"])
    d2=pd.read_parquet(second_delta,columns=["photo_id","observation_id","inat_taxon_id"])
    previous=pd.concat([prev,d1,d2],ignore_index=True)
    require(len(previous)==PRIOR_SAVED_PHOTO_ATTEMPTS and previous.photo_id.is_unique,
            "Old authentic image ledger not uniquely represented")
    require(not set(previous.photo_id.astype(str)) & set(frame.photo_id)
            and not (set(previous.observation_id.fillna("").astype(str))-{""})
                    & set(frame.observation_id),
            "Third-wave photos/observations overlap previously acquired originals")
    all_ids=pd.concat([previous,frame[["photo_id","observation_id","inat_taxon_id"]]],
                      ignore_index=True)
    options=pd.read_csv(global_options,compression="infer",dtype=str)
    queue=choose_third_wave(updated,all_ids,options,quota=quota)
    queue["measurement_id"]="MOIS-ZERO9-W4-"+queue.photo_id
    queue["image_filename"]=queue.measurement_id+".jpg"
    output.mkdir(parents=True)
    atlas_path=output/"mois_species_first_42111_after2304.csv.gz"
    delta_path=output/"mois_wave3_actual1024_strict9_photo_delta.parquet"
    queue_path=output/"mois_wave4_fresh_original1024_photo_queue.csv"
    updated.to_csv(atlas_path,index=False,compression="gzip")
    frame.to_parquet(delta_path,index=False)
    queue.to_csv(queue_path,index=False)
    gbif=updated[["inat_taxon_id","species","hue_status","dominant_hue9",
                  "observed_hue9_set"]].copy()
    gbif.to_csv(output/"mois_species_hue9_GBIF_source_after2304.csv.gz",
                index=False,compression="gzip")
    summary={"schema":SCHEMA,
            "previous_9ROI_1280_source_attempts":PRIOR_SAVED_PHOTO_ATTEMPTS,
            "new_original_image_attempts":len(frame),
            "new_strict_numeric_nineROI":int(frame.strict9.sum()),
            "new_clear_dominant_palette_hues":int(frame.dominant_label9.ne("").sum()),
            "new_technical_unknown":int((~frame.strict9).sum()),
            "cumulative_photo_attempts":PRIOR_SAVED_PHOTO_ATTEMPTS+1024,
            "cumulative_strict9_ROIs":PRIOR_STRICT_ROI9+int(frame.strict9.sum()),
            "source_nominal_taxa_with_any_strict9":int(updated.after_strict9.gt(0).sum()),
            "source_nominal_taxa_with_clear_hue":int(updated.dominant_hue9.fillna("").ne("").sum()),
            "source_nominal_taxa_without_strict9":int(updated.after_strict9.eq(0).sum()),
            "next_wave4_independent_new_source_photos":len(queue),
            "next_wave4_first_unmeasured_taxa":int(queue.acquisition_priority.eq("first_unmeasured").sum()),
            "next_wave4_retry":int(queue.acquisition_priority.eq("retry_unclear").sum()),
            "four_original_worker_receipts_sha256":proof,
            "new_source_photo_observation_overlap_with_previous":False,
            "GBIF_occurrences_collected_by_this_job":0,
            "species_hue_fixed_globally_proven":False,
            "output_sha256":{p.name:sha256(p) for p in
                              (atlas_path,delta_path,queue_path)}}
    (output/"receipt.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    print("MOIS_THIRD_WAVE_COMPACT_9PALETTE_CUMULATIVE_AND_WAVE4",
          json.dumps(summary,sort_keys=True),flush=True)
    return summary


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ("prior-atlas","source","actual","global-options","old-photos",
                "first-delta","second-delta","output"):
        p.add_argument("--"+arg,required=True)
    p.add_argument("--quota",type=int,default=1024)
    a=p.parse_args()
    run(a.prior_atlas,a.source,a.actual,a.global_options,a.old_photos,
        a.first_delta,a.second_delta,a.output,quota=a.quota)


if __name__=="__main__":
    main()
