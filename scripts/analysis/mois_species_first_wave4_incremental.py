#!/usr/bin/env python3
"""Append fourth-wave original-photo nine-colour ROI to MOIS 42,111-taxon atlas.

No model rerun of historical images, and no species-wide biological fixation
claim. Uses SHA-signed original worker output, excludes all old photo and
observation IDs, and selects fifth-wave first-photo targets by taxon breadth.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis.mois_species_first_palette9 import PALETTE, COLUMNS, class_labels, sha256
from scripts.analysis.mois_species_first_incremental_first256 import require
from scripts.analysis.mois_species_first_real1280_cumulative import choose_third_wave
from scripts.analysis.mois_species_first_wave3_incremental import add_to_catalogue

SCHEMA = "mois_species_first_fourth1024_original_roi9_incremental_v1"
SOURCE_ARM = "MOIS_SPECIES_FIRST_FOURTH_NEW_ZERO9_1024"
PRIOR_ATLAS_SHA = "4369818640dd557530274077610bd91637e8e86f3a6094dbeb12b5fd857aa697"
PRIOR_NUMERIC = 452458
PRIOR_ATTEMPTED = 1183450
PRIOR_STRICT_TAXA = 30857
SOURCE_BASE_SHA = "3fbd0f0ec7eb347b3bbc7ab509ec9312ec6d797ed932dcf7bb425959898c3835"
FIRST_DELTA_SHA = "4eb884b12322573143277858f8f774def3eeecc7c4ecf212d8d2a67a93d0199a"
SECOND_DELTA_SHA = "ee87a9046f630fc6941c49e01b435762a0785095b0fd1a1a60f943bdef8d570a"
THIRD_DELTA_SHA = "137ae9cc657f27dbf509f7d97f581efd1f9a355814f8bce0109d9630bf5309f7"
HUES = {"white", "yellow_orange", "red_pink", "blue_purple"}

def authentic_wave4(source_root, actual_root):
    source_root, actual_root = Path(source_root),Path(actual_root)
    manifest_path = source_root/"manifest.json"
    manifest = json.loads(manifest_path.read_text())
    require(manifest.get("schema")=="mois_full_source_observation_sampling_cohort_v1" and
            manifest.get("source_arm")==SOURCE_ARM and
            manifest.get("new_photo_attempts_planned")==1024 and
            manifest.get("source_previous_real2304_original_attempts")==PRIOR_ATTEMPTED and
            manifest.get("source_previous_strict9_realROIs")==PRIOR_NUMERIC and
            manifest.get("one_photo_per_distinct_species_observation") is True and
            manifest.get("worker_ids")==[0,1,2,3],
            "Wrong frozen fourth-wave authentic photographic source manifest")
    source,actual=[],[]
    proof={}
    for i in range(4):
        source_csv=source_root/"shards"/f"shard-{i:03d}.csv"
        terminal=actual_root/f"shard_{i:03d}"/"photos.csv"
        receipt=actual_root/f"shard_{i:03d}"/"receipt.json"
        require(source_csv.is_file() and terminal.is_file() and receipt.is_file(),
                "Original measured wave4 complete 256-photo terminal missing")
        require(sha256(source_csv)==manifest["photo_shard_sha256"][str(i)],
                "Source original photo-ID work shard SHA mismatch")
        r=json.loads(receipt.read_text())
        require(r.get("schema")=="mois_full_original93m_real_RGB_photo_unit_v1" and
                r.get("attempted")==256 and r.get("shard")==i and
                r.get("photos_sha256")==sha256(terminal) and
                r.get("source_sha256")==sha256(source_csv) and
                r.get("photo_colours_measured") is True,
                "Fourth-wave original pixel numeric terminal receipt SHA mismatch")
        a=pd.read_csv(source_csv,dtype=str,keep_default_na=False)
        b=pd.read_csv(terminal,dtype=str,keep_default_na=False)
        require(len(a)==len(b)==256 and a.photo_id.is_unique and
                b.photo_id.is_unique and set(a.photo_id)==set(b.photo_id),
                "Photo IDs lost between frozen source and measured worker outcome")
        for column in ("inat_taxon_id","observation_id"):
            aa=a.set_index("photo_id").sort_index()[column]
            bb=b.set_index("photo_id").sort_index()[column]
            require(aa.equals(bb),"Original wave4 taxon/observation ID changed after measurement")
        source.append(a)
        actual.append(b)
        proof[str(i)]=sha256(receipt)
    planned=pd.concat(source,ignore_index=True)
    raw=pd.concat(actual,ignore_index=True)
    require(len(planned)==len(raw)==1024 and
            planned.inat_taxon_id.is_unique and planned.observation_id.is_unique and
            raw.inat_taxon_id.is_unique and raw.observation_id.is_unique and
            raw.photo_id.is_unique and raw.photo_id.str.fullmatch(r"\d+").all(),
            "Fourth wave selected duplicate photo, plant taxon or source observation")
    matrix=raw[["flower_fraction_"+k for k in PALETTE]].apply(
        pd.to_numeric,errors="coerce").to_numpy(float)
    measured=(raw.morph.isin(HUES).to_numpy(bool) &
              np.isfinite(matrix).all(axis=1) & (matrix>=0).all(axis=1) &
              (matrix<=1.00001).all(axis=1) &
              (abs(np.nan_to_num(matrix,nan=0).sum(axis=1)-1)<=0.001))
    label=np.full(1024,"",dtype=object)
    label[measured]=class_labels(matrix[measured])
    delta=raw[["photo_id","observation_id","inat_taxon_id"]].copy()
    delta["strict9"]=measured
    delta["dominant_label9"]=label
    for i,k in enumerate(COLUMNS):
        delta[k]=np.where(measured,matrix[:,i],np.nan)
    return delta,proof


def add_and_queue(prior_atlas,wave4,prior_photo_ids,candidates,*,quota=1024,
                  expected_previous=PRIOR_ATTEMPTED):
    """Core species-first ledger update; explicit source count is a test seam."""
    require(len(prior_atlas)==42111 and prior_atlas.inat_taxon_id.is_unique,
            "42,111-taxon previous nine-colour photographic atlas missing")
    require(len(prior_photo_ids)==expected_previous and prior_photo_ids.photo_id.is_unique,
            "Previous original photo IDs incomplete or repeated")
    require(not (set(prior_photo_ids.photo_id.astype(str)) &
                 set(wave4.photo_id.astype(str))),
            "Actual fourth-wave image is a previously measured original photo")
    prior_obs=set(prior_photo_ids.observation_id.fillna("").astype(str))-{""}
    require(not prior_obs.intersection(set(wave4.observation_id.astype(str))),
            "Actual fourth-wave observation UUID was already measured previously")
    atlas=add_to_catalogue(prior_atlas,wave4,expected_photos=len(wave4))
    joined_ids=pd.concat([prior_photo_ids[["photo_id","observation_id","inat_taxon_id"]],
                          wave4[["photo_id","observation_id","inat_taxon_id"]]],
                         ignore_index=True)
    queue=choose_third_wave(atlas,joined_ids,candidates,quota=quota)
    queue["measurement_id"]="MOIS-ZERO9-W5-"+queue.photo_id
    queue["image_filename"]=queue.measurement_id+".jpg"
    require(len(queue)==quota and queue.inat_taxon_id.is_unique and
            queue.photo_id.is_unique and queue.observation_id.is_unique,
            "Fifth-wave source target missing or duplicates previous IDs")
    return atlas,queue


def run(prior_atlas,source,actual,old_photos,first_delta,second_delta,
        third_delta,options,out,*,quota=1024,verify_sha=True):
    out=Path(out)
    require(not out.exists(),"Refuse to overwrite existing fourth-wave genuine data receipt")
    expected={
        "prior":(prior_atlas,PRIOR_ATLAS_SHA),
        "old":(old_photos,SOURCE_BASE_SHA),
        "first":(first_delta,FIRST_DELTA_SHA),
        "second":(second_delta,SECOND_DELTA_SHA),
        "third":(third_delta,THIRD_DELTA_SHA),
    }
    if verify_sha:
        for name,(path,digest) in expected.items():
            require(sha256(path)==digest,"Prior genuine numeric/photo source SHA changed: "+name)
    atlas_before=pd.read_csv(prior_atlas,compression="infer",dtype={"inat_taxon_id":str})
    delta,signatures=authentic_wave4(source,actual)
    original=[]
    for filename in (old_photos,first_delta,second_delta,third_delta):
        original.append(pd.read_parquet(filename,
                    columns=["photo_id","observation_id","inat_taxon_id"]))
    previous=pd.concat(original,ignore_index=True)
    opts=pd.read_csv(options,compression="infer",dtype=str)
    updated,next_wave=add_and_queue(atlas_before,delta,previous,opts,quota=quota)
    require(int(atlas_before.after_strict9.gt(0).sum())==PRIOR_STRICT_TAXA,
            "Prior real taxon-level strict9 support drifted")
    now_good=int(delta.strict9.sum())
    require(int(updated.after_strict9.gt(0).sum())==PRIOR_STRICT_TAXA+now_good,
            "New strict9 original photo success must be a genuinely new taxon")
    out.mkdir(parents=True)
    saved={
        "source":out/"mois_wave4_genuine_original1024_strict9_delta.parquet",
        "atlas":out/"mois_42111_species_first_after3328.csv.gz",
        "next":out/"mois_wave5_next1024_new_taxa_source.csv",
        "gbif":out/"mois_wave4_42111_species9_GBIF_incremental_name_source.csv.gz",
    }
    delta.to_parquet(saved["source"],index=False)
    updated.to_csv(saved["atlas"],index=False,compression="gzip")
    next_wave.to_csv(saved["next"],index=False)
    gbif=updated[["inat_taxon_id","species","hue_status","dominant_hue9",
                  "observed_hue9_set"]].copy()
    gbif.to_csv(saved["gbif"],index=False,compression="gzip")
    report={
        "schema":SCHEMA,
        "original_cumulative_photo_attempts":PRIOR_ATTEMPTED+1024,
        "original_cumulative_numeric_strict9":PRIOR_NUMERIC+now_good,
        "source_taxa_in_frame":len(updated),
        "taxa_with_any_numeric_strict9":int(updated.after_strict9.gt(0).sum()),
        "taxa_without_numeric_strict9":int(updated.after_strict9.eq(0).sum()),
        "taxa_with_sampled_dominant_nine_palette":int(updated.dominant_hue9.fillna("").ne("").sum()),
        "new_original_photo_attempts":1024,
        "new_taxa_with_first_numeric_strict9":now_good,
        "new_taxa_with_clear_dominant_palette":int(delta.dominant_label9.ne("").sum()),
        "new_technical_ROI_unknowns":int((~delta.strict9).sum()),
        "next_wave5_new_taxa":int(next_wave.acquisition_priority.eq("first_unmeasured").sum()),
        "next_wave5_retry_unclear":int(next_wave.acquisition_priority.eq("retry_unclear").sum()),
        "original_worker_receipt_SHA256":signatures,
        "source_photo_observation_overlap":0,
        "GBIF_occurrences_added_in_this_step":0,
        "species_wide_color_fixation_proven":False,
        "cumulative_photographic_atlas_source_SHA256":sha256(prior_atlas),
        "output_SHA256":{k:sha256(v) for k,v in saved.items()}}
    (out/"receipt.json").write_text(json.dumps(report,sort_keys=True,indent=2)+"\n")
    print("MOIS_REAL_FOURTH_WAVE_SPECIES9_APPEND",json.dumps(report,sort_keys=True),flush=True)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("prior-atlas","source","actual","old-photos","first-delta",
                 "second-delta","third-delta","options","output"):
        p.add_argument("--"+name,required=True)
    p.add_argument("--quota",type=int,default=1024)
    z=p.parse_args()
    run(z.prior_atlas,z.source,z.actual,z.old_photos,z.first_delta,
        z.second_delta,z.third_delta,z.options,z.output,quota=z.quota)


if __name__=="__main__":
    main()
