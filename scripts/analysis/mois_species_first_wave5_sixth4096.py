#!/usr/bin/env python3
"""After genuine 5th ROI9 originals: increment 42,111 MOIS taxa, freeze 4096 fresh species.

The user goal is TAXON coverage, not 93,209,351 photograph enumeration.
All earlier original photo/observation IDs remain immutable; select first
new species photos in source metadata without looking at their hue.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
from scripts.analysis.mois_species_first_palette9 import PALETTE,COLUMNS,class_labels,sha256
from scripts.analysis.mois_species_first_incremental_first256 import require
from scripts.analysis.mois_species_first_real1280_cumulative import choose_third_wave
from scripts.analysis.mois_species_first_wave3_incremental import add_to_catalogue

SCHEMA="mois_species_first_fifth1024_real_roi9_append_sixth4096_v1"
PHOTO_ATLAS_SHA="4e22f1f7396c46b42a54f2095e8ad97da219f4d536abd96baf3e14e94fb9516c"
PHOTO_SHAS={
"baseline":"3fbd0f0ec7eb347b3bbc7ab509ec9312ec6d797ed932dcf7bb425959898c3835",
"first":"4eb884b12322573143277858f8f774def3eeecc7c4ecf212d8d2a67a93d0199a",
"second":"ee87a9046f630fc6941c49e01b435762a0785095b0fd1a1a60f943bdef8d570a",
"third":"137ae9cc657f27dbf509f7d97f581efd1f9a355814f8bce0109d9630bf5309f7",
"fourth":"e54c4990c6c1f4720e759c32c96e740d2b0023cedc700967387a224cb7a5fb60",
}
PRIOR_PHOTOS=1184474
PRIOR_NINE=452836
PRIOR_SUPPORTED=31235
PRIOR_LABELLED=28881
SOURCE_ARM="MOIS_SPECIES_FIRST_FIFTH_NEW_ZERO9_1024"
HUES={"white","yellow_orange","red_pink","blue_purple"}


def read_fifth(source,workers):
    source,workers=map(Path,(source,workers))
    manifest=json.loads((source/"manifest.json").read_text())
    require(manifest.get("schema")=="mois_full_source_observation_sampling_cohort_v1" and
            manifest.get("source_arm")==SOURCE_ARM and
            manifest.get("new_photo_attempts_planned")==1024 and
            manifest.get("source_previous_fourth_photo_attempts_plus_baseline")==PRIOR_PHOTOS and
            manifest.get("source_prior_no_reused_original_photo_or_observation") is True and
            manifest.get("worker_ids")==[0,1,2,3],
            "Wrong frozen 5th-wave new photo source")
    observed=[];source_rows=[];receipts={}
    for i in range(4):
        p=source/"shards"/f"shard-{i:03d}.csv"
        photo=workers/f"shard_{i:03d}"/"photos.csv"
        receipt=workers/f"shard_{i:03d}"/"receipt.json"
        require(p.is_file() and photo.is_file() and receipt.is_file(),
                "Missing genuine 5th-wave photo/receipt from FCP")
        require(sha256(p)==manifest["photo_shard_sha256"][str(i)],
                "Fifth source ID inventory SHA drift")
        r=json.loads(receipt.read_text())
        require(r.get("schema")=="mois_full_original93m_real_RGB_photo_unit_v1" and
                r.get("shard")==i and r.get("attempted")==256 and
                r.get("source_sha256")==sha256(p) and
                r.get("photos_sha256")==sha256(photo) and
                r.get("photo_colours_measured") is True,
                "Fifth genuine original photograph terminal/receipt mismatch")
        d=pd.read_csv(photo,dtype=str,keep_default_na=False)
        a=pd.read_csv(p,dtype=str,keep_default_na=False)
        require(len(d)==len(a)==256 and d.photo_id.is_unique and a.photo_id.is_unique and
                set(d.photo_id)==set(a.photo_id),
                "Source/original measured fifth photo ID contradiction")
        for col in ("inat_taxon_id","observation_id"):
            require(a.set_index("photo_id").sort_index()[col].equals(
                    d.set_index("photo_id").sort_index()[col]),
                    "Fifth original observation or source taxon changed")
        observed.append(d);source_rows.append(a)
        receipts[str(i)]=sha256(receipt)
    data=pd.concat(observed,ignore_index=True)
    src=pd.concat(source_rows,ignore_index=True)
    require(len(data)==1024 and
            len(data)==data.photo_id.nunique()==data.observation_id.nunique()==data.inat_taxon_id.nunique() and
            len(src)==src.photo_id.nunique()==src.observation_id.nunique(),
            "Fifth source photographs not one distinct photo per nominal species/observation")
    pixel_cols=["flower_fraction_"+p for p in PALETTE]
    numbers=data[pixel_cols].apply(pd.to_numeric,errors="coerce").to_numpy(float)
    valid=(data.morph.isin(HUES).to_numpy(bool) & np.isfinite(numbers).all(axis=1)&
           (numbers>=0).all(axis=1)&(numbers<=1.00001).all(axis=1)&
           (np.abs(np.nan_to_num(numbers).sum(axis=1)-1)<=0.001))
    label=np.full(1024,"",dtype=object)
    label[valid]=class_labels(numbers[valid])
    frame=data[["photo_id","observation_id","inat_taxon_id"]].copy()
    frame["strict9"]=valid;frame["dominant_label9"]=label
    for i,c in enumerate(COLUMNS):
        frame[c]=np.where(valid,numbers[:,i],np.nan)
    return frame,receipts


def append_fifth_and_choose_sixth(atlas,fifth,prior,options,*,quota=4096,
                                 expected_previous=PRIOR_PHOTOS):
    """Explicit old photo-ID and observation-based exclusion, preserve unknowns."""
    require(len(atlas)==42111 and atlas.inat_taxon_id.is_unique,
            "Genuine old 42111 photographic source-taxon atlas incomplete")
    require(int(atlas.after_strict9.gt(0).sum())==PRIOR_SUPPORTED and
            int(atlas.dominant_hue9.fillna("").ne("").sum())==PRIOR_LABELLED,
            "Signed previous numeric nine-colour species breadth changed")
    require(len(fifth)==1024 and fifth.photo_id.is_unique and fifth.inat_taxon_id.is_unique and
            fifth.observation_id.is_unique,
            "Fifth original photo/observation source duplicated")
    require(len(prior)==expected_previous and prior.photo_id.is_unique,
            "Old original photo-ID catalogue incomplete or duplicated")
    require(not set(prior.photo_id.astype(str)).intersection(set(fifth.photo_id.astype(str))),
            "A fifth-wave exact original photo ID was already previously attempted")
    old_obs=set(prior.observation_id.fillna("").astype(str))-{""}
    require(not old_obs.intersection(set(fifth.observation_id.astype(str))),
            "Fifth original source observation was already photographed")
    updated=add_to_catalogue(atlas,fifth,expected_photos=1024)
    ids=pd.concat([prior[["photo_id","observation_id","inat_taxon_id"]],
                   fifth[["photo_id","observation_id","inat_taxon_id"]]],ignore_index=True)
    queue=choose_third_wave(updated,ids,options,quota=quota)
    require(len(queue)==quota and queue.photo_id.is_unique and queue.observation_id.is_unique and
            queue.inat_taxon_id.is_unique and queue.acquisition_priority.eq("first_unmeasured").all(),
            "Sixth source must target fresh species with one independent photograph")
    queue["measurement_id"]="MOIS-SPECIES-FIRST-W6-"+queue.photo_id
    queue["image_filename"]=queue.measurement_id+".jpg"
    queue["shard"]=(np.arange(len(queue))//256).astype(str)
    return updated,queue


def run(atlas,source,workers,previous,options,out,*,quota=4096,verify_sha=True):
    out=Path(out)
    require(not out.exists(),"Refuse overwriting a previously signed fifth original ROI receipt")
    if verify_sha:
        require(sha256(atlas)==PHOTO_ATLAS_SHA,
                "Previous genuine fourth-wave 42,111-taxon atlas SHA drifted")
        for key,path in previous.items():
            require(sha256(path)==PHOTO_SHAS[key],
                    "Previously genuine photo IDs source byte SHA drift: "+key)
    d=pd.read_csv(atlas,compression="infer",dtype={"inat_taxon_id":str})
    delta,proof=read_fifth(source,workers)
    ids=pd.concat([pd.read_parquet(path,
               columns=["photo_id","observation_id","inat_taxon_id"])
                   for path in previous.values()],ignore_index=True)
    choices=pd.read_csv(options,compression="infer",dtype=str)
    updated,next_queue=append_fifth_and_choose_sixth(d,delta,ids,choices,quota=quota)
    success=int(delta.strict9.sum())
    require(int(updated.after_strict9.gt(0).sum())==PRIOR_SUPPORTED+success,
            "True fifth numeric9 success did not add same number of first-covered taxa")
    out.mkdir(parents=True)
    paths={
        "atlas":out/"mois_species_first_42111_after4352.csv.gz",
        "delta":out/"mois_fifth1024_real_original9_photo_delta.parquet",
        "next":out/"mois_sixth4096_first_new_species_photo_queue.csv",
    }
    updated.to_csv(paths["atlas"],index=False,compression="gzip")
    delta.to_parquet(paths["delta"],index=False)
    next_queue.to_csv(paths["next"],index=False)
    receipt={
        "schema":SCHEMA,"signed_previous_original_photo_attempts":PRIOR_PHOTOS,
        "previous_original_nineROI":PRIOR_NINE,
        "real_fifth_new_photo_attempts":len(delta),
        "real_fifth_new_numerical_nineROI":success,
        "real_fifth_dominant_palette_photo_labels":int(delta.dominant_label9.ne("").sum()),
        "real_fifth_technical_unknown":int((~delta.strict9).sum()),
        "cumulative_original_photo_attempts":PRIOR_PHOTOS+1024,
        "cumulative_genuine_nineROI":PRIOR_NINE+success,
        "source_taxa_total":len(updated),
        "source_taxa_with_numerical_nineROI":int(updated.after_strict9.gt(0).sum()),
        "source_taxa_with_dominant_nine_palette":int(updated.dominant_hue9.fillna("").ne("").sum()),
        "source_taxa_without_nineROI":int(updated.after_strict9.eq(0).sum()),
        "next_source_distinct_first_photo_taxa":len(next_queue),
        "next_source_retry_taxa":0,
        "fifth_original_worker_receipts_SHA256":proof,
        "prior_photos_and_observations_reused":False,
        "GBIF_species_occurrences_added":0,
        "species_wide_flower_colour_fixed_proven":False,
        "output_sha256":{key:sha256(path) for key,path in paths.items()},
    }
    (out/"receipt.json").write_text(json.dumps(receipt,sort_keys=True,indent=2)+"\n")
    print("MOIS_REAL_FIFTH_NINEHUE_APPEND_SIXTH4096_SOURCE",
          json.dumps(receipt,sort_keys=True),flush=True)
    return receipt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("atlas","source","workers","baseline","first","second","third",
              "fourth","global-options","output"):
        p.add_argument("--"+k,required=True)
    p.add_argument("--quota",type=int,default=4096)
    a=p.parse_args()
    run(a.atlas,a.source,a.workers,dict(baseline=a.baseline,first=a.first,
       second=a.second,third=a.third,fourth=a.fourth),a.global_options,
       a.output,quota=a.quota)


if __name__=="__main__":
    main()
