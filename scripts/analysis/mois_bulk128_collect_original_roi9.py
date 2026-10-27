#!/usr/bin/env python3
"""Data-only collector for 128 completed original photo ROI-v4 shards.

Never redownload a source photo or re-run YOLO11/EfficientSAM. Read original
numeric terminals and source SHA receipts; retain all unknown/failed images.
Reference and gap-targeted source arms remain separately identifiable.
"""
from __future__ import annotations
import argparse,collections,hashlib,json,math
from pathlib import Path
import numpy as np
import pandas as pd

from scripts.acquisition.mois_completed_pixel_checkpoint import collect_wave
from scripts.analysis.mois_reassemble_original_new4096_measurements import METHOD
from scripts.acquisition.mois_bulk128_tier10_source import (
    SCHEMA as SOURCE_SCHEMA,REFERENCE_ARM,GAP_ARM,NO_COLOUR_PRIORITY,
    THIN_PRIORITY,P5,
)

SCHEMA="mois_bulk128_real_original_photo_ROI9_complete_signed_v1"
PALETTE=("white","yellow","orange","red","pink","magenta","purple","blue","bronze")
FRACTIONS=tuple("flower_fraction_"+name for name in PALETTE)
HUES=("white","yellow_orange","red_pink","blue_purple")

def require(ok,msg):
    if not bool(ok):
        raise ValueError("MOIS_BULK128_COLLECT: "+msg)

def sha(file):
    h=hashlib.sha256()
    with Path(file).open("rb") as f:
        for chunk in iter(lambda:f.read(4<<20),b""):h.update(chunk)
    return h.hexdigest()

def validate_join(rows,root,manifest):
    d=pd.DataFrame(rows)
    require(len(d)==manifest["new_photo_attempts_planned"],
            "Real acquired source photos do not exhaust frozen 128-shard roster")
    for col in ("photo_id","observation_id","inat_taxon_id"):
        d[col]=d[col].fillna("").astype(str)
    require(d.photo_id.is_unique and d.observation_id.is_unique
            and not (d.photo_id=="").any() and not (d.observation_id=="").any(),
            "Original photo/observation identity repeated in ROI terminals")
    source=[]
    for part in range(128):
        f=Path(root)/"provenance"/f"part-{part:03d}.csv"
        r=json.loads((Path(root)/f"part-{part:03d}.receipt.json").read_text())
        require(sha(f)==r["source_photo_provenance_SHA256"]
                and sha(Path(root)/f"part-{part:03d}.receipt.json")==
                    manifest["per_128_shard_receipts_sha256"][str(part)],
                "Original per-shard source photo/taxon provenance SHA changed")
        source.append(pd.read_csv(
            f,dtype={"photo_id":str,"observation_id":str,"inat_taxon_id":str},
            low_memory=False))
    p=pd.concat(source,ignore_index=True)
    require(len(p)==len(d) and p.photo_id.is_unique
            and p.observation_id.astype(str).is_unique,
            "Lost or duplicated photos across independently original 128 source parts")
    p=p.rename(columns={
        "observation_id":"selected_observation_id",
        "inat_taxon_id":"selected_taxon_id",
    })
    joined=d.merge(p,on="photo_id",how="left",validate="one_to_one",sort=False)
    require(joined.selected_taxon_id.notna().all()
            and joined.inat_taxon_id.eq(joined.selected_taxon_id.astype(str)).all()
            and joined.observation_id.eq(joined.selected_observation_id.astype(str)).all()
            and ((joined.source_arm==GAP_ARM)==
                 joined.acquisition_design_cohort.eq("targeted")).all()
            and joined.source_arm.isin((REFERENCE_ARM,GAP_ARM)).all()
            and joined.original_source_priority.notna().all()
            and joined.loc[joined.acquisition_design_cohort.eq("targeted"),
                           "original_source_priority"].ne(P5).all(),
            "Frozen preimage sample design disagrees with real photo terminal IDs")
    return joined

def classify_nine(df):
    require(set(FRACTIONS)<=set(df),
            "Nine original numeric FCP biological palette fractions missing")
    values=df[list(FRACTIONS)].apply(pd.to_numeric,errors="coerce").to_numpy(float)
    finite=np.isfinite(values).all(axis=1)
    nonnegative=(np.where(np.isfinite(values),values,-1)>=0).all(axis=1)
    bounded=(np.where(np.isfinite(values),values,0)<=1+1e-6).all(axis=1)
    mass=np.where(np.isfinite(values),values,0).sum(axis=1)
    palette=finite&nonnegative&bounded&(mass>0)
    known=df.morph.isin(HUES).to_numpy(bool)
    failure=df.measurement_status.eq("image_acquisition_failed").to_numpy(bool)
    strict=palette&known
    require(not (strict&failure).any()
            and df.morph.isin((*HUES,"mixed_uncertain")).all()
            and known.sum()+sum(~known)==len(df),
            "Missing/technical original ROI wrongly relabelled as a flower hue")
    out=df.copy()
    out["strict9"]=strict
    out["original_four_hue_classified"]=known
    out["image_acquisition_failed"]=failure
    out["original_nine_component_available_but_not_strict"]=(palette&~known)
    normalized=np.full_like(values,np.nan)
    normalized[palette,:]=values[palette,:]/mass[palette,None]
    for i,name in enumerate(PALETTE):
        out["strict9_"+name]=np.where(strict,normalized[:,i],np.nan)
    arr=np.where(palette,normalized,0)
    entropy=-np.sum(np.where(arr>0,arr*np.log(np.maximum(arr,1e-12)),0),
                    axis=1)/math.log(len(PALETTE))
    out["strict9_palette_entropy"]=np.where(strict,entropy,np.nan)
    require(out.loc[~out.strict9,["strict9_"+x for x in PALETTE]].isna().all().all(),
            "Technical unknown pixels silently promoted to botanical trait")
    return out

def summary(df,group_cols):
    rows=[]
    for keys,g in df.groupby(group_cols,observed=True,dropna=False,sort=True):
        if not isinstance(keys,tuple):keys=(keys,)
        good=g[g.strict9]
        row={c:str(v) for c,v in zip(group_cols,keys)}
        row.update({
           "attempted_original_photos":len(g),
           "original_distinct_observations":int(g.observation_id.nunique()),
           "distinct_nominal_taxa_attempted":int(g.inat_taxon_id.nunique()),
           "strict9_ROI_colours":len(good),
           "automatic_four_hue_classified":int(g.original_four_hue_classified.sum()),
           "photo_acquisition_failed":int(g.image_acquisition_failed.sum()),
           "technical_or_ambiguous_unknown":int((~g.original_four_hue_classified).sum()),
           "strict9_newly_classified_zero_colour_taxa":int(
                good.loc[good.original_source_priority.isin(NO_COLOUR_PRIORITY),
                         "inat_taxon_id"].nunique()),
           "strict9_success_fraction":len(good)/len(g),
           "conditional_photo_white_fraction":(
               float(good.strict9_white.mean()) if len(good) else None),
        })
        rows.append(row)
    return rows

def collect(source,workers,output):
    source,workers,output=map(Path,(source,workers,output))
    require(not output.exists(),"Do not replace immutable measured 128-shard photo evidence")
    m=json.loads((source/"manifest.json").read_text())
    require(m.get("schema")=="mois_full_source_observation_sampling_cohort_v1"
            and m.get("sampling_schema")==SOURCE_SCHEMA
            and m.get("source_128_compute_shards")==128
            and m.get("work_units")==512
            and m.get("group_ids")==list(range(128))
            and m.get("P5_targeted_photo_attempts")==0
            and m.get("original_source_photo_IDs_and_observation_UUIDs_all_unique"),
            "Not the frozen authentic bulk128 original photo acquisition roster")
    completed,receipt=collect_wave(source,workers,METHOD)
    require(receipt["wave_complete"] and receipt["pending_photo_attempts"]==0
            and receipt["planned_photo_attempts"]==receipt["attempted"]
            ==m["new_photo_attempts_planned"]
            and len(receipt["completed_units"])==512,
            "Full 128 photo shards are not yet all genuinely processed")
    d=classify_nine(validate_join(completed,source,m))
    require(d.source_arm.isin((REFERENCE_ARM,GAP_ARM)).all()
            and len(d.loc[d.acquisition_design_cohort.eq("reference")])==
                     m["selected_reference_photo_attempts"]
            and len(d.loc[d.acquisition_design_cohort.eq("targeted")])==
                     m["selected_targeted_photo_attempts"]
            and int(d.image_acquisition_failed.sum())==receipt["acquisition_failed"]
            and int(d.original_four_hue_classified.sum())==receipt["classifiable"],
            "Original real-photo classifications and source design mismatch")
    total_good=int(d.strict9.sum())
    newsource=(
        d.loc[d.strict9&d.original_source_priority.isin(NO_COLOUR_PRIORITY),
              "inat_taxon_id"].nunique()
    )
    output.mkdir(parents=True)
    p=output/"bulk128_actual_original_photo_ROI9.parquet"
    d.to_parquet(p,index=False,compression="zstd")
    summaries={}
    for name,cols in (
        ("design",["acquisition_design_cohort"]),
        ("source_priority",["acquisition_design_cohort","original_source_priority"]),
        ("region_season",["acquisition_design_cohort","latitude_region",
                          "recorded_local_season"])):
        table=output/f"bulk128_actual_photo_ROI9_by_{name}.csv"
        pd.DataFrame(summary(d,cols)).to_csv(table,index=False)
        summaries[name]=sha(table)
    report={
       "schema":SCHEMA,
       "full_upstream_original_photo_IDs":93_209_351,
       "number_source_compute_shards":128,
       "number_completed_original_photo_model_units":512,
       "attempted_new_original_photos":len(d),
       "distinct_original_observations":int(d.observation_id.nunique()),
       "new_strict_original_nine_palette_ROIs":total_good,
       "new_original_four_hue_classifications":int(d.original_four_hue_classified.sum()),
       "new_original_images_acquisition_failed":int(d.image_acquisition_failed.sum()),
       "new_original_four_hue_unknown_images":int((~d.original_four_hue_classified).sum()),
       "newly_covered_source_taxa_with_previous_zero_strict9":int(newsource),
       "targeted_previous_P5_source_taxa":0,
       "reference_photo_attempts":m["selected_reference_photo_attempts"],
       "gap_targeted_photo_attempts":m["selected_targeted_photo_attempts"],
       "source_manifest_SHA256":sha(source/"manifest.json"),
       "original_ROI_photo_outcomes_Parquet_SHA256":sha(p),
       "summary_SHA256":summaries,
       "completed_real_photo_terminal_SHA_audit":receipt,
       "original_photos_downloaded_OR_colour_inferred_in_THIS_collector":0,
       "targeted_source_selection_not_unbiased_global_flora_sample":True,
       "botanically_verified_open_flower_wild_native_and_petal_status":False,
       "technical_unknown_preserved":True,
    }
    (output/"bulk128_actual_original_photo_ROI9_receipt.json").write_text(
         json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("MOIS_BULK128_ACTUAL_ORIGINAL_IMAGE_PHOTO_ROI9_SOURCE_AUDITED",
          json.dumps({k:report[k] for k in (
              "attempted_new_original_photos","new_strict_original_nine_palette_ROIs",
              "new_original_four_hue_classifications",
              "new_original_images_acquisition_failed",
              "newly_covered_source_taxa_with_previous_zero_strict9",
              "reference_photo_attempts","gap_targeted_photo_attempts",
              "targeted_previous_P5_source_taxa"
           )},sort_keys=True),flush=True)
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for x in ("source","workers","output"):
        p.add_argument("--"+x,required=True,type=Path)
    q=p.parse_args()
    collect(q.source,q.workers,q.output)

if __name__=="__main__":
    main()
