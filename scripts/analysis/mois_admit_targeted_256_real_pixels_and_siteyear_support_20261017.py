#!/usr/bin/env python3
"""Measure-then-admit ONLY first finished genuine 256 targeted source photo ROIs.

Scientific target is increased same-named-taxon/geographic-cell/season-year
photographic overlap, not global white-flower prevalence. This targeted source
CANNOT be pooled into genus-climate global moisture inference as unbiased data.
Every new RGB photo ROI is verified against exact original fourth-unmeasured
source reservation exclusion, the source-only selection SHA, original frozen
YOLO/EfficientSAM model attestation, 4x64 pixels/shard and immutable old309439
photo hue/date records. Unfinished 3840 candidates are NOT measured data.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import zipfile
from pathlib import Path
import numpy as np
import pandas as pd

from scripts.analysis.mois_reassemble_original_new4096_measurements import (
    BASE_COLS,HUES,METHOD,MEASURED_SCHEMA,
)
from scripts.analysis.mois_seasonal_267k_climate_normals import (
    labelled_calendar,require,
)
from scripts.analysis.mois_target_season_cell_original_photo_source_only_20261017 import (
    FREEZE_OUT,PHOTO_SRC,PHOTO_BASE_SHA,RESERVE_N,SHARD_SIZE,
    N_SHARDS,SELECT_ONLY_PIXEL_SHARD,sha,
)

ORIGINAL_SOURCE="mois_original309439_fullFCP256_actual_photo_hues_and_source_dates.csv.gz"
OUT_SOURCE="mois_targeted256_actual_ROI_309695_original_photo_and_dates.csv.gz"
OUT_RESULT="mois_targeted256_original_photo_repeated_season_siteyear_QC_receipt.json"
OUT_PHOTOS="mois_targeted256_original_photo_season_white_source.csv"
OUT_SUPPORT="mois_targeted256_original_photo_taxon_same025degree_year_summer_winter_support.csv"

def read_base(old_archive:Path):
    with zipfile.ZipFile(old_archive) as z:
        raw=z.read(ORIGINAL_SOURCE)
    require(hashlib.sha256(raw).hexdigest()==PHOTO_BASE_SHA,
            "Original 309439 photo/date source byte SHA drift")
    a=pd.read_csv(io.BytesIO(raw),compression="gzip",dtype=str,
                  keep_default_na=False,low_memory=False)
    require(len(a)==309439 and a.photo_id.is_unique and
            tuple(a.columns)==(*BASE_COLS,"qualified_date",
                               "original_source_date_provenance") and
            int(a.morph.isin(HUES).sum())==149662 and
            int((a.morph.isin(HUES)&a.qualified_date.str.len().eq(10)).sum())==149653,
            "Previous genuine original 309439 photo/date/four hues changed")
    return a

def read_plan(plan:Path):
    with zipfile.ZipFile(plan) as z:
        rec=json.loads(z.read(FREEZE_OUT))
        source=z.read(PHOTO_SRC)
        shardname=f"mois_source_only_shard_{SELECT_ONLY_PIXEL_SHARD:02d}_256.csv"
        one=z.read(shardname)
    require(rec["schema"]=="mois_3430502_targeted_opposite_season_after309439_and_reserved4096_v1"
            and rec["first_wave_original_source_photo_ids"]==RESERVE_N
            and rec["first_wave_original_shards"]==N_SHARDS
            and rec["photos_per_original_shard"]==SHARD_SIZE
            and rec["first_wave_original_4096_csv_sha256"]==
                hashlib.sha256(source).hexdigest()
            and rec["per_shard_sha256"][str(SELECT_ONLY_PIXEL_SHARD)]==
                hashlib.sha256(one).hexdigest()
            and rec["new_photo_pixels_actually_measured"] is False
            and rec["source_photo_pixels_stage_no_pixels_yet"] is True
            and rec["current_recruited_source_is_not_representative_global_photo_prevalence"] is True,
            "Original targeted photo-source selected IDs/provenance mutated")
    src=list(csv.DictReader(io.StringIO(one.decode())))
    require(len(src)==SHARD_SIZE and
            len({r["measurement_id"] for r in src})==SHARD_SIZE and
            all(int(r["shard"])==SELECT_ONLY_PIXEL_SHARD for r in src),
            "Actual first targeted source photo-ID shard not source-locked")
    return rec,src

def read_worker(root:Path,plan_rec:dict,src:list[dict]):
    k=SELECT_ONLY_PIXEL_SHARD
    one=list(root.rglob(f"mois_original_256_fresh_ROI_shard_{k:02d}.csv"))
    receipts=list(root.rglob(f"mois_original_256_ROI_shard_{k:02d}_receipt.json"))
    require(len(one)==len(receipts)==1,
            "At least one 256 original RGB photo worker terminal must REALLY finish")
    r=json.loads(receipts[0].read_text())
    require(r["schema"]==MEASURED_SCHEMA and int(r["shard"])==k and
            int(r["terminal_rows"])==int(r["original_source_ids"])==SHARD_SIZE and
            r["source_sha256"]==plan_rec["per_shard_sha256"][str(k)] and
            r["source_frozen_4096_selection_sha256"]==
               plan_rec["first_wave_original_4096_csv_sha256"] and
            r["source_model_attestation"]==METHOD and
            r["pixels_ephemeral_no_original_image_bytes_persisted"] is True and
            r["fcp_within_species_polymorphism_estimand_computed"] is False and
            r["original_new_photo_hue_terminal_SHA256"]==sha(one[0]),
            "Original photo RGB/pixel source/model SHA is not genuine")
    table=list(csv.DictReader(one[0].open(newline="",encoding="utf-8-sig")))
    source_id={v["measurement_id"]:v for v in src}
    seen=set();taken=[]
    for row in table:
        mid=row["measurement_id"]
        require(mid in source_id and mid not in seen,
                "Unselected or duplicate original photo outcome")
        seen.add(mid)
        old=source_id[mid]
        for name in ("photo_id","inat_taxon_id","species","observation_id","observer_id"):
            require(row[name]==old[name],
                    "Original source photo, named taxon or observation/observer identity changed")
        yes=row["morph"] in HUES
        require(row["morph"] in (*HUES,"mixed_uncertain")
                and (row["measurement_status"]=="classified_four_state_morph")==yes
                and row["classified_four_hue"]==("true" if yes else "false"),
                "Actual independent ROI-v4 photo outcome or technical-unknown state changed")
        vals={key:row.get(key,"") for key in BASE_COLS}
        for name in ("latitude","longitude","photo_license","observed_on"):
            vals[name]=old.get(name,"")
        vals["source_original_measured"]="mois_targeted_source_photo_season_actual_ROI_v4"
        vals["old_mois_or_new_fcp"]="mois_targeted_opposite_photo_season"
        vals["qualified_date"]=old["observed_on"]
        vals["original_source_date_provenance"]="20260927_FCP_numeric_source_photo_UUID_qualified_actual_day"
        taken.append(vals)
    require(len(seen)==SHARD_SIZE and
            int(r["four_hue_classified"])==sum(v["morph"] in HUES for v in taken) and
            int(r["unclassified"])==sum(v["morph"] not in HUES for v in taken),
            "Original terminal four-hue photograph counts do not reconcile")
    # Save records in the *pre-image source photo-ID order* rather than
    # worker/terminal ordering so metadata collection cannot reorder identity.
    identity={z["photo_id"]:z for z in taken}
    return [identity[z["photo_id"]] for z in src],r

def support(p:pd.DataFrame,geograin=.25):
    q=p[["photo_id","inat_taxon_id","observation_id","morph",
         "qualified_date","latitude","longitude"]].copy()
    q.rename(columns={"qualified_date":"observed_on"},inplace=True)
    for fld in ("latitude","longitude"):
        q[fld]=pd.to_numeric(q[fld],errors="coerce")
    seasonal=labelled_calendar(q)
    good=seasonal[seasonal.date.notna()&seasonal.classified&
          seasonal.zone.isin(("north_extratropical","south_extratropical"))&
          seasonal.season.isin(("summer","winter"))&
          seasonal.latitude.between(-90,90)&seasonal.longitude.between(-180,180)].copy()
    good["x"]=np.floor((good.longitude+180)/geograin).astype(int)
    good["y"]=np.floor((good.latitude+90)/geograin).astype(int)
    pair=good.groupby(["inat_taxon_id","zone","x","y","season_year","season"],
        observed=True).observation_id.nunique().unstack(fill_value=0)
    for season in ("summer","winter"):
        if season not in pair.columns:pair[season]=0
    summaries={}
    for hemi in ("north_extratropical","south_extratropical"):
        q=pair[pair.index.get_level_values("zone")==hemi]
        first=(q["summer"]>=1)&(q["winter"]>=1)
        strict=(q["summer"]>=2)&(q["winter"]>=2)
        summaries[hemi]={
          "same_original_taxon_025degree_cell_same_seasonyear_both_summer_winter_min1_obs":
             int(first.sum()),
          "same_original_taxon_025degree_cell_same_seasonyear_both_summer_winter_min2_obs":
             int(strict.sum())}
    return summaries,seasonal

def run(oldzip:Path,planzip:Path,workers:Path,output:Path):
    require(not output.exists() or not any(output.iterdir()),
            "Never overwrite a signed real original photo source")
    original=read_base(oldzip)
    plan,src=read_plan(planzip)
    accepted,worker=read_worker(workers,plan,src)
    new=pd.DataFrame(accepted,columns=original.columns)
    require(len(new)==SHARD_SIZE and new.photo_id.is_unique and
            new.photo_id.isin(original.photo_id).sum()==0,
            "Unfinished or previously seen genuine photo IDs cannot enter new cohort")
    appended=pd.concat([original,new],ignore_index=True)
    require(len(appended)==309695 and appended.photo_id.is_unique and
            appended.iloc[:309439].reset_index(drop=True).equals(original),
            "Original 309439 photo pixels/date/hue/ROI records were changed")
    nnew=int(new.morph.isin(HUES).sum())
    require(nnew==worker["four_hue_classified"] and nnew>0 and
            appended.morph.isin(HUES).sum()==149662+nnew and
            (appended.morph.isin(HUES)&appended.qualified_date.str.len().eq(10)).sum()==149653+nnew,
            "Actually measured source photographs or dates have not entered photo denominator")
    before_s,_=support(original)
    after_s,calendar=support(appended)
    output.mkdir(parents=True,exist_ok=True)
    path=output/OUT_SOURCE
    appended.to_csv(path,index=False,lineterminator="\n",
                    compression={"method":"gzip","mtime":0})
    stats=[]
    for zone in ("north_extratropical","south_extratropical"):
        q=calendar[calendar.zone.eq(zone)&calendar.date.notna()]
        for season in ("summer","winter"):
            v=q[q.season.eq(season)]
            v=v[v.classified]
            stats.append({
              "hemisphere":zone,"season":season,
              "original_classifiable_photo_ROIs":len(v),
              "white_original_classifiable_photo_ROIs":int(v.white.sum()),
              "white_fraction":round(float(v.white.mean()),8) if len(v) else None,
            })
    seasons=output/OUT_PHOTOS
    pd.DataFrame(stats).to_csv(seasons,index=False,lineterminator="\n")
    comp=output/OUT_SUPPORT
    with comp.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=[
            "hemisphere","before_same_taxon_siteyear_min1","after_same_taxon_siteyear_min1",
            "before_same_taxon_siteyear_min2","after_same_taxon_siteyear_min2"])
        w.writeheader()
        for hemi in ("north_extratropical","south_extratropical"):
            a=before_s[hemi];b=after_s[hemi]
            w.writerow({
                "hemisphere":hemi,
                "before_same_taxon_siteyear_min1":a[
                 "same_original_taxon_025degree_cell_same_seasonyear_both_summer_winter_min1_obs"],
                "after_same_taxon_siteyear_min1":b[
                 "same_original_taxon_025degree_cell_same_seasonyear_both_summer_winter_min1_obs"],
                "before_same_taxon_siteyear_min2":a[
                 "same_original_taxon_025degree_cell_same_seasonyear_both_summer_winter_min2_obs"],
                "after_same_taxon_siteyear_min2":b[
                 "same_original_taxon_025degree_cell_same_seasonyear_both_summer_winter_min2_obs"],
            })
    report={
        "schema":"mois_targeted256_actual_original_ROI_photo_season_support_v1",
        "source_3430502_photo_IDs_NOT_images":3430502,
        "source_reserved_targeted_photo_IDs_NOT_pixel_measurements":4096,
        "number_original_targeted_photo_pixels_actually_measured":SHARD_SIZE,
        "number_unmeasured_source_only_targeted_photo_IDs":RESERVE_N-SHARD_SIZE,
        "source_prior_309439_real_photo_SHA256":PHOTO_BASE_SHA,
        "targeted_original_photo_source_selection_SHA256":
           plan["first_wave_original_4096_csv_sha256"],
        "actual_original_pixel_ROI_model_frozen_SHA":METHOD,
        "frozen_targeted_worker_terminal_SHA256":
           worker["original_new_photo_hue_terminal_SHA256"],
        "target_source_selection_prior_hue_label_not_read":True,
        "original_photo_colour_model_source_attested":True,
        "before_309439_real_original_photo_IDs":309439,
        "before_149662_original_photo_ROIs_four_hue":149662,
        "before_149653_original_photo_ROIs_four_hue_plus_observation_day":149653,
        "additional_new_original_photos_with_actual_pixel_ROI":SHARD_SIZE,
        "additional_four_hue_photo_ROIs":nnew,
        "additional_original_photo_ROIs_technical_unknown":SHARD_SIZE-nnew,
        "additional_four_hue_photo_ROIs_with_original_observation_day":nnew,
        "expanded_real_original_photo_IDs":len(appended),
        "expanded_original_classified_four_hue_photo_ROIs":int(
            appended.morph.isin(HUES).sum()),
        "expanded_original_date_and_classified_four_hue_photo_ROIs":int(
            (appended.morph.isin(HUES)&appended.qualified_date.str.len().eq(10)).sum()),
        "before_named_taxon_cell_seasonyear_support":before_s,
        "after_named_taxon_cell_seasonyear_support":after_s,
        "expanded_actual_original_photo_source_SHA256":hashlib.sha256(
            path.read_bytes()).hexdigest(),
        "expanded_photo_white_source_season_csv_SHA256":sha(seasons),
        "expanded_same_taxon_siteyear_support_csv_SHA256":sha(comp),
        "strict_targeted_sampling_not_a_globally_representative_colour_census":True,
        "NO_original_global_moisture_genus_water_maxT_refit_or_causal_claim":True,
        "PHOTO_ROI_not_independently_confirmed_flower_open_petals":True,
        "FCP_intraspecific_D_W_or_50km_results_NOT_imported":True,
    }
    (output/OUT_RESULT).write_text(json.dumps(report,indent=2,
                                             ensure_ascii=False,sort_keys=True)+"\n")
    print("MOIS_REAL_TARGETED_ORIGINAL_PHOTO_PIXELS_AND_PHENOTIME_SUPPORT",
        json.dumps({k:report[k] for k in (
           "additional_new_original_photos_with_actual_pixel_ROI",
           "additional_four_hue_photo_ROIs",
           "expanded_real_original_photo_IDs",
           "expanded_original_date_and_classified_four_hue_photo_ROIs",
           "before_named_taxon_cell_seasonyear_support",
           "after_named_taxon_cell_seasonyear_support")},sort_keys=True),flush=True)
    return report


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--prior309439-original-source",required=True,type=Path)
    p.add_argument("--targeted4096-source-frozen-zip",required=True,type=Path)
    p.add_argument("--targeted256-real-pixel-workers",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    x=p.parse_args()
    run(x.prior309439_original_source,x.targeted4096_source_frozen_zip,
        x.targeted256_real_pixel_workers,x.output)
