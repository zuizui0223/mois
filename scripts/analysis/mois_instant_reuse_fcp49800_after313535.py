#!/usr/bin/env python3
"""MOIS speed-up: REUSE 49,800 already measured FCP original RGB/ROI-v4 photos.

Never re-download/reclassify an original image when FCP has already completed
the same frozen ROI-v4 method. Verify FCP's complete 256-terminal source,
preserve all previous MOIS 313535 photo colour IDs and original dates,
deduplicate by original photo_id, and exclude season + ongoing fifth4096
source-photo reservations to avoid future double counting.

This is an extra, nonprobability high-depth 498-taxon photographic sampling
arm. The 313535 source remains the PRIMARY worldwide descriptive denominator
until photo/source-propensity sensitivity is explicitly evaluated.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis import mois_global313535_water_light_heat_habitat_partition as prior
from scripts.analysis.mois_import_full_FCP_wave4_completed_original_photos import (
    read_fcp_wave4,
    EXPECTED_FCP_PHOTOS,
)
from scripts.analysis import mois_next4096_original_photo_ids_after313535 as wave

PHOTO_OUT="mois_primary313535_plus_FCP49800_nonoverlapping_verified_RGB_photo_labels.parquet"
NOVEL_OUT="mois_only_new_FCP_wave4_genuine_photo_ROIs_not_reserved.parquet"
AUDIT_OUT="mois_reused_FCP49800_photo_status_overlap_and_colour_quality.csv"
RECEIPT_OUT="mois_reused_complete_FCP49800_after313535_photoID_source_receipt.json"
HUES=prior.HUES
SOURCE_FCP="FCP_wave4_49800_original_ROI_v4_nonrandom_498taxa"
EXPECTED_ORIGINAL=prior.N_PHOTOS
DATE=re.compile(r"^\d{4}-\d{2}-\d{2}$")


def need(ok,message):
    if not ok:raise ValueError(message)


def digest(path:Path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(1<<20),b""):
            h.update(block)
    return h.hexdigest()


def _csv_reservation(root, filename, receipt_filename, expected_key, expected_n):
    plan=root/filename
    receiptfile=root/receipt_filename
    need(plan.is_file() and receiptfile.is_file(),
         "Missing original photo-ID source-only reservation evidence")
    doc=json.loads(receiptfile.read_text())
    need(doc.get("status")=="SOURCE_ONLY_4096_RESERVED_NO_NEW_PIXEL_HUES"
         and doc.get("new_photo_pixels_actually_measured") is False
         and doc.get(expected_key)==digest(plan)
         and doc.get("first_wave_original_source_photo_ids")==expected_n,
         "Source-only reservation mismatches signed original photo roster")
    t=pd.read_csv(plan,dtype=str,keep_default_na=False)
    need(len(t)==expected_n and t.photo_id.is_unique
         and t.photo_id.str.fullmatch(r"\d+").all(),
         "Original reserved 4096 photo IDs incomplete or duplicated")
    return set(t.photo_id), doc


def reservation_ids(fifth,season,original_ids):
    fifth_ids,fifth_rec=_csv_reservation(
        fifth,wave.SOURCE,wave.FREEZE,
        "first_wave_original_4096_csv_sha256",4096)
    need(fifth_rec.get("schema")==
         "mois_next4096_after_true313535_photo_ID_blind_source_v1"
         and fifth_rec.get("prior_real_original_photo_IDs")==EXPECTED_ORIGINAL,
         "Fifth original next4096 source reservation does not match 313535")
    season_ids,season_rec=_csv_reservation(
        season,"mois_targeted4096_photo_original_census_blind_selection.csv",
        "mois_first4096_original_source_only_freeze.json",
        "first_wave_original_4096_csv_sha256",4096)
    need(season_rec.get("schema")==
         "mois_3430502_targeted_opposite_season_after309439_and_reserved4096_v1"
         and season_rec.get("prior_measured_original_photo_ids")==309439,
         "Independent season-SI photo reservation source drift")
    need(fifth_ids.isdisjoint(season_ids)
         and fifth_ids.isdisjoint(original_ids)
         and season_ids.isdisjoint(original_ids),
         "Fifth/global and SI reservations overlap existing actual photo pixels")
    return fifth_ids,season_ids,fifth_rec,season_rec


def main_source(primary, fcp_root, fifth_root, season_root):
    p,manifest=prior.load_source(primary)
    # Explicit original PHOTO release; retain prior exact source cohort fields.
    original=pd.read_parquet(primary/prior.PHOTO)
    need(original.photo_id.is_unique and len(original)==EXPECTED_ORIGINAL
         and original.photo_id.tolist()==p.photo_id.tolist()
         and original.classified.astype(bool).sum()==prior.N_HUES
         and original.white.astype(bool).sum()==prior.WHITE_N,
         "Prior photographed phenotype colour source SHA identity drift")
    records,receipt=read_fcp_wave4(fcp_root)
    need(len(records)==EXPECTED_FCP_PHOTOS and
         receipt["measurement_original_fixed_colour_worker"] is True,
         "FCP 49800 real source/photo worker was not completed")
    original_ids=set(original.photo_id.astype(str))
    fifth,season,fifthrec,seasonrec=reservation_ids(
        fifth_root,season_root,original_ids)
    extra=[]
    counts=Counter()
    prior_label=original.set_index("photo_id")[["morph","inat_taxon_id"]]
    novel_photos=set()
    for row in records:
        pid=row["photo_id"]
        classified=row["morph"] in HUES
        counts["original_FCP_wave4_measured_photo_attempts"]+=1
        counts["original_FCP_wave4_four_hue_ROIs"]+=classified
        if pid in original_ids:
            counts["already_in_primary313535"]+=1
            prevrow=prior_label.loc[pid]
            counts["overlap_taxon_disagreement"]+=(
                str(prevrow["inat_taxon_id"])!=row["inat_taxon_id"])
            counts["overlap_both_classified"]+=(
                prevrow["morph"] in HUES and classified)
            counts["overlap_hue_disagreement"]+=(
                prevrow["morph"] in HUES and classified
                and prevrow["morph"]!=row["morph"])
            continue
        if pid in fifth:
            counts["deferred_to_running_fifth4096_wave"]+=1
            continue
        if pid in season:
            counts["excluded_SI_season_targeted4096"]+=1
            continue
        need(pid not in novel_photos,"FCP new photo source duplicate")
        novel_photos.add(pid)
        day=str(row.get("observed_on",""))
        if not DATE.fullmatch(day):
            day=""
        col=row["morph"] if classified else "mixed_uncertain"
        lon=float(row["longitude"])
        lat=float(row["latitude"])
        need(np.isfinite(lon) and np.isfinite(lat)
             and -180<=lon<=180 and -90<=lat<=90,
             "Previously ROI-verified FCP photograph has invalid source coordinates")
        new={
            "photo_id":pid,
            "inat_taxon_id":row["inat_taxon_id"],
            "species":row["species"],
            "morph":col,
            "observation_id":str(row["observation_id"]),
            "qualified_date":day,
            "old_mois_or_new_fcp":SOURCE_FCP,
            "latitude":lat,"longitude":lon,
            "classified":bool(classified),
            "white":bool(col=="white"),
        }
        need(set(new)==set(original.columns),
             "Historical published 313535 original photo source columns changed")
        extra.append(new)
        counts["actually_novel_ROI_processed_photo_IDs"]+=1
        counts["actually_novel_photo_four_hue_ROIs"]+=classified
        counts["actually_novel_photo_dated_four_hue_ROIs"]+=classified and bool(day)
        counts["novel_visible_white_ROIs"]+=col=="white"
        counts["novel_original_RBG_pixel_unknown"]+=not classified
    need(counts["original_FCP_wave4_measured_photo_attempts"]==
         counts["already_in_primary313535"]+
         counts["deferred_to_running_fifth4096_wave"]+
         counts["excluded_SI_season_targeted4096"]+
         counts["actually_novel_ROI_processed_photo_IDs"]==
         EXPECTED_FCP_PHOTOS,
         "FCP 49,800 photo-ID category accounting fails")
    # New model actual ratio is source-specific, not globally weighted.
    add=pd.DataFrame(extra,columns=list(original.columns))
    need(len(add)>0 and add.photo_id.is_unique
         and add.photo_id.isin(original.photo_id).sum()==0,
         "FCP source completed but no truly new photo pixels found")
    new=pd.concat((original,add),ignore_index=True)
    need(new.photo_id.is_unique
         and len(new)==EXPECTED_ORIGINAL+len(add)
         and new.iloc[:EXPECTED_ORIGINAL].reset_index(drop=True).equals(original),
         "Original MOIS313535 photo classification/date changed")
    return original,add,new,dict(counts),receipt,manifest,fifthrec,seasonrec


def run(mainrelease,fcp_source,fifth,season,output):
    need(not output.exists() or not any(output.iterdir()),
         "Refuse to overwrite immutable measured source-photo archive")
    original,extra,full,counts,fcp,base,fifthrec,seasonrec=main_source(
        mainrelease,fcp_source,fifth,season)
    output.mkdir(parents=True,exist_ok=True)
    full.to_parquet(output/PHOTO_OUT,index=False,compression="zstd")
    extra.to_parquet(output/NOVEL_OUT,index=False,compression="zstd")
    records=[]
    for key,value in sorted(counts.items()):
        records.append({"item":key,"real_original_photo_count":int(value)})
    pd.DataFrame(records).to_csv(output/AUDIT_OUT,index=False,
                                  lineterminator="\n")
    rates=extra.loc[extra.classified]
    report={
        "schema":"mois_FCP49800_real_RGB_photoID_reuse_after313535_v1",
        "original_MOIS_image_processed_photo_attempts":EXPECTED_ORIGINAL,
        "original_MOIS_classified_four_hue_ROIs":prior.N_HUES,
        "source_FCP_original_49800_actual_image_terminal_receipt_SHA256":
            digest(next(iter(fcp_source.rglob("terminal_result.json")))),
        "source_FCP_original_real_ROI_data_SHA256":fcp["measured_sha256"],
        "FCP_original_model_source_signature":fcp["source_metadata_sha256"],
        "verified_same_original_ROI_v4_method_no_new_RGB_inference":True,
        "FCP_original_photo_attempts_with_real_pixels":EXPECTED_FCP_PHOTOS,
        "reconciliation":counts,
        "new_original_photo_pixel_attempts_admitted":len(extra),
        "new_original_classified_four_hue_ROIs":int(extra.classified.sum()),
        "new_original_white_photo_ROIs":int(extra.white.sum()),
        "new_original_dated_four_hue_ROIs":int(
            (extra.classified&extra.qualified_date.str.len().eq(10)).sum()),
        "expanded_actual_original_photo_attempts":len(full),
        "expanded_original_classified_four_hue_ROIs":int(full.classified.sum()),
        "expanded_original_white_ROIs":int(full.white.sum()),
        "new_FCP_specific_photo_white_fraction":
            float(rates.white.mean()) if len(rates) else None,
        "base313535_original_photo_white_fraction":
            prior.WHITE_N/prior.N_HUES,
        "selection_recruitment":SOURCE_FCP,
        "new_498taxon_highdepth_arm_NOT_unweighted_global_flora_prevalence":True,
        "original_first313535_photo_hues_date_taxon_unchanged":True,
        "unclassified_FCP_photo_ROIs_not_nonwhite":True,
        "already_running_4096_real_pixel_reservation_excluded":True,
        "separate_SI_priority_reservation_excluded":True,
        "old_FCP_within_species_D_W_and_50km_inference_NOT_imported":True,
        "no_actual_open_flower_biochemical_UV_pigment_or_fitness_confirmation":True,
        "no_50_50_evolutionary_colour_equilibrium_assumed":True,
        "source_MOIS313535_original_photo_SHA256":
            base["files_sha256"][prior.PHOTO],
        "source_FCP_vs_fifth_photoID_reservation_SHA256":
            fifthrec["first_wave_original_4096_csv_sha256"],
        "source_SI_season_photoID_reservation_SHA256":
            seasonrec["first_wave_original_4096_csv_sha256"],
        "derived_file_SHA256":{p:digest(output/p) for p in (
            PHOTO_OUT,NOVEL_OUT,AUDIT_OUT)},
    }
    (output/RECEIPT_OUT).write_text(json.dumps(
        report,indent=2,sort_keys=True,ensure_ascii=False,
        allow_nan=False)+"\n")
    print("MOIS_INSTANT_REUSE_FCP_49800_ACTUAL_RGB_PHOTO_ROI",json.dumps({
        "already_processed_FCP_original_RGB":EXPECTED_FCP_PHOTOS,
        "new_original_photo_pixels_admitted":len(extra),
        "new_visible_four_hue_ROIs":int(extra.classified.sum()),
        "expanded_image_attempts":len(full),
        "old_photo_IDs_retained":EXPECTED_ORIGINAL,
        "FCP_source_specific_white_fraction":report[
            "new_FCP_specific_photo_white_fraction"],
        "overlap_audit":counts,
    },sort_keys=True),flush=True)
    return report


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--previous313535-Parquet-release",type=Path,required=True)
    p.add_argument("--FCP-49800-completed-original-photo-pixels",type=Path,required=True)
    p.add_argument("--ongoing-fifth4096-original-photo-reservation",type=Path,required=True)
    p.add_argument("--season4096-source-only-photo-reservation",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    run(args.previous313535_Parquet_release,
        args.FCP_49800_completed_original_photo_pixels,
        args.ongoing_fifth4096_original_photo_reservation,
        args.season4096_source_only_photo_reservation,args.output)
