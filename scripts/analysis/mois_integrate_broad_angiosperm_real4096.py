#!/usr/bin/env python3
"""Join ONLY real original ROI-v4 photographed RGB outcomes for global-breadth source.

This 4096-image wave was sampled to discover NEW angiosperm taxa/genera,
not deepen FCP's few spatially repeated species. Read all 16 original
256-photo SHA-signed image workers and never count source-only photo IDs
as visible petals. Preserve immutable 313535 old photo hue/date rows and
keep the new photographic source arm separate for bias diagnostics.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

import pandas as pd

from scripts.analysis.mois_reassemble_original_new4096_measurements import (
    HUES, METHOD, MEASURED_SCHEMA,
)
from scripts.acquisition import mois_global_angiosperm_breadth_source as broad

PHOTO="primary313535_plus_broad4096_genuine_original_hue_photos.parquet"
NOVEL="mois_first_broad4096_genuine_RGB_hue_photo_taxon_quality.csv"
RECEIPT="mois_global_angiosperm_broad4096_actual_photo_colour_taxa_receipt.json"
EXPECTED=313535
NEW=4096


def require(ok,why):
    if not ok:raise ValueError(why)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def verify_old(previous):
    r=json.loads((previous/broad.CORE_MANIFEST).read_text())
    src=previous/broad.CORE_PHOTO
    require(r["schema"]=="mois_313535_genuine_original_ROI_feature_lake_v2"
            and r["actual_original_photo_attempts"]==EXPECTED
            and r["four_hue_labels"]==151239
            and r["files_sha256"][broad.CORE_PHOTO]==sha(src),
            "Original primary global photograph/classifier baseline changed")
    old=pd.read_parquet(src)
    require(len(old)==EXPECTED and old.photo_id.is_unique
            and int(old.classified.sum())==151239
            and int(old.white.sum())==71427,
            "Immutable published source original-photo ID or 4-hue rows drifted")
    return old,r


def verified_blind_plan(plan):
    doc=json.loads((plan/broad.RECEIPT).read_text())
    filename=plan/broad.FIRST
    require(doc["schema"]=="mois_open_global_angiosperm_breadth_first4096_v1"
            and doc["status"]=="SOURCE_ONLY_4096_RESERVED_NO_NEW_PIXEL_HUES"
            and doc["first_wave_original_source_photo_ids"]==NEW
            and doc["first_wave_distinct_photo_ids"]==NEW
            and doc["new_photo_pixels_actually_measured"] is False
            and doc["first_wave_original_4096_csv_sha256"]==sha(filename)
            and doc["source_photo_IDs_disjoint_from_both_running_and_SI_reservations"] is True
            and doc["no_photo_colour_or_fitness_or_bloom_label_used"] is True,
            "Source-only global new-taxon photo plan not verified")
    p=pd.read_csv(filename,dtype=str,keep_default_na=False)
    require(len(p)==NEW and p.photo_id.is_unique
            and p.measurement_id.is_unique
            and p.observation_id.is_unique
            and p.inat_taxon_id.nunique()==doc["first_wave_named_taxa"]
            and p.genus.nunique()==doc["first_wave_named_genera"]
            and (p.source_snapshot_date==doc["source_snapshot_date"]).all(),
            "Worldwide angiosperm breadth/source photo identity drifted")
    return p,doc


def one_worker(k,expected,workers,doc):
    expected=expected.set_index("measurement_id",drop=False)
    candidates=list(workers.rglob(f"mois_original_256_fresh_ROI_shard_{k:02d}.csv"))
    receipts=list(workers.rglob(f"mois_original_256_ROI_shard_{k:02d}_receipt.json"))
    require(len(candidates)==len(receipts)==1,
            f"Unfinished source-true global-breadth 256-photo worker: {k}")
    receipt=json.loads(receipts[0].read_text())
    require(receipt["schema"]==MEASURED_SCHEMA
            and int(receipt["shard"])==k
            and receipt["terminal_rows"]==256
            and receipt["original_source_ids"]==256
            and receipt["source_sha256"]==
                doc["per_shard_sha256"][str(k)]
            and receipt["source_frozen_4096_selection_sha256"]==
                doc["first_wave_original_4096_csv_sha256"]
            and receipt["source_model_attestation"]==METHOD
            and receipt["pixels_ephemeral_no_original_image_bytes_persisted"] is True
            and receipt["original_new_photo_hue_terminal_SHA256"]==sha(candidates[0]),
            "Genuine original full-image ROI-v4 worker source/model terminal invalid")
    measured=pd.read_csv(candidates[0],dtype=str,keep_default_na=False)
    require(len(measured)==256 and measured.measurement_id.is_unique
            and set(measured.measurement_id)==set(expected.index),
            "Expected 256 original photo RGB terminal outcomes missing")
    status=Counter()
    for row in measured.to_dict("records"):
        v=expected.loc[row["measurement_id"]]
        for key in ("photo_id","inat_taxon_id","species","observer_id",
                    "latitude","longitude","photo_license"):
            require(str(row[key])==str(v[key]),
                    "Original RGB-image-photo taxon/observer/site/license differs")
        known=row["morph"] in HUES
        require(row["morph"] in (*HUES,"mixed_uncertain")
                and (row["measurement_status"]=="classified_four_state_morph")==known
                and row["classified_four_hue"]==str(known).lower(),
                "Original 4-hue classification status not observed in actual model")
        status["classified"]+=known
        status["unknown"]+=not known
    require(status["classified"]==receipt["four_hue_classified"]
            and status["unknown"]==receipt["unclassified"],
            "Original pixel/photo classification result differs from worker receipt")
    return measured,status,receipt


def join(prior,plan,workers,out):
    require(not out.exists() or not any(out.iterdir()),
            "Never overwrite a completed original image/photo source data archive")
    old,old_receipt=verify_old(prior)
    source,freeze=verified_blind_plan(plan)
    require(not source.photo_id.isin(old.photo_id).any(),
            "Globally broad photo ID source duplicates earlier actual source")
    all_results=[]
    receipts=[]
    for k in range(broad.N_SHARDS):
        shard=plan/f"{broad.PLAN_PREFIX}{k:02d}_256.csv"
        require(sha(shard)==freeze["per_shard_sha256"][str(k)],
                "Original broad-source shard tampered before pixel measurement")
        subset=source[source.shard.eq(str(k))]
        require(len(subset)==256,
                "Original taxonomic-breadth source worker shape invalid")
        actual,status,receipt=one_worker(k,subset,workers,freeze)
        all_results.append(actual)
        receipts.append({
            "shard":k,
            "original_photo_ROIs":256,
            "classified_visible_four_hue_ROIs":status["classified"],
            "unknown_original_photo_ROIs":status["unknown"],
            "terminal_SHA256":sha(next(iter(workers.rglob(
                f"mois_original_256_fresh_ROI_shard_{k:02d}.csv")))),
            "frozen_model_verified":True,
        })
    actual=pd.concat(all_results,ignore_index=True)
    require(len(actual)==NEW and actual.measurement_id.is_unique
            and actual.photo_id.is_unique
            and not actual.photo_id.isin(old.photo_id).any(),
            "Global novel source photo RGB outcomes duplicate old images")
    matched=source.merge(
        actual[["measurement_id","morph","measurement_status",
                "classified_four_hue"]],
        on="measurement_id",validate="one_to_one",how="inner")
    require(len(matched)==NEW,
            "New photographed original results lost frozen taxon/source rows")
    frame=pd.DataFrame({
        "photo_id":matched.photo_id,
        "inat_taxon_id":matched.inat_taxon_id,
        "species":matched.species,
        "morph":matched.morph,
        "observation_id":matched.observation_id,
        "qualified_date":matched.observed_on,
        "old_mois_or_new_fcp":matched.source_arm,
        "latitude":pd.to_numeric(matched.latitude,errors="coerce"),
        "longitude":pd.to_numeric(matched.longitude,errors="coerce"),
        "classified":matched.morph.isin(HUES),
        "white":matched.morph.eq("white"),
    })
    require(set(frame.columns)==set(old.columns),
            "Published primary photo sample schema differs from new worldwide source")
    frame=frame[old.columns]
    allphotos=pd.concat((old,frame),ignore_index=True)
    require(len(allphotos)==EXPECTED+NEW and allphotos.photo_id.is_unique
            and allphotos.iloc[:EXPECTED].reset_index(drop=True).equals(old),
            "Previous photos/hues/dates overwritten by global source")
    out.mkdir(parents=True,exist_ok=True)
    allphotos.to_parquet(out/PHOTO,index=False,compression="zstd")
    matched.to_csv(out/NOVEL,index=False,lineterminator="\n")
    oldtaxa=set(old.inat_taxon_id.astype(str))
    oldgenera=set(old.species.astype(str).str.split().str[0])
    genuinely_new_taxa=set(matched.loc[
        matched.morph.isin(HUES),"inat_taxon_id"])-oldtaxa
    genuinely_new_genera=set(matched.loc[
        matched.morph.isin(HUES),"genus"])-oldgenera
    classified=int(frame.classified.sum())
    r={
        "schema":"mois_global_breadth_first4096_ACTUAL_original_RGB_photo_ROIs_v1",
        "original_photo_RGB_attempts_baseline":EXPECTED,
        "genuine_new_breadth_photo_attempts":NEW,
        "expanded_actual_original_photo_attempts":len(allphotos),
        "new_original_four_hue_classifier_successes":classified,
        "new_original_four_hue_classifier_unknowns":NEW-classified,
        "new_original_white_visible_hue_ROIs":int(frame.white.sum()),
        "new_original_taxa_with_classified_visible_hue_not_in_core":
            len(genuinely_new_taxa),
        "new_original_genera_with_classified_visible_hue_not_in_core":
            len(genuinely_new_genera),
        "new_taxa_per_1000_real_photo_RGB_attempts":
            1000*len(genuinely_new_taxa)/NEW,
        "new_genera_per_1000_real_photo_RGB_attempts":
            1000*len(genuinely_new_genera)/NEW,
        "selected_original_source_named_taxa":freeze["first_wave_named_taxa"],
        "selected_original_source_named_genera":freeze["first_wave_named_genera"],
        "chosen_4hue_photo_classifier_success_rate":classified/NEW,
        "original_photo_IDs_hues_dates_unchanged":True,
        "separate_sampling_arm_tagged_and_NOT_unbiased_global_white_frequency":True,
        "all_16_original_photo_RGB_worker_receipts":receipts,
        "verified_original_source_metadata_first4096_photo_SHA256":
            freeze["first_wave_original_4096_csv_sha256"],
        "verified_previous313535_original_photo_parquet_SHA256":
            old_receipt["files_sha256"][broad.CORE_PHOTO],
        "new_combined_photo_parquet_SHA256":sha(out/PHOTO),
        "new_original_photo_taxon_quality_csv_SHA256":sha(out/NOVEL),
        "source_open_licensed_photo_but_botanical_bloom_not_verified":True,
        "new_photo_pixel_source_once_only":True,
        "no_unmeasured_metadata_used_as_hue":True,
        "no_photographic_50_50_equilibrium_pigment_cost_or_fitness_claim":True,
        "source_climate_not_yet_extracted_at_new_sites":True,
    }
    (out/RECEIPT).write_text(json.dumps(
        r,sort_keys=True,indent=2,ensure_ascii=False)+"\n")
    print("MOIS_REAL_NEW_ANGIOSPERM_SPECIES_FIRST_RGB_ROIS",json.dumps({
        k:r[k] for k in (
            "genuine_new_breadth_photo_attempts",
            "new_original_four_hue_classifier_successes",
            "new_original_taxa_with_classified_visible_hue_not_in_core",
            "new_original_genera_with_classified_visible_hue_not_in_core",
            "new_taxa_per_1000_real_photo_RGB_attempts",
            "new_genera_per_1000_real_photo_RGB_attempts")
    },sort_keys=True),flush=True)
    return r


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--prior313535-hue-Parquet-release",type=Path,required=True)
    p.add_argument("--source-plan-first4096",type=Path,required=True)
    p.add_argument("--actual16-broad-original-RGB-workers",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    join(a.prior313535_hue_Parquet_release,
         a.source_plan_first4096,a.actual16_broad_original_RGB_workers,a.output)
