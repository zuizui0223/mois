#!/usr/bin/env python3
"""Admit fifth 4096 image ROIs only after each of sixteen genuine 256 workers.

Frozen previous 313535 original source photo/day/label rows are immutable.
No missing worker is counted; no original licensed image bytes are stored.
A separately targeted 4096 source reservation is never mixed into this
source-balanced but not globally probability-representative new wave.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from pathlib import Path
import pandas as pd
from scripts.analysis.mois_reassemble_original_new4096_measurements import (
    BASE_COLS,HUES,METHOD,MEASURED_SCHEMA,
)
from scripts.analysis import mois_next4096_original_photo_ids_after313535 as src

OUT_SOURCE="mois_full317631_actual_original_photo_RGB_hues_and_dates.csv.gz"
OUT_RECEIPT="mois_fifth4096_real_original_photo_source_reconciled_receipt.json"
OUT_HUES="mois_fifth4096_real_new_photo_colour_count.csv"
EXPECTED_TOTAL=317631


def need(ok,msg):
    if not ok:raise ValueError(msg)


def file_sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def read_source(prior,prior_rec,plan):
    before,r=src.parse_old(prior,prior_rec)
    doc=json.loads((plan/src.FREEZE).read_text())
    expect=plan/src.SOURCE
    need(doc["schema"]=="mois_next4096_after_true313535_photo_ID_blind_source_v1"
         and doc["status"]=="SOURCE_ONLY_4096_RESERVED_NO_NEW_PIXEL_HUES"
         and doc["new_source_photo_IDs"]==4096
         and doc["new_photo_pixels_actually_measured"] is False
         and doc["first_wave_original_4096_csv_sha256"]==file_sha(expect)
         and doc["prior_full313535_source_SHA256"]==file_sha(prior)
         and doc["targeted4096_excluded_as_reserved_NOT_independent_world_sample"] is True,
         "Fifth new source photo plan not original source-verified")
    original=pd.read_csv(prior,compression="gzip",dtype=str,
                         keep_default_na=False,low_memory=False)
    need(tuple(original.columns)==(*BASE_COLS,"qualified_date",
                                   "original_source_date_provenance"),
         "Previous original ROI source schema drift")
    with expect.open(newline="",encoding="utf-8") as f:
        rows=list(csv.DictReader(f))
    need(len(rows)==4096 and len({r["photo_id"] for r in rows})==4096
         and len({r["measurement_id"] for r in rows})==4096
         and set(r["photo_id"] for r in rows).isdisjoint(before),
         "Source photo IDs reused from previously actually processed pixels")
    return original,rows,doc


def complete_worker(source,receipt_file,result_file,freeze,k):
    doc=json.loads(receipt_file.read_text())
    need(doc["schema"]==MEASURED_SCHEMA
         and doc["shard"]==k
         and doc["terminal_rows"]==256
         and doc["original_source_ids"]==256
         and doc["source_sha256"]==freeze["per_shard_sha256"][str(k)]
         and doc["source_frozen_4096_selection_sha256"]==
             freeze["first_wave_original_4096_csv_sha256"]
         and doc["source_model_attestation"]==METHOD
         and doc["pixels_ephemeral_no_original_image_bytes_persisted"] is True
         and doc["original_new_photo_hue_terminal_SHA256"]==file_sha(result_file),
         f"Actual original RGB/YOLO/EfficientSAM image worker {k} was not source-verified")
    with result_file.open(newline="",encoding="utf-8-sig") as f:
        records=list(csv.DictReader(f))
    expected={x["measurement_id"]:x for x in source}
    need(len(records)==256 and len({x["measurement_id"] for x in records})==256
         and set(x["measurement_id"] for x in records)==set(expected),
         "Completed original 256-photo terminal records not identical to source")
    classify=0
    for r in records:
        x=expected[r["measurement_id"]]
        need(all(r[key]==x[key] for key in (
             "photo_id","inat_taxon_id","species","observation_id","observer_id",
             "latitude","longitude","photo_license"))
             and r["morph"] in (*HUES,"mixed_uncertain")
             and (r["morph"] in HUES)==(
                r["measurement_status"]=="classified_four_state_morph")
             and r["classified_four_hue"]==(
                "true" if r["morph"] in HUES else "false"),
             "Incorrect photographic pixel, taxonomy, observer or classifier status")
        classify+=r["morph"] in HUES
    need(classify==doc["four_hue_classified"]
         and 256-classify==doc["unclassified"],
         "Finished worker original classifiable/missing accounting mismatch")
    return {r["measurement_id"]:r for r in records},doc


def run(prior,priorrec,plan,workers,out):
    need(not out.exists() or not any(out.iterdir()),
         "Do not overwrite verified previously measured source")
    old,planrows,freeze=read_source(prior,priorrec,plan)
    rows={}
    receipts=[]
    for k in range(16):
        selected=[x for x in planrows if int(x["shard"])==k]
        need(len(selected)==256,
             "Fifth original image source partition missing")
        target=plan/f"mois_source_only_shard_{k:02d}_256.csv"
        need(file_sha(target)==freeze["per_shard_sha256"][str(k)],
             "Original independent photographer source shard SHA changed")
        found=list(workers.rglob(f"mois_original_256_fresh_ROI_shard_{k:02d}.csv"))
        docs=list(workers.rglob(f"mois_original_256_ROI_shard_{k:02d}_receipt.json"))
        need(len(found)==len(docs)==1,
             "Not all 16 original image RGB ROI workers genuinely completed")
        data,r=complete_worker(selected,docs[0],found[0],freeze,k)
        need(not (set(rows)&set(data)),"Repeated original photo measurement ID across workers")
        rows.update(data);receipts.append({
            "shard":k,"classified":r["four_hue_classified"],
            "unclassified":r["unclassified"],
            "original_terminal_sha":r["original_new_photo_hue_terminal_SHA256"],
            "original_model_attested":True,
        })
    need(len(rows)==4096 and set(rows)=={x["measurement_id"] for x in planrows},
         "Full fourth/source fifth 4096 photo-RGB pixels not measured")
    new=[]
    for selected in planrows:
        pixel=rows[selected["measurement_id"]]
        r={key:"" for key in BASE_COLS}
        for key in BASE_COLS:
            r[key]=pixel.get(key,"") or selected.get(key,"") or ""
        r["source_original_measured"]="mois_fifth4096_new_original_ROI_v4"
        r["old_mois_or_new_fcp"]="mois_new_original_ROI"
        r["observed_on"]=selected["observed_on"]
        r["qualified_date"]=selected["observed_on"]
        r["original_source_date_provenance"]=(
            "fcp_full3430502_photo_observed_on_after_fifth_new4096_real_pixel_ROI")
        new.append(r)
    extra=pd.DataFrame(new,columns=list(old.columns)).fillna("")
    need(extra.qualified_date.str.fullmatch(r"\d{4}-\d{2}-\d{2}").all()
         and extra.photo_id.is_unique
         and int(extra.morph.isin(HUES).sum())==
                sum(z["classified"] for z in receipts),
         "New image source row/date/hue classification reconciliation failed")
    allrows=pd.concat((old,extra),ignore_index=True)
    need(len(allrows)==EXPECTED_TOTAL and allrows.photo_id.is_unique
         and old.equals(allrows.iloc[:len(old)].reset_index(drop=True))
         and int(allrows.morph.isin(HUES).sum())==
             src.N_CLASSIFIED+int(extra.morph.isin(HUES).sum()),
         "Genuine already measured source photos changed or duplicates counted")
    out.mkdir(parents=True,exist_ok=True)
    resultfile=out/OUT_SOURCE
    allrows.to_csv(resultfile,index=False,lineterminator="\n",
                   compression={"method":"gzip","mtime":0})
    hue_counts=extra.morph.value_counts().reindex(
        [*HUES,"mixed_uncertain"],fill_value=0)
    hue_counts.rename_axis("automatic_ROI_four_hue_or_unknown").reset_index(
        name="real_new_photo_RGB_attempts").to_csv(
            out/OUT_HUES,index=False,lineterminator="\n")
    full_class=int(allrows.morph.isin(HUES).sum())
    daterate=int((allrows.morph.isin(HUES)&
                  allrows.qualified_date.str.fullmatch(
                      r"\d{4}-\d{2}-\d{2}")).sum())
    report={
        "schema":"mois_fifth4096_original_real_ROI_photo_source_after313535_v1",
        "previous_real_original_photo_attempts":src.N_PRIOR,
        "new_original_RGB_image_attempts":src.N_NEW,
        "new_original_classifiable_ROIs":int(extra.morph.isin(HUES).sum()),
        "new_original_technical_unknown_ROIs":int((~extra.morph.isin(HUES)).sum()),
        "total_source_photo_pixels_actually_processed":len(allrows),
        "total_original_classified_four_hue":full_class,
        "total_original_date_plus_four_hue":daterate,
        "previous_313535_original_photo_rows_unchanged":True,
        "new_photo_source_priority_did_not_use_colour_or_climate":True,
        "source_original_photo_ID_and_observer_disjoint_within_wave":True,
        "source_new_photo_observations_not_world_flora_abundance_sample":True,
        "separately_prioritised_season4096_source_photo_IDs_excluded":True,
        "new_original_petal_anthesis_not_independently_verified":True,
        "frozen_model_method":METHOD,
        "worker_pixel_source_terminal_receipts":receipts,
        "source_photo_selection_SHA256":file_sha(plan/src.SOURCE),
        "previous_complete_original_source_SHA256":file_sha(prior),
        "expanded_complete_original_source_SHA256":file_sha(resultfile),
        "source_four_hue_counts_SHA256":file_sha(out/OUT_HUES),
        "not_photo_metadata_only_and_no_raw_image_files_persisted":True,
        "no_FCP_withinspecies_W_D_or_50km_estimand_imported":True,
    }
    (out/OUT_RECEIPT).write_text(
        json.dumps(report,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("MOIS_FIFTH4096_REAL_ORIGINAL_IMAGE_PIXELS_COMPLETE",
          json.dumps({k:report[k] for k in (
              "new_original_RGB_image_attempts",
              "new_original_classifiable_ROIs",
              "total_source_photo_pixels_actually_processed",
              "total_original_classified_four_hue",
              "total_original_date_plus_four_hue")},sort_keys=True),
          flush=True)
    return report


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--prior-full313535-real-photo-RGB",type=Path,required=True)
    p.add_argument("--prior-313535-source-receipt",type=Path,required=True)
    p.add_argument("--new-source-plan",type=Path,required=True)
    p.add_argument("--actual16-RGB-worker-terminals",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    run(a.prior_full313535_real_photo_RGB,
        a.prior_313535_source_receipt,a.new_source_plan,
        a.actual16_RGB_worker_terminals,a.output)
