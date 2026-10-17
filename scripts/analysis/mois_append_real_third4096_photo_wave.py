#!/usr/bin/env python3
"""Append real third 4096 original flower-hue ROI results to the actual 256504-source MOIS.

Photo reuse across FCP and MOIS is allowed; exact ID duplicates within MOIS
are forbidden. Require the prior two new-photo source receipts and ALL 16
third 256-image original ROI-v4 terminal SHA receipts. No unprocessed photo,
white-frequency imputation, or FCP within-species D/W/50km statistic.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
from pathlib import Path

from scripts.analysis.mois_reassemble_original_new4096_measurements import (
    load_real_wave,BASE_COLS,HUES,sha,require,
)

PRIOR_PHOTOS = 256504
PRIOR_NEW_PIXELS = 8192
THIRD_NEW_PIXELS = 4096
TOTAL_PHOTOS = 260600


def append(prior:Path,receipt:Path,third_plan:Path,third_roi:Path,out:Path)->dict:
    require(not out.exists() or not any(out.iterdir()),
            "Never overwrite an earlier completed original-photo source")
    base=json.loads(receipt.read_text(encoding="utf-8"))
    require(base.get("schema")=="mois_existing_248312_plus_completed_real_4096_photo_colour_waves_v1"
            and base.get("status")=="ACTUAL_COMPLETED_ORIGINAL_IMAGE_ROI_ADDITIONS_NOT_SOURCE_ONLY"
            and base.get("all_original_unique_photos")==PRIOR_PHOTOS
            and base.get("new_original_source_photos")==PRIOR_NEW_PIXELS
            and isinstance(base.get("all_four_hue_classified_photos"),int)
            and base.get("all_four_hue_classified_photos")>=120807
            and base.get("original_model_SHA_verified_for_every_photo") is True
            and set(base.get("source_waves",{}))=={"wave1","wave2"}
            and sha(prior)==base.get("all_real_source_photo_hue_csv_sha256"),
            "Prior full two-wave 256504 actual image-ROI source is not sealed")

    third,proof=load_real_wave(third_plan,third_roi,"wave3")
    require(len(third)==THIRD_NEW_PIXELS,"Third wave actual photo-ROI source unfinished")
    actual_new_classified=sum(r["morph"] in HUES for r in third)
    out.mkdir(parents=True,exist_ok=True)
    target=out/"mois_all_real_original_hue_measured_248312_plus_new_waves.csv.gz"
    unique=set()
    class_total=0
    with target.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw,mode="wb",mtime=0,filename="") as gz:
            with io.TextIOWrapper(gz,encoding="utf-8",newline="") as f:
                writer=csv.DictWriter(f,fieldnames=BASE_COLS,lineterminator="\n")
                writer.writeheader()
                with gzip.open(prior,"rt",newline="",encoding="utf-8") as b:
                    reader=csv.DictReader(b)
                    require(tuple(reader.fieldnames or ())==BASE_COLS,
                            "Original first+second photographed hue columns changed")
                    for row in reader:
                        pid=row["photo_id"]
                        require(pid and pid not in unique,
                                "Previous two waves duplicate original image ID")
                        unique.add(pid)
                        writer.writerow(row)
                        class_total+=row["morph"] in HUES
                require(len(unique)==PRIOR_PHOTOS and
                        class_total==base["all_four_hue_classified_photos"],
                        "Original 256504 source photo or colour denominator changed")
                for row in third:
                    pid=row["photo_id"]
                    require(pid and pid not in unique,
                            "Third original photographic photo ID overlaps previous source")
                    unique.add(pid)
                    require(set(row)==set(BASE_COLS),
                            "Third ROI source photo column names differ")
                    writer.writerow(row)
                    class_total+=row["morph"] in HUES
    require(len(unique)==TOTAL_PHOTOS and
            class_total==base["all_four_hue_classified_photos"]+actual_new_classified,
            "New actually measured third photo union incomplete")
    source_waves=dict(base["source_waves"])
    source_waves["wave3"]=proof
    receipt_out={
        "schema":"mois_existing_248312_plus_completed_real_4096_photo_colour_waves_v1",
        "status":"ACTUAL_COMPLETED_ORIGINAL_IMAGE_ROI_ADDITIONS_NOT_SOURCE_ONLY",
        "previous_source_original_photos":248312,
        "previous_source_photo_hues_classified":120807,
        "new_original_source_photos":PRIOR_NEW_PIXELS+THIRD_NEW_PIXELS,
        "new_four_hue_measured_photo_ids":base["new_four_hue_measured_photo_ids"]+actual_new_classified,
        "all_original_unique_photos":len(unique),
        "all_four_hue_classified_photos":class_total,
        "new_third_4096_source_four_hue_classified":actual_new_classified,
        "new_third_4096_source_unclassified":THIRD_NEW_PIXELS-actual_new_classified,
        "source_waves":source_waves,
        "original_historical_photo_source_sha256":base["original_historical_photo_source_sha256"],
        "prior_two_photo_wave_SHA256":sha(prior),
        "all_real_source_photo_hue_csv_sha256":sha(target),
        "genetic_polymorphism_spatial_D_W_or_50km_measured":False,
        "original_model_SHA_verified_for_every_photo":True,
        "new_source_photo_IDs_never_reclassified_as_confirmed_open_flower":True,
        "photo_based_genus_moisture_six_anchor_model_refitted":False,
        "original_RGB_photo_bytes_persisted":False,
    }
    (out/"mois_real_additional_photo_waves_receipt.json").write_text(
        json.dumps(receipt_out,ensure_ascii=False,indent=2,
                   sort_keys=True,allow_nan=False)+"\n",encoding="utf-8")
    print("MOIS_REAL_SOURCE_NEW_THIRD_4096_ORIGINAL_ROI_PHOTOS",
          json.dumps({"all_original_unique_photo_ids":len(unique),
                      "new_third_four_hues":actual_new_classified,
                      "total_four_hues":class_total},sort_keys=True),
          flush=True)
    return receipt_out


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--previous-256504-original-hue-gzip",type=Path,required=True)
    p.add_argument("--previous-two-waves-receipt",type=Path,required=True)
    p.add_argument("--third-4096-photo-source",type=Path,required=True)
    p.add_argument("--third-actual-roi-results",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    append(args.previous_256504_original_hue_gzip,
           args.previous_two_waves_receipt,
           args.third_4096_photo_source,
           args.third_actual_roi_results,
           args.output)
