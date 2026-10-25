#!/usr/bin/env python3
"""Restore original 9-colour ROI fractions lost ONLY in the flattened historic
911,767-photo checkpoint by exact photo-ID linking its authentic full NEW
source ROI rows. No photo image fetch/inference, no extrapolation of colours
from existing four-hue classes.

The original checkpoint stores a compact union and a separate
`new_completed_original_photo_ROIs.parquet` containing all frozen continuous
flower_fraction_* columns. Rejoin ONLY to that historic main source arm, and
reject photo/colour/measurement-status disagreements before restoration.
"""
from __future__ import annotations
import argparse,json,hashlib
from pathlib import Path
import numpy as np
import pandas as pd

from scripts.analysis.mois_real1043930_continuous9_source_lake import (
    OUTPUT, SCHEMA as PRIOR_SCHEMA, read_nine, PALETTE, PROVENANCE, sha, require,
)

SCHEMA="mois_latest1043930_recovered_history_original_checkpoint_nine_palette_v1"
OUTPUT_NEW="mois_real1043930_original_ROI9_with_recovered_historical_checkpoint.parquet"
HISTORY_FILE="new_completed_original_photo_ROIs.parquet"
HISTORY_SOURCE="historical_main"

def restore(old_lake,old_lake_receipt,history,history_receipt,out):
    old_lake,old_lake_receipt,history,history_receipt,out=map(
        Path,(old_lake,old_lake_receipt,history,history_receipt,out))
    require(not out.exists(),"Immutable continuous-palette source cannot be overwritten")
    old=json.loads(old_lake_receipt.read_text(encoding="utf-8"))
    evidence=json.loads(history_receipt.read_text(encoding="utf-8"))
    require(old.get("schema")==PRIOR_SCHEMA and
            old.get("source_real_original_photo_attempts")==1043930 and
            sha(old_lake)==old["continuous_trait_lake_parquet_sha256"],
            "Original 1,043,930-photo audit source has changed")
    require(evidence.get("schema")=="mois_sealed_actual_pixel_checkpoint_v1"
            and evidence.get("combined_unique_photo_attempts")==911767 and
            evidence.get("output_sha256",{}).get(HISTORY_FILE)==sha(history),
            "Historic original-photo complete terminal Parquet SHA differs")
    # Only 911k main source has missing 9-fraction columns.
    require(old["per_original_source_arm"]["historical_main"][
            "all_nine_original_fraction_columns_present"] is False,
            "Do not replace already-present frozen original continuous ROI features")
    original=pd.read_parquet(old_lake)
    require(len(original)==1043930 and original.photo_id.astype(str).is_unique,
            "Original real photo IDs not exact and unique")
    original["photo_id"]=original.photo_id.astype(str)
    photo,src=read_nine(history,source_arm="historical_checkpoint_new_full",
                       original_four_hues=int(
                          evidence["new_completed_four_hue_classified"])
                       if "new_completed_four_hue_classified" in evidence
                       else int(pd.read_parquet(history,columns=["morph"]).morph.isin(
                           ("white","yellow_orange","red_pink","blue_purple")).sum()))
    # Verify original image identity and forbid patching unrelated source arms.
    maps=original.set_index("photo_id")
    all_ids=photo.photo_id.tolist()
    require(len(photo)>0 and len(photo)<=911767 and
            set(all_ids).issubset(set(maps.index)),
            "Historical frozen new original ROI photo IDs outside measured world source")
    base=maps.loc[all_ids].reset_index()
    require(base.source_arm_original.eq(HISTORY_SOURCE).all(),
            "Historical new original ROI not from historical_main")
    require(base.original_four_hue.tolist()==photo.original_four_hue.tolist(),
            "Published exact original photo four-hue label differs between archives")
    status_old=base.original_measurement_status.fillna("").astype(str).to_numpy()
    status_new=photo.original_measurement_status.fillna("").astype(str).to_numpy()
    require(((status_old=="") | (status_new=="") |
             (status_old==status_new)).all(),
            "Original previously published ROI terminal status changed")
    require(not base.continuous9_strict.any() and
            not base.continuous9_exploratory_ambiguous.any(),
            "Existing measured colour vector must never be rewritten")

    strict=photo.continuous9_strict.to_numpy(bool)
    amb=photo.continuous9_exploratory_ambiguous.to_numpy(bool)
    selected=strict|amb
    keep_cols=[c for c in photo.columns if c.startswith("continuous9_")]
    keep_cols += ["original_nine_fraction_columns_present",
                  "nine_fraction_column_count","original_9palette_numeric_positive"]
    keep_cols=list(dict.fromkeys(keep_cols))
    loc=original.set_index("photo_id")
    selected_ids=photo.loc[selected,"photo_id"].tolist()
    for column in keep_cols:
        loc.loc[selected_ids,column]=photo.loc[selected,column].to_numpy()
    # Strict flags are encoded outside continuous9_*.
    loc.loc[selected_ids,"continuous9_strict"]=strict[selected]
    loc.loc[selected_ids,"continuous9_exploratory_ambiguous"]=amb[selected]
    loc.loc[selected_ids,"continuous9_status"]=photo.loc[selected,"continuous9_status"].to_numpy()
    # Restore raw status provenance but never mutate previously recorded four hue.
    new=loc.reset_index()
    require(len(new)==1043930 and new.photo_id.is_unique
            and new.morph.tolist()==original.morph.tolist(),
            "New full source photo ID rows or four-hue classifications changed")
    require(new.continuous9_strict.sum()==
                int(original.continuous9_strict.sum())+int(strict.sum()),
            "Restored 9-palette strict results did not reconcile")
    require(new.continuous9_exploratory_ambiguous.sum()==
                int(original.continuous9_exploratory_ambiguous.sum())+int(amb.sum()),
            "Exploratory ambiguous results did not reconcile")
    # Float values in the source are normalized from *original* ROI vectors,
    # never inferred from a colour name or environment.
    to_save=out/OUTPUT_NEW
    out.mkdir(parents=True)
    new.to_parquet(to_save,index=False,compression="zstd")
    report={
        "schema":SCHEMA,
        "original_signed_image_photo_IDs":1043930,
        "original_signed_four_hue_labels":459235,
        "historical_full_source_palette_partial_parquet_photo_IDs":len(photo),
        "historical_image_checkpoint_full_fraction_columns":src[
            "all_nine_original_fraction_columns_present"],
        "original_historical_full_roi_archive_sha256":sha(history),
        "original_historical_receipt_sha256":sha(history_receipt),
        "original_1043930_palette_lake_sha256":sha(old_lake),
        "new_strict_9composition_restored_from_historical_full_original_ROI":
            int(strict.sum()),
        "new_exploratory_ambiguous_restored":int(amb.sum()),
        "new_strict_plus_ambiguous_usable_original_palette_vectors":
            int(selected.sum()),
        "total_strict_original_9component_vectors":int(new.continuous9_strict.sum()),
        "total_ambiguous_original_9component_vectors":int(new.continuous9_exploratory_ambiguous.sum()),
        "total_real_9component_vectors_including_exploratory":int(
            (new.continuous9_strict|new.continuous9_exploratory_ambiguous).sum()),
        "remaining_historical_four_hues_without_complete9":int((
            new.source_arm_original.eq(HISTORY_SOURCE) &
            new.original_four_hue_classified &
            ~new.continuous9_strict).sum()),
        "original_photo_ROI_categories_unchanged":True,
        "technical_unknowns_imputed":0,
        "no_original_photo_bytes_downloaded_or_neural_inference":True,
        "output_full_lake_sha256":sha(to_save),
    }
    (out/"mois_1043930_historical_full9_restoration_receipt.json").write_text(
        json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("MOIS_1043930_ORIGINAL_HISTORICAL9_RESTORED_WITHOUT_IMAGE_FETCH",
          json.dumps(report,sort_keys=True),flush=True)
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("old-lake","old-receipt","historical-full","historical-receipt","output"):
        p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args()
    restore(a.old_lake,a.old_receipt,a.historical_full,a.historical_receipt,a.output)

if __name__=="__main__":
    main()
