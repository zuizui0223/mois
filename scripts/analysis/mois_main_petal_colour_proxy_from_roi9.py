#!/usr/bin/env python3
"""Dominant floral-region colour as the practical MOIS main-petal-colour PROXY.

No organ-specific petal, calyx, anther, bract or floret segmenter is required.
One original flower-detection ROI -> dominant (maximum fraction) of its nine
frozen colour-palette components. Main colour remains separate from fraction,
runner-up and quality: it is an inferred main PETAL proxy, not ground-truth
botanical petal pigment or independent validation of the original 4 classes.

All 1,043,930 source photo IDs are retained; strict source-verified nine-vector
measurements are primary and authentic ambiguous-palette rows sensitivity only.
No arbitrary 0.70 threshold deletes a mixed/graded floral image.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

SOURCE_SHA = "b20bbefe15020c21c3b5b0f7edf3f8e08148a4ee7dda8bf772a89bf6b126f447"
SOURCE_SCHEMA = "mois_real1043930_original9_with_FCP447500_independent_source_reuse_v1"
OUTPUT_SCHEMA = "mois_original_flower_ROI9_dominant_as_main_petal_colour_proxy_v1"
PALETTE = ("white","yellow","orange","red","pink","magenta",
           "purple","blue","bronze")
COLUMNS = tuple("continuous9_"+c for c in PALETTE)
GROUPED = ("white","yellow_orange","yellow_orange","red_pink",
           "red_pink","red_pink","blue_purple","blue_purple","yellow_orange")
N_PHOTOS = 1_043_930
N_STRICT = 390_541
N_AMBIGUOUS = 36_224

def check(condition, reason):
    if not condition:
        raise ValueError("MOIS_MAIN_PETAL_PROXY: "+reason)

def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda:stream.read(4<<20),b""):
            h.update(block)
    return h.hexdigest()

def assign_main_colour(frame):
    """Return photo-ID-conserving proxy, not colour-verified botanical petals."""
    needed={"photo_id","source_arm_original","original_four_hue",
            "continuous9_strict","continuous9_exploratory_ambiguous",*COLUMNS}
    check(needed.issubset(frame.columns),"Missing original nine-colour fraction columns")
    check(frame.photo_id.astype(str).is_unique,"Repeated original photo ID")
    strict=frame.continuous9_strict.to_numpy(dtype=bool)
    ambiguous=frame.continuous9_exploratory_ambiguous.to_numpy(dtype=bool)
    check(not np.any(strict&ambiguous),"Strict and ambiguous original ROI flags overlap")
    eligible=strict|ambiguous

    # Empty vectors are ABSENT. Do not treat NaN or all-zero vectors as black/white.
    fractions=frame[list(COLUMNS)].to_numpy(dtype=np.float64)
    x=fractions[eligible]
    check(np.isfinite(x).all() and (x>=0).all()
          and (x<=1+1e-6).all() and np.allclose(x.sum(axis=1),1,atol=2e-5),
          "Source has missing/invalid normalized nine-palette ROI fractions")
    check(np.isnan(fractions[~eligible]).all(),
          "An unaccepted source photo carries an unexpected colour estimate")
    ordered=np.argsort(-x,axis=1,kind="stable")
    first=ordered[:,0]
    second=ordered[:,1]
    top=np.take_along_axis(x,first[:,None],axis=1).reshape(-1)
    runner=np.take_along_axis(x,second[:,None],axis=1).reshape(-1)
    check(np.all(top>=runner),"Main ROI palette component is not the largest")
    categories=np.array(PALETTE,dtype=object)
    coarse=np.array(GROUPED,dtype=object)
    outcome=np.full(len(frame),"",dtype=object)
    source_state=np.full(len(frame),"no_authentic_original_ROI9",dtype=object)
    evidence_state=np.where(strict,"strict_original_ROI9",
                     np.where(ambiguous,"ambiguous_ROI9_sensitivity_only",
                              "no_authentic_original_ROI9"))
    outcome[eligible]=categories[first]
    first_source=np.full(len(frame),np.nan)
    second_source=np.full(len(frame),np.nan)
    second_color=np.full(len(frame),"",dtype=object)
    margin=np.full(len(frame),np.nan)
    num_two=np.zeros(len(frame),dtype=np.int8)
    main_coarse=np.full(len(frame),"",dtype=object)
    first_source[eligible]=top
    second_source[eligible]=runner
    second_color[eligible]=categories[second]
    margin[eligible]=top-runner
    num_two[eligible]=(x>=.10).sum(axis=1).astype(np.int8)
    main_coarse[eligible]=coarse[first]
    # The 70% signal is a continuous dominance descriptor only, not admission.
    quality=np.full(len(frame),"not_measured",dtype=object)
    quality[eligible]=np.where(top>=.7,"dominant_fraction_ge_0.70",
                                  "mixed_or_low_dominance_below_0.70")
    row=pd.DataFrame({
        "photo_id":frame.photo_id.astype(str).to_numpy(),
        "source_arm":frame.source_arm_original.astype(str).to_numpy(),
        "original_four_hue":frame.original_four_hue.fillna("").astype(str).to_numpy(),
        "main_petal_colour_proxy":outcome,
        "main_colour_is_a_botanically_verified_petal":np.full(len(frame),False),
        "main_colour_provenance":evidence_state,
        "main_colour_fraction":first_source,
        "second_colour_component":second_color,
        "second_colour_fraction":second_source,
        "dominant_minus_second_fraction":margin,
        "number_of_palette_components_ge_10percent":num_two,
        "main_colour_dominance_tag":quality,
        "main_colour_coarse_4group_for_legacy_only":main_coarse,
        "strict_main_colour_primary":strict,
        "ambiguous_main_colour_sensitivity_only":ambiguous,
    })
    check(len(row)==len(frame) and row.photo_id.is_unique
          and row.loc[~eligible,"main_petal_colour_proxy"].eq("").all()
          and row.loc[~eligible,"main_colour_fraction"].isna().all(),
          "Original source rows or colour unknowns changed")
    return row

def run(source,receipt,output):
    source,receipt,output=map(Path,(source,receipt,output))
    check(not output.exists(),"Refuse overwriting prior signed main-colour source")
    source_evidence=json.loads(receipt.read_text(encoding="utf-8"))
    check(source_evidence.get("schema")==SOURCE_SCHEMA
          and source_evidence.get("latest_real_photo_ID_attempts")==N_PHOTOS
          and source_evidence.get("strict_original9_vectors_now")==N_STRICT
          and source_evidence.get("ambiguous_original9_vectors_now")==N_AMBIGUOUS
          and source_evidence.get("output_photo_lake_SHA256")==SOURCE_SHA
          and digest(source)==SOURCE_SHA,
          "Original 1m genuinely measured floral ROI9 source SHA/count mismatched")
    fields=["photo_id","source_arm_original","original_four_hue",
            "continuous9_strict","continuous9_exploratory_ambiguous",*COLUMNS]
    original=pd.read_parquet(source,columns=fields)
    check(len(original)==N_PHOTOS and original.photo_id.astype(str).is_unique,
          "Expected exactly 1,043,930 genuinely attempted original photo IDs")
    rows=assign_main_colour(original)
    check(int(rows.strict_main_colour_primary.sum())==N_STRICT
          and int(rows.ambiguous_main_colour_sensitivity_only.sum())==N_AMBIGUOUS,
          "Actual source nine-ROI main colour counts did not reconcile")
    output.mkdir(parents=True)
    target=output/"mois_main_petal_colour_proxy_1043930_original_photo_IDs.parquet"
    rows.to_parquet(target,index=False,compression="zstd")
    strict=rows[rows.strict_main_colour_primary]
    prior=original.loc[original.continuous9_strict,"original_four_hue"].astype(str)
    record={
        "schema":OUTPUT_SCHEMA,
        "source_original_photo_attempts":N_PHOTOS,
        "strict_dominant_9colour_source_photo_ROIs":N_STRICT,
        "ambiguous_dominant_9colour_exploration_only":N_AMBIGUOUS,
        "source_images_without_accepted_ROI9_colour":N_PHOTOS-N_STRICT-N_AMBIGUOUS,
        "main_colour_proxy_palette_counts_strict":{
            v:int(strict.main_petal_colour_proxy.eq(v).sum()) for v in PALETTE},
        "strict_first_palette_at_least_70_percent":int(
            strict.main_colour_fraction.ge(.7).sum()),
        "strict_first_palette_below_70_percent_STILL_RETAINED":int(
            strict.main_colour_fraction.lt(.7).sum()),
        "strict_multiple_palette_components_at_least10pct":int(
            strict.number_of_palette_components_ge_10percent.ge(2).sum()),
        "strict_original_chromatic_source_with_dominant_white":int((
            strict.original_four_hue.ne("white")&
            strict.main_petal_colour_proxy.eq("white")).sum()),
        "source_nine_fraction_sha256":SOURCE_SHA,
        "source_original_receipt_sha256":digest(receipt),
        "main_colour_photo_table_SHA256":digest(target),
        "flower_roi_required":True,
        "petal_calyx_stamen_bract_instance_separation_required":False,
        "main_petal_colour_is_operational_ROI_dominance_proxy_not_anatomical_ground_truth":True,
        "strict_palette_fraction_not_discretized_as_a_4_hue_gate":True,
        "original_four_hue_labels_unchanged":True,
        "new_images_downloaded":0,
        "new_segmentation_or_colour_inference":0,
    }
    check(sum(record["main_colour_proxy_palette_counts_strict"].values())==N_STRICT,
          "Nine dominant-palette source categories did not sum to photo count")
    (output/"mois_main_petal_colour_proxy_1043930_receipt.json").write_text(
        json.dumps(record,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("MOIS_MAIN_PETAL_PROXY_ORIGINAL_1043930_SOURCE",
          json.dumps(record,sort_keys=True),flush=True)
    return record

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--receipt",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    run(a.source,a.receipt,a.output)

if __name__=="__main__":
    main()
