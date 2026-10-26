#!/usr/bin/env python3
"""Reality check of the proposed flower-box main colour against authentic
original EfficientSAM union-ROI nine palette fractions on EXACT photo IDs.

The 256 K1 calibration was constructed as 128 previously valid 4-state
ROIs and 128 ROI/flip-gate failures; it is NOT a globally representative
accuracy sample and neither ROI method is expert-labelled petals. This
diagnostic reveals whether a cheap YOLO BOX dominant CIELAB proxy can safely
replace costly segmentation in the 93.2m first pass, without inferring
unobserved image colour from categorical flower classes.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd

SOURCE_SHA="b20bbefe15020c21c3b5b0f7edf3f8e08148a4ee7dda8bf772a89bf6b126f447"
BOX_SHA="b4233c7818a56e67934fc7375ca3ddb110f79b7180d0f7e0adcb408b02d70c5c"
SCHEMA="mois_source_frozen_K1_same_photo_9palette_vs_500px_box_Lab_proxy_diagnostic_v1"
def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(4<<20),b""):h.update(b)
    return h.hexdigest()

def check(k,why):
    if not k:raise ValueError("MOIS_K1_BOX_VS_TRUE_ROI: "+why)

def summary(d):
    pass_roi=d[d.continuous9_strict]
    fail_roi=d[~d.continuous9_strict]
    box=d[d.screen_status.eq("candidate")].copy()
    valid=box.loc[box.continuous9_strict].copy()
    valid["white9"]=pd.to_numeric(valid.continuous9_white,errors="coerce")
    valid["box_white"]=pd.to_numeric(valid.floral_box_fraction_white_like,errors="coerce")
    valid["box_green"]=pd.to_numeric(valid.floral_box_fraction_green_like,errors="coerce")
    valid["box_brown"]=pd.to_numeric(valid.floral_box_fraction_brown_like,errors="coerce")
    valid["dominant9_white"]=valid.apply(lambda z:z.continuous9_white==
          max(float(z["continuous9_"+i]) for i in (
               "white","yellow","orange","red","pink","magenta","purple","blue","bronze")),axis=1)
    # Box white-like pixels are an arbitrary Lab threshold, not true petals.
    valid["abs_white9_vs_box_white"]=np.abs(valid.white9-valid.box_white)
    valid["box_mode_white_like"]=(valid.floral_box_main_L_star.ge(80)&
                                  valid.floral_box_main_C_star.le(15))
    corr=float(valid.white9.corr(valid.box_white,method="spearman")) if len(valid)>=5 else None
    pearson=float(valid.white9.corr(valid.box_white)) if len(valid)>=5 else None
    rmse=float(np.sqrt(np.mean((valid.white9-valid.box_white)**2))) if len(valid) else None
    result={
        "full_signed_original_K1_256_images":256,
        "known_4state_photo_ROI_source_sample":128,
        "previous_ROI_flip_gate_failed_source_sample":128,
        "box_photo_decoded":int(d.image_pixel_sha256.notna().sum()),
        "box_actual_colour_assays":int(len(box)),
        "known_ROI9_box_overlap":int(len(valid)),
        "known_source_success_but_no_box":int(len(pass_roi)-len(valid)),
        "prior_ROI_gate_failed_but_box_colour_available":int(len(
            fail_roi[fail_roi.screen_status.eq("candidate")])),
        "strict_ROI9_white_fraction_vs_BOX_lab_white_like_spearman":corr,
        "strict_ROI9_white_fraction_vs_BOX_lab_white_like_pearson":pearson,
        "strict_ROI9_white_fraction_vs_BOX_lab_white_like_MAE":
            float(valid.abs_white9_vs_box_white.mean()),
        "strict_ROI9_white_fraction_vs_BOX_lab_white_like_RMSE":rmse,
        "strict_ROI9_median_white_fraction":float(valid.white9.median()),
        "strict_box_median_lab_white_like_fraction":float(valid.box_white.median()),
        "strict_box_green_like_ge25pct":int(valid.box_green.ge(.25).sum()),
        "strict_box_brown_like_ge25pct":int(valid.box_brown.ge(.25).sum()),
        "strict_roi9_dominant_white_but_box_mode_NOT_white_like":int((
            valid.dominant9_white & ~valid.box_mode_white_like).sum()),
        "strict_roi9_dominant_nonwhite_but_box_mode_white_like":int((
            ~valid.dominant9_white & valid.box_mode_white_like).sum()),
        "strict_roi9_dominant_white_total":int(valid.dominant9_white.sum()),
        "strict_roi9_largest_white_palette_but_box_green_ge25pct":int((
            valid.dominant9_white&valid.box_green.ge(.25)).sum()),
        "all_box_median_dominant_lab_mode_fraction":float(
            box.floral_box_main_mode_fraction.median()),
        "all_box_median_central_border_deltaE":float(
            box.floral_box_central_vs_border_lab_delta_e.median()),
        "source_benchmark_not_blind_human_verified_petal_annotations":True,
        "roi9_and_500px_box_are_different_inference_methods_and_pixel_mixing":True,
        "five_pixel_category_thin_proxy_not_native_green_or_brown_flora":True,
    }
    return result,valid

def run(palette,palette_receipt,box,box_summary,output):
    palette,palette_receipt,box,box_summary,output=map(
        Path,(palette,palette_receipt,box,box_summary,output))
    check(not output.exists(),"Do not overwrite prior same-original-photo comparison")
    r=json.loads(palette_receipt.read_text(encoding="utf-8"))
    b=json.loads(box_summary.read_text(encoding="utf-8"))
    check(r.get("output_photo_lake_SHA256")==SOURCE_SHA==sha(palette)
          and r.get("latest_real_photo_ID_attempts")==1043930
          and r.get("strict_original9_vectors_now")==390541,
          "The original 1m-photo source nine-palette SHA is not verified")
    check(b.get("color_feature_mode")=="floral_box_lab_v2"
          and b.get("pilot_images_attempted")==256
          and b.get("flower_box_candidate_count")==251
          and b.get("stage1_pixel_packet_SHA256")==BOX_SHA==sha(box),
          "Not the actual 256 photo K1 box Lab image-source result")
    df=pd.read_csv(box,dtype={"photo_id":str})
    original=pd.read_parquet(palette,columns=[
        "photo_id","continuous9_strict","continuous9_white",
        *["continuous9_"+i for i in ("yellow","orange","red","pink",
           "magenta","purple","blue","bronze")],
        "original_four_hue","source_arm_original"])
    original.photo_id=original.photo_id.astype(str)
    d=df.merge(original,on="photo_id",how="left",sort=False,
               validate="one_to_one")
    check(len(d)==256 and d.photo_id.is_unique and
          d.source_arm_original.eq("K1").all()
          and int(d.continuous9_strict.sum())==128
          and int(d.screen_status.eq("candidate").sum())==251,
          "Source outcome-stratified K1 subset photo IDs or ROI statuses changed")
    for column in [
        "floral_box_fraction_white_like","floral_box_fraction_green_like",
        "floral_box_fraction_brown_like","floral_box_main_L_star",
        "floral_box_main_C_star","floral_box_main_mode_fraction",
        "floral_box_central_vs_border_lab_delta_e"]:
        d[column]=pd.to_numeric(d[column],errors="coerce")
    result,valid=summary(d)
    output.mkdir(parents=True)
    audit=output/"mois_K1_source_ROI9_vs_500px_flower_box_Lab_result.json"
    result["schema"]=SCHEMA
    result["nine_colour_original_ROI_lake_SHA256"]=sha(palette)
    result["K1_photo_box_real_RGB_pixel_source_SHA256"]=sha(box)
    result["photo_bytes_downloaded_this_analytic_join"]=0
    audit.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n")
    # Audit the largest continuous-white discrepancies for a possible blind
    # PHOTO/ROI check, NEVER require human petal annotations for the full 93m.
    flags=valid[["photo_id","original_four_hue","continuous9_white",
                 "floral_box_fraction_white_like","floral_box_fraction_green_like",
                 "floral_box_fraction_brown_like","floral_box_main_L_star",
                 "floral_box_main_C_star","floral_box_central_vs_border_lab_delta_e",
                 "abs_white9_vs_box_white"]].sort_values(
                    ["abs_white9_vs_box_white","photo_id"],
                     ascending=[False,True]).head(64)
    flags.to_csv(output/"mois_K1_photo_ID_high_ROI9_vs_box_white_difference_for_review.csv",
                  index=False)
    print("MOIS_K1_SAME_REAL_PHOTO_ROI9_VERSUS_FLORAL_BOX_LAB_VALIDITY",json.dumps(
        result,sort_keys=True,allow_nan=False),flush=True)
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("palette","palette-receipt","box","box-summary","output"):
        p.add_argument("--"+k,type=Path,required=True)
    z=p.parse_args()
    run(z.palette,z.palette_receipt,z.box,z.box_summary,z.output)
if __name__=="__main__":main()
