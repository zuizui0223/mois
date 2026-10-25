#!/usr/bin/env python3
"""REAL signed 256 K1 source images: observed 251 CIELAB floral BOX colours
by Arctic/South/season/moisture/dryness/elevation/island source gradients.

This is a deliberately balanced technical pilot: 128 previously valid ROIs
and 128 prior ROI failures, not a probability sample of flora, not 251
botanically verified petals or flowers. Raw green-like/brown-like box pixels
could be leaves, bracts, inflorescence structure or backgrounds.

No biology inferred for group sizes <10; hue is circular and unavailable when
colour is near achromatic. Never infer CIELAB from 9 palette fractions.
"""
from __future__ import annotations
import argparse,csv,gzip,json,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
from scripts.analysis.mois_real1043930_global_zonal_colour import (
    latitude_zone,season_label,elev_band,q5,
    CLIMATE_SHA,sha,check
)

SCHEMA="mois_real_K1_256_251_flower_box_CIELAB_all_ecological_axes_source_diagnostic_v1"
N=256
CANDIDATE=251
SOURCE_SHA="b4233c7818a6e67934fc7375ca3ddb110f79b7180d0f7e0adcb408b02d70c5c"
LAB_PROTOCOL="mois_flower_box_srgb_cielab_mode_v1"
LAB_FIELDS=(
    "floral_box_main_L_star",
    "floral_box_main_a_star","floral_box_main_b_star",
    "floral_box_main_C_star","floral_box_main_hue_degrees",
    "floral_box_main_mode_fraction","floral_box_second_mode_fraction",
    "floral_box_L50","floral_box_C50",
    "floral_box_fraction_white_like",
    "floral_box_fraction_green_like",
    "floral_box_fraction_brown_like",
    "floral_box_central_vs_border_lab_delta_e",
)

def circular_stats(angles,weights=None):
    z=np.asarray(angles,dtype=float)
    mask=np.isfinite(z)
    w=np.ones(len(z),dtype=float) if weights is None else np.asarray(weights,dtype=float)
    mask &=np.isfinite(w)&(w>0)
    if not mask.any():return {"hue_circular_mean_degrees":None,"hue_resultant_length":None}
    a=np.deg2rad(z[mask]);w=w[mask]
    re=np.sum(w*np.cos(a))/np.sum(w)
    im=np.sum(w*np.sin(a))/np.sum(w)
    length=float(np.hypot(re,im))
    return {"hue_circular_mean_degrees":
                float(np.degrees(np.arctan2(im,re))%360) if length>=1e-3 else None,
            "hue_resultant_length":length}

def med(x):
    z=np.asarray(x,dtype=float);z=z[np.isfinite(z)]
    return float(np.median(z)) if len(z) else None

def show(p,key,group):
    size=len(group)
    chosen=group[group.screen_status.eq("candidate")&
                 np.isfinite(group.floral_box_main_L_star)&
                 np.isfinite(group.floral_box_main_C_star)]
    nc=len(chosen)
    r={
      "axis":key,"stratum":str(p),
      "signed_original_K1_photo_ID_attempts_in_pilot":size,
      "colour_candidate_box_images":nc,
      "actual_Lab_original_pixel_assays":nc,
      "source_coordinate_data_provenance":"K1 original 256-photo deliberate old-ROI-half-success/half-failure",
      "N_10_for_minimum_descriptive_support":bool(nc>=10),
      "median_main_L_star":med(chosen.floral_box_main_L_star),
      "median_main_C_star":med(chosen.floral_box_main_C_star),
      "median_main_mode_fraction":med(chosen.floral_box_main_mode_fraction),
      "median_second_mode_fraction":med(chosen.floral_box_second_mode_fraction),
      "median_box_white_like_fraction":med(chosen.floral_box_fraction_white_like),
      "median_box_green_like_fraction":med(chosen.floral_box_fraction_green_like),
      "median_box_brown_like_fraction":med(chosen.floral_box_fraction_brown_like),
      "median_box_centre_boundary_deltaE":med(chosen.floral_box_central_vs_border_lab_delta_e),
      "n_box_green_fraction_ge_25pct":int(
          (chosen.floral_box_fraction_green_like>=.25).sum()),
      "n_box_brown_fraction_ge_25pct":int(
          (chosen.floral_box_fraction_brown_like>=.25).sum()),
      "n_box_white_fraction_ge_25pct":int(
          (chosen.floral_box_fraction_white_like>=.25).sum()),
      "n_main_hue_defined":int(chosen.floral_box_main_hue_degrees.notna().sum()),
      "box_colour_NOT_verified_petal":True,
    }
    r.update(circular_stats(chosen.floral_box_main_hue_degrees.to_numpy(),
                            chosen.floral_box_main_C_star.to_numpy()))
    return r

def run(pixel_dir,climate,geography,geo_receipt,output):
    pixel_dir,climate,geography,geo_receipt,output=map(
        Path,(pixel_dir,climate,geography,geo_receipt,output))
    check(not output.exists(),"Cannot rewrite prior source K1 real Lab pilot ecology")
    summary=json.loads((pixel_dir/"K1_medium_pixel_benchmark_summary.json").read_text())
    pkt=pixel_dir/"packet-0000000.csv.gz"
    proof=json.loads((pixel_dir/"packet-0000000.receipt.json").read_text())
    check(summary.get("schema")=="mois_K1_256_actually_decoded_medium_pixel_benchmark_v1"
          and summary.get("pilot_images_attempted")==N
          and summary.get("color_feature_mode")=="floral_box_lab_v2"
          and summary.get("flower_box_candidate_count")==CANDIDATE
          and summary.get("stage1_pixel_packet_SHA256")==SOURCE_SHA==sha(pkt)
          and proof["schema"]=="mois_all_original_photo_medium_pixel_stage1_CIELAB_dominant_v2"
          and proof["CIELAB_box_source_protocol"]==LAB_PROTOCOL
          and proof["output_sha256"]==sha(pkt),
          "Not already executed actual 256 source K1 photographs with true pixel Lab")
    frame=pd.read_csv(pkt,compression="gzip",dtype={"photo_id":str})
    check(len(frame)==N and frame.photo_id.is_unique
          and int(frame.screen_status.eq("candidate").sum())==CANDIDATE,
          "Not 256 signed original photos or 251 actual scanned flower boxes")
    for field in LAB_FIELDS:
        frame[field]=pd.to_numeric(frame[field],errors="coerce")
    original_c=frame[frame.screen_status.eq("candidate")]
    check(original_c.floral_box_main_L_star.notna().all()
          and original_c.floral_box_main_C_star.notna().all(),
          "At least one 256 original 500px source image lacks real Lab")
    check(int(original_c.floral_box_fraction_green_like.ge(.25).sum())==92
          and int(original_c.floral_box_fraction_brown_like.ge(.25).sum())==33
          and int(original_c.floral_box_fraction_white_like.ge(.25).sum())==37,
          "Independent source-pixel dominant Lab colour QC counts changed")
    p=pd.read_parquet(climate,columns=[
       "photo_id","wc2_bio12_10m","wc2_annual_vapr_10m","wc2_elev_2p5m"])
    g=pd.read_parquet(geography,columns=["photo_id","land_type"])
    receipt=json.loads(geo_receipt.read_text(encoding="utf-8"))
    check(len(p)==len(g)==1043930 and
          p.photo_id.astype(str).tolist()==g.photo_id.astype(str).tolist()
          and sha(climate)==CLIMATE_SHA
          and sha(geography)==receipt["output_geography_photo_ID_parquet_SHA256"],
          "Pilot true photo IDs not being joined to previously verified 1m climate/island")
    p.photo_id=p.photo_id.astype(str)
    g.photo_id=g.photo_id.astype(str)
    full=frame.merge(p,on="photo_id",how="left",validate="one_to_one",
                     sort=False).merge(g,on="photo_id",how="left",
                          validate="one_to_one",sort=False)
    check(len(full)==N and full.photo_id.tolist()==frame.photo_id.tolist()
          and full.land_type.notna().all(),
          "Authentic source K1 image climate/island photo IDs not joined")
    full["hemisphere"]=np.where(full.latitude>=0,"north",
                        np.where(full.latitude<0,"south","unknown"))
    full["latitude_zone"]=[latitude_zone(x) for x in full.latitude]
    month=pd.to_datetime(full.observed_on,errors="coerce").dt.month
    full["season"]=[season_label(a,b) for a,b in zip(
        full.latitude.to_numpy(float),month.to_numpy(float))]
    full["altitude_band"]=[elev_band(x) for x in full.wc2_elev_2p5m]
    full["vapour_quintile"],vapr_all_cut=q5(p.wc2_annual_vapr_10m.to_numpy(float))
    thresholds=np.asarray(vapr_all_cut,dtype=float)
    values=full.wc2_annual_vapr_10m.to_numpy(float)
    full["vapour_quintile"]=np.where(np.isfinite(values),
        np.asarray(["Q"+str(i) for i in np.searchsorted(thresholds,
                      np.where(np.isfinite(values),values,0),side="right")+1],
                   dtype=object),"unknown")
    _,rain_cut=q5(p.wc2_bio12_10m.to_numpy(float))
    values=full.wc2_bio12_10m.to_numpy(float)
    full["rainfall_quintile"]=np.where(np.isfinite(values),
        np.asarray(["Q"+str(i) for i in np.searchsorted(np.asarray(rain_cut),
                      np.where(np.isfinite(values),values,0),side="right")+1],
                   dtype=object),"unknown")
    groups=[
       ("all_actual_K1_256",[]),
       ("hemisphere",["hemisphere"]),
       ("latitude_zone",["latitude_zone"]),
       ("hemisphere_season",["hemisphere","season"]),
       ("latitude_zone_season",["latitude_zone","season"]),
       ("island_mainland",["land_type"]),
       ("elevation",["altitude_band"]),
       ("atmospheric_vapour_quintile",["vapour_quintile"]),
       ("annual_rainfall_quintile",["rainfall_quintile"]),
       ("hemisphere_island",["hemisphere","land_type"]),
       ("hemisphere_elevation",["hemisphere","altitude_band"]),
    ]
    summaries=[]
    for label,cols in groups:
        if not cols:
            summaries.append(show(full,label,"all"))
        else:
            for k,v in full.groupby(cols,dropna=False,observed=True):
                if not isinstance(k,tuple):k=(k,)
                summaries.append(show("|".join(map(str,k)),label,v))
    out={
       "schema":SCHEMA,
       "real_original_photo_attempts":N,
       "actual_CIELAB_flower_box_photos":CANDIDATE,
       "original_ROI_K1_outcome_balanced_128_success_128_failure":True,
       "sample_not_random_global_or_island_color_study":True,
       "source_lab_real_photo_packet_sha256":SOURCE_SHA,
       "annual_climate_photo_SHA256":sha(climate),
       "source_GSHHG_photo_sha256":sha(geography),
       "pilot_max5_argmax_photo_groups_not_a_taxon_equal_macro_estimate":True,
       "lab_grain":"500px YOLO detected flower BOX, possibly includes leaves or bracts",
       "lab_is_sRGB_D65_device_photo_colour_not_color_chart_calibrated":True,
       "lab_green_brown_are_box_pixel_traits_not_botanical_petal_species":True,
       "hue_circular_and_achromatic_hue_missing":True,
       "all_axes_descriptive_lab_pilot":summaries,
       "photo_ID_source_image_downloads_in_this_POSTmeasurement_analysis":0,
    }
    output.mkdir(parents=True)
    sfile=output/"mois_real256_flower_box_Lab_all_axes_descriptive.csv"
    pd.DataFrame(summaries).to_csv(sfile,index=False)
    res=output/"mois_real256_flower_box_Lab_all_axes_result.json"
    res.write_text(json.dumps(out,sort_keys=True,indent=2,allow_nan=False)+"\n")
    print("MOIS_K1_REAL_256_LAB_WHITE_GREEN_BROWN_GLOBAL_AXES_PILOT",json.dumps({
        "attempted":N,"Lab":CANDIDATE,"strata":len(summaries),
        "hemispheres":[v for v in summaries if v["axis"]=="hemisphere"],
        "island_mainland":[v for v in summaries if v["axis"]=="island_mainland"],
        "polar":[v for v in summaries if v["axis"]=="latitude_zone"
                 and "polar" in v["stratum"]],
        "source_file_SHA256":sha(sfile),
    },sort_keys=True,allow_nan=False),flush=True)
    return out

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("pixel-dir","climate","geography","geo-receipt","output"):
        p.add_argument("--"+key,required=True,type=Path)
    x=p.parse_args()
    run(x.pixel_dir,x.climate,x.geography,x.geo_receipt,x.output)
if __name__=="__main__":main()
