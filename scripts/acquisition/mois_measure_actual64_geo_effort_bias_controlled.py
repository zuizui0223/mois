#!/usr/bin/env python3
"""Exactly 64 genuine NEW original photo pixels: conservative bias-aware Stage1.

Only SHA-frozen original photo IDs selected in advance from independent original
93,209,351-photo source via source tier07, among taxa with >=100 source photos,
0 historical strict ROI9 values, and 1-4 genuine prior photo attempts.

NO biological petal truth is asserted: the warm YOLO first flower box gives
an uncalibrated device-sRGB CIELAB proxy (Lab, chroma and circular hue) plus
white/green/brown-like box pixel fractions. Every access/failure/screen-negative
is recorded and joined to predeclared north/south/tropical × prior-effort cells.
The 64 outcome-conditioned photos are NOT a population-unbiased global sample;
full-frame species, photo, observation or geography selection probabilities are
NOT inferred from this convenience/outcome-conditioned eligibility restriction.

Provider pace <=2 image requests/sec; any 403/429/503 aborts, not retried
indefinitely. No raw photo bytes are stored or published.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,json,math
from collections import Counter
from pathlib import Path
import pandas as pd
from scripts.acquisition.mois_full93m_medium_pixel_scan import (
    WarmFlowerBoxDetector, scan_source, download_medium, Pacer,
    COLUMNS_LAB,LAB_SCHEMA,sha_file
)
from scripts.acquisition.mois_select_new64_tier07_never_strict_100plus import (
    SCHEMA as SOURCE_SCHEMA,STRATA,EFFORT
)

SCHEMA="mois_20261026_genuine_new64_source_tier07_bias_controlled_photo_srgbLab_v1"
SOURCE_SHA="3debdd9353e423cf21b5f768f8aead9bb7587fcdc5f761ac7e381746d4099b0b"
PROVENANCE_SHA="c8b9be91bb27edbe66b52725e34fba0eab8dd5ec5617bca1b99bdaf893ab36ed"
SOURCE_PHOTOS=64
BASELINE_IMAGES=1043930
BASELINE_STRICT=400710
MODEL_SHA="f1aaeec4664fe2c178e5cf2bc1f508977bef3e4aa7b40613026cb8ae3de789d5"
MODEL_CLASS="item"
TIER=7

def require(ok,msg):
    if not ok:raise ValueError("MOIS_NEW64_REAL_PIXELS: "+msg)

def evaluate_inputs(source,receipt,provenance):
    source,receipt,provenance=map(Path,(source,receipt,provenance))
    doc=json.loads(receipt.read_text(encoding="utf-8"))
    require(doc.get("schema")==SOURCE_SCHEMA
            and doc.get("full_93m_parent_photo_ids")==93209351
            and doc.get("previous_9colour_success_source")==BASELINE_STRICT
            and doc.get("all_1043930_previously_attempted_original_photo_ids_excluded")==BASELINE_IMAGES
            and doc.get("actual_new_original_photo_IDs_SELECTED_not_yet_image_downloaded")==SOURCE_PHOTOS
            and doc.get("source_tier")==TIER
            and doc.get("all_16_original_tier07_source_tar_archives_SHA_verified") is True
            and doc.get("one_photo_one_source_observation_one_taxon_in_minimal_pilot") is True
            and doc.get("image_bytes_newly_downloaded")==0
            and doc.get("photo_source_SHA256")==SOURCE_SHA==sha_file(source)
            and doc.get("provenance_file_SHA256")==PROVENANCE_SHA==sha_file(provenance),
            "Source is not the predeclared 64 image/observation ID immutable selection")
    s=pd.read_csv(source,compression="gzip",dtype={"photo_id":str,
                     "inat_taxon_id":str,"observation_uuid":str})
    p=pd.read_csv(provenance,dtype={"photo_id":str,"taxon_id":str,
                     "observation_id":str})
    require(len(s)==len(p)==SOURCE_PHOTOS and
            s.photo_id.is_unique and s.observation_uuid.is_unique
            and s.inat_taxon_id.is_unique and p.photo_id.is_unique
            and p.observation_id.is_unique and p.taxon_id.is_unique
            and s.photo_id.tolist()==p.photo_id.tolist()
            and s.inat_taxon_id.tolist()==p.taxon_id.tolist()
            and s.observation_uuid.tolist()==p.observation_id.tolist(),
            "Exact original precolour source photo, taxon or observation IDs drifted")
    require(s.tier.astype(int).eq(TIER).all()
            and s.frozen_loader_supported.eq("true").all()
            and s.image_extension.isin(["jpg","jpeg","png"]).all()
            and s.photo_id.str.fullmatch(r"\d+").all(),
            "Original source image is unsupported, malformed or wrong tier")
    require(p.original_fullframe_taxon_photo_opportunities.astype(int).ge(100).all()
            and p.prior_valid9_colour_ROIs.astype(int).eq(0).all()
            and p.prior_photo_attempts.astype(int).between(1,4).all()
            and p.actual_original_photoID_already_attempted.eq(False).all()
            and p.original_observation_already_attempted.eq(False).all(),
            "Selected pilot not all novel from >=100-photographed zero-ROI taxa")
    crosses=p.groupby(["source_hemisphere_band","sampling_prior_attempt_class"]).size()
    expected={(geo,eff) for geo in STRATA for eff in EFFORT}
    require(set(crosses.index.to_list())==expected
            and crosses.eq(8).all(),
            "Original image pilot must retain eight outcome-blind selection slots per stratum")
    for col in ("source_tier_presentation_probability",
                "conditional_species_probability_given_eligible_source_tier",
                "conditional_one_observation_probability_given_eligible_taxon"):
        v=pd.to_numeric(p[col],errors="coerce")
        require(v.notna().all() and v.between(0,1,inclusive="right").all(),
                "Source strata conditional selection probability invalid: "+col)
    require(p.source_tier_presentation_probability.between(.012799,.012801).all(),
            "Original tier07 independent observation inclusion probability altered")
    return doc,s,p

def describe_scan(frame):
    out=[]
    groups=[("pilot_all",frame)]
    groups.extend(("%s|%s"%(x[0],x[1]),g) for x,g in frame.groupby(
        ["source_hemisphere_band","sampling_prior_attempt_class"],observed=True))
    groups.extend(("latitude:"+str(k),g) for k,g in frame.groupby(
        "source_hemisphere_band",observed=True))
    for label,data in groups:
        statuses=data.screen_status.fillna("unknown")
        usable=data[data.screen_status.eq("candidate") &
                    pd.to_numeric(data.floral_box_main_L_star,
                                  errors="coerce").notna()].copy()
        means={}
        for f in ("floral_box_main_L_star","floral_box_main_C_star",
                  "floral_box_main_mode_fraction","floral_box_second_mode_fraction",
                  "floral_box_fraction_white_like","floral_box_fraction_green_like",
                  "floral_box_fraction_brown_like","floral_box_central_vs_border_lab_delta_e"):
            v=pd.to_numeric(usable[f],errors="coerce")
            means["median_"+f]=float(v.median()) if v.notna().any() else None
        out.append({
            "design_stratum":label,
            "original_source_photo_IDs_selected":len(data),
            "photo_RGB_successfully_downloaded":int(data.image_pixel_sha256.notna().sum()),
            "flower_boxes_with_actual_Lab_proxy":len(usable),
            "YOLO_flower_box_screen_negatives":int(statuses.eq("negative").sum()),
            "RGB_access_or_model_errors":int(statuses.isin(["error","uncertain"]).sum()),
            "nominal_white_like_BOX_at_least25pct":int(
                pd.to_numeric(usable.floral_box_fraction_white_like,errors="coerce").ge(.25).sum()),
            "nominal_green_like_BOX_at_least25pct":int(
                pd.to_numeric(usable.floral_box_fraction_green_like,errors="coerce").ge(.25).sum()),
            "nominal_brown_like_BOX_at_least25pct":int(
                pd.to_numeric(usable.floral_box_fraction_brown_like,errors="coerce").ge(.25).sum()),
            "no_true_flower_petal_color_ground_truth":True,
            "not_global_photo_colour_prevalence":True,
            **means
        })
    return out

def run(source,receipt,provenance,weights,output,*,rps=2.0,workers=8):
    source,receipt,provenance,weights,output=map(Path,(
        source,receipt,provenance,weights,output))
    require(not output.exists(),"Original new64 real-pixel benchmark cannot overwrite previous receipt")
    require(0<rps<=2.0 and 1<=workers<=8,
            "Minimal photo acquisition paced to <=2 requests/sec and <=8 workers")
    proof,planned,design=evaluate_inputs(source,receipt,provenance)
    require(sha_file(weights)==MODEL_SHA,"Frozen official FCP YOLO11 original-image model drift")
    detector=WarmFlowerBoxDetector(weights,MODEL_SHA,0,MODEL_CLASS,
                                  device="cpu",conf=.02)
    pacer=Pacer(rps)
    pixels=output/"true_500px_photo_pixels"
    receipt_scan=scan_source(
        source,pixels,expected_source_sha256=SOURCE_SHA,
        expected_count=SOURCE_PHOTOS,part=0,tier=TIER,
        detector=detector,download=lambda url:download_medium(url,pacer),
        workers=workers,batch_size=SOURCE_PHOTOS,max_packets=1,
        source_role="bias_controlled_pilot",
        feature_mode="floral_box_lab_v2")
    require(receipt_scan["actual_unique_source_photo_ids_attempted"]==SOURCE_PHOTOS
            and receipt_scan["source_tier"]==TIER
            and receipt_scan["schema"]==LAB_SCHEMA
            and receipt_scan["source_gzip_sha256"]==SOURCE_SHA
            and receipt_scan["all_original_input_photos_attempted"] is False
            and receipt_scan["outcome_conditioned_supplement_not_population_unbiased"] is True,
            "Real measured 64 source photo pixels cannot be forged as full93m scan")
    packet=pixels/"packet-0000000.csv.gz"
    p=pd.read_csv(packet,compression="gzip",dtype={"photo_id":str})
    require(len(p)==SOURCE_PHOTOS and p.photo_id.is_unique
            and p.photo_id.tolist()==planned.photo_id.tolist()
            and p.source_photo_inclusion_probability.isna().all()
            and p.source_arm.eq("bias_controlled_tier07_missing9_nonrepresentative_pilot").all(),
            "Source photo IDs/order/status or false global propensity in pixel terminal")
    ndecoded=int(p.image_pixel_sha256.notna().sum())
    require(ndecoded==receipt_scan["medium_image_pixels_decoded"],
            "Real photo download count differs from authentic receipt")
    p=p.merge(design,on="photo_id",how="left",validate="one_to_one",
              sort=False,suffixes=("", "_source_design"))
    require(len(p)==SOURCE_PHOTOS and p.taxon_id.notna().all()
            and p.inat_taxon_id.astype(str).tolist()==p.taxon_id.astype(str).tolist()
            and p.observation_id.astype(str).tolist()==p.observation_id_source_design.astype(str).tolist(),
            "Selected source taxon or observed original UUID changed in image inference")
    # The original packet retains observation_id; precolour provenance uses
    # observation_id_source_design and must match every individual photo.
    p.to_csv(output/"mois_actual64_photoID_original_RGB_Lab_and_sampling.csv.gz",
             index=False,compression="gzip")
    stats=describe_scan(p)
    summaries=pd.DataFrame(stats)
    summaries.to_csv(output/"mois_actual64_photo_colour_by_bias_stratum.csv",index=False)
    counts={k:int(v) for k,v in p.screen_status.value_counts().items()}
    results={
        "schema":SCHEMA,
        "original_full_frame_93m_photo_opportunities":93209351,
        "source_42111_taxa_20260927":42111,
        "old_1043930_genuine_original_photo_IDs_all_excluded":True,
        "old_400710_real_roi9_colours_all_preserved":True,
        "new_original_source_photo_attempts":SOURCE_PHOTOS,
        "new_original_source_observations":int(p.observation_id.nunique()),
        "new_original_source_nominal_taxa":int(p.taxon_id.nunique()),
        "source_photo_IDs_ever_measured_before":0,
        "source_observation_IDs_ever_measured_before":0,
        "source_distinct_4geo_x_2history_bias_cells":8,
        "source_each_cell_photo_targets":8,
        "genuine_500px_medium_RGB_decoded":ndecoded,
        "screen_status_original_64_photos":counts,
        "real_flower_box_source_RGB_Lab_photo_count":int(
            (p.screen_status.eq("candidate")&
             pd.to_numeric(p.floral_box_main_L_star,errors="coerce").notna()).sum()),
        "genuine_original_nine_fraction_full_roi_reinferred":0,
        "no_unknown_photo_declared_negative_or_nonwhite":True,
        "technical_unknown_source_photo_ids_preserved":True,
        "source_and_selector_receipt_sha256":sha_file(receipt),
        "source_selection_provenance_sha256":PROVENANCE_SHA,
        "source_photo_ID_csv_gz_sha256":SOURCE_SHA,
        "FCP_frozen_flower_YOLO11_model_sha256":MODEL_SHA,
        "original_64_sRGB_pixels_packet_sha256":sha_file(packet),
        "new_photo_stratified_pixel_csv_sha256":
            sha_file(output/"mois_actual64_photoID_original_RGB_Lab_and_sampling.csv.gz"),
        "real_photo_CIELAB_L_star_C_star_hue_are_device_sRGB_BOX_proxy":True,
        "leaf_bract_calyx_background_contamination_may_be_large":True,
        "ROI_quality_sampling_not_global_flora_probability_sample":True,
        "precolour_design_conditioned_on_zero_previous_strict9_and_1to4_attempts":True,
        "source_support_N_within_each_geo_previous_effort_stratum":
            proof["geographic_and_previous_effort_strata"],
        "descriptive_strata_all_candidate_and_unknown_in_denominator":stats,
        "new_RGB_images_redistributed_as_bytes":0,
    }
    target=output/"mois_actual64_original_photo_pixel_bias_aware_receipt.json"
    target.write_text(json.dumps(results,sort_keys=True,indent=2,allow_nan=False)+"\n")
    print("MOIS_ACTUAL64_NEW_SOURCE_PIXELS_GEO_EFFORT_BIAS_CONTROLLED",json.dumps(
        {k:v for k,v in results.items() if k not in (
           "descriptive_strata_all_candidate_and_unknown_in_denominator",
           "source_support_N_within_each_geo_previous_effort_stratum")},
        sort_keys=True,allow_nan=False),flush=True)
    return results

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("source","receipt","provenance","weights","output"):
        p.add_argument("--"+name,required=True,type=Path)
    p.add_argument("--rps",default=2.0,type=float)
    p.add_argument("--workers",default=8,type=int)
    a=p.parse_args()
    run(a.source,a.receipt,a.provenance,a.weights,a.output,
        rps=a.rps,workers=a.workers)

if __name__=="__main__":
    main()
