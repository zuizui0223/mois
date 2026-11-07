#!/usr/bin/env python3
"""Integrate 128 sealed medium-500px flower-box and FULL original ROI-v4 photos.

Two different image-size acquisitions used the SAME immutable 131,072 photo
IDs. Count each original photo and observation ONCE, verify all 512 actual
FCP image terminal receipts with their source SHA and nine numeric fractions,
and join its medium YOLO-floral-box Lab status on three source identity keys.

The real full-resolution ROI9 outcome is scientifically different from the
500px box-Lab proxy. A 500px screen negative can STILL carry strict9 colour:
measure this discordance before considering any future high-recall exclusion.
No source flower-stage/native status is inferred. Neither adaptive arm is a
worldwide-flora white-frequency probability sample.
"""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib,json,math
from pathlib import Path
import numpy as np
import pandas as pd

from scripts.acquisition.mois_completed_pixel_checkpoint import collect_wave
from scripts.analysis.mois_reassemble_original_new4096_measurements import METHOD

SCHEMA="mois_bulk128_131072_one_photo_two_pixel_methods_strict9_mediumLab_integrated_v1"
SOURCE_SCHEMA="mois_original93m_tier10_bulk128_probability_reference_plus_gap_v1"
SOURCE_SHA="7e975e701281ba812451acd48be7e91a7d17f2bbe6f9f53ad15ef4f8db9e557e"
HUES=("white","yellow_orange","red_pink","blue_purple")
PALETTE=("white","yellow","orange","red","pink","magenta","purple","blue","bronze")
FRAC=tuple("flower_fraction_"+p for p in PALETTE)
REFERENCE="mois_bulk128_tier10_photo_reference_v1"
TARGETED="mois_bulk128_tier10_taxon_region_season_gap_v1"
MEDIUM_COLUMNS=("photo_id","inat_taxon_id","observation_id","design_cohort",
                "original_priority","latitude_region","local_season","screen_status",
                "medium_decoded","floral_box_main_L_star","floral_box_main_C_star",
                "floral_box_main_rgb_hex","unverified_floral_box_not_petal_ROI")
BASE_STRICT9=402859
BASE_SOURCE_TAXA_WITH_STRICT9=22787
BASE_PHOTO_ATTEMPTS=1045978+4096
ALL_SOURCE_TAXA=42111
N=131072

def require(condition,message):
    if not bool(condition):
        raise ValueError("MOIS_BULK128_DUAL_PIXEL_INTEGRATION: "+message)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for piece in iter(lambda:f.read(4<<20),b""):
            h.update(piece)
    return h.hexdigest()

def strict_nine(frame):
    """A full original ROI with nine nonnegative numeric masses AND fixed hue."""
    require({"morph","measurement_status"}|set(FRAC)<=set(frame.columns),
            "FCP original numeric palette or unknown-state columns absent")
    P=frame[list(FRAC)].apply(pd.to_numeric,errors="coerce").to_numpy(dtype=float)
    finite=np.isfinite(P).all(axis=1)
    nonneg=(np.where(np.isfinite(P),P,-1)>=0).all(axis=1)
    under=(np.where(np.isfinite(P),P,1e3)<=1+1e-6).all(axis=1)
    total=np.where(np.isfinite(P),P,0).sum(axis=1)
    clean=finite&nonneg&under&(total>0)
    fixed=frame.morph.isin(HUES).to_numpy(bool)
    failed=frame.measurement_status.eq("image_acquisition_failed").to_numpy(bool)
    good=clean&fixed&~failed
    require(not (fixed&failed).any()
            and not (good&~np.isfinite(total)).any(),
            "Source original photo acquisition failure coloured as real ROI9")
    mass=np.full(P.shape,np.nan)
    mass[clean]=P[clean]/total[clean,None]
    values=frame.copy()
    values["original_four_hue_known"]=fixed
    values["strict_original_roi9"]=good
    values["original_photo_download_failed"]=failed
    values["nine_raw_but_no_fixed_four_hue"]=clean&~fixed&~failed
    for i,pigment in enumerate(PALETTE):
        values["roi9_"+pigment]=np.where(good,mass[:,i],np.nan)
    entropy=-np.sum(np.where(np.isfinite(mass)&(mass>0),
                      mass*np.log(np.maximum(mass,1e-12)),0),axis=1)/math.log(9)
    values["roi9_entropy"]=np.where(good,entropy,np.nan)
    require(not values.loc[~values.strict_original_roi9,
                           ["roi9_"+x for x in PALETTE]].notna().any().any(),
            "Original technical no-ROI image given invented hue")
    return values

def integrate_frames(medium,original):
    """Three-key inner join requiring EXACT 131072 unique original source IDs."""
    n=len(original)
    require(n==len(medium)
            and {"photo_id","observation_id","inat_taxon_id","source_arm"}<=set(original)
            and set(MEDIUM_COLUMNS)<=set(medium),
            "Both original pixel methods must have same exact observed photo count/schema")
    for df in (medium,original):
        for k in ("photo_id","observation_id","inat_taxon_id"):
            df[k]=df[k].fillna("").astype(str)
        require(df.photo_id.is_unique and df.observation_id.is_unique
                and df.photo_id.str.fullmatch(r"\d+").all()
                and df.observation_id.ne("").all(),
                "Duplicate or invalid original image ID or observation UUID")
    m=medium.copy()
    require(m.design_cohort.isin(("reference","targeted")).all()
            and m.screen_status.isin(("candidate","negative","uncertain","error")).all()
            and m.unverified_floral_box_not_petal_ROI.astype(str).eq("True").all(),
            "500px image LAB labels are not source-verifiable flower-box proxies")
    # Preserve original large-resolution image SHA independently from medium500 SHA.
    raw=original[[
        "photo_id","observation_id","inat_taxon_id","source_arm","morph",
        "measurement_status","image_sha256",*FRAC
    ]].copy()
    original_ids=set(raw.photo_id); medium_ids=set(m.photo_id)
    require(original_ids==medium_ids,"Full-resolution FCP and medium500 use different original photo IDs")
    raw=raw.rename(columns={
        "morph":"original_four_hue",
        "measurement_status":"original_measurement_status",
        "image_sha256":"original_large_photo_image_sha256",
    })
    d=m.merge(raw,on=["photo_id","observation_id","inat_taxon_id"],
              how="inner",validate="one_to_one",sort=False)
    require(len(d)==n,"Same original photo ID had disagreeing observation or plant taxon between methods")
    expected=d.design_cohort.map({"reference":REFERENCE,"targeted":TARGETED})
    require(d.source_arm.eq(expected).all(),
            "Source reference/targeted role differs across signed actual pixel methods")
    result=strict_nine(
        d.rename(columns={"original_four_hue":"morph",
                          "original_measurement_status":"measurement_status"})
    )
    result=result.rename(columns={"morph":"original_four_hue",
                                 "measurement_status":"original_measurement_status"})
    result["medium_box_proxy_not_strict9"]=True
    result["medium_screen_negative_with_strict_original_roi9"]=(
        result.screen_status.eq("negative") & result.strict_original_roi9)
    result["medium_screen_candidate_but_no_strict_original_roi9"]=(
        result.screen_status.eq("candidate") & ~result.strict_original_roi9)
    result["medium_actual_RGB_decoded"]=result.medium_decoded.astype(str).eq("True")
    result["medium_Lab_box_measured"]=pd.to_numeric(
        result.floral_box_main_L_star,errors="coerce").notna()
    result["newly_colour_supported_taxon_in_frozen_P1P2_targeted"] = (
        result.design_cohort.eq("targeted") &
        result.original_priority.isin((
            "P1_source_rare_no_colour_1to9",
            "P2_source_has_photos_no_colour_10plus")) &
        result.strict_original_roi9)
    require(not result.loc[result.design_cohort.eq("targeted"),
                           "original_priority"].eq("P5_existing_10plus_strict9").any(),
            "P5 targeted original photo forbidden by precolour source")
    return result

def summarize_strata(d,keys):
    out=[]
    for values,g in d.groupby(keys,dropna=False,observed=True,sort=True):
        if not isinstance(values,tuple):values=(values,)
        base={key:str(val) for key,val in zip(keys,values)}
        true=g[g.strict_original_roi9]
        out.append({**base,
            "unique_original_photo_attempts":int(len(g)),
            "unique_original_observations":int(g.observation_id.nunique()),
            "nominal_source_taxa_in_stratum":int(g.inat_taxon_id.nunique()),
            "500px_images_decoded":int(g.medium_actual_RGB_decoded.sum()),
            "500px_flower_box_candidate":int(g.screen_status.eq("candidate").sum()),
            "500px_flower_box_Lab_measured":int(g.medium_Lab_box_measured.sum()),
            "full_FCP_original_fourhue_known":int(g.original_four_hue_known.sum()),
            "strict_original_roi9_real_images":int(g.strict_original_roi9.sum()),
            "large_image_acquisition_failures":int(g.original_photo_download_failed.sum()),
            "medium_screen_negative_but_strict_roi9":int(
                g.medium_screen_negative_with_strict_original_roi9.sum()),
            "medium_screen_candidate_but_no_strict_roi9":int(
                g.medium_screen_candidate_but_no_strict_original_roi9.sum()),
            "unique_P1P2_prior_no_colour_taxa_newly_strict_roi9":int(
                true.loc[true.newly_colour_supported_taxon_in_frozen_P1P2_targeted,
                         "inat_taxon_id"].nunique()),
        })
    return out

def run(source,roi_workers,medium_photo_file,medium_receipt,out,
        *,expected=N):
    source,roi_workers,medium_photo_file,medium_receipt,out=map(
        Path,(source,roi_workers,medium_photo_file,medium_receipt,out))
    require(not out.exists(),"Do not overwrite immutable original source data integration")
    src=source/"manifest.json"
    meta=json.loads(src.read_text())
    require(sha(src)==SOURCE_SHA
            and meta.get("schema")=="mois_full_source_observation_sampling_cohort_v1"
            and meta.get("sampling_schema")==SOURCE_SCHEMA
            and meta.get("new_photo_attempts_planned")==N
            and meta.get("work_units")==512
            and meta.get("worker_ids")==list(range(512))
            and meta.get("P5_targeted_photo_attempts")==0,
            "Frozen full source 128-part originals were not loaded")
    medium_doc=json.loads(medium_receipt.read_text())
    require(medium_doc.get("schema")==
            "mois_bulk128_true_original_medium500_picture_pixel_complete_audit_v1"
            and medium_doc.get("new_real_medium_picture_attempts")==expected
            and medium_doc.get("128_frozen_original_taxon_photo_partitions_authenticated")==128
            and medium_doc.get("new_original_strict9_color_ROI_count")==0
            and medium_doc.get("actual_500px_colour_status_photoID_file_SHA256")==
            sha(medium_photo_file),
            "Original photo medium500 source-id and signed colour proxy audit incompatible")
    from scripts.analysis.mois_reassemble_original_new4096_measurements import METHOD
    originals,seal=collect_wave(source,roi_workers,METHOD)
    require(seal["wave_complete"] and
            seal["attempted"]==seal["planned_photo_attempts"]==expected
            and seal["pending_photo_attempts"]==0
            and len(originals)==expected
            and len(seal["completed_units"])==512,
            "Missing SHA-signed original FCP image worker unit; never infer 9 colours from 500px")
    raw=pd.DataFrame(originals)
    original=pd.read_csv(medium_photo_file,compression="gzip",dtype=str,
                         keep_default_na=False,low_memory=False)
    d=integrate_frames(original,raw)
    require(len(d)==expected,"Joined original image pixel and box Lab photos were dropped")
    # Ecological interpretations must operate on typed strict ROI outcomes,
    # not on box proxy or technical missingness.
    source_balances=(d.groupby("design_cohort",observed=True).size().to_dict())
    require(source_balances=={"reference":32768,"targeted":98304},
            "One of two precoloured source inclusion arms mutated in image worker")
    missing=len(d)-int(d.strict_original_roi9.sum())
    gain_ids=d.loc[d.newly_colour_supported_taxon_in_frozen_P1P2_targeted,
                   "inat_taxon_id"].unique().tolist()
    out.mkdir(parents=True)
    output=out/"bulk128_original_photo_ID_500px_and_strict_roi9_ONE_ROW.parquet"
    # Drop raw worker fraction columns only AFTER retaining normalized 9 fractions.
    d=d.drop(columns=list(FRAC))
    d.to_parquet(output,index=False,compression="zstd")
    files={}
    files["by_design"]=out/"bulk128_photo_roi9_success_by_source_design.csv"
    pd.DataFrame(summarize_strata(d,["design_cohort"])).to_csv(
        files["by_design"],index=False)
    files["by_region_season"]=out/"bulk128_photo_roi9_success_by_region_and_local_season.csv"
    pd.DataFrame(summarize_strata(d,["design_cohort","latitude_region","local_season"])
                 ).to_csv(files["by_region_season"],index=False)
    files["by_priority"]=out/"bulk128_photo_roi9_by_original_taxon_gap_priority.csv"
    pd.DataFrame(summarize_strata(d,["design_cohort","original_priority"])
                 ).to_csv(files["by_priority"],index=False)
    files["medium_negative_cross_check"]=out/"bulk128_500px_screen_vs_real_original_roi9.csv"
    (pd.crosstab(d.screen_status,d.strict_original_roi9)
        .rename_axis(index="medium_screen_status",columns="strict_original9")
        .to_csv(files["medium_negative_cross_check"]))
    gfile=out/"bulk128_new_strict9_P1P2_nominal_taxa.csv"
    pd.DataFrame({"inat_taxon_id":sorted(gain_ids)}).to_csv(gfile,index=False)
    evidence={
        "schema":SCHEMA,
        "exact_93m_index_source_total":93209351,
        "true_original_photo_IDs_measured_once_by_both_methods":expected,
        "real_medium500_RGB_decoded":int(d.medium_actual_RGB_decoded.sum()),
        "real_medium500_Lab_flower_box_descriptors":int(d.medium_Lab_box_measured.sum()),
        "genuine_original_fullROI9_strict_vectors":int(d.strict_original_roi9.sum()),
        "genuine_original_fullROI9_four_hue_known":int(d.original_four_hue_known.sum()),
        "fullROI9_unknown_without_strict":missing,
        "full_original_photo_image_acquisition_failed":int(d.original_photo_download_failed.sum()),
        "500px_negative_but_fullROI9_strict":int(
            d.medium_screen_negative_with_strict_original_roi9.sum()),
        "500px_candidate_but_fullROI9_unknown":int(
            d.medium_screen_candidate_but_no_strict_original_roi9.sum()),
        "previously_no9_P1P2_targeted_taxa_with_new_strict9_from_bulk128":len(gain_ids),
        "reference_and_targeted_photo_counts":{k:int(v) for k,v in source_balances.items()},
        "reference_and_targeted_not_pooled_for_population_white_estimates":True,
        "all_strict9_from_original_full_size_not_500px_colour_proxy":True,
        "500px_box_colour_does_not_certify_actual_petals":True,
        "unclassified_and_failed_photos_kept_as_unknown":True,
        "prior_source_9palette_frames_NOT_modified":True,
        "independent_botanical_wild_native_flower_stage_verification":False,
        "new_image_downloads_or_model_inferences_by_integration":0,
        "source_metadata_manifest_SHA256":sha(src),
        "medium500_completion_receipt_SHA256":sha(medium_receipt),
        "original_512_work_unit_terminal_SHA_attestations":seal,
        "integrated_photo_ID_parquet_SHA256":sha(output),
        "source_balance_and_crossvalidation_files_SHA256":{
            p.name:sha(p) for p in sorted(files.values())
        },
        "newly_measured_prior0_taxa_list_SHA256":sha(gfile),
        # Baseline is source-frozen; other parallel main runs MAY have since
        # added valid original photos. Do not mislabel simple counts a current
        # merged main-branch global photo-ID census without a new crossarm join.
        "preexisting_signed_42111_source_taxa_with_strict9_at_bulk_selection":BASE_SOURCE_TAXA_WITH_STRICT9,
        "preexisting_signed_fullROI9_photos_at_bulk_selection":BASE_STRICT9,
        "main_branch_parallel_workers_not_reconciled_with_this_trial":True,
    }
    out_json=out/"bulk128_original_roi9_plus_medium500_integrated_receipt.json"
    out_json.write_text(json.dumps(evidence,sort_keys=True,indent=2,allow_nan=False)+"\n")
    print("MOIS_128_FULL_REAL_ORIGINAL_ROI9_AND_MEDIUM500_SOURCE_JOIN",json.dumps({
        "actual_photo_ids":expected,
        "strict_original9":evidence["genuine_original_fullROI9_strict_vectors"],
        "strict9_unknown":missing,
        "medium500_decoded":evidence["real_medium500_RGB_decoded"],
        "500px_screen_negatives_with_true_original_strict9":
            evidence["500px_negative_but_fullROI9_strict"],
        "frozen_prior0_P1P2_taxa_newly_measured":
            evidence["previously_no9_P1P2_targeted_taxa_with_new_strict9_from_bulk128"],
    },sort_keys=True),flush=True)
    return evidence

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("source","roi-workers","medium-photo-file","medium-receipt","output"):
        p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args()
    run(a.source,a.roi_workers,a.medium_photo_file,a.medium_receipt,a.output)

if __name__=="__main__":
    main()
