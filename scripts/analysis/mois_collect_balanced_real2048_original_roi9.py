#!/usr/bin/env python3
"""Collect ACTUALLY MEASURED 2048 independent original photo ROI9 outcomes.

The 256 reference and 1792 hemisphere/latitude/season targeted photos were
selected from complete signed source tier07 BEFORE RGB inference. Use frozen
FCP YOLO11+EfficientSAM ROI-v4 and retain source status for ALL photo attempts,
especially technical unknowns. Do NOT merge reference and targeted components
as though they were one global probability sample.

No re-inference or image transfer here: read ONLY original sealed photo pixels
receipt and frozen terminal photo fraction fields, source manifest and source
exact photo-by-observation sampling provenance. The old 1,043,930 photos and
their 400,710 authentic nine vectors are never changed by collection.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from pathlib import Path
import numpy as np
import pandas as pd

from scripts.acquisition.mois_completed_pixel_checkpoint import collect_wave
from scripts.analysis.mois_reassemble_original_new4096_measurements import METHOD
from scripts.acquisition.mois_balance1043930_tier07_2048_source import (
    SOURCE_ARMS,SCHEMA as SOURCE_SCHEMA,
    OLD_PHOTO_SHA,OLD_OBS_SHA,ZONE,SEASONS
)
SCHEMA="mois_20261026_actual2048_novel_geo_season_balance_original_FCP_ROI9_v1"
PALETTE=("white","yellow","orange","red","pink","magenta","purple","blue","bronze")
FRACTIONS=("flower_fraction_"+name for name in PALETTE)
FRACTIONS=tuple(FRACTIONS)
HUES=("white","yellow_orange","red_pink","blue_purple")
ATTEMPTS=2048
BASELINE_STRICT=400710

def require(ok,msg):
    if not ok:raise ValueError("MOIS_COLLECT_2048_BALANCED_ROI9: "+msg)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(4<<20),b""):h.update(b)
    return h.hexdigest()

def normalize_source_design_provenance(p):
    """Translate the signed tier07 `design_cohort` source field once.

    Pixel-worker rows do not contain this pre-RGB source allocation. Keep it
    coupled to the signed photo-ID provenance, never deduce cohort from hue.
    """
    require({"photo_id", "inat_taxon_id", "design_cohort"} <= set(p.columns)
            and "sampling_design_cohort" not in p.columns,
            "Signed original photo provenance design_cohort header changed")
    result = p.rename(columns={
        "inat_taxon_id": "taxon_id_provenance",
        "design_cohort": "sampling_design_cohort",
    }).copy()
    require(result.sampling_design_cohort.isin(("reference", "targeted")).all(),
            "Source provenance contains unrecognized or missing cohort labels")
    return result


def source_design_identity_matches_pixels(d):
    """Validate pre-colour source allocation against exact original image arm."""
    require({"source_arm", "sampling_design_cohort"} <= set(d.columns),
            "Signed sampling cohort lost during source-photo ID join")
    return bool(d.sampling_design_cohort.isin(("reference", "targeted")).all()
                and ((d.source_arm == SOURCE_ARMS["reference"]) ==
                     d.sampling_design_cohort.eq("reference")).all())


def summarize(d):
    results=[]
    for (role,zone,season),g in d.groupby(
         ["sampling_design_cohort","latitude_region","recorded_local_season"],
          observed=True,dropna=False,sort=True):
        good=g[g.continuous9_strict]
        white=good.continuous9_white
        n=len(g)
        results.append({
            "design_cohort":str(role),
            "latitude_zone":str(zone),
            "recorded_local_season":str(season),
            "genuinely_attempted_original_photos":n,
            "distinct_photo_taxa":int(g.inat_taxon_id.nunique()),
            "distinct_photo_observations":int(g.observation_id.nunique()),
            "original_four_hue_classified":int(g.original_four_hue.isin(HUES).sum()),
            "new_authentic_strict9_ROIs":len(good),
            "strict9_Roi_success_percent":100*len(good)/n if n else None,
            "strict_roi9_mean_white_fraction":float(white.mean()) if len(white) else None,
            "new_image_acquisition_failures":int(g.image_acquisition_failed.sum()),
            "new_nonclassifiable_image_Rois":int((
                ~g.original_four_hue.isin(HUES) & ~g.image_acquisition_failed).sum()),
            "all_photo_attempts_kept_even_when_ROI_fails":True,
        })
    return results

def collect(source,workers,out):
    source,workers,out=map(Path,(source,workers,out))
    require(not out.exists(),"No overwrite of new-source sealed original photo colours")
    meta=json.loads((source/"manifest.json").read_text())
    require(meta.get("schema")=="mois_full_source_observation_sampling_cohort_v1"
            and meta.get("sampling_schema")==SOURCE_SCHEMA
            and meta.get("existing_strict9_actual_original_roi_SHA256")==OLD_PHOTO_SHA
            and meta.get("all_previous_photo_attempts_sha256")==OLD_OBS_SHA
            and meta.get("new_photo_attempts_planned")==ATTEMPTS
            and meta.get("work_units")==8
            and meta.get("worker_ids")==list(range(8))
            and meta.get("reference_uniform_source_and_targeted_stratum_independently_separated")
            and meta.get("source_new_image_IDs_and_observations_unique_and_prior_disjoint"),
            "Not original 2048 signed source age/season/geography control experiment")
    rows,seal=collect_wave(source,workers,METHOD)
    require(seal["wave_complete"] and
            seal["attempted"]==seal["planned_photo_attempts"]==ATTEMPTS
            and seal["pending_photo_attempts"]==0
            and len(rows)==ATTEMPTS,
            "Original 2048 photo source terminal run not fully complete")
    d=pd.DataFrame(rows)
    d.photo_id=d.photo_id.astype(str)
    require(d.photo_id.is_unique
            and d.observation_id.astype(str).is_unique
            and len(d.source_arm.unique())==2,
            "Original-photo or biological original-observation source duplicate")
    expected=dict(meta["source_arm_counts"])
    require(d.source_arm.value_counts().to_dict()==expected,
            "Original source role of 2048 images was altered during photograph inference")
    provenance=source/"selected_original_photo_design_provenance.csv"
    require(meta["original_sampling_provenance_SHA256"]==sha(provenance),
            "Bias-balanced original observation selection probability proof mutated")
    p=pd.read_csv(provenance,dtype={"photo_id":str,
              "original_observation_id":str,"inat_taxon_id":str})
    require(len(p)==ATTEMPTS and p.photo_id.is_unique
            and p.original_observation_id.is_unique
            and set(p.photo_id)==set(d.photo_id),
            "Image inference lost or fabricated original source photo IDs")
    p=p[["photo_id","original_observation_id","inat_taxon_id",
         "design_cohort","latitude_region","recorded_local_season",
         "source_5degree_grid_cell",
         "source_original_tier07_probability_width",
         "conditional_stage2_observation_and_photo_probability"]].copy()
    p=normalize_source_design_provenance(p)
    d=d.merge(p,on="photo_id",validate="one_to_one",how="left",sort=False)
    require((d.inat_taxon_id.astype(str)==d.taxon_id_provenance.astype(str)).all()
            and (d.observation_id.astype(str)==
                  d.original_observation_id.astype(str)).all()
            and source_design_identity_matches_pixels(d),
            "Same exact original photo ID must carry exact sampling taxon/observation provenance")
    require(set(FRACTIONS).issubset(d.columns),
            "Source-original YOLO/EfficientSAM terminal did not record measured numeric nine fractions")
    P=d[list(FRACTIONS)].apply(pd.to_numeric,errors="coerce").to_numpy(dtype=float)
    # Technical empty/invalid ROI9 fractions are UNKNOWN, including cases
    # where source old four-hue classification exists without 9 raw fractions.
    finite=np.isfinite(P).all(axis=1)
    nonnegative=(np.nan_to_num(P,nan=-1)>=0).all(axis=1)
    positive=np.isfinite(P).all(axis=1)&(P.sum(axis=1)>0)
    raw_ok=finite&nonnegative&positive&(np.nanmax(
            np.where(np.isfinite(P),P,0),axis=1)<=1+1e-6)
    known=d.morph.isin(HUES).to_numpy(dtype=bool)
    strict=raw_ok&known
    normalized=np.full_like(P,np.nan)
    normalized[raw_ok]=P[raw_ok]/P[raw_ok].sum(axis=1,keepdims=True)
    d["continuous9_strict"]=strict
    # Ambiguous palette can be held for a sensitivity analysis ONLY; do not
    # assign a strict result to technical flip/flower failures.
    d["continuous9_exploratory_ambiguous"]=(raw_ok&~known&
        d.measurement_status.eq(
          "not_evaluable_ambiguous_palette_composition").to_numpy(bool))
    for j,name in enumerate(PALETTE):
        d["continuous9_"+name]=np.where(strict,normalized[:,j],np.nan)
    entropy=-np.sum(np.where(np.isfinite(normalized)&(normalized>0),
                 np.where(np.isfinite(normalized),normalized,0)*
                 np.log(np.maximum(np.where(np.isfinite(normalized),normalized,0),1e-12)),
                 0),axis=1)/math.log(9)
    d["continuous9_entropy"]=np.where(strict,entropy,np.nan)
    d["original_four_hue"]=d.morph.fillna("").astype(str)
    d["image_acquisition_failed"]=d.measurement_status.eq("image_acquisition_failed")
    require(int(d.continuous9_strict.sum())<=int(known.sum())
            and all((d.loc[d.continuous9_strict,FRACTIONS].apply(
                 pd.to_numeric,errors="coerce").to_numpy(dtype=float)>=0).flat)
            and not d.loc[~d.continuous9_strict,
                          ["continuous9_"+x for x in PALETTE]].notna().any().any(),
            "Strict original ROI9 outcome invalid or unknown photo assigned a colour")
    retain=[
      "photo_id","measurement_id","inat_taxon_id","species",
      "observation_id","observer_id","observed_on","latitude","longitude",
      "positional_accuracy","photo_license","source_arm",
      "selection_probability","image_sha256",
      "original_four_hue","measurement_status",
      "continuous9_strict","continuous9_exploratory_ambiguous",
      *(f"continuous9_{x}" for x in PALETTE),
      "continuous9_entropy",
      "sampling_design_cohort","latitude_region","recorded_local_season",
      "source_5degree_grid_cell",
      "source_original_tier07_probability_width",
      "conditional_stage2_observation_and_photo_probability",
      "image_acquisition_failed",
    ]
    out.mkdir(parents=True)
    result=out/"mois_true2048_biasbalanced_original_photo_ROI9_outcomes.parquet"
    d[retain].to_parquet(result,index=False,compression="zstd")
    summary=summarize(d)
    bystratum=out/"mois_true2048_biasbalanced_ROI9_success_by_region_season.csv"
    pd.DataFrame(summary).to_csv(bystratum,index=False)
    detail={
       "schema":SCHEMA,
       "old_original_photo_IDs_unchanged":1043930,
       "old_valid_original_9fraction_ROIs_unchanged":BASELINE_STRICT,
       "new_independent_original_photo_ids_attempted":ATTEMPTS,
       "new_original_observations_distinct":d.observation_id.nunique(),
       "new_separate_reference_original_photos":int(
           d.sampling_design_cohort.eq("reference").sum()),
       "new_separate_targeted_original_photos":int(
           d.sampling_design_cohort.eq("targeted").sum()),
       "new_genuine_original_photo_image_bytes_downloaded":int(
           seal["downloaded"]),
       "new_frozen_original_roi_four_hue_classified":int(known.sum()),
       "new_strict_original9_vectors":int(strict.sum()),
       "new_explicit_ambiguous9_sensitivity_only":int(
           d.continuous9_exploratory_ambiguous.sum()),
       "new_image_acquisition_failures":int(d.image_acquisition_failed.sum()),
       "new_original_photo_unknown_hue_without_strict9":int((~strict).sum()),
       "source_original_roi_method_SHA_attested":METHOD,
       "old_prior_400710_plus_new_strict9_simple_source_asset_total_NOT_equal_global_prevalence":
              int(BASELINE_STRICT+strict.sum()),
       "selection_depends_on_OLD_roi_success_probability_but_not_new_hue":True,
       "two_reference_and_targeted_design_arms_not_naively_poolable":True,
       "new_photo_dominance_without_petal_specific_ground_truth":True,
       "new_photo_inclusion_probabilities_are_conditional_not_unconditional_full93m":True,
       "unsuccessful_original_photo_ROIs_preserved_as_unknown":True,
       "all_new_photos_were_original_image_RGB_measured":True,
       "original_photo_source_manifest_SHA256":sha(source/"manifest.json"),
       "original_photo_design_provenance_SHA256":sha(provenance),
       "new_actual_2048_outcomes_SHA256":sha(result),
       "per_region_season_ROI9_evidence_SHA256":sha(bystratum),
       "all_worker_method_receipts":seal,
       "new_photo_IDs_do_not_overlap_prior_1043930_by_source_ID_attestation":True,
       "new_any_original_photo_bytes_redistributed":False,
    }
    (out/"mois_true2048_biasbalanced_original_pixel_collection_receipt.json").write_text(
        json.dumps(detail,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("MOIS_GENUINE2048_REGIONAL_SEASONAL_BIAS_REPAIR_ROI9_PIXEL_RESULT",
          json.dumps({
             "new_photo_attempts":ATTEMPTS,
             "new_original_four_hue":int(known.sum()),
             "new_valid_original9":int(strict.sum()),
             "new_acquisition_failures":int(d.image_acquisition_failed.sum()),
             "reference_valid9":int(d.loc[d.sampling_design_cohort.eq("reference"),
                                       "continuous9_strict"].sum()),
             "targeted_valid9":int(d.loc[d.sampling_design_cohort.eq("targeted"),
                                      "continuous9_strict"].sum()),
             "regions":[{"cohort":str(role),"zone":str(z),
                         "n":int(len(sub)),"strict9":int(sub.continuous9_strict.sum())}
                       for (role,z),sub in d.groupby(
                           ["sampling_design_cohort","latitude_region"],observed=True)],
             "actual_pixel_source_SHA256":sha(result),
          },sort_keys=True,allow_nan=False),flush=True)
    return detail

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("source","workers","output"):
        p.add_argument("--"+k,required=True,type=Path)
    x=p.parse_args()
    collect(x.source,x.workers,x.output)

if __name__=="__main__":
    main()
