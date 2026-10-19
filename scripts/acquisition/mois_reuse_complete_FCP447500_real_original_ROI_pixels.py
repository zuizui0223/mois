#!/usr/bin/env python3
"""Reuse all 447500 source-verified FCP original image-pixel ROI results in MOIS.

Acquisition/data engineering only. FCP's full 447500 original photo outcomes
were already measured and published; do not re-download/licence-copy images,
re-run classifier, re-use taxon summary D/W inference or re-measure overlap.
MOIS' prior 350399 original image outcomes always take precedence by photo ID.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis.mois_actual_photo_cohort_season_audit_20261019 import (
    HUES, assemble, require, sha256
)

EXPECTED_FCP=447500
EXPECTED_FCP_TAXA=4475
MOIS_BASE=350399
FCP_SCHEMA="fcp_all4730_permanent_measured_original_photo_results_v1"
FCP_COMBINE_SCHEMA="fcp_all4730_original_149900_plus_all_new_pixel_measurements_v1"
OUT_NEW="mois_FCP447500_actual_original_pixel_photo_IDs_not_previously_in350399.parquet"
OUT_ALL="mois_actual350399_plus_FCP447500_photo_ID_deduplicated_original_hues.parquet"
OUT_AUDIT="mois_FCP447500_signed_original_photo_reuse_on350399_receipt.json"
OUT_COUNTS="mois_FCP447500_photo_source_cohort_measured_colour_yields.csv"


def read_fcp_genuine(path, release_receipt, combine_receipt,
                     *, expected_rows=EXPECTED_FCP, expected_taxa=EXPECTED_FCP_TAXA):
    release=json.loads(Path(release_receipt).read_text())
    combined=json.loads(Path(combine_receipt).read_text())
    require(release.get("schema")==FCP_SCHEMA and
            release.get("genuinely_measured_photo_rows")==expected_rows and
            release.get("genuinely_measured_100photo_species")==expected_taxa and
            release.get("old_149900_original_photos_unchanged") is True and
            release.get("no_new_ecological_statistic_or_classifier") is True,
            "FCP release is not the complete authenticated actual 447500-photo ROI output")
    require(combined.get("schema")==FCP_COMBINE_SCHEMA
            and combined.get("total_original_pixel_terminal_photo_rows")==expected_rows
            and combined.get("total_original_100photo_measured_species")==expected_taxa
            and combined.get("original_149900_photo_rows_unchanged") is True
            and combined.get("new_U100_source_reuse_excluded_taxa")==3,
            "FCP complete original photo terminal 5-cohort evidence not verified")
    require(sha256(path)==release.get("source_archive_sha256")==
            combined.get("merged_photo_table_sha256") and
            sha256(combine_receipt)==release.get("source_join_receipt_sha256"),
            "The published FCP authentic original PHOTO PIXEL outcomes were altered")
    q=pd.read_csv(path,dtype={
        "photo_id":"string","inat_taxon_id":"string","observation_id":"string",
        "species":"string"},low_memory=False)
    required={"photo_id","inat_taxon_id","species","morph",
              "global_classifiable","fcp_source_cohort"}
    require(required <=set(q.columns),
            "FCP genuine source missing original photo/colour fields: "+
            repr(sorted(required-set(q.columns))))
    for key in ("photo_id","inat_taxon_id","species"):
        q[key]=q[key].fillna("").astype(str).str.strip()
    require(len(q)==expected_rows and q.photo_id.is_unique and
            q.photo_id.str.fullmatch(r"\d+").all() and
            q.inat_taxon_id.str.fullmatch(r"\d+").all() and
            q.species.ne("").all() and q.inat_taxon_id.nunique()==expected_taxa,
            "FCP independent original measured-photo IDs/taxa incomplete")
    require(q.groupby("inat_taxon_id",observed=True).size().eq(100).all(),
            "FCP original every-species 100 photo measurement contract changed")
    q["fcp_classifiable"]=q.global_classifiable.fillna(False).astype(str).str.strip().str.lower().isin(
        ("true","1","yes"))
    src=q.morph.fillna("").astype(str)
    require((~q.fcp_classifiable | src.isin(HUES)).all(),
            "FCP classified source original pixel outcome is not a recognised four-hue result")
    require(int(q.fcp_classifiable.sum())==
            release.get("genuinely_four_state_classified_photos"),
            "FCP authentically classifiable photo count mismatches signed public release")
    require(set(q.fcp_source_cohort.dropna().astype(str))==
            {"discovery","validation","prospective","wave4","new_U100"},
            "FCP original five source cohorts missing")
    return q,release,combined


def canonical_fcp(q):
    # Source CSV and the verified read_fcp_genuine frame must have identical
    # behaviour: use original FCP global_classifiable whenever a test/source
    # frame has not precomputed the verified helper boolean.
    require("global_classifiable" in q,
            "Original FCP machine ROI classifiability source field missing")
    from_source=q.global_classifiable.fillna(False).astype(str).str.strip().str.lower().isin(
        ("true","1","yes"))
    valid=q.fcp_classifiable.astype(bool) if "fcp_classifiable" in q else from_source
    require(valid.eq(from_source).all(),
            "FCP classifiability helper diverged from the original raw ROI source")
    x=pd.DataFrame()
    x["photo_id"]=q.photo_id.astype(str).to_numpy()
    x["inat_taxon_id"]=q.inat_taxon_id.astype(str).to_numpy()
    x["species"]=q.species.astype(str).to_numpy()
    x["observation_id"]=(q.observation_id.fillna("").astype(str).to_numpy()
                         if "observation_id" in q else [""]*len(q))
    for field in ("latitude","longitude"):
        x[field]=(pd.to_numeric(q[field],errors="coerce").to_numpy(float)
                  if field in q else np.full(len(q),np.nan))
    x["observed_on"]=(q.observed_on.fillna("").astype(str).to_numpy()
                      if "observed_on" in q else [""]*len(q))
    x["source_original_fcp_morph"]=q.morph.fillna("").astype(str).to_numpy()
    # Explicit source original global_classifiable, not a guess based only on
    # descriptive morphology. Retain the exact raw morph separately even if
    # invalid due to technical QC; never treat those as coloured flowers.
    x["morph"]=np.where(valid,
                        q.morph.fillna("").astype(str),
                        "mixed_uncertain")
    x["classified"]=valid.to_numpy(bool)
    x["white"]=x.morph.eq("white")
    x["arm"]="FCP_final_447500_actual_original_ROI"
    x["fcp_source_cohort"]=q.fcp_source_cohort.fillna("").astype(str).to_numpy()
    x["source_fcp_global_classifiable"]=valid.to_numpy(bool)
    x["photo_colour_provenance"]="authentic_FCP_original_ROI_v4_same_iNaturalist_source"
    x["measurement_status"]=(q.measurement_status.fillna("").astype(str).to_numpy()
                              if "measurement_status" in q else ["not_reported_by_FCP"]*len(q))
    x["technical_failure"]=x.measurement_status.eq("image_acquisition_failed")
    require(x.photo_id.is_unique and
            x.classified.eq(x.morph.isin(HUES)).all(),
            "FCP photo source normalization inferred or duplicated an actual pixel hue")
    return x


def canonical_old(old):
    s=old[["photo_id","inat_taxon_id","species","observation_id",
           "latitude","longitude","morph","classified","white","arm"]].copy()
    s["observed_on"]=old.observed_on.fillna("").astype(str) if "observed_on" in old else ""
    s["measurement_status"]=(old.measurement_status.fillna("").astype(str)
                             if "measurement_status" in old else "")
    s["technical_failure"]=s.measurement_status.eq("image_acquisition_failed")
    s["source_original_fcp_morph"]=""
    s["source_fcp_global_classifiable"]=False
    s["fcp_source_cohort"]=""
    s["photo_colour_provenance"]="original_MOIS_source_ROI_v4_prior350399"
    return s


def join(old, fcp, *, expected_old=MOIS_BASE):
    require(len(old)==expected_old and old.photo_id.is_unique,
            "MOIS previously measured original photo IDs were not frozen")
    require(fcp.photo_id.is_unique and len(fcp)>0,"FCP imported source missing real photo IDs")
    overlapping=fcp.photo_id.isin(set(old.photo_id))
    identical=old.set_index("photo_id")[["inat_taxon_id","species","morph","classified"]].join(
        fcp.loc[overlapping,["photo_id","inat_taxon_id","species","morph","classified"]].
            set_index("photo_id"), how="inner",lsuffix="_old",rsuffix="_FCP")
    inconsistent_taxon=int(identical.inat_taxon_id_old.ne(
        identical.inat_taxon_id_FCP).sum())
    inconsistent_hue=int((identical.classified_old &
                          identical.classified_FCP &
                          identical.morph_old.ne(identical.morph_FCP)).sum())
    new=fcp.loc[~overlapping].copy()
    old_base=canonical_old(old)
    combined=pd.concat([old_base,new],ignore_index=True,sort=False)
    require(len(combined)==expected_old+len(new)
            and combined.photo_id.is_unique
            and combined.iloc[:expected_old][["photo_id","morph","classified","white"]].
                reset_index(drop=True).equals(
                    old[["photo_id","morph","classified","white"]].
                        reset_index(drop=True)),
            "Existing MOIS authentic original photo RGB measurements changed")
    require(combined.classified.eq(combined.morph.isin(HUES)).all()
            and combined.white.eq(combined.morph.eq("white")).all(),
            "Reused FCP genuine ROI classifiability/morph not preserved")
    old_taxa=set(old.loc[old.classified,"inat_taxon_id"])-{""}
    added_taxa=set(new.loc[new.classified,"inat_taxon_id"])-{""}
    return combined,new,{
        "full_original_FCP_real_photo_attempts":len(fcp),
        "MOIS_previous_real_photo_attempts":len(old),
        "original_FCP_photo_IDs_already_in_MOIS":int(overlapping.sum()),
        "new_genuinely_measured_FCP_photo_IDs":len(new),
        "new_FCP_four_hue_classified":int(new.classified.sum()),
        "new_FCP_unclassifiable_original_ROIs":int((~new.classified).sum()),
        "already_measured_photo_ID_overlap_source_taxon_disagreements":inconsistent_taxon,
        "already_measured_photo_ID_overlap_both_classified_hue_disagreements":inconsistent_hue,
        "added_classifiable_nominal_taxa_not_in_old_MOIS":len(added_taxa-old_taxa),
        "union_total_genuinely_original_photo_IDs":len(combined),
        "union_total_four_hue_classified":int(combined.classified.sum()),
        "old_source_image_ROI_results_preserved":True,
        "source_photo_IDs_unambiguous":True,
        "this_is_not_a_new_provider_or_new_independent_photo_source":True,
        "no_original_images_downloaded_or_redistributed":True,
        "no_FCP_D_W_50km_or_ecology_results_imported":True,
        "unclassifiable_FCP_original_images_not_counted_as_nonwhite":True,
        "WCVP_native_or_environmental_analysis_run":False,
    }


def run(core,broad,wave2,cached8192,tail16384,
        fcp_measured,fcp_receipt,fcp_combine_receipt,out):
    out=Path(out)
    require(not out.exists() or not any(out.iterdir()),
            "Original source-pixel reuse must never overwrite a published result")
    old,groups=assemble(core,broad,wave2,cached8192,tail16384,strict=True)
    fcp,rc,original=read_fcp_genuine(fcp_measured,fcp_receipt,fcp_combine_receipt)
    combined,new,report=join(old,canonical_fcp(fcp))
    out.mkdir(parents=True,exist_ok=True)
    combined.to_parquet(out/OUT_ALL,index=False,compression="zstd")
    new.to_parquet(out/OUT_NEW,index=False,compression="zstd")
    rows=[]
    for name,s in combined.groupby("arm",observed=True):
        known=s.loc[s.classified]
        rows.append({
            "arm":name,"attempted_original_photo_ROIs":len(s),
            "four_hue_classified":len(known),
            "technically_unclassifiable_original_photo_ROIs":int((~s.classified).sum()),
            **{f"hue_{h}":int(known.morph.eq(h).sum()) for h in HUES},
            "original_nominal_taxa_with_four_hue":int(known.inat_taxon_id.nunique())
        })
    pd.DataFrame(rows).to_csv(out/OUT_COUNTS,index=False)
    report.update({
        "schema":"mois_import_source_signed_FCP447500_true_original_ROIs_on350399_v1",
        "original_FCP_full_release_schema":rc["schema"],
        "original_FCP_archive_source_photo_SHA256":sha256(fcp_measured),
        "original_FCP_public_release_receipt_SHA256":sha256(fcp_receipt),
        "original_FCP_original_cohort_join_receipt_SHA256":sha256(fcp_combine_receipt),
        "original_MOIS_source_photo_SHA256":{
            k:sha256(p) for k,p in [
                ("core",core),("broad",broad),("wave2",wave2),
                ("cached8192",cached8192),("tail16384",tail16384)]},
        "output_SHA256":{p:sha256(out/p) for p in (OUT_ALL,OUT_NEW,OUT_COUNTS)}
    })
    (out/OUT_AUDIT).write_text(json.dumps(report,sort_keys=True,
                                          ensure_ascii=False,indent=2,
                                          allow_nan=False)+"\n")
    print("MOIS_REUSE_FULL_FCP447500_REAL_ORIGINAL_ROI_PIXELS",json.dumps({
        k:report[k] for k in (
            "original_FCP_photo_IDs_already_in_MOIS",
            "new_genuinely_measured_FCP_photo_IDs",
            "new_FCP_four_hue_classified",
            "added_classifiable_nominal_taxa_not_in_old_MOIS",
            "union_total_genuinely_original_photo_IDs",
            "union_total_four_hue_classified"
        )},sort_keys=True),flush=True)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("core","broad","wave2","cached8192","tail16384",
                 "fcp-measured","fcp-receipt","fcp-combine-receipt","output"):
        p.add_argument("--"+name,type=Path,required=True)
    a=p.parse_args()
    run(a.core,a.broad,a.wave2,a.cached8192,a.tail16384,
        a.fcp_measured,a.fcp_receipt,a.fcp_combine_receipt,a.output)

if __name__=="__main__":
    main()
