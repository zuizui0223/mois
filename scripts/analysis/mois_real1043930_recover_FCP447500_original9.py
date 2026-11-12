#!/usr/bin/env python3
"""No-image additional continuous9 restoration from 447,500 original FCP ROIs.

The 1,043,930-photo MOIS four-hue source includes ALL original FCP photo IDs,
but earlier import intentionally collapsed its palette fractions into four
labels. The authentic FCP 447,500-photo original CSV still stores ALL nine
continuous photo ROI fractions, source SHA and global-classifiable status.
Restore ONLY previously missing continuous vectors. Strict source photo-ID
AND hue agreement is required; conflicts and ROI unknowns remain unresolved.
"""
from __future__ import annotations
import argparse,json,gzip,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
from scripts.analysis.mois_k1_nine_palette_cluster_audit import PALETTE,FRACTIONS,AMBIGUOUS
from scripts.analysis.mois_real1043930_continuous9_source_lake import PROVENANCE,sha,require

SCHEMA="mois_real1043930_original9_with_FCP447500_independent_source_reuse_v1"
OUTPUT="mois_real1043930_photo_continuous9_with_FCP_original_reuse.parquet"
FCP_EXPECTED=447500
FCP_CLASSIFIED=227107

def restore(prior_lake,prior_receipt,fcp_csv,fcp_receipt,out):
    prior_lake,prior_receipt,fcp_csv,fcp_receipt,out=map(Path,
        (prior_lake,prior_receipt,fcp_csv,fcp_receipt,out))
    require(not out.exists(),"Cannot replace old restored original photo trait data")
    doc=json.loads(prior_receipt.read_text(encoding="utf-8"))
    signed=json.loads(fcp_receipt.read_text(encoding="utf-8"))
    require(doc.get("schema")=="mois_latest1043930_recovered_history_original_checkpoint_nine_palette_v1"
            and doc.get("original_signed_image_photo_IDs")==1043930
            and doc.get("original_signed_four_hue_labels")==459235
            and doc.get("output_full_lake_sha256")==sha(prior_lake),
            "The recovered historical source 1043930 nine-palette trait lake changed")
    require(signed.get("schema")=="fcp_all4730_permanent_measured_original_photo_results_v1"
            and signed.get("genuinely_measured_photo_rows")==FCP_EXPECTED
            and signed.get("genuinely_four_state_classified_photos")==FCP_CLASSIFIED
            and signed.get("source_archive_sha256")==sha(fcp_csv),
            "Not the authentic full 447500 FCP original ROI source CSV")
    original=pd.read_parquet(prior_lake)
    original.photo_id=original.photo_id.astype(str)
    require(len(original)==1043930 and original.photo_id.is_unique
            and int(original.continuous9_strict.sum())==
                doc["total_strict_original_9component_vectors"],
            "Historical already-restored 1m-photo nine palette does not align")
    src_fields=set(pd.read_csv(fcp_csv,compression="gzip",nrows=0).columns)
    required={"photo_id","morph","global_classifiable",*FRACTIONS}
    require(required<=src_fields,
            "FCP original archived ROIs are missing true nine channels")
    cols=["photo_id","morph","global_classifiable",*FRACTIONS]
    if "measurement_status" in src_fields:cols.append("measurement_status")
    fcp=pd.read_csv(fcp_csv,compression="gzip",usecols=cols,
                    dtype={"photo_id":str},low_memory=False)
    require(len(fcp)==FCP_EXPECTED and fcp.photo_id.is_unique
            and fcp.photo_id.str.fullmatch(r"\d+").all(),
            "The real FCP original photo ID source is incomplete or repeated")
    measured=fcp.global_classifiable.astype(str).str.strip().str.lower().isin(
        ("true","1","yes"))
    require(int(measured.sum())==FCP_CLASSIFIED,"Frozen FCP global classifiability changed")
    current=original.set_index("photo_id")
    require(set(fcp.photo_id).issubset(set(current.index)),
            "Original FCP ROI photo ID absent from MOIS globally measured archive")
    exact=current.loc[fcp.photo_id.tolist()]
    require(exact.source_arm_original.eq("historical_main").all(),
            "Original FCP real pixel result appeared in a new/double-counted source arm")
    raw=np.column_stack([pd.to_numeric(fcp[c],errors="coerce").to_numpy(dtype=float)
                         for c in FRACTIONS])
    good=(np.isfinite(raw).all(axis=1)&(raw>=0).all(axis=1)&
          (raw<=1+1e-6).all(axis=1)&(raw.sum(axis=1)>0))
    require(np.all(good[measured.to_numpy(bool)]),
            "Published FCP globally classifiable photo has no real nine-part palette")
    ref_hue=exact.original_four_hue.to_numpy(str)
    original_hue=fcp.morph.fillna("").astype(str).to_numpy(str)
    fcp_success=measured.to_numpy(bool)
    label_agreement=(fcp_success & (ref_hue==original_hue) &
                     (exact.original_four_hue_classified.to_numpy(bool)))
    already=exact.continuous9_strict.to_numpy(bool)
    patch=label_agreement & ~already
    fcp_conflicts=fcp_success & (ref_hue!=original_hue)
    # A separate source provides independently computed ROI feature vectors,
    # not proof they are the same numerical mask. Do not overwrite previously
    # stored ROI fractions for any source with the same photo ID.
    status=(fcp["measurement_status"].fillna("").astype(str).to_numpy(str)
            if "measurement_status" in fcp else np.repeat("",len(fcp)))
    exploratory=((~fcp_success)&(ref_hue=="mixed_uncertain")&
                 (status==AMBIGUOUS)&good&
                 ~exact.continuous9_exploratory_ambiguous.to_numpy(bool)&
                 ~already)
    eligible=patch|exploratory
    ids=fcp.loc[eligible,"photo_id"].tolist()
    ratios=raw[eligible]/raw[eligible].sum(axis=1,keepdims=True)
    nearest=pd.Index(ids)
    for j,k in enumerate(PALETTE):
        current.loc[nearest,"continuous9_"+k]=ratios[:,j].astype("float32")
    entropy=(-(np.where(ratios>0,ratios*np.log(np.maximum(ratios,1e-15)),0)
               ).sum(axis=1)/np.log(len(PALETTE)))
    gini=1-(ratios*ratios).sum(axis=1)
    current.loc[nearest,"continuous9_entropy"]=entropy
    current.loc[nearest,"continuous9_gini_simpson"]=gini
    current.loc[nearest,"continuous9_dominant_fraction"]=ratios.max(axis=1)
    current.loc[nearest,"continuous9_strict"]=patch[eligible]
    current.loc[nearest,"continuous9_exploratory_ambiguous"]=exploratory[eligible]
    current.loc[nearest,"continuous9_status"]=np.where(
        patch[eligible],PROVENANCE["strict"],PROVENANCE["ambiguous"])
    current.loc[nearest,"original_nine_fraction_columns_present"]=True
    current.loc[nearest,"nine_fraction_column_count"]=9
    current.loc[nearest,"original_9palette_numeric_positive"]=True
    result=current.reset_index()
    require(len(result)==1043930 and result.photo_id.is_unique
            and result.original_four_hue.tolist()==original.original_four_hue.tolist()
            and int(result.continuous9_strict.sum())==
                int(original.continuous9_strict.sum())+int(patch.sum()),
            "Never modify the frozen original four-color source or 9-palette continuity")
    require(int(result.continuous9_exploratory_ambiguous.sum())==
                int(original.continuous9_exploratory_ambiguous.sum())+
                int(exploratory.sum()),"No ambiguous palette candidates may be fabricated")
    out.mkdir(parents=True)
    target=out/OUTPUT
    result.to_parquet(target,compression="zstd",index=False)
    summary={
        "schema":SCHEMA,"latest_real_photo_ID_attempts":1043930,
        "four_hue_original_labels_unchanged":459235,
        "FCP_original_CSV_sha256":sha(fcp_csv),
        "FCP_signed_original_archive_receipt_sha256":sha(fcp_receipt),
        "FCP_original_ROI_photo_IDs":len(fcp),
        "FCP_source_globally_four_hue_classifiable":int(fcp_success.sum()),
        "FCP_same_photo_id_four_hue_label_agreement":int(label_agreement.sum()),
        "FCP_same_photo_id_four_hue_label_conflict":int(fcp_conflicts.sum()),
        "FCP_photo_vectors_already_present_preserved":int((label_agreement&already).sum()),
        "FCP_new_strict_original_photo9_vectors_recovered":int(patch.sum()),
        "FCP_new_ambiguous_nonclassifiable_photo9_vectors":int(exploratory.sum()),
        "strict_original9_vectors_now":int(result.continuous9_strict.sum()),
        "ambiguous_original9_vectors_now":int(result.continuous9_exploratory_ambiguous.sum()),
        "strict_and_exploratory_total_photo9_now":int(
            (result.continuous9_strict|result.continuous9_exploratory_ambiguous).sum()),
        "all_photo_IDs_retained_including_unknowns":True,
        "same_photo_classification_conflicts_do_not_become_colour_features":True,
        "unusable_ROI_gate_photo_numbers_not_relabelled":True,
        "new_original_RGB_images_downloaded":0,
        "new_flower_detection_or_EfficientSAM_inference":0,
        "all_original_source_photo9_from_immutable_previous_FCP_raw_pixel_Roi":True,
        "output_photo_lake_SHA256":sha(target)
    }
    (out/"mois1043930_FCP447500_nine_palette_restoration_receipt.json").write_text(
        json.dumps(summary,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("MOIS_1043930_FCP447500_ORIGINAL9_RECOVERY",json.dumps(summary,sort_keys=True),flush=True)
    return summary

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for f in ("prior-lake","prior-receipt","fcp-csv","fcp-receipt","output"):
        p.add_argument("--"+f,required=True,type=Path)
    a=p.parse_args()
    restore(a.prior_lake,a.prior_receipt,a.fcp_csv,a.fcp_receipt,a.output)
if __name__=="__main__":
    main()
