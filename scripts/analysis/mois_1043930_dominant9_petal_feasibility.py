#!/usr/bin/env python3
"""Assess whether a single dominant original nine-palette ROI component suffices.

All input 1,043,930 source photographs, strict 390,541 measured nine-vectors,
and the separate 36,224 ambiguous-palette vectors are kept distinct. No new
images, guessed petal-only masks, supervised cluster names, re-inference or
four-hue class changes. Fractions describe the original *flower-detection
union ROI*, not independently checked petals.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

PALETTE=("white","yellow","orange","red","pink","magenta",
         "purple","blue","bronze")
COARSE={
    "white":"white",
    "yellow":"yellow_orange","orange":"yellow_orange","bronze":"yellow_orange",
    "red":"red_pink","pink":"red_pink","magenta":"red_pink",
    "purple":"blue_purple","blue":"blue_purple",
}
HUES=("white","yellow_orange","red_pink","blue_purple")
PINS={
    "photo_rows":1043930,
    "strict":390541,
    "ambiguous":36224,
    "known_four":459235,
    "source_sha256":"b20bbefe15020c21c3b5b0f7edf3f8e08148a4ee7dda8bf772a89bf6b126f447",
}
SCHEMA="mois_20261025_originalROI9_dominant_palette_diagnostic_v1"
FIELDS=("photo_id","original_four_hue","source_arm_original","continuous9_strict",
        "continuous9_exploratory_ambiguous",*(f"continuous9_{x}" for x in PALETTE))

def require(ok,reason):
    if not ok:
        raise ValueError("MOIS_DOMINANT_ORIGINAL_ROI9: "+reason)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(4<<20),b""):h.update(chunk)
    return h.hexdigest()

def describe(x):
    if len(x)==0:
        return None
    v=np.asarray(x,float)
    return {"n":len(v),"q10":float(np.quantile(v,.1)),
            "q25":float(np.quantile(v,.25)),
            "median":float(np.median(v)),
            "q75":float(np.quantile(v,.75)),
            "q90":float(np.quantile(v,.9)),
            "mean":float(v.mean())}

def source_summary(data,label):
    """Output aggregate ONLY: no source IDs, no observer or coordinates."""
    selected=data[data[label]].copy()
    P=selected[[f"continuous9_{x}" for x in PALETTE]].to_numpy(dtype=np.float64)
    require(P.shape==(len(selected),9)
            and np.isfinite(P).all() and (P>=0).all()
            and np.allclose(P.sum(axis=1),1,atol=2e-5,rtol=0),
            "Photo ROI9 fractions were not actual normalized source measurements")
    ordered=np.sort(P,axis=1)
    top=ordered[:,-1]
    second=ordered[:,-2]
    peak_idx=P.argmax(axis=1)
    most=np.asarray(PALETTE,dtype=object)[peak_idx]
    original=selected.original_four_hue.fillna("").astype(str).to_numpy()
    strict=label=="continuous9_strict"
    if strict:
        require(np.isin(original,HUES).all(),
                "Strict original ROI9 no longer agrees with source known-hue status")
    stats={
        "photo_count":len(selected),
        "source_mask_is_petal_only_semantic_classification":False,
        "top_one_palette_fraction_distribution":describe(top),
        "gap_largest_minus_second":describe(top-second),
        "percent_at_least_one_palette_fraction":{
            f"{threshold:.2f}":float((top>=threshold).mean()*100)
            for threshold in (.5,.6,.7,.8,.9,.95)},
        "at_least_two_palette_fractions_10pct":
            float(((P>=.1).sum(axis=1)>=2).mean()*100),
        "at_least_three_palette_fractions_10pct":
            float(((P>=.1).sum(axis=1)>=3).mean()*100),
        "dominant_palette_counts":{
            name:int(np.count_nonzero(most==name)) for name in PALETTE},
        "dominant_palette_share_percent":{
            name:float(np.mean(most==name)*100) for name in PALETTE},
        "quality_not_authenticated_as_petals":True,
    }
    if strict:
        coarse=np.asarray([COARSE[x] for x in most],dtype=object)
        agreed=coarse==original
        # Agreement is only against an existing algorithmic label, not
        # botanical accuracy or an independently heldout flower ground truth.
        stats["argmax_9_to_four_original_machine_label_agreement_percent"]=float(
            agreed.mean()*100)
        stats["original_chromatic_with_white_largest_palette_count"]=int(
            np.count_nonzero((original!="white")&(most=="white")))
        stats["original_white_with_nonwhite_largest_palette_count"]=int(
            np.count_nonzero((original=="white")&(most!="white")))
        stats["coarse_cross_tab_original_4hue_by_largest_9palette"]={
            key:{name:int(np.count_nonzero((original==key)&(most==name)))
                 for name in PALETTE}
            for key in HUES}
        stats["dominance_by_original_coarse_class"]={
            key:{"original_photo_count":int(np.sum(original==key)),
                 "dominant_fraction_q50":float(np.median(top[original==key])),
                 "at_least_0_70_dominant_percent":float(
                     (top[original==key]>=.7).mean()*100),
                 "largest_white_fraction_percent":float(
                     (most[original==key]=="white").mean()*100)}
            for key in HUES}
    return stats

def run(source,source_receipt,out):
    source,source_receipt,out=map(Path,(source,source_receipt,out))
    require(not out.exists(),"Do not overwrite previous source-only diagnostic")
    signed=json.loads(source_receipt.read_text(encoding="utf-8"))
    require(signed.get("schema")==
            "mois_real1043930_original9_with_FCP447500_independent_source_reuse_v1"
            and signed.get("latest_real_photo_ID_attempts")==PINS["photo_rows"]
            and signed.get("four_hue_original_labels_unchanged")==PINS["known_four"]
            and signed.get("strict_original9_vectors_now")==PINS["strict"]
            and signed.get("ambiguous_original9_vectors_now")==PINS["ambiguous"]
            and signed.get("output_photo_lake_SHA256")==PINS["source_sha256"]
            and sha(source)==PINS["source_sha256"],
            "Original 1,043,930-photo continuous palette source changed")
    d=pd.read_parquet(source,columns=list(FIELDS))
    require(len(d)==PINS["photo_rows"]
            and d.photo_id.astype(str).is_unique
            and int(d.continuous9_strict.sum())==PINS["strict"]
            and int(d.continuous9_exploratory_ambiguous.sum())==PINS["ambiguous"]
            and not (d.continuous9_strict&d.continuous9_exploratory_ambiguous).any(),
            "Photo ROI9 strict vs ambiguous evidence became inconsistent")
    strict=source_summary(d,"continuous9_strict")
    ambiguous=source_summary(d,"continuous9_exploratory_ambiguous")
    result={
        "schema":SCHEMA,
        "original_photo_attempts":PINS["photo_rows"],
        "original_4hue_classifiable_source":PINS["known_four"],
        "strict_nine_fraction_ROI_photos":PINS["strict"],
        "ambiguous_nine_fraction_ROI_photos_exploratory_only":PINS["ambiguous"],
        "strict_source_measurements":strict,
        "ambiguity_sensitivity_only":ambiguous,
        "site_or_environment_predictors_used_to_select_source":False,
        "no_photo_downloads":True,
        "no_image_ROI_inferences":True,
        "major_palette_is_not_an_independently_verified_major_petal_colour":True,
        "photographic_colour_response_should_be_reestimated_with_petal_masks":True,
        "original_photo_sha256":sha(source),
    }
    out.mkdir(parents=True)
    dest=out/"real1043930_originalROI9_dominant_palette_diagnostic.json"
    dest.write_text(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+"\n")
    rows=[]
    for key,n in strict["dominant_palette_counts"].items():
        rows.append({"component":key,"n_dominant_original_ROI9":n,
                     "fraction_percent":strict["dominant_palette_share_percent"][key]})
    pd.DataFrame(rows).to_csv(
        out/"real1043930_originalROI9_dominant_component_counts.csv",index=False)
    print("MOIS_REAL1043930_ORIGINALROI9_DOMINANCE_RESULT",json.dumps({
        "strict_photos":PINS["strict"],
        "ambiguous_photos":PINS["ambiguous"],
        "strict_dominant_distribution":strict["top_one_palette_fraction_distribution"],
        "pct_strict_dominant":strict["percent_at_least_one_palette_fraction"],
        "pct_strict_multiple_components_10pct":strict[
            "at_least_two_palette_fractions_10pct"],
        "original_chromatic_with_white_largest_palette_count":strict[
            "original_chromatic_with_white_largest_palette_count"],
        "original_white_with_nonwhite_largest_palette_count":strict[
            "original_white_with_nonwhite_largest_palette_count"],
        "original_4hue_machine_agreement_only":strict[
            "argmax_9_to_four_original_machine_label_agreement_percent"],
        "dominant_palette_counts":strict["dominant_palette_counts"],
    },sort_keys=True),flush=True)
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",required=True,type=Path)
    p.add_argument("--source-receipt",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    run(a.source,a.source_receipt,a.output)
if __name__=="__main__":main()
