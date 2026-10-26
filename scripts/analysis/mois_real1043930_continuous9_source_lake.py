#!/usr/bin/env python3
"""Recover measured 9-component colour compositions across ALL 1,043,930
MOIS original image-photo outcomes without downloading or re-inferring photos.

Follow the original six historical/K1/U100 source arms and 35 immutable FCP5
groups. The *four-hue* union is the signed photo-ID/accounting spine only:
the union deliberately discarded photo-level continuous fractions. Re-read
the SHA-attested source ROI Parquets, not that four-hue projection.

The 9 frozen "flower_fraction_*" values are palette proportions, NOT petal
CIELAB, UV reflectance, chemical anthocyanin or verified open-flower pixels.
Only source-verified four-hue-classified originals get strict continuous
composition; explicit ambiguous-palette failures may be included separately
in an exploratory sensitivity dataset. All other technical unknowns remain
unknown, even if a numeric palette vector happens to survive in a failed ROI.
"""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from scripts.analysis.mois_latest1043930_photo_union import (
    ARMS, HUES, SOURCE_N, BASE_HUES, DEPTH5_N, BASE_N, load_one, sha,
    CUM_SCHEMA, GROUP_SCHEMA,
)
from scripts.analysis.mois_k1_nine_palette_cluster_audit import (
    PALETTE, FRACTIONS, CLASSIFIED, AMBIGUOUS,
)

SCHEMA="mois_latest1043930_source_verified_photo_roi_continuous9_v1"
OUTPUT="mois_real1043930_original_photo_continuous9_trait_lake.parquet"
RECEIPT="mois_real1043930_continuous9_coverage_receipt.json"
K1_ID="K1"
EXPECTED_K1_KNOWN=27689
EXPECTED_K1_AMBIGUOUS_VALID=2510
PROVENANCE={"strict":"four_hue_original_ROI_with_valid_9composition",
            "ambiguous":"explicit_ambiguous_palette_original_ROI_unverified",
            "other":"no_usable_original_flower_colour_ROI"}

def require(cond, msg):
    if not cond:
        raise ValueError("MOIS_CONTINUOUS9_SOURCE: "+msg)

def read_nine(path, *, source_arm, original_four_hues):
    """Read only original photo identity/status and 9 fraction columns.

    A genuine archive may lack these old columns. Treat this as missing
    *source measurements*, never fabricate colour from four class labels.
    """
    source=Path(path)
    pf=pq.ParquetFile(source)
    fields=set(pf.schema.names)
    cols=["photo_id","morph"]
    optional=["measurement_status","image_sha256","flower_effective_pixels",
              "nuisance_pixel_fraction"]
    cols += [s for s in optional if s in fields]
    available=[s for s in FRACTIONS if s in fields]
    cols+=available
    require({"photo_id","morph"} <= fields,
            source_arm+": source no longer has original photo ID/hue")
    chunks=[]
    for arrow in pf.iter_batches(batch_size=65536,columns=cols):
        d=arrow.to_pandas()
        d["photo_id"]=d.photo_id.fillna("").astype(str)
        d["morph"]=d.morph.fillna("").astype(str)
        require(d.photo_id.str.fullmatch(r"\d+").all(),
                source_arm+": invalid original numeric photo ID")
        known=d.morph.isin(HUES).to_numpy(bool)
        states=(d.measurement_status.fillna("").astype(str).to_numpy()
                if "measurement_status" in d else
                np.repeat("",len(d)))
        recorded=(states!="")
        require(np.array_equal((states[recorded]==CLASSIFIED),known[recorded]),
                source_arm+": original explicit status contradicts 4-hue label")
        if len(available)==len(FRACTIONS):
            raw=np.column_stack(
                [pd.to_numeric(d[f],errors="coerce").to_numpy(dtype=float)
                 for f in FRACTIONS])
            finite=np.isfinite(raw).all(axis=1)
            nonneg=(raw>=0).all(axis=1)
            bounded=(raw<=1+1e-6).all(axis=1)
            mass=np.where(finite[:,None],raw,0).sum(axis=1)
            valid=finite&nonneg&bounded&(mass>0)
            palette=np.full((len(d),9),np.nan,dtype=np.float32)
            palette[valid]=(raw[valid]/mass[valid,None]).astype(np.float32)
        else:
            valid=np.zeros(len(d),dtype=bool)
            palette=np.full((len(d),9),np.nan,dtype=np.float32)
        # Do not turn an image-access failure into a coloured flower.
        technical_failure=states=="image_acquisition_failed"
        require(not np.any(valid&technical_failure),
                source_arm+": image acquisition failure carries forged ROI palette")
        strict=known&valid
        ambiguous=(~known)&(states==AMBIGUOUS)&valid
        provenance=np.full(len(d),PROVENANCE["other"],dtype=object)
        provenance[strict]=PROVENANCE["strict"]
        provenance[ambiguous]=PROVENANCE["ambiguous"]
        permitted=strict|ambiguous
        palette[~permitted,:]=np.nan
        denom=np.nansum(palette,axis=1)
        # The empty/missing denominator never becomes an outcome.
        entropy=np.full(len(d),np.nan,dtype=np.float32)
        complexity=np.full(len(d),np.nan,dtype=np.float32)
        largest=np.full(len(d),np.nan,dtype=np.float32)
        if permitted.any():
            mix=palette[permitted].astype(float)
            entropy[permitted]=(-np.where(
                mix>0,mix*np.log(np.maximum(mix,1e-15)),0).sum(axis=1)
                /np.log(len(PALETTE))).astype(np.float32)
            complexity[permitted]=(1-np.square(mix).sum(axis=1)).astype(np.float32)
            largest[permitted]=mix.max(axis=1).astype(np.float32)
        result=pd.DataFrame({
            "photo_id":d.photo_id,
            "source_arm_original":source_arm,
            "original_four_hue":d.morph.where(known,"mixed_uncertain"),
            "original_four_hue_classified":known,
            "original_measurement_status":states,
            "original_nine_fraction_columns_present":len(available)==9,
            "nine_fraction_column_count":len(available),
            "original_9palette_numeric_positive":valid,
            "continuous9_status":provenance,
            "continuous9_strict":strict,
            "continuous9_exploratory_ambiguous":ambiguous,
            "continuous9_entropy":entropy,
            "continuous9_gini_simpson":complexity,
            "continuous9_dominant_fraction":largest,
        })
        for i,name in enumerate(PALETTE):
            result["continuous9_"+name]=palette[:,i]
        chunks.append(result)
    data=pd.concat(chunks,ignore_index=True)
    require(data.photo_id.is_unique and
            int(data.original_four_hue_classified.sum())==original_four_hues,
            source_arm+": original classified hue and photo-ID conservation failed")
    return data,{
        "source_arm":source_arm,
        "original_photo_attempts":len(data),
        "source_four_hue_labels":int(data.original_four_hue_classified.sum()),
        "all_nine_original_fraction_columns_present":len(available)==9,
        "available_fraction_columns":available,
        "strict_original_4hue_ROIs_with_all_9_fractions":int(data.continuous9_strict.sum()),
        "ambiguous_explicit_palette_usable_exploratory":int(data.continuous9_exploratory_ambiguous.sum()),
        "original_photos_with_positive_all9_numerics":int(data.original_9palette_numeric_positive.sum()),
        "missing_strict_9vectors_even_if_4hue_known":int(
            (data.original_four_hue_classified&~data.continuous9_strict).sum()),
        "technical_9fraction_invalid_or_untrusted_photo_attempts":int(
            (data.continuous9_status==PROVENANCE["other"]).sum()),
        "source_parquet_sha256":sha(source),
    }

def create(source, latest, latest_receipt, output, *,
           verify_latest=True):
    source,latest,latest_receipt,output=map(
        Path,(source,latest,latest_receipt,output))
    require(not output.exists(),"Refuse rewriting immutable trait lake")
    master=json.loads(latest_receipt.read_text(encoding="utf-8"))
    require(master.get("source_pixel_attempts")==SOURCE_N
            and master.get("source_four_hue_classified")==459235
            and master.get("photo_parquet_sha256")==sha(latest),
            "Latest full original-photo union is not source-authorized")
    original=pd.read_parquet(latest,columns=["photo_id","morph","arm","date",
                                             "inat_taxon_id","species","latitude","longitude"])
    original.photo_id=original.photo_id.astype(str)
    require(len(original)==SOURCE_N and original.photo_id.is_unique,
            "Latest measured photo-ID universe not exactly complete")
    # Authenticate each original source independently, not just the flattened
    # latest four-hue projection.
    frame=[];proof=[]
    for arm in ARMS:
        ref,receipt,_,hue=ARMS[arm]
        expected,stamp=load_one(source,arm)
        part,audit=read_nine(source/ref,source_arm=arm,original_four_hues=hue)
        require(len(part)==len(expected) and
                part.photo_id.tolist()==expected.photo_id.tolist() and
                part.original_four_hue.tolist()==expected.morph.tolist(),
                "Continuous ROI source not photo aligned with SHA-verified "+arm)
        frame.append(part);proof.append({**audit,"receipt_SHA256":stamp["receipt_SHA256"]})
    cumulative=source/"fcp5"/"cumulative"/"mois_FCP5_partial_genuine_original_pixel_cumulative.json"
    cum=json.loads(cumulative.read_text(encoding="utf-8"))
    require(cum["schema"]==CUM_SCHEMA and
            cum["new_actual_finished_original_photo_groups"]==35 and
            cum["new_genuine_original_photo_pixel_attempts"]==DEPTH5_N and
            cum["previous_actual_independently_verified_original_photo_attempts"]==BASE_N,
            "The entire 35-group original RGB source did not finish")
    groups={int(r["group"]):r for r in cum["group_receipts"]}
    require(set(groups)==set(range(35)),"Missing frozen true-photo FCP5 group")
    for group in range(35):
        arm=f"FCP5_group_{group:02d}"
        folder=source/"fcp5"/("group-%02d"%group)
        path=folder/"original_photo_ROIs.parquet"
        receipt_path=folder/"original_group_receipt.json"
        record=json.loads(receipt_path.read_text(encoding="utf-8"))
        require(record["schema"]==GROUP_SCHEMA and
                int(record["original_worker_group"])==group and
                record["saved_file_SHA256"]["original_photo_ROIs.parquet"]==sha(path) and
                sha(receipt_path)==groups[group]["source_group_receipt_SHA256"],
                "FCP5 source photo pixels/receipt SHA no longer match")
        part,audit=read_nine(path,source_arm=arm,
                            original_four_hues=int(record["new_actual_four_hue_classified"]))
        require(len(part)==int(record["new_actual_original_photo_ID_attempts"]),
                "Missing sealed original FCP5 photo outcomes")
        frame.append(part)
        proof.append({**audit,"receipt_SHA256":sha(receipt_path)})
    all_rows=pd.concat(frame,ignore_index=True)
    require(len(all_rows)==SOURCE_N and all_rows.photo_id.is_unique
            and int(all_rows.original_four_hue_classified.sum())==459235,
            "Continuous trait source does not include all exact measured original photos")
    output.parent.mkdir(parents=True,exist_ok=True)
    full=original[["photo_id","morph","arm","date","inat_taxon_id","species",
                    "latitude","longitude"]].merge(
        all_rows,on="photo_id",how="left",sort=False,validate="one_to_one")
    require(len(full)==SOURCE_N and full.photo_id.tolist()==original.photo_id.tolist()
            and full.original_four_hue.notna().all(),
            "Not photo-ID complete continuous trait lake")
    require(full.morph.tolist()==full.original_four_hue.tolist(),
            "Original frozen four-hue source conflict; not safe to merge")
    require(int(full.continuous9_strict.sum())<=459235,
            "Invalid strict 9-palette source classification")
    require((full.loc[full.continuous9_strict,"continuous9_white"].between(0,1)).all(),
            "Strict white component must remain a continuous proportion")
    unknown=~full.original_four_hue_classified
    require((~full.loc[unknown & ~full.continuous9_exploratory_ambiguous,
                       "continuous9_white"].notna()).all(),
            "Technical unknown was silently promoted to flower-colour data")
    parquet=output/OUTPUT
    output.mkdir()
    full.to_parquet(parquet,index=False,compression="zstd")
    by_arm={}
    for v in proof:
        by_arm[v["source_arm"]]={k:val for k,val in v.items() if k!="source_arm"}
    report={
        "schema":SCHEMA,"source_real_original_photo_attempts":SOURCE_N,
        "source_original_four_hue_classified":459235,
        "source_technical_four_hue_unknown":SOURCE_N-459235,
        "old_original_pixel_photo_cumulative_SHA256":sha(latest),
        "original_FCP5_35group_cumulative_sha256":sha(cumulative),
        "strict_9component_real_colour_vectors":int(full.continuous9_strict.sum()),
        "ambiguous_9component_sensitivity_only":int(full.continuous9_exploratory_ambiguous.sum()),
        "strict_plus_ambiguous_nine_component_vectors":int(
            (full.continuous9_strict|full.continuous9_exploratory_ambiguous).sum()),
        "all_photo_IDs_kept_even_if_no_valid_roi":True,
        "all_9composition_values_nan_when_non_evaluable":True,
        "nine_palette_is_NOT_physical_lab_or_true_petals":True,
        "new_image_downloads":0,"new_roi_inferences":0,
        "per_original_source_arm":by_arm,
        "continuous_trait_lake_parquet_sha256":sha(parquet),
    }
    (output/RECEIPT).write_text(json.dumps(
        report,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("MOIS_1043930_ORIGINAL_CONTINUOUS9_REAL_PHOTO_TRAIT_LAKE",
          json.dumps({k:v for k,v in report.items() if k!="per_original_source_arm"},
                     sort_keys=True),flush=True)
    print("MOIS_1043930_CONTINUOUS9_SOURCE_ARMS",
          json.dumps({k:{key:v[key] for key in (
            "original_photo_attempts","source_four_hue_labels",
            "all_nine_original_fraction_columns_present",
            "strict_original_4hue_ROIs_with_all_9_fractions",
            "ambiguous_explicit_palette_usable_exploratory")}
            for k,v in by_arm.items()},sort_keys=True),flush=True)
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--latest",type=Path,required=True)
    p.add_argument("--latest-receipt",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    create(a.source,a.latest,a.latest_receipt,a.output)
if __name__=="__main__":
    main()
