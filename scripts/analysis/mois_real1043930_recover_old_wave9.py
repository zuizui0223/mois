#!/usr/bin/env python3
"""Recover 10,169 genuine existing flower ROI9 image measures from pre-350399.

The original image-based wave2/cached8192/tail16384 source Parquets retained
all nine measured flower_fraction_* columns. The old 350399->historical-main
flattening dropped them, leaving 10,169 four-hue-classifiable images with
NO 9-fraction vector in the current genuine 1,043,930-photo atlas. Source
photo IDs, original 4-colour labels, original source SHA, and terminal statuses
MUST match; never regenerate nine fractions from a four-state label.

Strictly conservative append to NEW immutable continuous9 Parquet. 390,541
previous authentic vectors + 10,169 original source-recovered => 400,710
authentic strict ROI9 vectors, 36,224 ambiguous are unchanged. No original
images downloaded, no new ROI inference, no colour-petal biology asserted.
"""
from __future__ import annotations

import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd

SOURCE_SHA="b20bbefe15020c21c3b5b0f7edf3f8e08148a4ee7dda8bf772a89bf6b126f447"
OUTPUT_SCHEMA="mois_real1043930_oldwaves_10169_original_roi9_restoration_v1"
PALETTE=("white","yellow","orange","red","pink","magenta",
         "purple","blue","bronze")
FRACTIONS=tuple("flower_fraction_"+x for x in PALETTE)
HUES=("white","yellow_orange","red_pink","blue_purple")
WAVES=(
    ("wave2",4096,1436,"61dcc82a810b02c92d28f579ef393d5b4caaba858083023a46b34158cc45c6c2"),
    ("cached8192",8192,2780,"11be4908cbd92221bf5cc3df4c1fc6d9b9fa412ccd2702b90762b9b1abfc22fe"),
    ("tail16384",16384,5953,"56dca89f67b43eca35ebebee1b59bf7a305a4a07f180bd7a6fb80557b898c307"),
)
PREVIOUS_STRICT=390541
NEW_STRICT=400710
AMBIGUOUS=36224
NEW_ORIGINAL_ROI9=10169
N=1043930

def require(ok,msg):
    if not ok:
        raise ValueError("MOIS_OLD_WAVE9_RECOVERY: "+msg)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for c in iter(lambda:f.read(4<<20),b""):h.update(c)
    return h.hexdigest()

def read_source(path,receipt,arm,expected,expected_sha):
    path,receipt=Path(path),Path(receipt)
    r=json.loads(receipt.read_text())
    require(r.get("photo_data_sha256")==expected_sha==sha(path)
            and r.get("attempted")==expected,
            "Frozen original-photo worker is not signed or complete: "+arm)
    data=pd.read_parquet(path,columns=["photo_id","morph",*FRACTIONS,
             *(["measurement_status"] if "measurement_status" in
                   __import__("pyarrow.parquet",fromlist=["ParquetFile"]).
                     ParquetFile(path).schema.names else [])])
    data.photo_id=data.photo_id.astype(str)
    require(len(data)==expected and data.photo_id.is_unique
            and data.photo_id.str.fullmatch(r"\d+").all(),
            "Original wave photo identity count is altered: "+arm)
    return data

def candidate_patch(atlas,source,arm):
    """Require an exact original source photo ID and terminal colour agreement.

    Return a fully observed source patch of ORIGINAL numeric fractions, never
    infer from the categorical morphology. A source classified photo without
    its original nine fractions is explicitly NOT recoverable.
    """
    require(set(("photo_id","original_four_hue","continuous9_strict",
                 "source_arm_original")).issubset(atlas.columns),
            "Expected exact source atlas columns")
    require(set(("photo_id","morph",*FRACTIONS)).issubset(source.columns),
            "Old real image source lacks full numerical 9component colour")
    ids=source.photo_id.astype(str)
    require(ids.is_unique and atlas.index.is_unique,"Input photo IDs are repeated")
    rows=atlas.reindex(ids.tolist())
    require(rows.original_four_hue.notna().all(),
            "Old source original photo not in true 1043930-photo source")
    require(rows.source_arm_original.eq("historical_main").all(),
            "One original image appears in wrong source arm; no rewriting")
    source_hue=source.morph.fillna("").astype(str).to_numpy()
    prior_hue=rows.original_four_hue.fillna("").astype(str).to_numpy()
    known=np.isin(source_hue,HUES)
    prior_success=np.isin(prior_hue,HUES)
    valid_raw=source[list(FRACTIONS)].apply(
        pd.to_numeric,errors="coerce").to_numpy(dtype=float)
    positive=(np.isfinite(valid_raw).all(axis=1)&
              (valid_raw>=0).all(axis=1)&
              (valid_raw<=1+1e-6).all(axis=1)&
              (valid_raw.sum(axis=1)>0))
    if "measurement_status" in source:
        statuses=source.measurement_status.fillna("").astype(str).to_numpy()
        explicit=(statuses!="")
        require(np.array_equal(
            (statuses[explicit]=="classified_four_state_morph"),
            known[explicit]),
            "Original wave terminal status/hue contradiction: "+arm)
    agreed=known&prior_success&(source_hue==prior_hue)
    already=rows.continuous9_strict.to_numpy(bool)
    candidate=agreed&positive&~already
    output=source.loc[candidate,["photo_id","morph",*FRACTIONS]].copy()
    output.photo_id=output.photo_id.astype(str)
    vals=valid_raw[candidate]
    output[list(FRACTIONS)]=(vals/vals.sum(axis=1,keepdims=True))
    require(output.photo_id.is_unique,
            "A true original photo appears twice in same recovery wave")
    audit={
       "original_wave":arm,"source_photo_attempts":len(source),
       "source_hue_success":int(known.sum()),
       "source_label_agreements_with_latest":int(agreed.sum()),
       "source_positive_9fraction_vectors":int(positive.sum()),
       "already_retained_strict9":int(already.sum()),
       "strict_new_9vectors_photoID_hue_and_source_verified":len(output),
       "source_disagreement_known_known":int((known&prior_success&
                                               (source_hue!=prior_hue)).sum()),
    }
    return output,audit

def restore(trait,trait_receipt,paths,out):
    trait,trait_receipt,out=map(Path,(trait,trait_receipt,out))
    require(not out.exists(),"Immutable ROI9 recovery output already exists")
    r=json.loads(trait_receipt.read_text())
    require(r.get("schema")==
            "mois_real1043930_original9_with_FCP447500_independent_source_reuse_v1"
            and r.get("latest_real_photo_ID_attempts")==N
            and r.get("strict_original9_vectors_now")==PREVIOUS_STRICT
            and r.get("ambiguous_original9_vectors_now")==AMBIGUOUS
            and r.get("output_photo_lake_SHA256")==SOURCE_SHA==sha(trait),
            "Unrecognized previously signed true photo9 atlas")
    frame=pd.read_parquet(trait)
    frame.photo_id=frame.photo_id.astype(str)
    require(len(frame)==N and frame.photo_id.is_unique
            and int(frame.continuous9_strict.sum())==PREVIOUS_STRICT
            and int(frame.continuous9_exploratory_ambiguous.sum())==AMBIGUOUS,
            "Original full measured photo ID universe not conserved")
    baseline_hues=frame.original_four_hue.astype(str).tolist()
    baseline_classified=frame.original_four_hue.isin(HUES).sum()
    atlas=frame.set_index("photo_id",verify_integrity=True)
    proof=[];collected=[]
    for name,n,expected,expected_sha in WAVES:
        sourcefile,receiptfile=paths[name]
        original=read_source(sourcefile,receiptfile,name,n,expected_sha)
        recovered,audit=candidate_patch(atlas,original,name)
        require(len(recovered)==expected,
                "Previously audited photoID-specific 9vectors disappeared: "+name)
        collected.append(recovered)
        proof.append({**audit,"source_parquet_sha256":sha(sourcefile),
                      "source_receipt_sha256":sha(receiptfile)})
    all_recover=pd.concat(collected,ignore_index=True)
    require(all_recover.photo_id.is_unique
            and len(all_recover)==NEW_ORIGINAL_ROI9
            and atlas.loc[all_recover.photo_id,"continuous9_strict"].eq(False).all()
            and atlas.loc[all_recover.photo_id,
                    "continuous9_exploratory_ambiguous"].eq(False).all(),
            "The original three disjoint new wave9 data are not novel authentic vectors")
    raw=all_recover[list(FRACTIONS)].to_numpy(dtype=float)
    require(np.isfinite(raw).all() and (raw>=0).all()
            and np.allclose(raw.sum(axis=1),1,atol=1e-6),
            "Original nine-source hue fractions cannot be normalized")
    ids=all_recover.photo_id.tolist()
    for name in PALETTE:
        atlas.loc[ids,"continuous9_"+name]=all_recover[
            "flower_fraction_"+name].to_numpy(dtype="float32")
    E=-np.sum(np.where(raw>0,raw*np.log(np.maximum(raw,1e-15)),0),axis=1)/np.log(9)
    G=1-np.square(raw).sum(axis=1)
    atlas.loc[ids,"continuous9_entropy"]=E.astype("float32")
    atlas.loc[ids,"continuous9_gini_simpson"]=G.astype("float32")
    atlas.loc[ids,"continuous9_dominant_fraction"]=raw.max(axis=1).astype("float32")
    atlas.loc[ids,"continuous9_strict"]=True
    atlas.loc[ids,"continuous9_status"]="four_hue_original_ROI_with_valid_9composition"
    atlas.loc[ids,"original_nine_fraction_columns_present"]=True
    atlas.loc[ids,"nine_fraction_column_count"]=9
    atlas.loc[ids,"original_9palette_numeric_positive"]=True
    updated=atlas.reset_index()
    require(len(updated)==N and updated.photo_id.is_unique
            and updated.photo_id.tolist()==frame.photo_id.tolist()
            and updated.original_four_hue.astype(str).tolist()==baseline_hues
            and int(updated.original_four_hue.isin(HUES).sum())==baseline_classified==459235
            and int(updated.continuous9_strict.sum())==NEW_STRICT
            and int(updated.continuous9_exploratory_ambiguous.sum())==AMBIGUOUS,
            "Recovered colour trait lake changed source original IDs/labels/unknowns")
    require(updated.loc[ids,"continuous9_white"].notna().all()
             if ids and all(i in updated.index for i in ids) else
             updated.loc[updated.photo_id.isin(ids),"continuous9_white"].notna().all(),
             "Recovered original nine-colour white fraction missing")
    out.mkdir(parents=True)
    target=out/"mois_real1043930_continuous9_original_source_recovered_400710_strict.parquet"
    updated.to_parquet(target,index=False,compression="zstd")
    result={
        "schema":OUTPUT_SCHEMA,
        "full_original_photo_attempts":N,
        "old_9component_strict_photo_count":PREVIOUS_STRICT,
        "source_photo9_restored_original_wave2_cached8k_tail16k":NEW_ORIGINAL_ROI9,
        "new_original_strict9_photo_count":NEW_STRICT,
        "ambiguous9_photo_count_unchanged":AMBIGUOUS,
        "four_hue_classified_unchanged":459235,
        "four_hue_known_but_still_no_9fraction":459235-NEW_STRICT,
        "remaining_technical_no9_source_photos":N-NEW_STRICT-AMBIGUOUS,
        "original_oldwave_genuine_source_checks":proof,
        "input_1043930_frozen_photo9_SHA256":SOURCE_SHA,
        "new_source_photo9_lake_SHA256":sha(target),
        "all_exact_original_photo_IDs_preserved":True,
        "original_hue_labels_preserved":True,
        "older_390541_photo9_records_not_reinferred_or_modified":True,
        "new_original_image_downloads":0,
        "new_neural_ROI_inferences":0,
        "true_anatomical_flower_petals_and_CIELAB_not_claimed":True,
    }
    (out/"mois_real1043930_recovered_400710_original9_source_receipt.json").write_text(
        json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("MOIS_REAL1043930_OLDWAVE_EXISTING_ROI9_RESTORED",json.dumps(
        result,sort_keys=True),flush=True)
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ("trait","trait-receipt","wave2","wave2-receipt",
                "cached8192","cached8192-receipt",
                "tail16384","tail16384-receipt","output"):
        p.add_argument("--"+arg,required=True,type=Path)
    z=p.parse_args()
    inputs={
        "wave2":(z.wave2,z.wave2_receipt),
        "cached8192":(z.cached8192,z.cached8192_receipt),
        "tail16384":(z.tail16384,z.tail16384_receipt)}
    restore(z.trait,z.trait_receipt,inputs,z.output)

if __name__=="__main__":main()
