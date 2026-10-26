#!/usr/bin/env python3
"""Assemble the ACTUAL latest 1,043,930 distinct MOIS four-hue photo attempts.

Photo-level evidence only: four previously signed historical pixel sources,
fully processed K1, FCP U100 depth1 new-only, and ALL 35 signed FCP5 original
image groups. Each original image ID is used once, including unknown statuses.
Does not assign colours to metadata-only candidate photos, and does not require
the 93m cheap-stage YOLO to be independently botanically correct.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import pandas as pd
import numpy as np
import pyarrow.parquet as pq

SOURCE_N=1043930
BASE_N=1009032
DEPTH5_N=34898
BASE_HUES=445805
HUES=("white","yellow_orange","red_pink","blue_purple")
ARMS={
 "historical_main":("historical/main/combined_actual_photo_checkpoint.parquet",
                   "historical/main/checkpoint_receipt.json",911767,410323),
 "historical_group00":("historical/group00/measured_photos.parquet",
                       "historical/group00/receipt.json",4138,1573),
 "historical_R1":("historical/r1/measured_photos.parquet",
                  "historical/r1/receipt.json",4468,1760),
 "historical_R2":("historical/r2/measured_photos.parquet",
                  "historical/r2/receipt.json",8765,3259),
 "K1":("k1/measured_photos.parquet","k1/receipt.json",76753,27689),
 "FCP_U100_depth1":("fcp/new_FCP_depth1_original_photo_ROIs.parquet",
                   "fcp/receipt.json",3141,1201),
}
GROUP_SCHEMA="mois_FCP_U100_nested5_new34898_genuine_original_photo_completed_group_v1"
CUM_SCHEMA="mois_FCP_nested5_genuine_partial_original_photo_roi_global_cumulative_v1"

def require(condition,why):
    if not condition:raise ValueError("MOIS_LATEST_REAL_PHOTOS: "+why)
def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(4<<20),b""):h.update(chunk)
    return h.hexdigest()

def canonical(frame,arm):
    """One uniform photo observation grain, no hue inference from technical unknown."""
    d=frame.copy()
    for col in ("photo_id","inat_taxon_id","species","observation_id",
                "morph","measurement_status","latitude","longitude","observed_on",
                "qualified_date"):
        if col not in d:d[col]=""
    for col in ("photo_id","inat_taxon_id","species","observation_id",
                "morph","measurement_status"):
        d[col]=d[col].fillna("").astype(str).str.strip()
        d.loc[d[col].str.lower().isin(("nan","none","null")) ,col]=""
    d["photo_id"]=d.photo_id.astype(str)
    require(d.photo_id.str.fullmatch(r"\d+").all()
            and d.photo_id.is_unique,arm+" photo IDs are not unique numeric originals")
    d["observed_on"]=d["qualified_date"].fillna("").astype(str).replace("nan","").where(
        d["qualified_date"].notna() & d["qualified_date"].astype(str).ne(""),d["observed_on"].fillna("").astype(str))
    d["date"]=pd.to_datetime(d.observed_on,errors="coerce")
    for col in ("latitude","longitude"):
        d[col]=pd.to_numeric(d[col],errors="coerce")
    d["morph"]=d.morph.where(d.morph.isin(HUES),"mixed_uncertain")
    d["classified"]=d.morph.isin(HUES)
    status=d.measurement_status
    strict=status.ne("")
    require((status[strict].eq("classified_four_state_morph").to_numpy()
             ==d.loc[strict,"classified"].to_numpy()).all(),
            arm+" terminal classifiable status conflicts with original hue")
    d["arm"]=arm
    d["white"]=d.morph.eq("white").astype("int8")
    d["known_image_roi"]=d.classified.astype("int8")
    d["date"]=d.date.dt.strftime("%Y-%m-%d").fillna("")
    return d[["photo_id","inat_taxon_id","species","observation_id","morph",
              "measurement_status","latitude","longitude","date","arm","white","known_image_roi"]]

def load_one(root,arm):
    path,receipt,n,n_hues=ARMS[arm]
    photo=root/path
    rec=json.loads((root/receipt).read_text())
    if arm=="historical_main":
        require(rec["schema"]=="mois_sealed_actual_pixel_checkpoint_v1"
                and rec["combined_unique_photo_attempts"]==n
                and rec["combined_four_hue_classified"]==n_hues,
                arm+" signed baseline count drift")
        pin=rec["output_sha256"][photo.name]
    elif arm.startswith("historical") or arm=="K1":
        require(rec["schema"]=="mois_global93m_observation_baseline_real_RGB_photo_results_v1"
                and rec["attempted"]==n and rec["classifiable"]==n_hues,
                arm+" signed genuine photo count drift")
        pin=rec["measured_photo_parquet_sha256"]
    else:
        require(rec["newly_true_image_measured_original_photo_attempts"]==n
                and rec["newly_true_image_FCP_four_hue_classified"]==n_hues,
                arm+" genuine FCP U100 photo count drift")
        pin=rec["actual_new_image_parquet_SHA256"]
    require(sha(photo)==pin,arm+" original photo parquet source SHA changed")
    d=canonical(pd.read_parquet(photo),arm)
    require(len(d)==n and int(d.known_image_roi.sum())==n_hues,
            arm+" source photo ROI counts do not match signed receipt")
    return d,{"arm":arm,"original_photo_ids":n,"classified":n_hues,
              "parquet_SHA256":pin,"receipt_SHA256":sha(root/receipt)}

def assemble(source,output):
    source,output=Path(source),Path(output)
    require(not output.exists(),"Refuse to rewrite previous complete photo pixel cohort")
    cumulative=source/"fcp5"/"cumulative"/"mois_FCP5_partial_genuine_original_pixel_cumulative.json"
    cur=json.loads(cumulative.read_text(encoding="utf-8"))
    require(cur.get("schema")==CUM_SCHEMA
            and cur.get("new_actual_finished_original_photo_groups")==35
            and cur.get("new_finished_original_photo_group_IDs")==list(range(35))
            and cur.get("new_genuine_original_photo_pixel_attempts")==DEPTH5_N
            and cur.get("previous_actual_independently_verified_original_photo_attempts")==BASE_N
            and cur.get("previous_actual_independently_verified_four_hue_ROIs")==BASE_HUES
            and cur.get("exact_precolour_disjoint_source_global_original_photo_attempts")==SOURCE_N,
            "Not full 35/35 signed FCP5 genuine original photo group census")
    groups={r["group"]:r for r in cur["group_receipts"]}
    require(set(groups)==set(range(35)),"35 all signed original FCP groups not provided")
    frames=[]
    audits=[]
    for arm in ARMS:
        data,proof=load_one(source,arm)
        frames.append(data)
        audits.append(proof)
    count=0
    hues=0
    for g in range(35):
        files=source/"fcp5"/("group-%02d"%g)
        path=files/"original_photo_ROIs.parquet"
        receipt_path=files/"original_group_receipt.json"
        rec=json.loads(receipt_path.read_text(encoding="utf-8"))
        require(rec.get("schema")==GROUP_SCHEMA
                and rec.get("original_worker_group")==g
                and sha(receipt_path)==groups[g]["source_group_receipt_SHA256"],
                "FCP5 cumulative original receipt SHA or group changed")
        digest=sha(path)
        require(digest==rec["saved_file_SHA256"][path.name],
                "Published FCP5 original photo pixels source SHA changed")
        n=int(rec["new_actual_original_photo_ID_attempts"])
        n_hues=int(rec["new_actual_four_hue_classified"])
        require(n>0 and n<=1024 and
                sum(rec["new_original_hue_counts"][h] for h in HUES)==n_hues,
                "FCP5 group original source hue counts malformed")
        df=canonical(pd.read_parquet(path), "FCP5_group_"+str(g).zfill(2))
        require(len(df)==n and int(df.known_image_roi.sum())==n_hues,
                "FCP5 original photo group not actually measured")
        count+=n;hues+=n_hues;frames.append(df)
        audits.append({"arm":str(g),"original_photo_ids":n,
                       "classified":n_hues,"parquet_SHA256":digest,
                       "receipt_SHA256":sha(receipt_path)})
    require(count==DEPTH5_N
            and hues==cur["new_genuine_automatic_four_hue_ROIs"]
            and sum(r["original_photo_ids"] for r in audits)==SOURCE_N
            and sum(r["classified"] for r in audits)==BASE_HUES+hues
            and BASE_HUES+hues==cur["exact_precolour_disjoint_source_global_four_hue_ROIs"],
            "Latest 35 group attempts and four-hue source census differ")
    d=pd.concat(frames,ignore_index=True)
    require(len(d)==SOURCE_N and d.photo_id.is_unique
            and int(d.known_image_roi.sum())==BASE_HUES+hues,
            "Global latest actual photo-ID union is not exactly disjoint")
    output.mkdir(parents=True)
    parquet=output/"mois_latest_1043930_true_original_photo_hue_observation.parquet"
    d.to_parquet(parquet,compression="zstd",index=False)
    dates=int((d.date!="").sum())
    recognized=int(((d.latitude.between(-90,90)) &
                   (d.longitude.between(-180,180))).sum())
    result={
        "schema":"mois_fully_signed_latest_original_photo_colour_1043930_v1",
        "source_pixel_attempts":len(d),
        "source_original_unique_photo_ids":int(d.photo_id.nunique()),
        "source_four_hue_classified":int(d.known_image_roi.sum()),
        "original_technical_colour_unknowns":int((1-d.known_image_roi).sum()),
        "four_hue_counts":{h:int(d.morph.eq(h).sum()) for h in HUES},
        "source_observation_day_present":dates,
        "source_coordinates_valid":recognized,
        "source_original_photo_arm_counts":audits,
        "fcp5_global_cumulative_receipt_sha256":sha(cumulative),
        "latest_fcp5_genuine_complete_groups":35,
        "new_photo_image_downloads":0,
        "new_ROI_photo_inferences":0,
        "photo_parquet_sha256":sha(parquet),
        "not_93m_already_measured":True,
    }
    (output/"mois_latest_1043930_photo_union_receipt.json").write_text(
        json.dumps(result,sort_keys=True,indent=2,ensure_ascii=False)+"\n",
        encoding="utf-8")
    print("MOIS_LATEST_1043930_REAL_SOURCE_HUE_CENSUS",
          json.dumps({k:v for k,v in result.items() if k!="source_original_photo_arm_counts"},
                     sort_keys=True),flush=True)
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    assemble(a.source,a.output)
if __name__=="__main__":
    main()
