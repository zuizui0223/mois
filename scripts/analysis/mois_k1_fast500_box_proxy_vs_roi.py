#!/usr/bin/env python3
"""Assess cheap 500px YOLO box-HSV fractions against frozen K1 ROI-v4 colours.

This is a DATA-ONLY calibration on genuine source photo IDs; no new image
downloads, no flower stage certification and NO recolouring old unknowns.
Use the full 76,753 photo K1 original source and the previously completed
256-medium-image test packet (128 ROI-labelled / 128 ROI-gate unknown).
Only the 128 original four-class labels can serve as reference, and even
those are machine labels rather than human-validated flower petals.
"""
from __future__ import annotations
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

HUES=("white","yellow_orange","red_pink","blue_purple")
SCHEMA="mois_K1_256_500px_original_yolo_box_proxy_vs_frozen_4hue_calibration_v1"
K1_SHA="0d3de3374a06b1e7b9ad0281e2ab7844632ed139db871488524061191c8f34c1"
K1_STATUS_COLOUR="classified_four_state_morph"
K1_STATUS_UNKNOWN="not_evaluable_roi_or_flip_gate"
N=256
NUM_EACH=128
MEASUREMENT_STATUSES=(K1_STATUS_COLOUR,K1_STATUS_UNKNOWN)

def require(cond,reason):
    if not cond:
        raise ValueError("MOIS_K1_FAST_PIXEL_PROXY: "+reason)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(4<<20),b""):
            h.update(b)
    return h.hexdigest()

def analyze(k1_parquet, stage1_packet, output, *,
            expected_k1_sha=K1_SHA, expected_total=N,
            expected_group_size=NUM_EACH, min_margin=.10):
    import pandas as pd
    from scripts.acquisition.mois_full93m_medium_pixel_scan import COLUMNS_OUT
    k1_parquet,stage1_packet,output=map(Path,(k1_parquet,stage1_packet,output))
    require(not output.exists(),"Do not overwrite published photographic calibration results")
    require(sha(k1_parquet)==expected_k1_sha,"Original K1 photo ROI source bytes changed")
    require(math.isfinite(min_margin) and 0<=min_margin<=1,"Invalid prespecified proxy margin")
    frame=pd.read_parquet(k1_parquet,columns=["photo_id","measurement_status","morph"])
    require(frame.photo_id.astype(str).is_unique,"Original K1 photo IDs repeated")
    labels={str(pid):(str(status),str(hue))
            for pid,status,hue in
            zip(frame.photo_id,frame.measurement_status,frame.morph)}
    read=[]
    with gzip.open(stage1_packet,"rt",encoding="utf-8",newline="") as file:
        rd=csv.DictReader(file)
        require(tuple(rd.fieldnames or [])==tuple(COLUMNS_OUT),
                "Medium 500px actual-photo packet fields changed")
        read=list(rd)
    photo_ids=[r["photo_id"] for r in read]
    require(len(read)==expected_total and len(set(photo_ids))==expected_total,
            "Every original photo ID must appear exactly once in 256-image pilot")
    require(set(photo_ids)<=set(labels),"Medium image photo IDs not from frozen K1 source")
    status_counts=Counter(labels[pid][0] for pid in photo_ids)
    require(all(status_counts[s]==expected_group_size for s in MEASUREMENT_STATUSES)
            and sum(status_counts.values())==expected_total,
            "K1 ROI-success/ROI-gate calibration strata are not unchanged")
    confusion={k:{v:0 for v in HUES} for k in HUES}
    source_class_counts=Counter()
    total_status=Counter()
    coverage=Counter()
    exact=Counter()
    considered=Counter()
    safe=Counter()
    unknown_box=Counter()
    invalid=0
    for row in read:
        pid=row["photo_id"]
        status,hue=labels[pid]
        detected=row["screen_status"]=="candidate"
        total_status[(status,row["screen_status"])]+=1
        if status==K1_STATUS_COLOUR:
            require(hue in HUES,"Frozen four-hue labels changed")
            source_class_counts[hue]+=1
        if not detected:
            continue
        try:
            p=[float(row["box_proxy_"+key]) for key in HUES]
            require(all(math.isfinite(x) and 0<=x<=1 for x in p),
                    "Invalid simple box colour fraction")
            sorted_four=sorted(enumerate(p),key=lambda item:(-item[1],item[0]))
            best=sorted_four[0][0]
            strongest,second=p[best],sorted_four[1][1]
            predicted=HUES[best]
            margin=strongest-second
        except (ValueError,TypeError,KeyError):
            invalid+=1
            continue
        coverage[status]+=1
        if status==K1_STATUS_COLOUR:
            confusion[hue][predicted]+=1
            exact["all"]+=int(predicted==hue)
            considered["all"]+=1
            if margin>=min_margin:
                exact["confident"]+=int(predicted==hue)
                considered["confident"]+=1
        else:
            # This is an image-box proxy present, NOT a rescued verified petal.
            unknown_box["candidate_proxy"]+=1
            unknown_box["margin_ge_prespecified"]+=int(margin>=min_margin)
        if margin>=min_margin:
            safe[status]+=1

    per_class={
        hue:{
            "machine_labelled_K1_reference_photos":source_class_counts[hue],
            "box_reference_agreement_when_candidate":
                (confusion[hue][hue]/sum(confusion[hue].values())
                 if sum(confusion[hue].values()) else None),
            "predicted_box_proxy_counts":confusion[hue],
        } for hue in HUES
    }
    r={
        "schema":SCHEMA,
        "K1_real_original_source_parquet_sha256":sha(k1_parquet),
        "medium_genuine_pixel_packet_sha256":sha(stage1_packet),
        "original_photos_evaluated":expected_total,
        "original_ROI_success_reference":status_counts[K1_STATUS_COLOUR],
        "original_ROI_gate_unknown_without_reference":
            status_counts[K1_STATUS_UNKNOWN],
        "original_ROI_status_by_medium_screen":{
            str(status):{state:n for (s,state),n in total_status.items()
                         if s==status}
            for status in MEASUREMENT_STATUSES},
        "prespecified_four_proxy_margin":min_margin,
        "roi_success_four_hue_comparable_photo_count":considered["all"],
        "original_roi_success_box_proxy_argmax_exact_agreement":
            (exact["all"]/considered["all"] if considered["all"] else None),
        "confident_proxy_coverage_in_roi_success":
            considered["confident"]/status_counts[K1_STATUS_COLOUR],
        "confident_proxy_exact_agreement":
            (exact["confident"]/considered["confident"]
             if considered["confident"] else None),
        "reference_by_original_hue":per_class,
        "roi_gate_unknown_box_proxy_present_NOT_rescued":unknown_box["candidate_proxy"],
        "roi_gate_unknown_margin_ge_threshold_NOT_rescued":
            unknown_box["margin_ge_prespecified"],
        "candidate_boxes_with_unusable_colour_proxy":invalid,
        "automatically_human_certified_open_flower_petals":0,
        "unknown_original_roi_labels_promoted_to_colours":0,
        "new_image_bytes_fetched":0,
        "stage1_box_proxy_is_not_a_verified_petal_mask":True,
        "confidence_threshold_or_model_tuned_on_this_same_sample":False,
    }
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(r,ensure_ascii=False,sort_keys=True,indent=2)+"\n",
                      encoding="utf-8")
    return r

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--k1-parquet",type=Path,required=True)
    p.add_argument("--stage1-packet",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    print("MOIS_K1_REAL500_PIXEL_FAST_PROXY_FOUR_HUE_AGREEMENT",
          json.dumps(analyze(a.k1_parquet,a.stage1_packet,a.output),sort_keys=True))

if __name__=="__main__":
    main()
