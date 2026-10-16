#!/usr/bin/env python3
"""Phenology-FIRST execution order for the SAME Mois photo sample.

This never selects, drops or changes any source observation/photo ID. It only
reorders the already SHA-sealed, 128-partition Mois sample ROI packets so
original observations annotated with iNaturalist Flowers (term 12/value 13)
are evaluated before unannotated, unverified or explicit no-flower evidence.

The full selected sample must STILL reach terminal completeness. A partial
Flowers-priority processing wave is not a global biological result, and no
filter is allowed to change sampling probabilities or species/genera weights.
The four-table original FCP 2026-09-27 metadata did NOT contain annotations.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import shutil
import sys
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))
from scripts.acquisition.mois_stratified_photo_selection_20261016 import sha256_file
from scripts.acquisition.mois_build_sample_roi_handoff_20261016 import (
    WORKER, ACQUISITION, N_PARTITIONS,
)
from scripts.analysis.mois_phenology_observer_bias_audit_20261016 import load_overlay

PRI={
    "FLOWERS":0,
    "NO_ANNOTATION":1,
    "NOT_VERIFIED":2,
    "BUDS_ONLY":3,
    "FRUITS_ONLY":4,
    "NO_FLOWERS_FRUITS":5,
}


def reorder(source:Path, overlay:Path, overlay_sha256:str, destination:Path)->dict:
    stage=destination.with_name(destination.name+".INCOMPLETE")
    if destination.exists() or stage.exists():
        raise ValueError("Refuse already-existing or interrupted flowering-priority packet handoff")
    root_receipt=source/"mois_roi_packet_receipt.json"
    evidence=json.loads(root_receipt.read_text(encoding="utf8"))
    if (evidence.get("schema")!="mois_sample_scoped_original_roi_v4_blind_packets_v1"
        or evidence.get("status")!="SAMPLE_ONLY_BLIND_METADATA_SEALED_PIXELS_NOT_MEASURED"
        or evidence.get("only_mois_sample_not_all_fcp_source_photos") is not True
        or evidence.get("original_pixels_measured") is not False):
        raise ValueError("No admitted complete original Mois-only photo selection")
    annotation,origin=load_overlay(overlay,overlay_sha256)
    # Every selected observation UUID must be matched exactly once.
    priority={}
    metadata=source/"private_metadata_join.csv"
    if sha256_file(metadata)!=evidence["file_sha256"].get("private_metadata_join.csv"):
        raise ValueError("Original observation/photo join key changed after selection")
    stage.mkdir(parents=True)
    with metadata.open(newline="",encoding="utf8") as f:
        reader=csv.DictReader(f)
        if not {"measurement_id","observation_id","photo_id","photo_design_weight"}.issubset(reader.fieldnames or []):
            raise ValueError("Original protected metadata and inclusion probabilities missing")
        for row in reader:
            ident,uuid=row["measurement_id"],row["observation_id"]
            if not ident or ident in priority or uuid not in annotation:
                raise ValueError("Missing/duplicated selected ROI ID or unverified original observation key")
            priority[ident]=PRI[annotation[uuid]]
    if len(priority)!=evidence["selected_photo_ids"]:
        raise ValueError("Annotation priority mapping lost or gained selected original photographs")
    # One sample photo per original observation is guaranteed by the sealed
    # source selector; exact overlay coverage prevents erasure of unannotated.
    if set(annotation)!={r["observation_id"] for r in csv.DictReader(metadata.open(newline="",encoding="utf8"))}:
        raise ValueError("Annotation overlay must cover the selected original observation universe exactly")
    counts=Counter()
    copy_hashes={}
    batch_order=[]
    for p in range(N_PARTITIONS):
        name_w=f"worker_partition_{p:03}.csv"
        name_a=f"acquisition_partition_{p:03}.csv"
        with (source/name_w).open(newline="",encoding="utf8") as wf, \
             (source/name_a).open(newline="",encoding="utf8") as af:
            w=csv.DictReader(wf);a=csv.DictReader(af)
            if w.fieldnames!=WORKER or a.fieldnames!=ACQUISITION:
                raise ValueError("Invalid protected FCP-blind worker packet interface")
            workers=list(w);urls=list(a)
        if (sha256_file(source/name_w)!=evidence["file_sha256"][name_w]
            or sha256_file(source/name_a)!=evidence["file_sha256"][name_a]
            or len(workers)!=len(urls)
            or len(workers)!=evidence["partition_photo_counts"][p]):
            raise ValueError("Changed photo source partition or missing original photo IDs")
        if len({row["measurement_id"] for row in workers})!=len(workers):
            raise ValueError("Duplicate photo source ID in ROI worker packet")
        paired=[]
        for ordinal,(w,a) in enumerate(zip(workers,urls)):
            ident=w["measurement_id"]
            if ident not in priority or any(w[k]!=a[k] for k in WORKER):
                raise ValueError("Original worker/acquisition identity mismatch")
            paired.append((priority[ident],ordinal,w,a))
            counts[priority[ident]]+=1
        paired.sort(key=lambda item:(item[0],item[1]))
        for b in range((len(paired)+63)//64):
            chunk=paired[b*64:(b+1)*64]
            batch_order.append({
                "partition":p,"batch":b,"batch_size":64,
                "first_priority":chunk[0][0],
                "n_flowers_annotated":sum(z[0]==0 for z in chunk),
                "n_selected":len(chunk),
            })
        for name,cols,position in ((name_w,WORKER,2),(name_a,ACQUISITION,3)):
            with (stage/name).open("w",newline="",encoding="utf8") as out:
                writer=csv.DictWriter(out,fieldnames=cols)
                writer.writeheader()
                for entry in paired:
                    writer.writerow(entry[position])
            copy_hashes[name]=sha256_file(stage/name)
    if sum(counts.values())!=evidence["selected_photo_ids"]:
        raise ValueError("Flowering-priority reordering did not preserve photo sample")
    with (stage/"private_phenology_roi_work_order.csv").open("w",newline="",encoding="utf8") as out:
        writer=csv.DictWriter(out,fieldnames=[
            "partition","batch","batch_size","first_priority",
            "n_flowers_annotated","n_selected"
        ])
        writer.writeheader()
        writer.writerows(sorted(batch_order,key=lambda z:(
            z["first_priority"],-z["n_flowers_annotated"],z["partition"],z["batch"]
        )))
    copy_hashes["private_phenology_roi_work_order.csv"]=sha256_file(
        stage/"private_phenology_roi_work_order.csv"
    )
    shutil.copyfile(metadata,stage/"private_metadata_join.csv")
    copy_hashes["private_metadata_join.csv"]=sha256_file(stage/"private_metadata_join.csv")
    result={**evidence,
        "status":"SAMPLE_ONLY_BLIND_METADATA_SEALED_PIXELS_NOT_MEASURED",
        "file_sha256":copy_hashes,
        "phenology_prioritized_processing_order":True,
        "phenology_priority":list(PRI),
        "phenology_annotation_overlay_sha256":overlay_sha256,
        "phenology_source_date":origin["annotation_snapshot_date"],
        "phenology_priority_counts":{label:counts[value] for label,value in PRI.items()},
        "phenology_first_batch_schedule_preserves_all_sample_units":True,
        "parent_unmodified_sample_packet_receipt_sha256":sha256_file(root_receipt),
        "photographs_dropped_or_added":0,
        "sampling_inclusion_probabilities_changed":False,
        "flowers_only_early_results_are_not_worldwide_inference":True,
        "original_pixels_measured":False,
    }
    (stage/"mois_roi_packet_receipt.json").write_text(json.dumps(result,sort_keys=True,indent=2)+"\n",encoding="utf8")
    stage.replace(destination)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original-sample-roi-packets",required=True,type=Path)
    p.add_argument("--private-observation-annotation-overlay",required=True,type=Path)
    p.add_argument("--annotation-overlay-sha256",required=True)
    p.add_argument("--out-private-prioritized-packets",required=True,type=Path)
    a=p.parse_args()
    result=reorder(a.original_sample_roi_packets,a.private_observation_annotation_overlay,
                   a.annotation_overlay_sha256,a.out_private_prioritized_packets)
    print(json.dumps({k:result[k] for k in (
        "selected_photo_ids","phenology_priority_counts",
        "photographs_dropped_or_added","sampling_inclusion_probabilities_changed")},sort_keys=True))

if __name__=="__main__":
    main()
