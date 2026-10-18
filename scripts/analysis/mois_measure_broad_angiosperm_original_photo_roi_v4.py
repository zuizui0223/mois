#!/usr/bin/env python3
"""Actually measure a 256-photo taxonomic-breadth MOIS angiosperm cohort.

Source is worldwide distinct species/genera from open iNaturalist metadata,
not 1000 high-depth FCP species. Photo-id source has no new hue until genuine
RGB YOLO11/EfficientSAM ROI-v4 model terminals are completed and SHA checked.
The frozen original model and detector hashes remain identical to core MOIS.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

HUES=("white","yellow_orange","red_pink","blue_purple")
BIO=("white","yellow","orange","red","pink","magenta","purple","blue","bronze")
EXTRA=("measurement_id","photo_id","inat_taxon_id","species","observation_id",
       "observer_id","latitude","longitude","photo_license","photo_url_large",
       "source_snapshot_date","original_photo_source","morph","measurement_status",
       "classified_four_hue","method_source_commit",
       *(f"flower_fraction_{h}" for h in BIO))


def sha(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for v in iter(lambda:f.read(1024*1024),b""):
            h.update(v)
    return h.hexdigest()


def validate(shard:Path,source_receipt:Path,index:int):
    receipt=json.loads(source_receipt.read_text())
    if (receipt.get("schema")!="mois_open_global_angiosperm_breadth_first4096_v1"
        or receipt["status"]!="SOURCE_ONLY_4096_RESERVED_NO_NEW_PIXEL_HUES"
        or int(receipt["first_wave_original_source_photo_ids"])!=4096
        or receipt["per_shard_sha256"][str(index)]!=sha(shard)
        or receipt["no_photo_colour_or_fitness_or_bloom_label_used"] is not True):
        raise ValueError("MOIS frozen source photo IDs or eligibility changed")
    with shard.open(newline="",encoding="utf-8") as f:
        rows=list(csv.DictReader(f))
    if len(rows)!=256 or len({r["measurement_id"] for r in rows})!=256 or len(
            {r["photo_id"] for r in rows})!=256:
        raise ValueError("Source 256 photo shard incomplete")
    if any(int(r["shard"])!=index or r["source_snapshot_date"]!=receipt["source_snapshot_date"]
           or not r["photo_url_large"].startswith("https://inaturalist-open-data.s3.amazonaws.com/photos/")
           or "morph" in r or "hue" in r for r in rows):
        raise ValueError("Original scientific source photo blinded fields drift")
    return rows,receipt


class PhotoSource:
    def __init__(self,rows:list[dict],shard:int):
        self.rows=rows
        self.shard=shard

    def get_batch(self,part:int,batch:int,size:int):
        from pandas import DataFrame
        if part!=self.shard or not 0<=batch<4 or size!=64:
            raise ValueError("Original 64-photo fixed compute batch index invalid")
        segment=self.rows[batch*64:(batch+1)*64]
        if len(segment)!=64:raise ValueError("Absent original photo IDs")
        worker=DataFrame([{
            "measurement_id":r["measurement_id"],
            "image_filename":r["image_filename"],
            "photo_license":r["photo_license"]} for r in segment])
        acquisition=DataFrame([{
            "measurement_id":r["measurement_id"],
            "image_filename":r["image_filename"],
            "photo_license":r["photo_license"],
            "photo_url_large":r["photo_url_large"]} for r in segment])
        return worker,acquisition,4


def run(shard:Path,source_receipt:Path,shard_id:int,fcp_root:Path,
        sam:Path,output:Path):
    rows,freeze=validate(shard,source_receipt,shard_id)
    from scripts.acquisition.run_fcp_fresh_image_partition import run_batch
    original=fcp_root/"scripts/analysis/measure_random_photo_first_blind_partition.py"
    detector=fcp_root/"data/atlas/qualification/roi_v4_training/jrc_yolo11n_last_v4.pt"
    index=PhotoSource(rows,shard_id)
    batches=output/"sealed_original_roi_v4_batches"
    allmeasured={}
    results=[]
    fingerprint=None
    for batch in range(4):
        result=run_batch(Path("."),shard_id,batch,64,batches,
                         source=original,detector=detector,sam=sam,
                         indexed_source=index)
        if (result["status"]!="COMPLETE_WORK_UNIT" or
            result["n_source_photo_ids"]!=64 or result["n_terminal_photo_ids"]!=64 or
            result["reused_historical_measurement_csv"] is not False or
            result["image_bytes_persisted"] is not False):
            raise ValueError("Original ROI-v4 worker returned false source completion")
        if fingerprint is None:
            fingerprint=result["method"]
        elif fingerprint!=result["method"]:
            raise ValueError("Original frozen ROI model changed across four 64-photo batches")
        table=batches/f"terminal_p{shard_id:03}_b{batch:07}.csv"
        if sha(table)!=result["terminal_csv_sha256"]:
            raise ValueError("Original ROI-v4 resulting source CSV changed")
        with table.open(newline="",encoding="utf-8") as f:
            for x in csv.DictReader(f):
                if x["measurement_id"] in allmeasured:
                    raise ValueError("ROI duplicated original image measurement")
                allmeasured[x["measurement_id"]]=x
        print("MOIS_256_ORIGINAL_PIXEL_ROI_BATCH",json.dumps({
            "shard":shard_id,"batch":batch,
            "classified":result["n_classifiable"],
            "acquisition_failed":result["n_acquisition_failed"]}),flush=True)
    if len(allmeasured)!=256 or set(allmeasured)!={r["measurement_id"] for r in rows}:
        raise ValueError("Some original source photos lost without terminal colour-missingness")
    target=output/f"mois_original_256_fresh_ROI_shard_{shard_id:02d}.csv"
    target.parent.mkdir(parents=True,exist_ok=True)
    statuses=Counter()
    with target.open("x",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=EXTRA,lineterminator="\n")
        writer.writeheader()
        for r in rows:
            measured=allmeasured[r["measurement_id"]]
            status=measured["measurement_status"]
            hue=measured["morph"]
            classified=status=="classified_four_state_morph"
            if classified!=(hue in HUES):
                raise ValueError("Original four-hue model status not consistent with hue")
            o={key:r.get(key,"") for key in EXTRA}
            o.update({key:measured.get(key,"") for key in ("morph","measurement_status",
                        *(f"flower_fraction_{h}" for h in BIO))})
            o["original_photo_source"]="MOIS_worldwide_angiosperm_species_first_open_S3"
            o["classified_four_hue"]=str(classified).lower()
            o["method_source_commit"]=fingerprint["source_commit"]
            writer.writerow(o)
            statuses["four_hue_classified" if classified else "unclassified"]+=1
            statuses[status]+=1
    rec={
        "schema":"mois_256_source_original_roi_v4_four_hue_terminal_v1",
        "source_sha256":sha(shard),"source_freeze_sha256":sha(source_receipt),
        "source_frozen_4096_selection_sha256":freeze["first_wave_original_4096_csv_sha256"],
        "shard":shard_id,"original_source_ids":256,
        "distinct_original_observers":len({r["observer_id"] for r in rows}),
        "terminal_rows":sum(statuses[x] for x in ("four_hue_classified","unclassified")),
        "four_hue_classified":statuses["four_hue_classified"],
        "unclassified":statuses["unclassified"],
        "source_model_attestation":fingerprint,
        "original_new_photo_hue_terminal_SHA256":sha(target),
        "fcp_within_species_polymorphism_estimand_computed":False,
        "pixels_ephemeral_no_original_image_bytes_persisted":True,
        "original_source_photos_not_true_flowering_stage_verified":True,
    }
    (output/f"mois_original_256_ROI_shard_{shard_id:02d}_receipt.json").write_text(
        json.dumps(rec,indent=2,sort_keys=True)+"\n")
    print("MOIS_256_REAL_SOURCE_IMAGE_ROI_SHARD_PASS",
          json.dumps({"shard":shard_id,"classified":rec["four_hue_classified"],
                      "unclassified":rec["unclassified"]}),flush=True)
    return rec


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--photo-shard",type=Path,required=True)
    p.add_argument("--precolour-receipt",type=Path,required=True)
    p.add_argument("--shard-id",type=int,required=True)
    p.add_argument("--fcp-runtime-root",type=Path,required=True)
    p.add_argument("--original-sam",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    a=p.parse_args()
    run(a.photo_shard,a.precolour_receipt,a.shard_id,
        a.fcp_runtime_root,a.original_sam,a.outdir)
