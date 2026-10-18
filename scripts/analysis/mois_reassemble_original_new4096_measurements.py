#!/usr/bin/env python3
"""MOIS: materialise REAL newly measured original flower-hue photos.

Reconcile SHA-frozen source-ID selections with four original ROI-v4 batch
receipts for every 256-photo work shard. Missing shards are fatal, never
implicitly measured. No FCP polymorphism D/W/geographic estimand is fitted.
Original photos and model byte hashes are not modified.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
from collections import Counter
from pathlib import Path

BASE_N=248312
BASE_CLASSIFIED=120807
PER_WAVE=4096
SHARDS=16
PER_SHARD=256
HUES=("white","yellow_orange","red_pink","blue_purple")
MEASURED_SCHEMA="mois_256_source_original_roi_v4_four_hue_terminal_v1"
METHOD={
 "source_commit":"9fae6ccdf684a46026f72ba12e98de2c5c54bf2a",
 "measurement_git_blob":"b1d7091db44682559bbe4b48a9f3b8db648501a6",
 "detector_sha256":"f1aaeec4664fe2c178e5cf2bc1f508977bef3e4aa7b40613026cb8ae3de789d5",
 "efficient_sam_encoder_sha256":"84ed466ffcc5c1f8d08409bc34a23bb364ab2c15e402cb12d4335a42be0e0951",
 "efficient_sam_decoder_sha256":"a62f8fa5ea080447c0689418d69e58f1e83e0b7adf9c142e2bd9bcc8045c0b11",
 "location_blind":True,"no_old_measurement_csv_inputs":True,
 "persisted_image_pixels":False
}
BASE_COLS=("photo_id","inat_taxon_id","species","morph","source_original_measured",
 "old_mois_or_new_fcp","measurement_status","latitude","longitude",
 "photo_license","observation_id","observed_on","time_observed_at",
 "wc_elevation_m","soil_pH","soil_SOC","soil_N","soil_clay",
 "soil_available_water_proxy","wc_bio1","wc_bio5","wc_bio12","wc_bio15")


def sha(path:Path)->str:
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for part in iter(lambda:f.read(1<<20),b""): h.update(part)
    return h.hexdigest()


def require(ok:bool,reason:str):
    if not ok:raise ValueError(reason)


def load_plan(root:Path,source_name:str)->tuple[list[dict],dict]:
    receipts=list(root.rglob("mois_first4096_original_source_only_freeze.json"))
    require(len(receipts)==1,f"Missing or duplicate original {source_name} source freeze")
    freeze=json.loads(receipts[0].read_text(encoding="utf-8"))
    require(freeze["status"]=="SOURCE_ONLY_4096_RESERVED_NO_NEW_PIXEL_HUES"
            and freeze["first_wave_original_source_photo_ids"]==PER_WAVE
            and freeze["first_wave_original_shards"]==SHARDS
            and freeze["photos_per_original_shard"]==PER_SHARD
            and freeze["first_wave_distinct_source_observers"]==PER_WAVE
            and freeze["no_photo_colour_or_fitness_or_bloom_label_used"] is True,
            "Not an authenticated original outcome-blind 4096-photo sample")
    source_file={
        "wave1":"mois_source_only_4096_original_photo_selection.csv",
        "wave2":"mois_original_full_3430502_outside_first_screen_WAVE2_4096.csv",
        "wave3":"mois_original_full_3430502_outside_first_screen_WAVE3_4096.csv",
        "wave5":"mois_next4096_after313535_photo_ID_source_blind_wave5.csv",
    }
    require(source_name in source_file,"Unknown original-wave ID selection")
    allfiles=list(root.rglob(source_file[source_name]))
    require(len(allfiles)==1 and sha(allfiles[0])==freeze["first_wave_original_4096_csv_sha256"],
            "Original photo-first source ID frozen CSV changed")
    with allfiles[0].open(newline="",encoding="utf-8") as f:
        records=list(csv.DictReader(f))
    require(len(records)==PER_WAVE and
            len({r["photo_id"] for r in records})==PER_WAVE and
            len({r["measurement_id"] for r in records})==PER_WAVE and
            len({r["observer_id"] for r in records})==PER_WAVE,
            "Photo/measurement/observer source identities are not disjoint")
    for k in range(SHARDS):
        f=list(root.rglob(f"mois_source_only_shard_{k:02d}_256.csv"))
        require(len(f)==1 and sha(f[0])==freeze["per_shard_sha256"][str(k)],
                f"Original {source_name} source shard {k} SHA changed")
    return records,freeze


def load_real_wave(plan_root:Path,artifact_root:Path,
                   source_name:str)->tuple[list[dict],dict]:
    plan,freeze=load_plan(plan_root,source_name)
    expected={r["measurement_id"]:r for r in plan}
    results={}
    status=Counter()
    for k in range(SHARDS):
        receipts=list(artifact_root.rglob(
            f"mois_original_256_ROI_shard_{k:02d}_receipt.json"))
        terminals=list(artifact_root.rglob(
            f"mois_original_256_fresh_ROI_shard_{k:02d}.csv"))
        require(len(receipts)==len(terminals)==1,
                f"Original {source_name} shard {k} not complete and SHA sealed")
        r=json.loads(receipts[0].read_text(encoding="utf-8"))
        require(r["schema"]==MEASURED_SCHEMA and r["shard"]==k and
                r["terminal_rows"]==r["original_source_ids"]==PER_SHARD and
                r["source_sha256"]==freeze["per_shard_sha256"][str(k)] and
                r["source_frozen_4096_selection_sha256"]==
                freeze["first_wave_original_4096_csv_sha256"] and
                r["source_model_attestation"]==METHOD and
                r["pixels_ephemeral_no_original_image_bytes_persisted"] is True and
                r["fcp_within_species_polymorphism_estimand_computed"] is False and
                sha(terminals[0])==r["original_new_photo_hue_terminal_SHA256"],
                f"Unsealed/mutated original photo or model output {source_name}:{k}")
        records=[]
        with terminals[0].open(newline="",encoding="utf-8-sig") as f:
            for row in csv.DictReader(f):
                records.append(row)
                mid=row["measurement_id"]
                require(mid in expected and int(expected[mid]["shard"])==k,
                        f"ROI worker {source_name} returned an unselected/source-foreign photo")
                require(mid not in results and row["photo_id"]==expected[mid]["photo_id"] and
                        row["inat_taxon_id"]==expected[mid]["inat_taxon_id"] and
                        row["species"]==expected[mid]["species"] and
                        row["observation_id"]==expected[mid]["observation_id"] and
                        row["observer_id"]==expected[mid]["observer_id"],
                        f"Post-blind original source identity drift {source_name}:{k}")
                require(row["morph"] in (*HUES,"mixed_uncertain") and
                        (row["measurement_status"]=="classified_four_state_morph")==
                        (row["morph"] in HUES) and
                        row["classified_four_hue"]==
                        ("true" if row["morph"] in HUES else "false"),
                        "Original 100543 method coarse-hue status contract broken")
                results[mid]=row
                status[row["morph"]]+=1
        require(len(records)==PER_SHARD and
                sum(row["morph"] in HUES for row in records)==r["four_hue_classified"] and
                sum(row["morph"] not in HUES for row in records)==r["unclassified"],
                f"Unreported unclassifiable original flowers in {source_name}:{k}")
    require(len(results)==PER_WAVE and set(expected)==set(results),
            f"Original {source_name} 4096 source photo models not fully finished")
    output=[]
    for source in plan:
        q=results[source["measurement_id"]]
        outcome={k:"" for k in BASE_COLS}
        for k in BASE_COLS:
            if k in q and q[k] is not None:
                outcome[k]=q[k]
        outcome["source_original_measured"]="mois_"+source_name+"_new_original_ROI_v4"
        outcome["old_mois_or_new_fcp"]="mois_new_original_ROI"
        # First wave was originally deliberately date-blind; second wave
        # includes the real observed_on field from the full original snapshot.
        outcome["observed_on"]=source.get("observed_on","")
        output.append(outcome)
    return output,{
       "status":"REAL_4096_ORIGINAL_ROI_PHOTO_COLOURS_FULLY_SEALED",
       "wave":source_name,"n_real_original_photo_IDs":len(output),
       "n_four_hue_classified":sum(o["morph"] in HUES for o in output),
       "source_freeze_sha256":sha(next(iter(plan_root.rglob(
           "mois_first4096_original_source_only_freeze.json")))),
       "source_original_selection_sha256":freeze["first_wave_original_4096_csv_sha256"],
       "model_attestation":METHOD,
       "no_imported_FCP_intraspecific_D_W_or_50km_measures":True,
       "original_photo_pixels_not_persisted":True,
    }


def assemble(base:Path,base_receipt:Path,wave_specs:list[tuple[str,Path,Path]],
             output_dir:Path):
    require(not output_dir.exists() or not any(output_dir.iterdir()),
            "Cannot overwrite a previous actual enlarged photo source")
    freeze=json.loads(base_receipt.read_text(encoding="utf-8"))
    require(freeze["final_photo_ids"]==BASE_N and
            freeze["final_classified_hue_ids"]==BASE_CLASSIFIED and
            sha(base)==freeze["result_source_photo_csv_sha256"],
            "Frozen first 248312 unique source photos differ from main result")
    add=[]
    metadocs={}
    for label,pr,ar in wave_specs:
        new,receipt=load_real_wave(pr,ar,label)
        add.extend(new)
        metadocs[label]=receipt
    priorids=set()
    allclass=0
    output_dir.mkdir(parents=True,exist_ok=True)
    dest=output_dir/"mois_all_real_original_hue_measured_248312_plus_new_waves.csv.gz"
    with dest.open("wb") as outfile:
        with gzip.GzipFile(fileobj=outfile,mode="wb",mtime=0,filename="") as compressed:
            with io.TextIOWrapper(compressed,encoding="utf-8",newline="") as outf:
                writer=csv.DictWriter(outf,fieldnames=BASE_COLS,lineterminator="\n")
                writer.writeheader()
                with gzip.open(base,"rt",newline="",encoding="utf-8") as f:
                    r=csv.DictReader(f)
                    require(tuple(r.fieldnames or ())==BASE_COLS,
                            "Old photographic source columns or provenance changed")
                    for row in r:
                        pid=row["photo_id"]
                        require(pid and pid not in priorids,
                                "Original 248312 photo source is not ID-unique")
                        priorids.add(pid)
                        allclass+=row["morph"] in HUES
                        writer.writerow(row)
                require(len(priorids)==BASE_N and allclass==BASE_CLASSIFIED,
                        "Original historical+FCP photo source totals drifted")
                for row in add:
                    pid=row["photo_id"]
                    require(pid and pid not in priorids,
                            "New source original photo ID overlaps historical or first/second wave")
                    priorids.add(pid)
                    allclass+=row["morph"] in HUES
                    writer.writerow(row)
    require(len(priorids)==BASE_N+PER_WAVE*len(wave_specs),
            "Missing or duplicate real photo outcomes in appended source waves")
    receipt={
      "schema":"mois_existing_248312_plus_completed_real_4096_photo_colour_waves_v1",
      "status":"ACTUAL_COMPLETED_ORIGINAL_IMAGE_ROI_ADDITIONS_NOT_SOURCE_ONLY",
      "previous_source_original_photos":BASE_N,
      "previous_source_photo_hues_classified":BASE_CLASSIFIED,
      "new_original_source_photos":len(add),
      "new_four_hue_measured_photo_ids":sum(r["morph"] in HUES for r in add),
      "all_original_unique_photos":len(priorids),
      "all_four_hue_classified_photos":allclass,
      "source_waves":metadocs,
      "original_historical_photo_source_sha256":sha(base),
      "all_real_source_photo_hue_csv_sha256":sha(dest),
      "genetic_polymorphism_spatial_D_W_or_50km_measured":False,
      "original_model_SHA_verified_for_every_photo":True,
      "new_source_photo_IDs_never_reclassified_as_confirmed_open_flower":True,
      "photo_based_genus_moisture_six_anchor_model_refitted":False,
      "original_RGB_photo_bytes_persisted":False,
    }
    (output_dir/"mois_real_additional_photo_waves_receipt.json").write_text(
        json.dumps(receipt,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("MOIS_REAL_APPENDED_ACTUAL_PHOTO_HUE_COUNTS",
          json.dumps({k:receipt[k] for k in (
              "new_original_source_photos","new_four_hue_measured_photo_ids",
              "all_original_unique_photos","all_four_hue_classified_photos")},sort_keys=True),
          flush=True)
    return receipt


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--prior-248312",type=Path,required=True)
    p.add_argument("--prior-receipt",type=Path,required=True)
    p.add_argument("--first-plan",type=Path,required=True)
    p.add_argument("--first-real-results",type=Path,required=True)
    p.add_argument("--second-plan",type=Path)
    p.add_argument("--second-real-results",type=Path)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    require((a.second_plan is None)==(a.second_real_results is None),
            "Second photo wave requires both original source and original ROI outcomes")
    groups=[("wave1",a.first_plan,a.first_real_results)]
    if a.second_plan:
        groups.append(("wave2",a.second_plan,a.second_real_results))
    assemble(a.prior_248312,a.prior_receipt,groups,a.output)
