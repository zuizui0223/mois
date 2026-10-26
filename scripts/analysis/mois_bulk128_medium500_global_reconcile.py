#!/usr/bin/env python3
"""Reconcile 128 genuine bulk photo 500px/box-Lab source-selected shards.

Only original image-pixel packets whose signed photo/observation identities
match the frozen 128-shard source are aggregated. No new images are fetched,
no EfficientSAM is run, and no photo is promoted from an uncertain detection
to botanical white/pigmented or a strict original nine-palette outcome.
"""
from __future__ import annotations

import argparse,collections,csv,gzip,hashlib,json
from pathlib import Path

from scripts.acquisition.mois_bulk128_medium500_photo_worker import (
    SCHEMA as PART_SCHEMA,WEIGHT_SHA,SHARD_MANIFEST_SHA
)
from scripts.acquisition.mois_bulk128_tier10_source import (
    SCHEMA as SOURCE_SCHEMA
)

SCHEMA="mois_bulk128_true_original_medium500_picture_pixel_complete_audit_v1"
COLUMNS=(
    "photo_id","inat_taxon_id","observation_id","design_cohort","original_priority",
    "latitude_region","local_season","screen_status","medium_decoded",
    "floral_box_main_L_star","floral_box_main_C_star",
    "floral_box_main_rgb_hex","unverified_floral_box_not_petal_ROI",
)
SOURCE_SUMMARY_SHA="efe19bb147c439d14c9342b5a70a4ada6d4cb00da0494a8882653ec30520425d"
SOURCE_MANIFEST_SHA=SHARD_MANIFEST_SHA

def require(cond,message):
    if not bool(cond):
        raise ValueError("MOIS_BULK128_500PX_AGGREGATE: "+message)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(4<<20),b""):
            h.update(chunk)
    return h.hexdigest()

def audit(source_manifest,source_summary,workers,output,*,expected_parts=128,
          verify_sha=True):
    source_manifest,source_summary,workers,output=map(
        Path,(source_manifest,source_summary,workers,output))
    require(not output.exists(),"Do not replace previously signed source pixel audit")
    if verify_sha:
        require(sha(source_manifest)==SOURCE_MANIFEST_SHA
                and sha(source_summary)==SOURCE_SUMMARY_SHA,
                "Bulk 128 source manifest or per-part photo targets changed")
    m=json.loads(source_manifest.read_text())
    require(m.get("schema")=="mois_full_source_observation_sampling_cohort_v1"
            and m.get("sampling_schema")==SOURCE_SCHEMA
            and m.get("source_128_compute_shards")==128
            and m.get("P5_targeted_photo_attempts")==0
            and m.get("source_tier")==10,
            "Wrong original source population or biased sampling schema")
    with source_summary.open(newline="",encoding="utf-8") as f:
        source_rows=list(csv.DictReader(f))
    source={int(x["part"]):x for x in source_rows}
    require(len(source_rows)==len(source)==128
            and sum(int(x["photo_targets"]) for x in source_rows)==
            m["new_photo_attempts_planned"],
            "Original photo-source per-part target census does not reconcile")
    require(1<=expected_parts<=128 and
            (expected_parts==128 or not verify_sha),
            "Partial synthetic audit must explicitly disable immutable full-run check")
    seen_photos=set();seen_observations=set()
    parts=[];aggregate=collections.Counter();by_design=collections.defaultdict(collections.Counter)
    cross=collections.defaultdict(collections.Counter)
    output.mkdir(parents=True)
    all_file=output/"bulk128_true_medium_boxLab_photo_IDs.csv.gz"
    with gzip.open(all_file,"wt",encoding="utf-8",newline="") as dest:
        writer=csv.DictWriter(dest,fieldnames=COLUMNS,lineterminator="\n")
        writer.writeheader()
        for part in range(expected_parts):
            directory=workers/f"mois-bulk128-med500-part-{part:03d}"
            receipt=directory/"bulk128_medium_photo_completion_receipt.json"
            csv_file=directory/"photo_ID_matched_reference_vs_targeted_boxLab.csv.gz"
            require(receipt.is_file() and csv_file.is_file(),
                    "Incomplete actual original medium photo shard: "+str(part))
            r=json.loads(receipt.read_text())
            require(r.get("schema")==PART_SCHEMA
                    and r.get("original_source_part")==part
                    and r.get("original_model_detector_weight_sha256")==WEIGHT_SHA
                    and r.get("original_source_metadata_snapshot_SHA256")==
                    SOURCE_MANIFEST_SHA
                    and r.get("actual_pixel_terminal_photo_and_cohort_SHA256")==sha(csv_file)
                    and r.get("expected_original_photos")==int(source[part]["photo_targets"])
                    and r.get("actually_medium_photo_IDs_attempted")==int(source[part]["photo_targets"])
                    and r.get("reference_and_adaptive_targeted_are_distinct")
                    and r.get("real_strict9_flower_pixel_vectors_generated")==0
                    and r.get("original_FCP_four_hue_full_image_ROI_inferred") is False
                    and r.get("raw_licensed_image_bytes_archived") is False,
                    "One attempted original photo shard is inconsistent with signed ROIs/source")
            n=0;decoded=0;roles=collections.Counter()
            with gzip.open(csv_file,"rt",encoding="utf-8",newline="") as f:
                rd=csv.DictReader(f)
                require(tuple(rd.fieldnames or ())==COLUMNS,
                        "Actual medium pixel source output fieldnames changed")
                for row in rd:
                    pid,oid=row["photo_id"],row["observation_id"]
                    require(pid.isdecimal() and oid and
                            pid not in seen_photos and oid not in seen_observations,
                            "One original photo/observation ID duplicated across 128 shards")
                    seen_photos.add(pid);seen_observations.add(oid)
                    role=row["design_cohort"]
                    require(role in ("reference","targeted")
                            and row["screen_status"] in (
                                "candidate","negative","uncertain","error")
                            and row["unverified_floral_box_not_petal_ROI"]=="True",
                            "Photo detection or petal verification field was silently converted")
                    if role=="targeted":
                        require(row["original_priority"]!="P5_existing_10plus_strict9"
                                and row["original_priority"] not in ("","unknown"),
                                "Previously deeply measured taxa selected as gap source")
                    writer.writerow(row)
                    n+=1;decoded+=(row["medium_decoded"]=="True")
                    roles[role]+=1
                    by_design[role][row["screen_status"]]+=1
                    if row["medium_decoded"]=="True":
                        by_design[role]["medium_decoded"]+=1
                    if row["floral_box_main_L_star"]:
                        by_design[role]["noncalibrated_Lab_modal"]+=1
                    cross[(role,row["latitude_region"],row["local_season"])][row["screen_status"]]+=1
            require(n==r["expected_original_photos"] and decoded==
                    r["real_medium_photo_pixels_decoded"]
                    and all(roles[k]==r["by_design"][k]["photo_attempts"]
                            for k in roles),
                    "Authenticated per-part original photo counts did not reconcile")
            aggregate["photo_attempts"]+=n
            aggregate["medium_decoded"]+=decoded
            aggregate["real_flower_candidate_box"]+=sum(
                r["by_design"].get(role,{}).get("flower_boxes",0)
                for role in ("reference","targeted"))
            parts.append({
                "part":part,
                "original_photo_attempts":n,
                "actual_medium_decoded":decoded,
                "reference_photo_attempts":roles["reference"],
                "targeted_photo_attempts":roles["targeted"],
                "scan_receipt_SHA256":sha(receipt),
                "photo_colour_status_SHA256":sha(csv_file),
            })
    if expected_parts==128:
        require(aggregate["photo_attempts"]==131072==
                m["new_photo_attempts_planned"]
                and len(seen_photos)==len(seen_observations)==131072
                and sum(by_design[k][x] for k in by_design for x in (
                    "candidate","negative","uncertain","error"))==131072,
                "128 source picture pixels were not all retrieved or marked unknown")
    with (output/"bulk128_picture_success_by_source_part.csv").open(
            "w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(parts[0]),lineterminator="\n")
        w.writeheader();w.writerows(parts)
    rows=[]
    for (role,zone,season),counts in sorted(cross.items()):
        rows.append({
            "design_cohort":role,
            "latitude_region":zone,
            "recorded_local_season":season,
            "original_photo_attempts":sum(counts.values()),
            **{status:counts[status] for status in (
                "candidate","negative","uncertain","error")},
        })
    with (output/"bulk128_true_photo_regions_and_seasons.csv").open(
            "w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n")
        w.writeheader();w.writerows(rows)
    receipt={
        "schema":SCHEMA,
        "original_full_source_photo_opportunities":93_209_351,
        "128_frozen_original_taxon_photo_partitions_authenticated":expected_parts,
        "new_real_medium_picture_attempts":aggregate["photo_attempts"],
        "real_medium_RGB_downloaded_decoded":aggregate["medium_decoded"],
        "actual_original_photo_flower_candidate_boxes":aggregate["real_flower_candidate_box"],
        "reference_and_gap_targeted_original_photo_states":{
            key:dict(value) for key,value in sorted(by_design.items())
        },
        "prior_P5_gap_targeted_images":0,
        "new_original_strict9_color_ROI_count":0,
        "unverified_noncalibrated_photo_box_colour_only":True,
        "source_design_is_not_unbiased_global_flora_prevalence":True,
        "all_unclassifiable_photo_IDs_preserved_as_unknown":True,
        "new_image_downloads_or_ROI_inferences_by_collector":0,
        "no_raw_licensed_photo_bytes_redistributed":True,
        "original_bulk128_precolour_manifest_SHA256":sha(source_manifest),
        "source_exact_shard_photo_counts_SHA256":sha(source_summary),
        "actual_500px_colour_status_photoID_file_SHA256":sha(all_file),
    }
    (output/"bulk128_real_medium_boxLab_completion_receipt.json").write_text(
        json.dumps(receipt,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("MOIS_BULK128_REAL_IMAGE_MEDIUM_PHOTO_CENSUS_AGGREGATED",
          json.dumps({k:receipt[k] for k in (
              "128_frozen_original_taxon_photo_partitions_authenticated",
              "new_real_medium_picture_attempts",
              "real_medium_RGB_downloaded_decoded",
              "actual_original_photo_flower_candidate_boxes")},sort_keys=True),flush=True)
    return receipt

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("source-manifest","source-summary","workers","output"):
        p.add_argument("--"+k,required=True,type=Path)
    a=p.parse_args()
    audit(a.source_manifest,a.source_summary,a.workers,a.output)

if __name__=="__main__":
    main()
