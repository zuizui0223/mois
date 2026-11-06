#!/usr/bin/env python3
"""One WARM YOLO detector scans eight signed tier×photo partitions once.

Deduplicates previously fully scanned original photo packets by matching
photo/source/model receipts; not biological preselection. Source original tier
archive, original part gzips and original photo ID counts are SHA verified.
Actual pixels scanned in one bounded provider-paced process. The source group
is only one 8-part piece of all 93,209,351 photographs, not 93m completion.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import tarfile
import tempfile

from scripts.acquisition.mois_full93m_medium_pixel_scan import (
    WarmFlowerBoxDetector, Pacer, download_medium, scan_source, sha_file,
    SCHEMA as SCAN_SCHEMA,
)
from scripts.acquisition.mois_onepass_full93m_target_tiers import (
    RESULT_SCHEMA, COLUMNS, EXPECTED_PHOTOS,
)

SCHEMA = "mois_source_authenticated_warm_8part_medium_image_group_v1"
REGISTRY_SCHEMA = "mois_complete_full93m_original_photo_target_registry_v1"

def require(condition,why):
    if not condition:
        raise ValueError("MOIS_FULL93M_WARM_GROUP: "+why)

def run(registry, group_receipt, archive, output, weights, weight_sha, *,
        tier=0, group=0, max_photo_ids=1000, rps=2.0, workers=8,
        reusable_part0=None, device="cpu", conf=0.02):
    registry,group_receipt,archive,output,weights=map(
        Path,(registry,group_receipt,archive,output,weights))
    require(not output.exists(), "Never overwrite a sealed group source")
    reg=json.loads(registry.read_text(encoding="utf-8"))
    proof=json.loads(group_receipt.read_text(encoding="utf-8"))
    require(reg.get("schema")==REGISTRY_SCHEMA and
            reg.get("exact_original_photo_IDs_in_all_disjoint_tiers")==EXPECTED_PHOTOS
            and reg.get("native_vs_introduced_not_hard_filtered") is True
            and reg.get("source_only_NO_photo_pixels_or_hue_measured") is True,
            "Not the complete authentic all93m source registry")
    require(proof.get("schema")==RESULT_SCHEMA
            and proof.get("source_snapshot_SHA256")==reg["full_source_manifest_sha256"],
            "Source group receipt does not match full93m registry")
    require(0<=tier<14 and 0<=group<16 and
            archive.name==f"tier-{tier:02d}-group-{group:02d}.tar" and
            reg["published_tier_archive_SHA256"].get(archive.name)==sha_file(archive),
            "Unsigned or wrong immutable original-photo tier archive")
    parts={int(x["part"]):x for x in proof["parts"]}
    require(set(parts)==set(range(group*8,group*8+8)),
            "Eight complete original photo subparts must be present")
    total_expected=sum(int(part["counts_by_disjoint_tier"][str(tier)])
                       for part in parts.values())
    require(total_expected==int(proof["counts_by_disjoint_tier"][str(tier)])
            and total_expected <= max_photo_ids and total_expected > 0,
            "Original photo group exceeds explicit bounded real pixel-work gate")
    require(0<rps<=4 and 1<=workers<=16,
            "Unbounded provider rate or concurrency is forbidden")
    detector=WarmFlowerBoxDetector(
        weights,weight_sha,0,"item",device=device,conf=conf)
    pacer=Pacer(rps)
    output.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="mois-all93m-warm-source-",dir=output.parent) as tmp:
        working=Path(tmp)
        source=working/"source"
        stage=working/"outputs"
        source.mkdir();stage.mkdir()
        expected_names={
            f"tier-{tier:02d}/part-{i:03d}.csv.gz" for i in parts
        }
        with tarfile.open(archive) as archive_obj:
            members={m.name:m for m in archive_obj if m.isfile()}
            require(set(members)==expected_names,
                    "Source archive must contain exactly 8 original part gzips")
            for number,part in sorted(parts.items()):
                name=f"tier-{tier:02d}/part-{number:03d}.csv.gz"
                compressed=archive_obj.extractfile(members[name]).read()
                in_file=source/f"part-{number:03d}.csv.gz"
                in_file.write_bytes(compressed)
                require(sha_file(in_file)==part["ranked_target_csv_sha256"][str(tier)],
                        "Original source-part compressed photo IDs SHA changed")
        decoded=candidates=negatives=uncertains=errors=total=0
        new_attempts=0
        receipts=[]
        for number,part in sorted(parts.items()):
            expected=int(part["counts_by_disjoint_tier"][str(tier)])
            sub=stage/f"part-{number:03d}"
            reused=0
            if number==0 and reusable_part0 is not None:
                original=Path(reusable_part0)
                require(original.is_dir(), "Prior scanned original part missing")
                proof_prior=json.loads((original/"complete_receipt.json").read_text())
                require(proof_prior.get("source_role")=="full93m"
                        and proof_prior.get("source_part")==0
                        and proof_prior.get("source_tier")==tier
                        and proof_prior.get("expected_photo_ids")==expected
                        and proof_prior.get("source_gzip_sha256")
                        ==part["ranked_target_csv_sha256"][str(tier)]
                        and proof_prior.get("detector_identity")
                        =="yolo11:"+weight_sha+":item:"+str(detector.conf),
                        "Previous original pixels cannot be reused under different model/source")
                shutil.copytree(original,sub)
                reused=expected
            result=scan_source(
                source/f"part-{number:03d}.csv.gz",sub,
                expected_source_sha256=part["ranked_target_csv_sha256"][str(tier)],
                expected_count=expected,part=number,tier=tier,
                detector=detector,download=lambda url:download_medium(url,pacer),
                batch_size=512,workers=workers,source_role="full93m")
            require(result["all_original_input_photos_attempted"] is True
                    and result["actual_unique_source_photo_ids_attempted"]==expected,
                    "Partial original photo input cannot be counted complete")
            decoded+=result["medium_image_pixels_decoded"]
            candidates+=result["screen_counts"].get("candidate",0)
            negatives+=result["screen_counts"].get("negative",0)
            uncertains+=result["screen_counts"].get("uncertain",0)
            errors+=result["screen_counts"].get("error",0)
            total+=expected;new_attempts+=expected-reused
            receipts.append({
                "source_part":number,"expected_photos":expected,
                "reused_prior_source_images":reused,
                "medium_images_decoded":result["medium_image_pixels_decoded"],
                "original_source_complete_receipt_sha256":sha_file(sub/"complete_receipt.json"),
            })
        require(total==total_expected
                and 0<=decoded<=total
                and candidates+negatives+uncertains+errors==total,
                "Original complete group status or pixel totals changed")
        summary={
            "schema":SCHEMA,
            "original_full_source_93m_metadata_photo_ids":EXPECTED_PHOTOS,
            "source_complete_93m":False,
            "source_group_tier":tier,"source_group":group,
            "actual_real_photo_stage1_completed_in_group":total,
            "previous_same_full_photo_model_photo_ids_reused_without_redownload":
                total-new_attempts,
            "new_original_medium_pixel_attempts_in_this_run":new_attempts,
            "genuinely_decoded_medium_photos_in_group":decoded,
            "stage1_candidate_photos":candidates,
            "stage1_negative_photos_NOT_removed":negatives,
            "stage1_uncertain":uncertains,
            "stage1_error":errors,
            "frozen_yolo_sha256":weight_sha,
            "frozen_yolo_one_warm_model_initialization":1,
            "medium_S3_per_worker_request_rate_ceiling":rps,
            "source_registry_sha256":sha_file(registry),
            "source_group_receipt_sha256":sha_file(group_receipt),
            "source_tier_tar_sha256":sha_file(archive),
            "completed_eight_part_shas":receipts,
            "no_hue_label_imputed_or_source_rois_relabelled":True,
            "no_raw_image_bytes_published":True,
        }
        (stage/"whole_group_receipt.json").write_text(
            json.dumps(summary,sort_keys=True,indent=2)+"\n",encoding="utf-8")
        stage.rename(output)
        return summary

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("registry","group-receipt","archive","weights","output"):
        p.add_argument("--"+name,required=True,type=Path)
    p.add_argument("--weight-sha256",required=True)
    p.add_argument("--tier",type=int,default=0)
    p.add_argument("--group",type=int,default=0)
    p.add_argument("--max-photo-ids",type=int,default=1000)
    p.add_argument("--requests-per-second",type=float,default=2.0)
    p.add_argument("--workers",type=int,default=8)
    p.add_argument("--reusable-part0",type=Path)
    p.add_argument("--device",default="cpu")
    x=p.parse_args()
    d=run(x.registry,x.group_receipt,x.archive,x.output,x.weights,
          x.weight_sha256,tier=x.tier,group=x.group,
          max_photo_ids=x.max_photo_ids,rps=x.requests_per_second,
          workers=x.workers,reusable_part0=x.reusable_part0,device=x.device)
    print("MOIS_REAL_FULL93M_STAGE1_WARM_8PART_COMPLETE",json.dumps(d,sort_keys=True))

if __name__=="__main__":
    main()
