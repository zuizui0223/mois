#!/usr/bin/env python3
"""Admit genuine, FINISHED FCP wave4 photo-pixel partitions without waiting for all 49800.

The original SHA-locked FCP worker assignment and blinded measurement IDs were
fixed before image outcome inspection. Every ingested photo has a completed ROI
terminal, source-matching original metadata and an original classifier state.
This is a PRELIMINARY executable source expansion, not a claim that unfinished
source partitions or FCP D/W within-species polymorphism were measured.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import re
from collections import Counter
from pathlib import Path

EXPECTED_SOURCE_PHOTOS = 49800
EXPECTED_SOURCE_SPECIES = 498
EXPECTED_MODEL_PARTITIONS = 256
MIN_COMPLETE_PARTITIONS = 16
BASE_PHOTOS = 252408
BASE_HUE_LABELS = 122383
HUES = ("white","yellow_orange","red_pink","blue_purple")
FIREWALL_METADATA_SHA = "3e270153134830496357f537f0dee0544647d492dcf751e856988d36ce1e502a"
FCP_AUTHORITATIVE_SOURCE_SHA = "8d98b14856dd2539b54e1db287aa08adea9c191e112223ff380cb49ac0f8a554"
FCP_RUN_ID = 38992543115
FCP_COMMIT = "ef448c9514e0c702e7887e00b288ee966728b0d5"
PATTERN = re.compile(r"b([01])_partition_s(\d\d)_p(\d\d)\.(csv|json)$")
COLS = (
    "photo_id","inat_taxon_id","species","morph","source_original_measured",
    "old_mois_or_new_fcp","measurement_status","latitude","longitude",
    "photo_license","observation_id","observed_on","time_observed_at",
    "wc_elevation_m","soil_pH","soil_SOC","soil_N","soil_clay",
    "soil_available_water_proxy","wc_bio1","wc_bio5","wc_bio12","wc_bio15"
)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for payload in iter(lambda:stream.read(1<<20),b""):
            h.update(payload)
    return h.hexdigest()


def parse_csv(path):
    with Path(path).open(newline="",encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def read_firewall(root):
    receipt_path=root/"measurement_firewall_manifest.json"
    m=json.loads(receipt_path.read_text(encoding="utf-8"))
    require(
        m["schema"]=="fcp_same_149900_method_wave4_firewall_v1"
        and m["new_photo_rows"]==EXPECTED_SOURCE_PHOTOS
        and m["new_species"]==EXPECTED_SOURCE_SPECIES
        and m["partitions"]==EXPECTED_MODEL_PARTITIONS
        and m["authorized_metadata_sha256"]==FCP_AUTHORITATIVE_SOURCE_SHA
        and m["metadata_join_sha256"]==FIREWALL_METADATA_SHA
        and m["blind_worker_columns"]==["measurement_id","image_filename","photo_license"]
        and m["original_frozen_results_changed"] is False
        and m["photo_pixels_opened"] is False
        and m["source_coordinates_in_worker"] is False
        and m["source_photo_urls_in_worker"] is False,
        "FCP original wave4 source / no-outcome worker firewall changed",
    )
    original=root/"sealed_keys/metadata_join_key.csv"
    require(digest(original)==FIREWALL_METADATA_SHA,
            "FCP wave4 photographic source ID metadata checksum changed")
    meta=parse_csv(original)
    require(len(meta)==EXPECTED_SOURCE_PHOTOS,"FCP original selected photograph source not complete")
    keys={}
    pids=set()
    for row in meta:
        mid,pid=row["measurement_id"],row["photo_id"]
        require(mid and mid not in keys and pid.isdecimal() and pid not in pids,
                "Duplicate or missing original FCP wave4 photo/measurement ID")
        require(row["inat_taxon_id"].isdecimal() and row["photo_license"]
                and row["species"] and row["observation_id"],"Invalid original FCP photo identity/license")
        keys[mid]=row
        pids.add(pid)
    require(len({r["inat_taxon_id"] for r in meta})==EXPECTED_SOURCE_SPECIES,
            "FCP original 498 source taxa changed")
    assignments=parse_csv(root/"partition_assignments.csv")
    require(len(assignments)==EXPECTED_SOURCE_PHOTOS,
            "FCP photo partition lottery/assignment incomplete")
    assigned={}
    for r in assignments:
        mid=r["measurement_id"]
        require(mid in keys and mid not in assigned,"Original worker matrix includes wrong or duplicate photo")
        b,s,p=(int(r[x]) for x in ("batch","semantic_shard","compute_partition"))
        require(0<=b<=1 and 0<=s<32 and 0<=p<4,
                "Original FCP pixel measurement matrix index invalid")
        require(int(keys[mid]["measurement_batch"])==b,
                "Original source-bound randomized batch ID changed")
        assigned[mid]=(b,s,p)
    require(set(keys)==set(assigned),"Unassigned original wave4 photograph IDs")
    return keys,assigned,m


def sealed_terminals(root,keys,assignments):
    found={}
    for f in root.rglob("*"):
        if not f.is_file():continue
        match=PATTERN.fullmatch(f.name)
        if not match:
            require("partition_s" not in f.name,
                    "Unrecognized original FCP pixel classifier terminal file")
            continue
        b,s,p,fmt=int(match[1]),int(match[2]),int(match[3]),match[4]
        part=(b,s,p)
        require(part not in found or fmt not in found[part],
                f"Duplicated original FCP result partition {part}:{fmt}")
        found.setdefault(part,{})[fmt]=f
    require(len(found)>=MIN_COMPLETE_PARTITIONS,
            "Less than sixteen actual FCP source pixel partitions have finished")
    require(len(found)<=EXPECTED_MODEL_PARTITIONS,
            "More than 256 source photo worker partitions returned")
    for part,files in found.items():
        require(set(files)=={"csv","json"},
                f"Source ROI partition files incomplete {part}")
    original=[]
    report=[]
    used=set()
    for part,files in sorted(found.items()):
        b,s,p=part
        q=json.loads(files["json"].read_text())
        require(
            q.get("status")=="complete_random_photo_first_terminal_partition"
            and q.get("semantic_shard")==s and q.get("compute_partition")==p
            and q.get("source_urls_present") is False
            and q.get("species_present") is False
            and q.get("coordinates_present") is False,
            f"FCP original image pixel ROI source terminal not genuine {part}"
        )
        table=parse_csv(files["csv"])
        require(q["terminal_rows"]==q["selected_rows"]==len(table) and
                q["classified_rows"]+q["mixed_uncertain_rows"]==len(table),
                f"Actual FCP source photo count is inconsistent {part}")
        for row in table:
            mid=row["measurement_id"]
            require(mid not in used and mid in keys and assignments[mid]==part,
                    f"An original photo ID was not selected for its completed terminal {part}")
            used.add(mid)
            morph,status=row["morph"],row["measurement_status"]
            classified=morph in HUES
            require(morph in (*HUES,"mixed_uncertain") and
                    classified==(status=="classified_four_state_morph"),
                    "Original ROI-v4 four-state photo status/morph diverged")
            if status=="classified_four_state_morph":
                require(row["image_sha256"].strip()!="",
                        "A classifiable ROI has no genuine original photo RGB SHA256")
            source=keys[mid]
            out={k:source.get(k,"") for k in COLS}
            out["morph"]=morph
            out["measurement_status"]=status
            out["source_original_measured"]="fcp_wave4_partial_original_ROI_v4"
            out["old_mois_or_new_fcp"]="fcp_additional_ROI"
            original.append(out)
        report.append({
            "batch":b,"semantic_shard":s,"partition":p,
            "actual_source_original_photo_rows":len(table),
            "actual_classifier_four_hue_rows":q["classified_rows"],
            "unclassifiable_original_photo_rows":q["mixed_uncertain_rows"],
            "terminal_csv_sha256":digest(files["csv"]),
            "terminal_receipt_sha256":digest(files["json"]),
        })
    require(len(used)==len(original)==sum(q["actual_source_original_photo_rows"] for q in report),
            "Actual completed-photo terminal count unreconciled")
    return original,report


def merge(old:Path,old_receipt:Path,firewall:Path,
          partitions:Path,fcp_run_metadata:Path,output:Path):
    require(not output.exists() or not any(output.iterdir()),
            "No replacement of a previous original-photo measured record")
    prior=json.loads(old_receipt.read_text())
    require(prior["all_original_unique_photos"]==BASE_PHOTOS
            and prior["all_four_hue_classified_photos"]==BASE_HUE_LABELS
            and prior["original_model_SHA_verified_for_every_photo"] is True
            and digest(old)==prior["all_real_source_photo_hue_csv_sha256"],
            "Previously completed MOIS 252408-photo source SHA or labels changed")
    run=json.loads(fcp_run_metadata.read_text())
    require(run["id"]==FCP_RUN_ID and run["head_sha"]==FCP_COMMIT
            and run["head_repository"]["full_name"]=="zuizui0223/fcp",
            "FCP pixel-photo source branch/version has changed")
    src,assigned,m=read_firewall(firewall)
    genuine,partition_receipts=sealed_terminals(partitions,src,assigned)
    output.mkdir(parents=True,exist_ok=True)
    allids=set()
    baseline={}
    old_count=old_hue=0
    newly=Counter()
    joined=output/"mois_fcp_wave4_PARTIAL_completed_pixel_original_photo_union.csv.gz"
    with joined.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw,mode="wb",mtime=0,filename="") as gz:
            with io.TextIOWrapper(gz,encoding="utf-8",newline="") as stream:
                writer=csv.DictWriter(stream,fieldnames=COLS,lineterminator="\n")
                writer.writeheader()
                with gzip.open(old,"rt",encoding="utf-8",newline="") as original:
                    reader=csv.DictReader(original)
                    require(tuple(reader.fieldnames or ())==COLS,
                            "Original 252408 photographic source columns changed")
                    for row in reader:
                        pid=row["photo_id"]
                        require(pid not in allids and pid,"Original MOIS photo ID duplicated")
                        allids.add(pid);baseline[pid]=(row["morph"],row["inat_taxon_id"])
                        old_count+=1;old_hue+=row["morph"] in HUES
                        writer.writerow(row)
                require((old_count,old_hue)==(BASE_PHOTOS,BASE_HUE_LABELS),
                        "252408 prior original photographic source changed")
                for row in genuine:
                    pid=row["photo_id"]
                    if pid in allids:
                        newly["photo_ID_overlaps_with_prior_MOIS"]+=1
                        original_hue,original_tax=baseline[pid]
                        newly["classifiable_original_overlap_disagreement"]+=int(
                            original_hue in HUES and row["morph"] in HUES
                            and original_hue!=row["morph"])
                        newly["taxon_ID_overlap_disagreement"]+=int(
                            original_tax!=row["inat_taxon_id"])
                        continue
                    allids.add(pid)
                    newly["original_new_photo_IDs_admitted"]+=1
                    newly["new_classified_four_hue_original_photos"]+=row["morph"] in HUES
                    writer.writerow(row)
    require(len(allids)==BASE_PHOTOS+newly["original_new_photo_IDs_admitted"],
            "Original FCP photo ID overlap/union count does not reconcile")
    require(newly["original_new_photo_IDs_admitted"]+
            newly["photo_ID_overlaps_with_prior_MOIS"]==len(genuine),
            "Lost a real original FCP photo terminal at the join")
    require(newly["original_new_photo_IDs_admitted"]>0,
            "FCP partial photo source contains no new MOIS original photo IDs")
    receipt={
      "schema":"mois_authentic_completed_original_FCP_wave4_ROI_PARTITIONS_source_v1",
      "status":("COMPLETE_ACTUAL_PIXELS" if len(partition_receipts)==256
                else "PARTIAL_ACTUAL_PIXELS_DO_NOT_CLAIM_FULL_FCP_49800"),
      "FCP_live_original_actions_run":FCP_RUN_ID,
      "FCP_original_run_head_sha":FCP_COMMIT,
      "FCP_run_status_when_snapshot":run["status"],
      "FCP_model_source_frozen":True,
      "FCP_photo_source_firewall_metadata_sha256":FIREWALL_METADATA_SHA,
      "FCP_source_original_49800_photographs":EXPECTED_SOURCE_PHOTOS,
      "FCP_completed_ROI_partitions":len(partition_receipts),
      "FCP_total_ROI_partitions":EXPECTED_MODEL_PARTITIONS,
      "FCP_real_completed_original_image_pixel_rows":len(genuine),
      "FCP_real_completed_original_four_hue_classifications":sum(r["morph"] in HUES for r in genuine),
      "FCP_real_terminal_partitions":partition_receipts,
      "original_MOIS_prior_photo_IDs":BASE_PHOTOS,
      "original_MOIS_prior_classified_hues":BASE_HUE_LABELS,
      "overlap_photographic_identifiers":newly["photo_ID_overlaps_with_prior_MOIS"],
      "additional_new_original_photo_IDs":newly["original_new_photo_IDs_admitted"],
      "additional_new_real_four_hue_classified_photo_IDs":newly["new_classified_four_hue_original_photos"],
      "original_FCP_white_or_other_four_hue_overlap_disagreement":newly["classifiable_original_overlap_disagreement"],
      "original_FCP_taxon_ID_overlap_disagreement":newly["taxon_ID_overlap_disagreement"],
      "combined_unique_photo_IDs":len(allids),
      "combined_four_hue_classified_photo_IDs":BASE_HUE_LABELS+newly["new_classified_four_hue_original_photos"],
      "source_photo_level_csv_sha256":digest(joined),
      "source_old_MOIS_photo_sha256":digest(old),
      "source_original_FCP_firewall_receipt_sha256":digest(firewall/"measurement_firewall_manifest.json"),
      "new_FCP_pixel_source_photo_classifications_NOT_independent_provider_replication":True,
      "FCP_genetic_or_within_species_D_W_50km_polymorphism_statistics_imported":False,
      "original_MOIS_100543_moisture_999genus_maxT_NOT_YET_refitted_on_this_partial_source":True,
      "missing_FCP_wave4_partitions_are_not_hues":True,
      "partial_FCP_pixel_completion_order_may_not_be_geographically_or_taxonomically_representative":True,
    }
    (output/"mois_fcp_wave4_PARTIAL_real_completed_original_ROI_receipt.json").write_text(
        json.dumps(receipt,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("MOIS_FCP_WAVE4_REAL_COMPLETED_PIXEL_PARTITIONS",json.dumps({
        k:receipt[k] for k in (
           "FCP_completed_ROI_partitions","FCP_real_completed_original_image_pixel_rows",
           "additional_new_original_photo_IDs",
           "additional_new_real_four_hue_classified_photo_IDs",
           "combined_unique_photo_IDs","combined_four_hue_classified_photo_IDs")},
        sort_keys=True),flush=True)
    return receipt


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--prior-252408-original-photos",type=Path,required=True)
    ap.add_argument("--prior-252408-photo-source-receipt",type=Path,required=True)
    ap.add_argument("--fcp-wave4-original-blind-firewall",type=Path,required=True)
    ap.add_argument("--fcp-wave4-completed-original-pixel-partitions",type=Path,required=True)
    ap.add_argument("--fcp-wave4-original-run-snapshot",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args()
    merge(a.prior_252408_original_photos,a.prior_252408_photo_source_receipt,
          a.fcp_wave4_original_blind_firewall,
          a.fcp_wave4_completed_original_pixel_partitions,
          a.fcp_wave4_original_run_snapshot,a.output)
