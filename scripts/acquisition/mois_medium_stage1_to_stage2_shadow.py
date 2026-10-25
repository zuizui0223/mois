#!/usr/bin/env python3
"""Convert verified actual medium-pixel Stage1 packets into original-photo Stage2.

SQLite disk-backed sort by numeric photo ID, not 93m in-memory sets.
No images fetched, hue relabelling, or negative-screen exclusions.
Feeds existing MOIS signed Stage2 probability gate in SHADOW mode only.
"""
from __future__ import annotations
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import sqlite3
import tempfile

from scripts.acquisition.mois_93m_stage2_probability_queue import (
    plan as make_stage2, sha256, SCHEMA as STAGE2_SCHEMA,
)
from scripts.acquisition.mois_full93m_medium_pixel_scan import (
    SCHEMA as STAGE1_SCHEMA, SOURCE_ARM, COLUMNS_OUT,
)

INPUT = ("photo_id", "observation_id", "source_arm",
         "selection_probability", "audit_stratum", "screen_status")
STATES = ("candidate", "negative", "uncertain", "error")
SCHEMA = "mois_verified_medium_pixel_stage1_to_existing_stage2_shadow_bridge_v1"

def require(condition, why):
    if not condition:
        raise ValueError("MOIS_STAGE1_STAGE2_BRIDGE: " + why)

def bridge(source_dir, output, *, expected_count=None):
    source_dir = Path(source_dir); output = Path(output)
    require(not output.exists(), "Refuse overwrite")
    receipt = json.loads((source_dir / "complete_receipt.json").read_text())
    require(receipt.get("schema") == STAGE1_SCHEMA
            and receipt.get("source_role") == "full93m"
            and receipt.get("all_original_input_photos_attempted") is True
            and receipt.get("calibration_subset_not_global_census") is False,
            "Only completed real full93m source shard may enter this bridge")
    count = int(receipt["expected_photo_ids"])
    require(expected_count is None or count == expected_count,
            "Expected exact photo count differs from signed complete source")
    require(count > 0 and receipt["actual_unique_source_photo_ids_attempted"] == count,
            "Incomplete original photo attempts")
    root = output.parent
    root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="mois-stage1-stage2-bridge-",dir=root) as tmp:
        temp = Path(tmp)
        database = sqlite3.connect(temp / "photo.sqlite")
        database.execute("CREATE TABLE photos (photo_id INTEGER PRIMARY KEY, "
            "observation_id TEXT NOT NULL, source_arm TEXT NOT NULL,"
            "selection_probability TEXT NOT NULL, audit_stratum TEXT NOT NULL,"
            "screen_status TEXT NOT NULL)")
        packets = sorted(source_dir.glob("packet-*.csv.gz"))
        require(len(packets) == receipt["packets_sealed"],
                "Missing or repeated signed terminal packets")
        stats=Counter(); scanned=0
        for packet in packets:
            proof = json.loads(packet.with_suffix("").with_suffix(".receipt.json").read_text())
            require(proof.get("schema") == STAGE1_SCHEMA
                    and proof.get("source_role") == "full93m"
                    and proof.get("source_gzip_sha256") == receipt["source_gzip_sha256"]
                    and proof.get("detector_identity") == receipt["detector_identity"]
                    and sha256(packet) == proof.get("output_sha256"),
                    "Mismatch in verified Stage1 packet or model")
            within=0
            with gzip.open(packet,"rt",encoding="utf-8",newline="") as stream:
                reader=csv.DictReader(stream)
                require(tuple(reader.fieldnames or ()) == COLUMNS_OUT,
                        "Stage1 output columns were changed")
                for row in reader:
                    pid=str(row["photo_id"])
                    require(pid.isdecimal() and 0 < int(pid) < 2**63,
                            "Malformed full-source photo ID")
                    require(row["source_arm"] == SOURCE_ARM
                            and row["source_photo_inclusion_probability"] == "1"
                            and row["screen_status"] in STATES
                            and row["observation_id"]
                            and row["audit_stratum"],
                            "Stage1 source identity or candidate status changed")
                    try:
                        database.execute("INSERT INTO photos VALUES (?,?,?,?,?,?)",
                            (int(pid),row["observation_id"], row["source_arm"],
                             "1",row["audit_stratum"],row["screen_status"]))
                    except sqlite3.IntegrityError as ex:
                        raise ValueError("Duplicate photo ID across source packets: "+pid) from ex
                    stats[row["screen_status"]]+=1
                    within+=1
            require(within == proof["photo_ids"],
                    "Stage1 archived packet row count changed")
            scanned+=within
        require(scanned==count and dict(stats)==receipt["screen_counts"],
                "Stage1 all statuses and source photo IDs not reconciled")
        target_csv = temp / "source_photo_medium_stage1_complete_sorted.csv"
        with target_csv.open("w",newline="",encoding="utf-8") as stream:
            writer=csv.writer(stream,lineterminator="\n")
            writer.writerow(INPUT)
            for row in database.execute(
                "SELECT CAST(photo_id AS TEXT),observation_id,source_arm,"
                "selection_probability,audit_stratum,screen_status FROM photos ORDER BY photo_id"
            ):
                writer.writerow(row)
        database.close()
        # Full source shadow keeps negatives; 1.0 is never a production
        # non-flower exclusion rate and requires no retrospective review waiver.
        actual_strata=set()
        with target_csv.open(newline="",encoding="utf-8") as stream:
            for row in csv.DictReader(stream): actual_strata.add(row["audit_stratum"])
        policy=temp/"negative_sampling_SHADOW_only.json"
        policy.write_text(json.dumps({
            "schema":"mois_93m_prespecified_negative_audit_policy_v1",
            "negative_audit_rates":{k:1.0 for k in sorted(actual_strata)},
        },sort_keys=True,indent=2)+"\n",encoding="utf-8")
        stage2_dir=temp/"stage2_queue"
        stage2=make_stage2(target_csv,policy,stage2_dir,
                           mode="shadow",expected_input_sha256=sha256(target_csv))
        require(stage2["schema"] == STAGE2_SCHEMA
                and stage2["counts"]["all_screened"] == count
                and stage2["counts"]["queued_full_roi"] == count
                and stage2["counts"]["negative_not_sampled"] == 0,
                "Shadow Stage2 unexpectedly deleted original images")
        info={
            "schema":SCHEMA,
            "source_original_photo_count":count,
            "source_stage1_complete_receipt_sha256":sha256(source_dir/"complete_receipt.json"),
            "source_medium_photo_attempts":scanned,
            "original_stage1_states":dict(stats),
            "original_photo_identity_duplicates":0,
            "source_original_photo_subset_is_not_complete_93m":True,
            "stage2_mode":"shadow",
            "stage2_full_ROI_candidates_including_screen_negatives":count,
            "screen_negatives_excluded":0,
            "original_4hue_or_petal_phenotypes_assigned":0,
            "original_image_downloads_by_bridge":0,
            "sorted_stage1_sha256":sha256(target_csv),
            "stage2_receipt_sha256":sha256(stage2_dir/"stage2_assignment_receipt.json"),
        }
        (temp/"bridge_receipt.json").write_text(json.dumps(info,sort_keys=True,indent=2)+"\n")
        temp.rename(output)
        return info

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source",required=True,type=Path)
    ap.add_argument("--output",required=True,type=Path)
    ap.add_argument("--expected-count",type=int)
    a=ap.parse_args()
    report=bridge(a.source,a.output,expected_count=a.expected_count)
    print("MOIS_ORIGINAL_93M_STAGE1_REAL_PIXELS_TO_STAGE2_SHADOW",json.dumps(report,sort_keys=True))

if __name__=="__main__":
    main()
