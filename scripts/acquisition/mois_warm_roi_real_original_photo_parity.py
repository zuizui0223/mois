#!/usr/bin/env python3
"""Cold/warm ROI-v4 parity on genuinely NEW signed tier06 original source photos.

One deterministic, outcome-blind list of 32 source observations is reserved
from the already-published 93,209,351-photo one-pass registry. Original photos
are downloaded ONCE, kept in an ephemeral temp directory, then passed to the
unchanged frozen FCP main() via two cold subprocess packets and two packets
with one warm estimator. No raw images or synthetic/rescored floral hue labels
enter a global ecological result. Refuse source/identity/method drift.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
import time

import pandas as pd

from scripts.acquisition.mois_onepass_full93m_target_tiers import (
    COLUMNS, EXPECTED_PHOTOS, EXPECTED_TAXA, RESULT_SCHEMA,
    obs_draw, normalized_extension, sha256,
)
from scripts.acquisition.mois_warm_original_roi import WarmFrozenROI, FROZEN_SOURCE_BLOB, git_blob

TIER = 6
CANDIDATES = 32
SUCCESSFUL = 16
BENCHMARK_SEED = "mois-warm-real-technical-signed-tier06-v1"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def select_unseen_candidate_observations(registry, group_receipt, archive,
                                         previous_parquet, previous_receipt,
                                         temp_directory, *,
                                         expected_photo_total=EXPECTED_PHOTOS,
                                         expected_taxa=EXPECTED_TAXA,
                                         max_candidates=CANDIDATES):
    """Photo/observation identities and exact source bytes precede all pixels."""
    registry = Path(registry)
    group_receipt = Path(group_receipt)
    archive = Path(archive)
    previous_parquet = Path(previous_parquet)
    receipt = json.loads(Path(previous_receipt).read_text())
    universe = json.loads(registry.read_text())
    source = json.loads(group_receipt.read_text())
    require(universe.get("schema") == "mois_complete_full93m_original_photo_target_registry_v1"
            and universe.get("exact_original_photo_IDs_in_all_disjoint_tiers") == expected_photo_total
            and universe.get("original_named_taxa_from_full_source") == expected_taxa
            and universe.get("source_only_NO_photo_pixels_or_hue_measured") is True,
            "Not a complete authentic unmeasured original-photo registry")
    require(source.get("schema") == RESULT_SCHEMA
            and source.get("source_snapshot_SHA256") == universe["full_source_manifest_sha256"]
            and source.get("images_downloaded") == 0 and source.get("colour_pixels_measured") == 0,
            "Partial source-tier receipt is unverified")
    require(archive.name == "tier-06-group-00.tar"
            and sha256(archive) == universe["published_tier_archive_SHA256"][archive.name],
            "Target tier06 archive has changed")
    require(receipt.get("schema") == "mois_sealed_actual_pixel_checkpoint_v1"
            and receipt.get("combined_unique_photo_attempts", 0) >= 700000,
            "Existing prior source-reconciled actual-photo checkpoint missing")
    require(sha256(previous_parquet) ==
            receipt["output_sha256"]["combined_actual_photo_checkpoint.parquet"],
            "Prior genuine photo-ID table SHA mismatch")
    prior = pd.read_parquet(previous_parquet, columns=["photo_id"]).photo_id.astype(str)
    require(len(prior) == receipt["combined_unique_photo_attempts"]
            and prior.is_unique and prior.str.fullmatch(r"\d+").all(),
            "Previously attempted original photo IDs not uniquely certified")
    previous = set(prior)
    parts = {int(x["part"]):x for x in source["parts"]}
    require(set(parts) == set(range(8)), "Technical tier source must be signed group-00 parts 0-7")

    with tarfile.open(archive) as t:
        members = t.getmembers()
        names = {m.name for m in members if m.isfile()}
        expected = {f"tier-06/part-{i:03d}.csv.gz" for i in range(8)}
        require(names == expected and len(members) <= 16
                and all(m.isfile() or m.isdir() for m in members),
                "Unexpected raw source tar members")
        t.extractall(temp_directory, filter="data")

    options = []
    seen_ids = set()
    source_rows = 0
    skipped_prior = 0
    for part_number, part in sorted(parts.items()):
        path = Path(temp_directory) / f"tier-06/part-{part_number:03d}.csv.gz"
        require(sha256(path) == part["ranked_target_csv_sha256"][str(TIER)],
                "Frozen source partition SHA mismatch")
        count = 0
        with gzip.open(path, "rt", encoding="utf-8", newline="") as f:
            reader = csv.DictReader(f)
            require(tuple(reader.fieldnames or ()) == COLUMNS,
                    "Original source photo columns changed")
            for row in reader:
                count += 1
                pid, obs = row["photo_id"], row["observation_uuid"]
                require(pid.isdecimal() and obs and pid not in seen_ids,
                        "Source photo/observation identity missing or duplicated")
                seen_ids.add(pid)
                require(row["tier"] == str(TIER)
                        and float(row["observation_inclusion_score"]) == obs_draw(obs),
                        "Original observation selection hash changed")
                if pid in previous:
                    skipped_prior += 1
                    continue
                ext = normalized_extension(row["image_extension"])
                if ext is None:
                    continue
                rank = hashlib.sha256(
                    (BENCHMARK_SEED + "|" + obs + "|" + pid).encode()
                ).hexdigest()
                options.append((rank, row, ext))
        require(count == int(part["counts_by_disjoint_tier"][str(TIER)]),
                "Actual original photo IDs differ from signed tier partition")
        source_rows += count
    require(source_rows == int(source["counts_by_disjoint_tier"][str(TIER)]),
            "Signed group00 tier06 total is wrong")
    candidates = []
    obs_used = set()
    for _, row, ext in sorted(options, key=lambda x:x[0]):
        obs = row["observation_uuid"]
        if obs in obs_used:
            continue
        obs_used.add(obs)
        candidates.append({
            "photo_id":row["photo_id"],
            "observation_id":obs,
            "photo_license":row["photo_license"],
            "image_extension":ext,
            "photo_url_large":f"https://inaturalist-open-data.s3.amazonaws.com/photos/{row['photo_id']}/large.{ext}",
        })
        if len(candidates) == max_candidates:
            break
    require(len(candidates) == max_candidates, "Fewer than 32 unmeasured true source observations")
    provenance = {
        "schema":"mois_new_original_tier06_warm_benchmark_precolour_reservation_v1",
        "full_upstream_photo_IDs":expected_photo_total,
        "source_archive_sha256":sha256(archive),
        "source_registry_sha256":sha256(registry),
        "source_group_receipt_sha256":sha256(group_receipt),
        "prior_checkpoint_receipt_sha256":sha256(previous_receipt),
        "prior_checkpoint_parquet_sha256":sha256(previous_parquet),
        "signed_source_tier":TIER,
        "signed_tier_group_00_original_photos":source_rows,
        "previously_attempted_photo_IDs_excluded":skipped_prior,
        "observation_and_photo_disjoint_candidates":len(candidates),
        "precolour_candidate_photo_IDs":[x["photo_id"] for x in candidates],
        "all_selected_source_photo_ids_unseen_at_baseline":True,
        "selection_based_on_photo_hue":False,
        "technical_source_only_not_ecological_sampling":True,
    }
    return candidates, provenance


def benchmark(runtime, sam, registry, group_receipt, archive,
              previous_parquet, previous_receipt, output):
    from scripts.acquisition.run_fcp_fresh_image_partition import (
        image_loader, PermanentOriginalPhotoError,
    )
    runtime = Path(runtime)
    sam = Path(sam)
    output = Path(output)
    require(not output.exists(), "Refuse overwrite of source-authenticated technical benchmark")
    source_code = runtime / "scripts/analysis/measure_random_photo_first_blind_partition.py"
    detector = runtime / "data/atlas/qualification/roi_v4_training/jrc_yolo11n_last_v4.pt"
    require(git_blob(source_code) == FROZEN_SOURCE_BLOB and detector.is_file() and sam.is_dir(),
            "Frozen reference YOLO/EfficientSAM model not restored")
    with tempfile.TemporaryDirectory(prefix="mois_warm_new_original_pixels_") as tmp:
        root = Path(tmp)
        candidates, provenance = select_unseen_candidate_observations(
            registry, group_receipt, archive, previous_parquet, previous_receipt,
            root,
        )
        (root/"images").mkdir()
        acquired = []
        attempts = []
        for option in candidates:
            photo_id = option["photo_id"]
            try:
                # Exactly one paced original-photo network attempt, never once per trial.
                payload = image_loader(option["photo_url_large"])
            except PermanentOriginalPhotoError as exc:
                attempts.append({"photo_id":photo_id, "observation_id":option["observation_id"],
                                 "acquired":False, "status":"permanently_unavailable"})
                continue
            filename = "PHOTO_" + photo_id + ".jpg"
            (root/"images"/filename).write_bytes(payload)
            acquired.append({**option,
                             "filename":filename,
                             "photo_sha256":hashlib.sha256(payload).hexdigest(),
                             "measurement_id":"MOISWARM-"+hashlib.sha256(photo_id.encode()).hexdigest()[:24].upper()})
            attempts.append({"photo_id":photo_id, "observation_id":option["observation_id"],
                             "acquired":True, "status":"technical_reference_acquired"})
            if len(acquired) == SUCCESSFUL:
                break
        require(len(acquired) == SUCCESSFUL,
                "Insufficient retrieved source originals for fixed 16-photo paired validation")
        packets = []
        for part in range(2):
            path = root / f"blind_{part}.csv"
            with path.open("w",newline="") as f:
                writer = csv.DictWriter(f, fieldnames=["measurement_id","image_filename","photo_license"],
                                        lineterminator="\n")
                writer.writeheader()
                for p in acquired[part*8:(part+1)*8]:
                    writer.writerow({"measurement_id":p["measurement_id"],
                                     "image_filename":p["filename"],
                                     "photo_license":p["photo_license"]})
            packets.append(path)

        env = os.environ.copy()
        env["PYTHONPATH"] = os.pathsep.join(x for x in (
            str(runtime), str(Path(__file__).resolve().parents[2]),
            env.get("PYTHONPATH","")
        ) if x)
        cold, warmed, baseline_seconds, warm_seconds = [], [], 0., 0.
        for i, packet in enumerate(packets):
            dest = root/f"cold_{i}.csv"
            started = time.perf_counter()
            subprocess.run([
                sys.executable,str(source_code),
                "--worker-manifest",str(packet),
                "--images-dir",str(root/"images"),
                "--detector-weight",str(detector),
                "--efficient-sam-dir",str(sam),
                "--output-csv",str(dest),
            ],check=True,env=env)
            baseline_seconds += time.perf_counter() - started
            cold.append(pd.read_csv(dest,dtype={"measurement_id":str}))
        warm = WarmFrozenROI()
        for i, packet in enumerate(packets):
            dest = root/f"warm_{i}.csv"
            started = time.perf_counter()
            warmed.append(warm.evaluate(
                packet,root/"images",source_code,detector,sam,dest
            ))
            warm_seconds += time.perf_counter() - started
        mismatches = []
        for i,(x,y) in enumerate(zip(cold,warmed)):
            try:
                pd.testing.assert_frame_equal(x,y,check_exact=True,check_dtype=True)
            except AssertionError as exc:
                mismatches.append({"blind_packet":i,"description":str(exc)[:300]})
        require(warm.calls == 2 and warm.estimator_initializations == 1,
                "Original model sessions not reused once")
        # Every attempted new original ID remains a technical source reservation;
        # no model-labelled hue is copied to the ecological MOIS data lake.
        output.mkdir(parents=True)
        reserve_path=output/"technical_original_photo_id_reservations.csv"
        with reserve_path.open("w",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=["photo_id","observation_id","acquired","status"],
                                  lineterminator="\n")
            writer.writeheader();writer.writerows(attempts)
        report = {
            **provenance,
            "schema":"mois_real_original_tier06_frozen_cold_warm_parity_v1",
            "technical_reserved_attempted_photo_ids":len(attempts),
            "original_photos_downloaded_once":len(acquired),
            "original_image_sha256_by_photo_ID":{x["photo_id"]:x["photo_sha256"] for x in acquired},
            "technical_photo_id_reservations_SHA256":sha256(reserve_path),
            "frozen_FCP_evaluator_git_blob":FROZEN_SOURCE_BLOB,
            "old_cold_model_initializations":2,
            "warm_model_initializations":warm.estimator_initializations,
            "original_cold_seconds":baseline_seconds,
            "warm_seconds":warm_seconds,
            "observed_16_real_photo_speed_ratio":baseline_seconds/warm_seconds,
            "cold_warm_all_original_photo_terminal_fields_exactly_equal":not mismatches,
            "mismatch_details":mismatches,
            "new_photo_colour_labels_published_to_ecological_cohort":0,
            "future_tier06_workers_must_deduplicate_these_technical_photo_ids":True,
            "flowering_stage_independently_validated":False,
            "small_16_original_technical_photo_sample_not_93m_throughput":True,
        }
        (output/"real_original_cold_warm_parity_receipt.json").write_text(
            json.dumps(report,indent=2,sort_keys=True)+"\n")
        print("MOIS_REAL_SOURCE_TIER06_FROZEN_WARM_PARITY",
              json.dumps({k:report[k] for k in (
                  "technical_reserved_attempted_photo_ids","original_photos_downloaded_once",
                  "original_cold_seconds","warm_seconds","observed_16_real_photo_speed_ratio",
                  "cold_warm_all_original_photo_terminal_fields_exactly_equal"
              )},sort_keys=True),flush=True)
        require(not mismatches,"Cold/warm original photo outcomes changed: disable warm option")
        return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("runtime","sam","registry","group-receipt","archive",
                "previous-parquet","previous-receipt","output"):
        p.add_argument("--"+key,type=Path,required=True)
    x=p.parse_args()
    benchmark(x.runtime,x.sam,x.registry,x.group_receipt,x.archive,
              x.previous_parquet,x.previous_receipt,x.output)


if __name__=="__main__":
    main()
