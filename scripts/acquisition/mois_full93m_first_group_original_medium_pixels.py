#!/usr/bin/env python3
"""Scan eight signed MOIS source tier parts with ONE paced medium-image worker.

The first 136-photo part can be reused from SHA-signed previous actual
completion. No source photos are redownloaded for an already completed part.
No image stage is inferred from metadata; all other parts read REAL 500px
pixels, run SHA-pinned YOLO only, and emit immutable receipt per photo.
"""
from __future__ import annotations
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import tarfile

from scripts.acquisition.mois_full93m_medium_pixel_scan import (
    Pacer, WarmFlowerBoxDetector, download_medium, scan_source, sha_file
)

SCHEMA = "mois_full93m_one_true_source_group_eight_part_500px_yolo_v1"
SOURCE_TAG = "mois-complete-93209351-original-photo-target-tiers-v1"
PLAN_SHA = "6b607a4b7b0863d65a0eb514d4158434453e7ae4d909f4cf1e3372b424eb9279"
SOURCE_PARTS = 8

def require(ok, reason):
    if not ok:
        raise ValueError("MOIS_93M_GROUP_REAL_PIXEL_SCAN: " + reason)

def scan(tar, registry, source_plan, result_dir, *,
         weights, weights_sha, tier=0, group=0, rps=2.0,
         workers=8, use_previous_part0=True,
         required_plan_sha=PLAN_SHA):
    tar, registry, source_plan, result_dir = map(
        Path, (tar, registry, source_plan, result_dir))
    require(result_dir.is_dir() and rps > 0 and rps <= 3
            and 1 <= workers <= 12, "Bounded single-worker group only")
    require(tier == 0 and group == 0,
            "This first verified live run only uses Tier 0 / Group 0")
    original = json.loads(registry.read_text(encoding="utf-8"))
    require(original["schema"] == "mois_complete_full93m_original_photo_target_registry_v1"
            and original["exact_original_photo_IDs_in_all_disjoint_tiers"] == 93209351,
            "Invalid original 93m source receipt")
    require(tar.name == f"tier-{tier:02d}-group-{group:02d}.tar"
            and sha_file(tar) ==
              original["published_tier_archive_SHA256"][tar.name],
            "Original 93m source tier-group TAR SHA changed")
    require(sha_file(source_plan) == required_plan_sha,
            "The original full-93m / 1792-job photo plan changed")
    expected = {}
    with source_plan.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["tier"] == str(tier) and 0 <= int(row["part"]) < SOURCE_PARTS:
                part = int(row["part"])
                require(part not in expected, "Duplicate original plan photo part")
                expected[part] = row
    require(set(expected) == set(range(SOURCE_PARTS)),
            "Missing signed group source parts")
    folder = result_dir / "signed_source"
    folder.mkdir(exist_ok=True)
    with tarfile.open(tar) as archive:
        files = {m.name for m in archive.getmembers() if m.isfile()}
        wanted = {f"tier-{tier:02d}/part-{p:03d}.csv.gz"
                  for p in range(SOURCE_PARTS)}
        require(files == wanted,
                "Source group TAR includes unknown or missing photo files")
        for part in range(SOURCE_PARTS):
            name = f"tier-{tier:02d}/part-{part:03d}.csv.gz"
            member = archive.getmember(name)
            payload = archive.extractfile(member).read()
            out = folder / f"part-{part:03d}.csv.gz"
            if out.exists():
                require(out.read_bytes() == payload,
                        "Previously staged signed source bytes changed")
            else:
                out.write_bytes(payload)
            require(sha_file(out) == expected[part]["source_gzip_sha256"],
                    "Original part source SHA mismatch")

    det = WarmFlowerBoxDetector(
        weights, weights_sha, 0, "item", device="cpu", conf=0.02)
    pacer = Pacer(rps)
    download = lambda url: download_medium(url, pacer)
    counts = Counter()
    decomp = 0
    new_count = 0
    manifest = []
    for part in range(SOURCE_PARTS):
        old_part = part == 0 and use_previous_part0
        where = result_dir / f"part-{part:03d}"
        require(where.is_dir() if old_part else not where.exists(),
                "Expected immutable prior part or fresh output location missing")
        planned = expected[part]
        receipt = scan_source(
            folder / f"part-{part:03d}.csv.gz", where,
            expected_source_sha256=planned["source_gzip_sha256"],
            expected_count=int(planned["expected_photo_rows"]),
            tier=tier, part=part,
            detector=det, download=download,
            workers=workers, batch_size=512)
        require(receipt["all_original_input_photos_attempted"] is True
                and receipt["expected_photo_ids"] ==
                    int(planned["expected_photo_rows"])
                and receipt["actual_unique_source_photo_ids_attempted"] ==
                    receipt["expected_photo_ids"]
                and receipt["source_part"] == part
                and receipt["source_tier"] == tier
                and receipt["detector_identity"].startswith(
                    "yolo11:" + weights_sha),
                "Not an original photo complete source part")
        manifest.append({
            "source_part": part, "source_tier": tier,
            "actual_original_photo_attempts":
                receipt["actual_unique_source_photo_ids_attempted"],
            "actual_medium_image_pixels_decoded":
                receipt["medium_image_pixels_decoded"],
            "new_medium_image_source": not old_part,
            "source_complete_receipt_SHA256": sha_file(where / "complete_receipt.json"),
        })
        counts.update(receipt["screen_counts"])
        decomp += receipt["medium_image_pixels_decoded"]
        if not old_part:
            new_count += receipt["actual_unique_source_photo_ids_attempted"]
    total = sum(x["actual_original_photo_attempts"] for x in manifest)
    require(total == sum(int(x["expected_photo_rows"]) for x in expected.values())
            and total == sum(counts.values()) and decomp <= total,
            "Actual image outcomes do not reconcile to the eight original source photo parts")
    report = {
        "schema": SCHEMA,
        "original_full_photo_ID_catalogue": 93209351,
        "original_source_registry_sha256": sha_file(registry),
        "original_tier_8part_tar_sha256": sha_file(tar),
        "original_1792_job_source_manifest_sha256": sha_file(source_plan),
        "tier": tier, "group": group,
        "parts_total": len(manifest), "full_group_original_photo_IDs": total,
        "previously_completed_actual_photo_IDs_reused_without_reinference":
            int(expected[0]["expected_photo_rows"]) if use_previous_part0 else 0,
        "new_500px_photo_ID_attempts": new_count,
        "medium_image_pixels_decoded": decomp,
        "colour_hue_labels_promoted_to_botanical_petals": 0,
        "all_source_photo_ids_in_group_attempted": True,
        "incomplete_93m_pixel_source_census": True,
        "screen_counts": dict(counts),
        "parts": manifest,
    }
    summary = result_dir / "tier00_group00_complete_receipt.json"
    if summary.exists():
        require(json.loads(summary.read_text()) == report,
                "Previously published group receipt changed")
    else:
        summary.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for flag in ("tar","registry","source-plan","result-dir",
                 "weights","weights-sha"):
        parser.add_argument("--"+flag, required=True)
    parser.add_argument("--tier", type=int, default=0)
    parser.add_argument("--group", type=int, default=0)
    parser.add_argument("--requests-per-second", type=float, default=2.0)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--do-not-reuse-part0", action="store_true")
    a=parser.parse_args()
    r=scan(a.tar, a.registry, a.source_plan, a.result_dir,
           weights=a.weights, weights_sha=a.weights_sha,
           tier=a.tier, group=a.group, rps=a.requests_per_second,
           workers=a.workers, use_previous_part0=not a.do_not_reuse_part0)
    print("MOIS_TIER00_TRUE_996_PHOTO_PIXELS_FINISHED",
          json.dumps(r,sort_keys=True))

if __name__=="__main__":
    main()
