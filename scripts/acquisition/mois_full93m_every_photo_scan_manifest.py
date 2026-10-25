#!/usr/bin/env python3
"""Plan EVERY source-photo pixel scan job from MOIS's existing signed 93m tiers.

Source stage only: no image fetch, model call, floral colour classification,
native filter, species cap, or resampling. All 14 disjoint observation-hash
tiers x all 128 original source partitions are retained (1792 jobs).
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import tempfile

SCHEMA = "mois_full93m_every_original_photo_medium_pixel_scan_jobs_v1"
EXPECTED_PHOTOS = 93_209_351
EXPECTED_PARTS = 128
EXPECTED_TIERS = 14
REGISTRY_SCHEMA = "mois_complete_full93m_original_photo_target_registry_v1"
GROUP_SCHEMA = "mois_onepass_full93m_disjoint_photo_target_tiers_v1"
FIELDS = ("part", "tier", "relative_source_gzip_path", "source_gzip_sha256",
          "expected_photo_rows", "photo_level_sampling_probability")

def require(ok, why):
    if not ok:
        raise ValueError("MOIS_FULL93M_SCAN_PLAN: " + why)

def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4 << 20), b""):
            h.update(block)
    return h.hexdigest()

def build_jobs(registry_path, receipts, output_dir, *,
               expected_photos=EXPECTED_PHOTOS,
               expected_parts=EXPECTED_PARTS,
               expected_tiers=EXPECTED_TIERS,
               verify_files_root=None):
    """Reconcile complete photo universe once; never load 93m IDs in memory."""
    registry_path = Path(registry_path)
    output_dir = Path(output_dir)
    require(not output_dir.exists(), "Existing immutable manifest must not be overwritten")
    registry = json.loads(registry_path.read_text(encoding="utf-8"))
    require(registry.get("schema") == REGISTRY_SCHEMA, "Wrong full-source registry schema")
    require(registry.get("exact_original_photo_IDs_in_all_disjoint_tiers") == expected_photos,
            "Not the complete 93m source population")
    require(registry.get("source_only_NO_photo_pixels_or_hue_measured") is True,
            "The registry must describe source IDs, not labels")
    require(registry.get("native_vs_introduced_not_hard_filtered") is True,
            "Source was prefiltered by nativeness")
    require(len(registry.get("new_target_counts_by_tier", {})) == expected_tiers,
            "Incomplete cumulative source tier counts")
    receipts = sorted(map(Path, receipts))
    require(bool(receipts), "No source-authenticated group receipts provided")
    seen = {}
    for path in receipts:
        doc = json.loads(path.read_text(encoding="utf-8"))
        require(doc.get("schema") == GROUP_SCHEMA
                and doc.get("source_snapshot_SHA256") ==
                   registry.get("full_source_manifest_sha256"),
                "Source manifest mismatch: " + str(path))
        for part in doc.get("parts", []):
            p = int(part["part"])
            require(0 <= p < expected_parts and p not in seen,
                    "Overlapping or invalid source partitions")
            counts = part.get("counts_by_disjoint_tier", {})
            hashes = part.get("ranked_target_csv_sha256", {})
            require(set(counts) == set(map(str, range(expected_tiers)))
                    and set(hashes) == set(map(str, range(expected_tiers))),
                    "One or more disjoint photo tiers missing")
            require(all(type(counts[str(t)]) is int and counts[str(t)] >= 0
                        for t in range(expected_tiers)), "Invalid tier source counts")
            require(all(len(hashes[str(t)]) == 64
                        and all(ch in "0123456789abcdef" for ch in hashes[str(t)])
                        for t in range(expected_tiers)), "Invalid source sha256")
            require(sum(counts.values()) == part["original_photo_IDs_in_part"],
                    "Photo count does not reconcile inside source partition")
            seen[p] = part
    require(set(seen) == set(range(expected_parts)),
            "Not all 128 original source partitions are present")
    totals = {str(t): sum(seen[p]["counts_by_disjoint_tier"][str(t)]
                          for p in range(expected_parts))
              for t in range(expected_tiers)}
    require(totals == registry["new_target_counts_by_tier"]
            and sum(totals.values()) == expected_photos,
            "Tier/source total is not exactly the complete 93m photo population")

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir()
    try:
        manifest = output_dir / "all_original_photos_stage1_jobs.csv"
        with manifest.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
            writer.writeheader()
            for part in range(expected_parts):
                entry = seen[part]
                for tier in range(expected_tiers):
                    name = "tier-%02d/part-%03d.csv.gz" % (tier, part)
                    expected_hash = entry["ranked_target_csv_sha256"][str(tier)]
                    if verify_files_root is not None:
                        src = Path(verify_files_root) / name
                        require(src.is_file() and sha256(src) == expected_hash,
                                "Source compressed shard is missing or modified: " + name)
                    writer.writerow({
                        "part": part, "tier": tier,
                        "relative_source_gzip_path": name,
                        "source_gzip_sha256": expected_hash,
                        "expected_photo_rows": entry["counts_by_disjoint_tier"][str(tier)],
                        "photo_level_sampling_probability": "1",
                    })
        report = {
            "schema": SCHEMA,
            "signed_full_source_registry_sha256": sha256(registry_path),
            "signed_group_receipts_sha256": [
                {"path": str(p), "sha256": sha256(p)} for p in receipts
            ],
            "source_original_photo_ID_count": expected_photos,
            "source_original_partitions": expected_parts,
            "source_disjoint_tiers": expected_tiers,
            "all_disjoint_scan_jobs_including_empty": expected_parts * expected_tiers,
            "counts_by_tier": totals,
            "manifest_sha256": sha256(manifest),
            "all_photo_ids_in_scope_without_colour_or_native_preselection": True,
            "photo_level_scan_intended_inclusion_probability": 1.0,
            "actual_images_downloaded_or_scanned": 0,
            "all_photo_source_bytes_checked": verify_files_root is not None,
            "scan_complete": False,
        }
        with (output_dir / "plan_receipt.json").open("w", encoding="utf-8") as handle:
            json.dump(report, handle, ensure_ascii=False, sort_keys=True, indent=2)
            handle.write("\n")
        return report
    except BaseException:
        for file in output_dir.iterdir():
            file.unlink()
        output_dir.rmdir()
        raise

def main():
    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument("--registry", required=True, type=Path)
    a.add_argument("--group-receipt", required=True, nargs="+", type=Path)
    a.add_argument("--output", required=True, type=Path)
    a.add_argument("--verify-files-root", type=Path)
    x = a.parse_args()
    result = build_jobs(x.registry, x.group_receipt, x.output,
                        verify_files_root=x.verify_files_root)
    print("MOIS_FULL93M_SOURCE_ALL_PHOTO_STAGE1_PLAN",
          json.dumps(result, sort_keys=True))

if __name__ == "__main__":
    main()
