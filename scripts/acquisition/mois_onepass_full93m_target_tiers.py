#!/usr/bin/env python3
"""Materialize ALL photo-ID targets once, by nested observation-probability tier.

Source input is the signed 93,209,351 original source PHOTO/OBSERVATION index.
Do NOT filter out introduced/unresolved taxa: regional WCVP is a separate,
non-deleting species x coarse region annotation; it does NOT identify an
individual wild/native image or a specifically native island population.
No image is fetched, no flower is certified, and no colour is inferred.
Each source archive group of 8 original taxon-modulo shards is processed once.
Subsequent acquisition waves can read ONLY the new disjoint tier, not 93m IDs.
"""
from __future__ import annotations

import argparse
import bisect
from collections import Counter
from contextlib import ExitStack
import csv
from datetime import date
import gzip
import io
import hashlib
import json
import math
from pathlib import Path

SOURCE_SCHEMA = "mois_full_42111_93209351_uncapped_original_photo_observation_index_v1"
RESULT_SCHEMA = "mois_onepass_full93m_disjoint_photo_target_tiers_v1"
EXPECTED_PHOTOS = 93_209_351
EXPECTED_TAXA = 42_111
EXPECTED_SHARDS = 128
SOURCE_SEED = "mois-global-observation-v1"
# These boundaries are cumulative probabilities. Each bin is a disjoint new wave.
TIER_ENDS = (.0002, .0004, .0008, .0016, .0032, .0064, .0128,
             .0256, .0512, .1024, .2048, .4096, .8192, 1.0)
SOURCE_PHOTO_COLUMNS = ("photo_id", "observation_uuid", "inat_taxon_id",
                        "photo_license", "image_extension", "photographer_id")
SOURCE_OBS_COLUMNS = ("observation_uuid", "inat_taxon_id", "observer_id",
                      "observed_on", "latitude", "longitude", "positional_accuracy")
COLUMNS = ("photo_id", "observation_uuid", "inat_taxon_id", "observer_id",
           "photographer_id", "observed_on", "latitude", "longitude",
           "positional_accuracy", "photo_license", "image_extension",
           "frozen_loader_supported", "valid_observed_date", "valid_geopoint",
           "observation_inclusion_score", "tier")
IMAGE_EXTENSIONS = frozenset({"jpg", "jpeg", "png"})


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(4 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def obs_draw(uuid):
    value = hashlib.sha256((SOURCE_SEED + "|" + str(uuid)).encode()).digest()[:8]
    return int.from_bytes(value, "big") / 2**64


def tier_of_score(score):
    require(0 <= score < 1, "An observation inclusion hash must be in [0,1)")
    return bisect.bisect_right(TIER_ENDS, score)


def usable_date(value):
    try:
        if len(value) != 10:
            return False
        return date.fromisoformat(value).isoformat() == value
    except (ValueError, TypeError):
        return False


def usable_point(lat, lon):
    try:
        x, y = float(lat), float(lon)
        return math.isfinite(x) and math.isfinite(y) and -90 <= x <= 90 and -180 <= y <= 180
    except (ValueError, TypeError):
        return False


def normalized_extension(value):
    value = str(value or "").strip().lower().lstrip(".")
    return value if value in IMAGE_EXTENSIONS else None


def csv_gz(path, fields):
    with gzip.open(path, mode="rt", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        require(tuple(reader.fieldnames or ()) == tuple(fields),
                f"Unexpected original source columns in {path}")
        yield from reader


def materialize_parts(index, output, parts, *, expected_total=EXPECTED_PHOTOS,
                      expected_taxa=EXPECTED_TAXA, expected_shards=EXPECTED_SHARDS):
    """Materialize specified original shard IDs (typically 8) and source receipt.

    Each source photo is visited exactly once. Observation hash is computed once
    per observation and reused by its photos. Exact 93m coverage is checked by
    aggregating disjoint group manifests, not by storing 93m Python ID strings.
    """
    index, output = Path(index), Path(output)
    require(not output.exists() or not any(output.iterdir()),
            "Refuse to overwrite any previously published pre-colour target roster")
    manifest_path = index / "manifest.json"
    source = json.loads(manifest_path.read_text(encoding="utf-8"))
    require(source.get("schema") == SOURCE_SCHEMA
            and source.get("photos_emitted") == expected_total
            and source.get("taxa_in_frame") == expected_taxa
            and len(source.get("shards", [])) == expected_shards
            and source.get("photo_caps") is None
            and source.get("original_photos_are_globally_unique_by_source_SQLite_primary_key") is True
            and source.get("photo_colours_measured") is False,
            "Must use the complete frozen unique-photo source, never a capped/photo-hue cohort")
    parts = tuple(sorted(set(int(x) for x in parts)))
    require(bool(parts) and all(0 <= n < expected_shards for n in parts),
            "Source partition index outside the immutable full source")
    root = output / "targets"
    group_counts = Counter()
    part_results = []
    for n in parts:
        entry = source["shards"][n]
        require(int(entry["part"]) == n,
                "Original source partition order changed")
        photo_path = index / "photos" / f"part-{n:03d}.csv.gz"
        obs_path = index / "observations" / f"part-{n:03d}.csv.gz"
        require(sha256(photo_path) == entry["photos_sha256"]
                and sha256(obs_path) == entry["observations_sha256"],
                "Original source PHOTO or OBSERVATION gzip SHA changed")
        # This map exists only for one original partition, not for 93m photos.
        observations = {}
        observed_rows = 0
        for o in csv_gz(obs_path, SOURCE_OBS_COLUMNS):
            uuid = o["observation_uuid"]
            require(uuid and uuid not in observations,
                    f"Duplicate or empty original observation UUID in partition {n}")
            score = obs_draw(uuid)
            observations[uuid] = (
                o["inat_taxon_id"], o["observer_id"], o["observed_on"],
                o["latitude"], o["longitude"], o["positional_accuracy"],
                score, tier_of_score(score),
                usable_date(o["observed_on"]),
                usable_point(o["latitude"], o["longitude"]),
            )
            observed_rows += 1
        require(observed_rows == entry["observation_count"],
                f"Original observation count mismatch in partition {n}")
        ids = set()
        counts = Counter()
        valid_date = valid_site = decoder_supported = 0
        output_paths = {}
        with ExitStack() as stack:
            writers = []
            for tier in range(len(TIER_ENDS)):
                folder = root / f"tier-{tier:02d}"
                folder.mkdir(parents=True, exist_ok=True)
                dest = folder / f"part-{n:03d}.csv.gz"
                raw = stack.enter_context(dest.open("wb"))
                zipped = stack.enter_context(gzip.GzipFile(
                    filename="", mode="wb", fileobj=raw, compresslevel=2, mtime=0))
                f = stack.enter_context(io.TextIOWrapper(
                    zipped, encoding="utf-8", newline=""))
                writer = csv.DictWriter(f, fieldnames=COLUMNS, lineterminator="\n")
                writer.writeheader()
                writers.append(writer)
                output_paths[tier] = dest
            for p in csv_gz(photo_path, SOURCE_PHOTO_COLUMNS):
                pid = p["photo_id"]
                require(pid.isdecimal() and pid not in ids,
                        f"Invalid or duplicate source photo ID in partition {n}")
                ids.add(pid)
                uuid = p["observation_uuid"]
                require(uuid in observations,
                        f"Source photo lacks its original observation in partition {n}")
                taxon, observer, observed_on, lat, lon, accuracy, score, tier, has_day, has_site = observations[uuid]
                require(p["inat_taxon_id"] == taxon,
                        "Photo and observation IDs point to different plant taxa")
                ext = normalized_extension(p["image_extension"])
                supported = ext is not None
                writers[tier].writerow({
                    "photo_id":pid,
                    "observation_uuid":uuid,
                    "inat_taxon_id":taxon,
                    "observer_id":observer,
                    "photographer_id":p["photographer_id"],
                    "observed_on":observed_on,
                    "latitude":lat, "longitude":lon,
                    "positional_accuracy":accuracy,
                    "photo_license":p["photo_license"],
                    # Keep unsupported extension in the source-only ledger.
                    "image_extension":ext if supported else p["image_extension"],
                    "frozen_loader_supported":"true" if supported else "false",
                    "valid_observed_date":"true" if has_day else "false",
                    "valid_geopoint":"true" if has_site else "false",
                    "observation_inclusion_score":format(score, ".17g"),
                    "tier":str(tier),
                })
                counts[tier] += 1
                valid_date += int(has_day)
                valid_site += int(has_site)
                decoder_supported += int(supported)
        require(len(ids) == entry["photo_count"] and sum(counts.values()) == entry["photo_count"],
                f"Part {n} photo-ID accounting failed after splitting full frame into tiers")
        # Never propagate 93m ID set outside this original partition.
        del ids
        del observations
        checksums = {str(t):sha256(file) for t,file in output_paths.items()}
        group_counts.update(counts)
        part_results.append({
            "part":n,
            "original_photo_IDs_in_part":int(sum(counts.values())),
            "original_observations_in_part":observed_rows,
            "counts_by_disjoint_tier":{str(k):int(counts[k]) for k in range(len(TIER_ENDS))},
            "source_photo_SHA256":entry["photos_sha256"],
            "source_observation_SHA256":entry["observations_sha256"],
            "ranked_target_csv_sha256":checksums,
            "valid_date_photo_candidates":valid_date,
            "valid_coordinate_photo_candidates":valid_site,
            "supported_frozen_loader_photo_candidates":decoder_supported,
        })
        print("MOIS_ONEPASS_FULL_SOURCE_TARGET_PART", json.dumps({
            "part":n, "photos":sum(counts.values()),
            "tier0":counts[0], "tier1":counts[1],
            "image_supported":decoder_supported
        },sort_keys=True),flush=True)
    output.mkdir(parents=True, exist_ok=True)
    report = {
        "schema":RESULT_SCHEMA,
        "source_snapshot_SHA256":sha256(manifest_path),
        "full_source_photo_IDs_expected":expected_total,
        "full_source_taxa_expected":expected_taxa,
        "source_shard_count":expected_shards,
        "processed_parts":list(parts),
        "original_photo_IDs_in_this_group":sum(group_counts.values()),
        "disjoint_tier_probability_upper_bounds":list(TIER_ENDS),
        "counts_by_disjoint_tier":{str(k):int(group_counts[k]) for k in range(len(TIER_ENDS))},
        "parts":part_results,
        "inclusion_hash_source_exact":SOURCE_SEED,
        "same_observation_uses_same_selection_score":True,
        "source_original_photo_IDs_uniqueness_proved_by_upstream_primary_key":True,
        "all_source_photo_IDs_retained_in_one_and_only_one_tier":True,
        "taxon_genus_observation_photo_caps":None,
        "images_downloaded":0,
        "flower_anthesis_certified":False,
        "colour_pixels_measured":0,
        "native_TDWG3_individual_or_island_status_inferred":False,
        "native_status_policy":"annotate_later_broad_TDWG3_never_hard_filter_photo_target",
        "status":"source_only_metadata_candidates_not_measured_hues",
    }
    (output / "receipt.json").write_text(json.dumps(
        report, ensure_ascii=False, sort_keys=True, indent=2) + "\n",encoding="utf-8")
    return report


def verify_all_groups(index_manifest, group_receipts, archives_dir=None):
    """Prove that 16 disjoint 8-part manifests account for all 93m source photos."""
    m = json.loads(Path(index_manifest).read_text())
    source_sha = sha256(index_manifest)
    expected = {int(x["part"]):int(x["photo_count"]) for x in m["shards"]}
    require(len(expected) == EXPECTED_SHARDS and m["photos_emitted"] == EXPECTED_PHOTOS,
            "A signed complete 93m source index is mandatory")
    seen = {}
    counts = Counter()
    archived = {}
    for path in group_receipts:
        g = json.loads(Path(path).read_text())
        require(g.get("schema") == RESULT_SCHEMA and g.get("source_snapshot_SHA256") == source_sha,
                "Target group receipt does not match the frozen 93m source")
        require(g.get("images_downloaded") == 0 and g.get("colour_pixels_measured") == 0,
                "Source target materialization cannot invent observed flower RGB")
        if archives_dir is not None:
            files = g.get("archive_sha256_by_tier")
            require(isinstance(files, dict) and len(files) == len(TIER_ENDS),
                    "Each source group must archive all 14 disjoint photo target tiers")
            for name, expected_hash in files.items():
                require(name not in archived and name.startswith("tier-")
                        and name.endswith(".tar"),
                        "Duplicate or incorrectly named source target tier archive")
                file = Path(archives_dir) / name
                require(file.is_file() and sha256(file) == expected_hash,
                        "Archived original 93m photo target bytes changed: " + name)
                archived[name] = expected_hash
        for part in g["parts"]:
            i = int(part["part"])
            require(i in expected and i not in seen, "Target shard included twice/missing from source")
            require(part["original_photo_IDs_in_part"] == expected[i],
                    "Original photo counts lost between index and disjoint target tiers")
            require(sum(part["counts_by_disjoint_tier"].values()) == expected[i],
                    "Original source photo missing from a one-time target tier")
            seen[i] = part
        counts.update({int(k):int(v) for k,v in g["counts_by_disjoint_tier"].items()})
    require(set(seen) == set(expected) and sum(counts.values()) == EXPECTED_PHOTOS,
            "Full 93m original photo IDs not covered by the complete tier registry")
    return {
        "schema":"mois_complete_full93m_original_photo_target_registry_v1",
        "exact_original_photo_IDs_in_all_disjoint_tiers":sum(counts.values()),
        "original_named_taxa_from_full_source":m["taxa_in_frame"],
        "full_source_manifest_sha256":source_sha,
        "disjoint_tier_probability_upper_bounds":list(TIER_ENDS),
        "new_target_counts_by_tier":{str(k):int(counts[k]) for k in range(len(TIER_ENDS))},
        "source_only_NO_photo_pixels_or_hue_measured":True,
        "native_vs_introduced_not_hard_filtered":True,
        "original_photo_pixels_or_images_redistributed":False,
        "published_tier_archive_SHA256":dict(sorted(archived.items())),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="cmd",required=True)
    a = sub.add_parser("group")
    a.add_argument("--index",type=Path,required=True)
    a.add_argument("--output",type=Path,required=True)
    a.add_argument("--parts",type=int,nargs="+",required=True)
    a = sub.add_parser("verify")
    a.add_argument("--index-manifest",type=Path,required=True)
    a.add_argument("--group-receipts",nargs="+",type=Path,required=True)
    a.add_argument("--output",type=Path,required=True)
    a.add_argument("--tier-archives",type=Path)
    args = p.parse_args()
    if args.cmd == "group":
        materialize_parts(args.index, args.output, args.parts)
    else:
        report=verify_all_groups(args.index_manifest,args.group_receipts,args.tier_archives)
        require(not args.output.exists(),"Do not overwrite already signed target manifest")
        args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
        print("MOIS_ALL_93209351_IMAGE_TARGET_SOURCE_IDS_REGISTERED_NOT_MEASURED",
              json.dumps({"target_photo_ids":report["exact_original_photo_IDs_in_all_disjoint_tiers"],
                          "tiers":report["new_target_counts_by_tier"]},sort_keys=True),flush=True)


if __name__ == "__main__":
    main()