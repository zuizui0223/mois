#!/usr/bin/env python3
"""Direct, outcome-blind Mois sampling from admitted FCP source shards.

Avoids a redundant 93-million-row climate-annotated intermediate CSV. Reads
and authenticates each original photo_manifest.ndjson exactly once and writes
only a <=40-observation-per-taxon *private* photo selection. Nothing downloads
image pixels, opens old hue outcomes, or estimates ecology.

Production requires completed FCP photo-ID admission, all 128 immutable source
shards, a byte-pinned WGS84 WorldClim BIO12 GeoTIFF and independent historical
BIO12 site calibration. Original site coordinates remain in restricted outputs.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from scripts.acquisition.mois_stratified_photo_selection_20261016 import (
    SNAPSHOT, SOURCE_SPECIES, SOURCE_PHOTOS, CAP,
    classify_site, read_verified_counts, sample_species, PRIVATE_FIELDS,
    sha256_file,
)

FRAME_SHA = "5d871bd1f1190c68fd513bb527e6c93de1606dc35ff8cfff5f6187850dc382cc"
MANIFEST_COLUMNS = {
    "inat_taxon_id", "species", "observation_id", "photo_id",
    "latitude", "longitude", "photo_license", "photo_url_large",
    "attribution", "observer_id", "observed_on",
}
PHOTO_URL_PREFIX = "https://inaturalist-open-data.s3.amazonaws.com/photos/"
HASH_RE = re.compile(r"^[0-9a-f]{64}$")


class PinnedBio12:
    """Read a relatively small, original-site WGS84 BIO12 raster once."""

    def __init__(self, tif: Path, expected_sha256: str):
        if not HASH_RE.fullmatch(expected_sha256):
            raise ValueError("Missing 64-character independent BIO12 raster SHA256 pin")
        if sha256_file(tif) != expected_sha256:
            raise ValueError("WorldClim BIO12 raster byte identity does not match the pin")
        try:
            import rasterio
        except ImportError as exc:
            raise RuntimeError("Install rasterio on the trusted compute host for the pinned GeoTIFF") from exc
        with rasterio.open(tif) as source:
            t = source.transform
            if (source.count != 1 or source.crs is None or source.crs.to_epsg() != 4326
                or abs(t.b) > 1e-10 or abs(t.d) > 1e-10 or t.a <= 0 or t.e >= 0
                or source.width * source.height > 12_000_000
                or source.bounds.left > -179.9 or source.bounds.right < 179.9):
                raise ValueError("Expected a small, north-up, global WGS84 BIO12 source raster")
            self.band = source.read(1, masked=True)
            self.transform = t
            self.width = source.width
            self.height = source.height
            self.nodata = source.nodata
            self.provenance = {
                "bio12_tif_sha256": expected_sha256,
                "bio12_crs": "EPSG:4326",
                "bio12_raster_width": source.width,
                "bio12_raster_height": source.height,
                "bio12_raster_nodata": self.nodata,
                "bio12_sample_method": "north-up raster cell containing original observation site",
            }

    def value(self, lat: str, lon: str) -> str:
        if not lat and not lon:
            return ""
        if not lat or not lon:
            raise ValueError("Half-empty source coordinates cannot be raster sampled")
        la, lo = float(lat), float(lon)
        if not math.isfinite(la) or not math.isfinite(lo) or la < -90 or la > 90 or lo < -180 or lo > 180:
            raise ValueError("Invalid original-site coordinate; no climate interpolation")
        col = math.floor((lo - self.transform.c) / self.transform.a)
        row = math.floor((la - self.transform.f) / self.transform.e)
        if row < 0 or col < 0 or row >= self.height or col >= self.width:
            return ""
        cell = self.band[row, col]
        if getattr(cell, "mask", False) is True or bool(getattr(cell, "mask", False)):
            return ""
        val = float(cell)
        if not math.isfinite(val) or val < 0 or val > 30000:
            return ""
        return f"{val:.12g}"


def calibrate(lookup, site_csv: Path, expected_sha256: str, *, min_sites: int = 10) -> dict:
    """Require original-photo-site historical BIO12 agreement BEFORE real sampling."""
    if not HASH_RE.fullmatch(expected_sha256) or sha256_file(site_csv) != expected_sha256:
        raise ValueError("Historical independent-site climate calibration SHA256 mismatch")
    with site_csv.open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames != ["latitude", "longitude", "bio12_mm"]:
            raise ValueError("Expected archived-site BIO12 calibration fields")
        n = 0
        regions = set()
        for row in reader:
            if not all(row.values()):
                raise ValueError("Calibration rows must be complete original-site observations")
            old = float(row["bio12_mm"])
            new = lookup.value(row["latitude"], row["longitude"])
            if not new or abs(float(new) - old) > 1.0:
                raise ValueError("BIO12 original-site extraction differs from archived FCP source")
            regions.add(classify_site(row)[0])
            n += 1
    if n < min_sites or len(regions) < (1 if min_sites < 10 else 3):
        raise ValueError("Too few geographically spread archived BIO12 calibration sites")
    return {"calibration_sites_matching": n,
            "calibration_geographic_groups": len(regions),
            "private_calibration_sha256": expected_sha256,
            "bio12_match_tolerance_mm": 1.0}


def index_shards(root: Path, n: int) -> dict[int, Path]:
    manifests: dict[int, Path] = {}
    for path in root.glob("shard_*_of_*/species/*/photo_manifest.ndjson"):
        if (not path.parent.name.isdecimal()
            or not re.fullmatch(r"shard_\d{3}_of_128", path.parents[2].name)):
            raise ValueError(f"Unexpected original FCP manifest layout: {path}")
        tax = int(path.parent.name)
        if tax < 1 or tax in manifests:
            raise ValueError("Duplicate or invalid taxon across FCP frozen source shards")
        manifests[tax] = path
    if len(manifests) != n:
        raise ValueError(f"Missing frozen per-taxon manifests: {len(manifests)} versus {n}")
    return manifests


def checked_source(root: Path, n_taxa: int, n_photos: int, synthetic: bool) -> None:
    bulk = json.loads((root / "bulk_source_receipt.json").read_text(encoding="utf-8"))
    if (bulk.get("schema") != "fcp_open_data_bulk_complete_source_v1"
        or bulk.get("status") != "COMPLETE_SNAPSHOT_METADATA_ALL_NOMINATED_TAXA"
        or bulk.get("source_snapshot_date") != SNAPSHOT
        or bulk.get("taxa_in_initial_frame") != n_taxa
        or bulk.get("n_target_taxon_licensed_photos") != n_photos
        or bulk.get("number_independently_sealed_source_shards") != 128
        or bulk.get("photo_cap", "INVALID") is not None
        or bulk.get("old_photo_measurements_read") is not False
        or bulk.get("image_pixels_opened") is not False):
        raise ValueError("FCP source shard package is not a completed original snapshot")
    if not synthetic and bulk.get("original_species_frame_sha256") != FRAME_SHA:
        raise ValueError("Original 42,111 taxon frame drift")


def original_rows(path: Path, tax: int, species: str, digest):
    """Authenticate bytes in the exact single pass used to select observations."""
    with path.open("rb") as source:
        for raw in source:
            digest.update(raw)
            if not raw.strip():
                raise ValueError("Unexpected empty original photo manifest row")
            try:
                obj = json.loads(raw)
            except Exception as exc:
                raise ValueError("Unreadable original FCP source photo manifest") from exc
            if (not isinstance(obj, dict) or not MANIFEST_COLUMNS.issubset(obj)
                or str(obj["inat_taxon_id"]) != str(tax)
                or str(obj["species"]) != species
                or str(obj.get("snapshot_date")) != SNAPSHOT):
                raise ValueError("Original source taxon/species/snapshot changed")
            photo = str(obj["photo_id"])
            url = str(obj["photo_url_large"])
            if (not photo.isdecimal() or not url.startswith(PHOTO_URL_PREFIX + photo + "/large.")
                or not str(obj["photo_license"]).strip()
                or not str(obj["observation_id"]).strip()):
                raise ValueError("Invalid original photo/licence/observation identity")
            lat, lon = obj["latitude"], obj["longitude"]
            yield {
                "inat_taxon_id": str(tax), "genus": species.split()[0],
                "photo_id": photo, "observation_id": str(obj["observation_id"]),
                "latitude": "" if lat is None else str(lat),
                "longitude": "" if lon is None else str(lon),
                "photo_license": str(obj["photo_license"]),
                "photo_url_large": url,
                "attribution": str(obj.get("attribution") or ""),
                "observer_id": str(obj.get("observer_id") or ""),
                "observed_on": str(obj.get("observed_on") or ""),
                "positional_accuracy_m": str(obj.get("positional_accuracy_m") or ""),
            }


def selected_rows_with_bio12(raw_rows, lookup):
    """No cross-species climate imputation; identical observation sites cached."""
    sites: dict[str, tuple[str, str, str]] = {}
    for row in raw_rows:
        oid = row["observation_id"]
        pair = (row["latitude"], row["longitude"])
        if oid not in sites:
            rain = lookup.value(*pair)
            sites[oid] = (pair[0], pair[1], rain)
        elif sites[oid][:2] != pair:
            raise ValueError("Inconsistent coordinates for the same original observation")
        row["bio12_mm"] = sites[oid][2]
        yield row


def run(root: Path, complete_counts: Path, admission: Path, bio12_lookup,
        output_csv: Path, output_receipt: Path, *,
        calibration: dict, synthetic: bool = False,
        species_sampler=None, selection_protocol: str | None = None) -> dict:
    if (species_sampler is None) != (selection_protocol is None):
        raise ValueError("Alternate sampling protocol must be named and have a dedicated sampler")
    if selection_protocol not in (None,
            "OUTCOME_BLIND_NESTED_LOCAL_SEASON_WITHIN_ORIGINAL_REGION_BIO12_V2"):
        raise ValueError("Unknown/unreviewed original-photo probability-sampling design")
    sampler = species_sampler or sample_species
    if output_csv.exists() or output_receipt.exists():
        raise ValueError("Existing photo selection or receipt cannot be overwritten")
    counts, count_sha = read_verified_counts(complete_counts, admission, synthetic=synthetic)
    manifests = index_shards(root, len(counts))
    checked_source(root, len(counts), sum(n for _, n in counts.values()), synthetic)
    # Byte-verified per-species manifest SHA values come from the *admitted*
    # full-source count ledger; totals alone cannot establish photo identity.
    manifests_sha = {}
    with complete_counts.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            tax = int(row["inat_taxon_id"])
            token = row["source_manifest_sha256"]
            if not HASH_RE.fullmatch(token):
                raise ValueError("Unsealed original photo manifest digest")
            manifests_sha[tax] = token
    temp = output_csv.with_name(output_csv.name + ".partial")
    if temp.exists():
        raise ValueError("Unresolved earlier interrupted private photo selection")
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    output_receipt.parent.mkdir(parents=True, exist_ok=True)
    photos = observations = selected = pilot = unknown_site = unknown_bio12 = 0
    sum_local_season_nested = sum_local_season_fallback = 0
    signature = hashlib.sha256()
    try:
        with temp.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=PRIVATE_FIELDS)
            writer.writeheader()
            for tax in sorted(counts):
                genus, expected = counts[tax]
                source_path = manifests[tax]
                checkpoint = json.loads((source_path.parent / "checkpoint.json").read_text(encoding="utf-8"))
                expected_sha = manifests_sha[tax]
                if (checkpoint.get("status") != "complete"
                    or checkpoint.get("taxon_id") != tax
                    or checkpoint.get("source_photo_rows") != expected
                    or checkpoint.get("metadata_sha256") != expected_sha
                    or checkpoint.get("photo_cap", "INVALID") is not None
                    or checkpoint.get("colour_pixels_opened") is not False):
                    raise ValueError(f"Taxon {tax} has unsealed FCP source checkpoint")
                digest = hashlib.sha256()
                # Preserve the exact source species name in the checkpoint,
                # rather than reconstructing the species from its genus.
                scientific_name = str(checkpoint.get("species", ""))
                if not scientific_name or scientific_name.split()[0] != genus:
                    raise ValueError("Source genus identity changed")
                if expected:
                    src = original_rows(source_path, tax, scientific_name, digest)
                    items = selected_rows_with_bio12(src, bio12_lookup)
                    n, o, k, p, bins = sampler(tax, genus, items, writer)
                else:
                    n = o = k = p = 0
                    bins = {"__unknown_site__": 0, "__unknown_climate__": 0}
                    # The 16 true zero-photo taxa still require verified
                    # zero-byte manifests, not omitted files.
                    with source_path.open("rb") as inp:
                        for raw in inp:
                            digest.update(raw)
                            raise ValueError("Zero-photo taxon has a nonempty manifest")
                if n != expected or digest.hexdigest() != expected_sha:
                    raise ValueError(f"Original manifest photo count/SHA mismatch for taxon {tax}")
                signature.update(f"{tax}:{expected_sha}\n".encode())
                photos += n
                observations += o
                selected += k
                pilot += p
                unknown_site += bins["__unknown_site__"]
                unknown_bio12 += bins["__unknown_climate__"]
                sum_local_season_nested += bins.get("__nested_observations__", 0)
                sum_local_season_fallback += bins.get("__fallback_observations__", 0)
        if photos != sum(n for _, n in counts.values()):
            raise ValueError("Whole-source photo opportunity denominator changed")
        selection_sha = sha256_file(temp)
        summary = {
            "schema": ("mois_direct_fcp_source_stratified_photo_selection_v2"
                       if selection_protocol else
                       "mois_direct_fcp_source_stratified_photo_selection_v1"),
            "sampling_protocol": (selection_protocol or "ORIGINAL_REGION_BIO12_V1"),
            "status": "SOURCE_ADMITTED_SELECTED_PRIVATE_IDS_ONLY_NO_PHOTO_PIXELS_MEASURED",
            "source_snapshot": SNAPSHOT,
            "source_species": len(counts),
            "source_species_with_photos": sum(n > 0 for _, n in counts.values()),
            "source_photo_ids_authenticated": photos,
            "source_observations": observations,
            "selected_original_photo_ids": selected,
            "technical_pilot_selected_photo_ids": pilot,
            "selection_photo_cap_per_taxon": CAP,
            "complete_source_counts_sha256": count_sha,
            "full_source_admission_sha256": sha256_file(admission),
            "per_taxon_manifest_digest_chain_sha256": signature.hexdigest(),
            "selected_private_photo_csv_sha256": selection_sha,
            "selected_missing_coordinates": unknown_site,
            "selected_missing_bio12": unknown_bio12,
            "selected_local_season_nested_observations": (
                sum_local_season_nested if selection_protocol else None),
            "selected_local_season_fallback_observations": (
                sum_local_season_fallback if selection_protocol else None),
            "site_bio12_provenance": getattr(bio12_lookup, "provenance", {}),
            "site_bio12_calibration": calibration,
            "historical_colour_outcomes_seen": False,
            "pixels_opened_or_colours_measured": False,
            "original_13var_environment_join_complete": False,
            "note": "Private selected image URLs/attribution and coordinates must not be committed to public GitHub.",
        }
        # Commit both files only after every original source taxon and manifest
        # has passed its SHA/checkpoint integrity tests.
        temp.replace(output_csv)
        output_receipt.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return summary
    except Exception:
        if temp.exists():
            temp.unlink()
        raise


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fcp-source-root", type=Path, required=True)
    ap.add_argument("--complete-counts", type=Path, required=True)
    ap.add_argument("--source-admission", type=Path, required=True)
    ap.add_argument("--bio12-geotiff", type=Path, required=True)
    ap.add_argument("--bio12-sha256", type=str, required=True)
    ap.add_argument("--archived-site-calibration", type=Path, required=True)
    ap.add_argument("--archived-site-calibration-sha256", type=str, required=True)
    ap.add_argument("--out-private-photo-csv", type=Path, required=True)
    ap.add_argument("--out-aggregate-receipt", type=Path, required=True)
    args = ap.parse_args()
    lookup = PinnedBio12(args.bio12_geotiff, args.bio12_sha256)
    calibration = calibrate(lookup, args.archived_site_calibration,
                            args.archived_site_calibration_sha256)
    out = run(args.fcp_source_root, args.complete_counts, args.source_admission,
              lookup, args.out_private_photo_csv, args.out_aggregate_receipt,
              calibration=calibration)
    print(json.dumps({key: out[key] for key in
                      ("status", "source_species", "source_photo_ids_authenticated",
                       "selected_original_photo_ids", "technical_pilot_selected_photo_ids")},
                     sort_keys=True))


if __name__ == "__main__":
    main()
