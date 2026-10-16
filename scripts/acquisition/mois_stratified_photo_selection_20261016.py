#!/usr/bin/env python3
"""Outcome-blind, taxon-wide, climate/geography-stratified Mois photo sampler.

Consumes *admitted*, globally deduplicated FCP photo IDs after the all-ID freeze.
The candidate CSV is sorted by taxon ID and annotated with nominal genus and
site BIO12 from the original coordinates (never with any flower-colour outcome).
The output is PRIVATE: it includes third-party photo and location identifiers.

One eligible photo is chosen uniformly (hash-priority) per observation, then
up to 40 independent observations per species are sampled across fixed strata.
The design preserves all nonempty strata, and reports observation/photo
selection probabilities for survey-weighted ecological analysis.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path

SNAPSHOT = "2026-09-27"
COUNT_SHA = "ef26c16252dd785f2d08bde3b18938c7624d73a54386fa3d3974dceb09ed97a9"
SOURCE_SPECIES = 42111
SOURCE_PHOTOS = 93209351
CAP = 40
SEED = "MOIS_20261016_GLOBAL_WHITE_13ENV_BLIND_V1"
REQUIRED = {"inat_taxon_id", "genus", "photo_id", "observation_id",
            "latitude", "longitude", "bio12_mm", "photo_license"}
COUNT_COLUMNS = ["inat_taxon_id", "species", "n_unique_target_photo_ids",
                 "source_status", "source_manifest_sha256"]
PRIVATE_FIELDS = ["inat_taxon_id", "genus", "observation_id", "photo_id",
                  "photo_license", "photo_url_large", "attribution",
                  "observer_id", "observed_on", "positional_accuracy_m",
                  "latitude", "longitude", "bio12_mm",
                  "region_stratum", "precip_stratum", "source_photos_in_observation",
                  "available_observations_in_stratum", "sampled_observations_in_stratum",
                  "observation_inclusion_probability", "photo_inclusion_probability",
                  "photo_design_weight", "technical_pilot_flag"]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def priority(kind: str, identifier: str) -> int:
    value = f"{SEED}\x1f{kind}\x1f{identifier}".encode("utf-8")
    return int.from_bytes(hashlib.sha256(value).digest()[:16], "big")


def classify_site(row: dict[str, str]) -> tuple[str, str]:
    lat, lon, rain = (row[key].strip() for key in ("latitude", "longitude", "bio12_mm"))
    if bool(lat) != bool(lon):
        raise ValueError("Half-missing coordinates: do not invent a geographic stratum")
    region = "G_UNKNOWN"
    if lat:
        la, lo = float(lat), float(lon)
        if not (math.isfinite(la) and math.isfinite(lo) and -90 <= la <= 90 and -180 <= lo <= 180):
            raise ValueError("Invalid photographed-site coordinates")
        belt = "S" if la < -23.5 else ("N" if la > 23.5 else "T")
        sector = min(2, int((lo + 180) // 120))
        region = f"G_{belt}{sector}"
    moisture = "P_UNKNOWN"
    if rain:
        mm = float(rain)
        if not math.isfinite(mm) or mm < 0 or mm > 30000:
            raise ValueError("Invalid annual precipitation BIO12")
        moisture = "P_DRY" if mm < 500 else ("P_MEDIUM" if mm < 1500 else "P_WET")
    return region, moisture


def read_verified_counts(counts_path: Path, admission_path: Path, *, synthetic: bool = False):
    counts_sha = sha256_file(counts_path)
    admission = json.loads(admission_path.read_text(encoding="utf-8"))
    if (admission.get("schema") != "fcp_complete_photo_ID_source_vs_independent_provider_audit_v1"
        or admission.get("status") != "REAL_ALL_42111_ORIGINAL_PHOTO_ID_COUNTS_RECONCILED_NO_PIXELS"
        or admission.get("snapshot_date") != SNAPSHOT
        or admission.get("fully_frozen_count_ledger_sha256") != counts_sha
        or admission.get("original_photo_IDs_covered_without_count_loss") is not True
        or admission.get("full_study_original_pixels_classified") is not False
        or admission.get("flower_colour_H2_W_or_environment_results") is not False):
        raise ValueError("A count-only scan is NOT an admitted full photo-ID census")
    if not synthetic and (admission.get("source_count_ledger_sha256") != COUNT_SHA
                          or admission.get("species") != SOURCE_SPECIES
                          or admission.get("distinct_original_photo_ids") != SOURCE_PHOTOS):
        raise ValueError("Original 42,111-taxon photo-ID population fingerprint changed")
    counts = {}
    with counts_path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames != COUNT_COLUMNS:
            raise ValueError("Unexpected complete-source count schema")
        for row in reader:
            tax, n = row["inat_taxon_id"], row["n_unique_target_photo_ids"]
            if (not tax.isdecimal() or not n.isdecimal() or int(tax) < 1
                or int(tax) in counts or row["source_status"] != "complete"
                or not row["species"].strip()):
                raise ValueError("Unfrozen, malformed, or duplicate source count taxon")
            counts[int(tax)] = (row["species"].strip().split()[0], int(n))
    if (len(counts) != admission.get("species")
        or sum(n for _, n in counts.values()) != admission.get("distinct_original_photo_ids")):
        raise ValueError("Completed count ledger disagrees with admitted full source")
    return counts, counts_sha


def quotas(n_by_stratum: dict[tuple[str, str], int]) -> dict[tuple[str, str], int]:
    """Outcome-independent sqrt(N) allocation with >=1 in every occupied cell."""
    if not n_by_stratum or len(n_by_stratum) > CAP:
        raise ValueError("Unsupported/nonpositive stratification frame")
    target = min(CAP, sum(n_by_stratum.values()))
    q = {s: 1 for s in n_by_stratum}
    for _ in range(target - len(q)):
        eligible = [s for s in q if q[s] < n_by_stratum[s]]
        if not eligible:
            raise AssertionError("Quota allocation capacity mismatch")
        chosen = min(eligible, key=lambda s: (q[s] / math.sqrt(n_by_stratum[s]), s))
        q[chosen] += 1
    assert sum(q.values()) == target
    return q


def sample_species(taxon: int, genus: str, rows, writer) -> tuple[int, int, int, int, dict[str, int]]:
    obs = {}
    seen_photo_ids = set()
    n_photos = 0
    for row in rows:
        n_photos += 1
        photo, oid = row["photo_id"].strip(), row["observation_id"].strip()
        if not photo or not oid or photo in seen_photo_ids:
            raise ValueError("Missing or duplicate original photo/observation identity")
        if not row["photo_license"].strip() or row["genus"].strip() != genus:
            raise ValueError("Original photo licence or species-to-genus mapping drift")
        seen_photo_ids.add(photo)
        group = classify_site(row)
        old = obs.get(oid)
        if old is not None:
            if old["group"] != group:
                raise ValueError("A single observation has inconsistent site/precipitation bins")
            old["count"] += 1
        p = priority("representative-photo", photo)
        if old is None or p < old["score"]:
            obs[oid] = {"group": group, "score": p, "count": 1 if old is None else old["count"],
                        "row": row}
    if not obs:
        raise ValueError("A nonzero original source taxon has no valid photographs")
    by_stratum = {}
    for oid, value in obs.items():
        by_stratum.setdefault(value["group"], []).append((priority("observation", f"{taxon}/{oid}"), oid, value))
    quota = quotas({s: len(v) for s, v in by_stratum.items()})
    selected = 0
    unknown_site = 0
    unknown_climate = 0
    bin_counts = {}
    # Technical pilot is only a preassigned processing stage; no colour feedback.
    pilot = int(priority("taxon-pilot", str(taxon)) % 10 == 0)
    for group in sorted(by_stratum):
        candidates = by_stratum[group]
        n, k = len(candidates), quota[group]
        for _, oid, representative in sorted(candidates, key=lambda x: (x[0], x[1]))[:k]:
            row = representative["row"]
            # One photo of m for a sampled observation; no substitutions after seeing colour.
            obs_pi = k / n
            photo_pi = obs_pi / representative["count"]
            writer.writerow({**{key: row[key] for key in ("inat_taxon_id", "genus", "photo_id",
                                                        "observation_id", "photo_license", "latitude",
                                                        "longitude", "bio12_mm")},
                             **{key: str(row.get(key) or "") for key in
                                ("photo_url_large", "attribution", "observer_id",
                                 "observed_on", "positional_accuracy_m")},
                             "region_stratum": group[0], "precip_stratum": group[1],
                             "source_photos_in_observation": representative["count"],
                             "available_observations_in_stratum": n,
                             "sampled_observations_in_stratum": k,
                             "observation_inclusion_probability": f"{obs_pi:.12g}",
                             "photo_inclusion_probability": f"{photo_pi:.12g}",
                             "photo_design_weight": f"{1/photo_pi:.12g}",
                             "technical_pilot_flag": pilot})
            selected += 1
            unknown_site += (group[0] == "G_UNKNOWN")
            unknown_climate += (group[1] == "P_UNKNOWN")
            bin_counts[f"{group[0]}:{group[1]}"] = bin_counts.get(f"{group[0]}:{group[1]}", 0) + 1
    return n_photos, len(obs), selected, pilot * selected, {**bin_counts,
                "__unknown_site__": unknown_site, "__unknown_climate__": unknown_climate}


def select(candidate_csv: Path, counts_csv: Path, admission_json: Path,
           out_csv: Path, out_receipt: Path, *, _synthetic: bool = False) -> dict:
    if out_csv.exists() or out_receipt.exists():
        raise ValueError("Refuse to overwrite a preexisting selection or source receipt")
    counts, count_sha = read_verified_counts(counts_csv, admission_json, synthetic=_synthetic)
    seen_taxa = set()
    summed_photos = summed_obs = summed_selected = summed_pilot = 0
    selected_bins = {}
    temp_csv = out_csv.with_name(out_csv.name + ".partial")
    if temp_csv.exists():
        raise ValueError("Unresolved interrupted photo selection .partial file")
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    out_receipt.parent.mkdir(parents=True, exist_ok=True)
    try:
        with candidate_csv.open(newline="", encoding="utf-8") as inf, temp_csv.open("w", newline="", encoding="utf-8") as outf:
            reader = csv.DictReader(inf)
            columns = set(reader.fieldnames or ())
            if not REQUIRED.issubset(columns):
                raise ValueError("Candidate frame lacks required original identities/sites/BIO12/licence")
            if any(any(t in c.lower() for t in ("colour", "color", "white", "hue", "pigment", "roi_", "flower_class"))
                   for c in columns):
                raise ValueError("Colour outcomes are forbidden in the outcome-blind sampling frame")
            writer = csv.DictWriter(outf, fieldnames=PRIVATE_FIELDS)
            writer.writeheader()
            last_taxon = 0
            for raw_taxon, group in itertools.groupby(reader, key=lambda row: row["inat_taxon_id"].strip()):
                if not raw_taxon.isdecimal() or int(raw_taxon) <= last_taxon:
                    raise ValueError("Candidate frame must be strictly sorted into taxon blocks")
                taxon = int(raw_taxon)
                if taxon not in counts or counts[taxon][1] == 0:
                    raise ValueError("Candidate taxon absent from admitted nonzero source counts")
                genus, expected = counts[taxon]
                n, o, k, pilot_k, bins = sample_species(taxon, genus, group, writer)
                if n != expected:
                    raise ValueError(f"Photo source loss or addition for taxon {taxon}: {n} vs {expected}")
                seen_taxa.add(taxon)
                last_taxon = taxon
                summed_photos += n
                summed_obs += o
                summed_selected += k
                summed_pilot += pilot_k
                for key, value in bins.items():
                    selected_bins[key] = selected_bins.get(key, 0) + value
        if seen_taxa != {tax for tax, (_, n) in counts.items() if n > 0}:
            raise ValueError("Not all originally photographed species are present in the candidate manifest")
        if summed_photos != sum(n for _, n in counts.values()):
            raise ValueError("Input photo denominator is not the admitted full 42,111-taxon source")
        receipt = {
            "schema": "mois_global_white_13env_outcome_blind_photo_sample_v1",
            "status": "PHOTO_IDS_SELECTED_NO_IMAGE_PIXELS_OR_COLOURS_MEASURED",
            "source_snapshot": SNAPSHOT,
            "original_full_source_admission_sha256": sha256_file(admission_json),
            "full_source_taxon_count_ledger_sha256": count_sha,
            "annotated_candidate_frame_sha256": sha256_file(candidate_csv),
            "private_selected_photo_csv_sha256": sha256_file(temp_csv),
            "source_taxa": len(counts), "source_taxa_with_photos": len(seen_taxa),
            "source_photo_opportunities_verified": summed_photos,
            "distinct_source_observations": summed_obs,
            "selected_original_photo_ids": summed_selected,
            "technical_pilot_photo_ids": summed_pilot,
            "max_selected_per_taxon": CAP,
            "pilot_taxon_probability": 0.1,
            "fixed_strata": "3 latitude belts x 3 longitude sectors + unknown; BIO12 <500, 500-1499, >=1500 mm + unknown",
            "allocation": "at_least_one_per_nonempty_stratum_then_sqrt_observation_count",
            "sampling_unit": "one_original_photo_per_distinct_observation",
            "observed_photo_outcomes_used_for_selection": False,
            "missing_classification_results_are_not_replaced": True,
            "private_selected_geo_missing": selected_bins.pop("__unknown_site__", 0),
            "private_selected_precip_missing": selected_bins.pop("__unknown_climate__", 0),
            "selected_stratum_counts": dict(sorted(selected_bins.items())),
            "ready_for_13_environment_association": False,
            "notes": "Image retrieval, ROI-v4 quality assessment, photo-site 13-variable joining and cohort-specific inference have NOT run.",
        }
        temp_csv.replace(out_csv)
        out_receipt.write_text(json.dumps(receipt, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")
        return receipt
    except Exception:
        if temp_csv.exists():
            temp_csv.unlink()
        raise


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--candidate-csv", type=Path, required=True)
    p.add_argument("--complete-counts", type=Path, required=True)
    p.add_argument("--admission", type=Path, required=True)
    p.add_argument("--out-csv", type=Path, required=True)
    p.add_argument("--out-receipt", type=Path, required=True)
    a = p.parse_args()
    result = select(a.candidate_csv, a.complete_counts, a.admission, a.out_csv, a.out_receipt)
    print(json.dumps({k: result[k] for k in ("status", "source_taxa", "source_photo_opportunities_verified",
                                                   "selected_original_photo_ids", "technical_pilot_photo_ids")}, sort_keys=True))


if __name__ == "__main__":
    main()
