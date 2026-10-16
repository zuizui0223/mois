#!/usr/bin/env python3
"""Build sample-scoped, location-blind FCP ROI-v4 input packets from Mois sample.

This is NOT the FCP uncapped 93M-photo census firewall. Original FCP ROI-v4
source and weights are reused only as a measurement method; Mois has its own
complete, sampled and source-attested photo universe with fixed probabilities.
Never pass these subset packets to FCP's full-census source admission or
pretend they represent all 93,209,351 original photo opportunities.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sqlite3
from contextlib import ExitStack
from pathlib import Path

from scripts.acquisition.mois_stratified_photo_selection_20261016 import (
    PRIVATE_FIELDS, sha256_file, CAP, SOURCE_SPECIES, SOURCE_PHOTOS,
)
from scripts.acquisition.mois_sample_direct_from_fcp_source_20261016 import PHOTO_URL_PREFIX

PHOTO_SALT = "FCP_ALL_OBSERVABLE_PHOTO_BLIND_20261014_V1"
N_PARTITIONS = 128
WORKER = ["measurement_id", "image_filename", "photo_license"]
ACQUISITION = ["measurement_id", "image_filename", "photo_url_large", "photo_license"]
JOIN = ["measurement_id", *PRIVATE_FIELDS]
URL_RE = re.compile(r"^https://inaturalist-open-data\.s3\.amazonaws\.com/photos/([0-9]+)/large\.(jpg|jpeg|png)$", re.I)


def blind_id(photo_id: str) -> str:
    original = f"{PHOTO_SALT}\x1fphoto\x1f{photo_id}".encode()
    return "FCPR-" + hashlib.sha256(original).hexdigest().upper()[:24]


def worker_partition(measurement_id: str) -> int:
    return int(hashlib.sha256(measurement_id.encode()).hexdigest()[:16], 16) % N_PARTITIONS


def verify_sample(sample: Path, receipt_path: Path, *, synthetic: bool = False) -> dict:
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    if (receipt.get("schema") != "mois_direct_fcp_source_stratified_photo_selection_v1"
        or receipt.get("status") != "SOURCE_ADMITTED_SELECTED_PRIVATE_IDS_ONLY_NO_PHOTO_PIXELS_MEASURED"
        or receipt.get("selected_private_photo_csv_sha256") != sha256_file(sample)
        or receipt.get("historical_colour_outcomes_seen") is not False
        or receipt.get("pixels_opened_or_colours_measured") is not False
        or receipt.get("original_13var_environment_join_complete") is not False
        or receipt.get("selection_photo_cap_per_taxon") != CAP):
        raise ValueError("No genuine frozen Mois photo-ID sample to convert into ROI packets")
    if not synthetic and (receipt.get("source_species") != SOURCE_SPECIES
                          or receipt.get("source_photo_ids_authenticated") != SOURCE_PHOTOS):
        raise ValueError("Mois sample not bound to full original 42,111/93M population")
    return receipt


def assemble(sample: Path, receipt_path: Path, outdir: Path, *, synthetic: bool = False) -> dict:
    if outdir.exists() or outdir.with_name(outdir.name + ".incomplete").exists():
        raise ValueError("Refuse overwriting a sealed or unresolved private ROI handoff")
    receipt = verify_sample(sample, receipt_path, synthetic=synthetic)
    stage = outdir.with_name(outdir.name + ".incomplete")
    stage.mkdir(parents=True)
    ids_sqlite = stage / "_temporary_sample_identity_check.sqlite"
    total = pilot = 0
    by_species: dict[str, int] = {}
    part_counts = [0] * N_PARTITIONS
    try:
        with sqlite3.connect(ids_sqlite) as db, ExitStack() as stack:
            db.execute("CREATE TABLE ids (photo_id TEXT PRIMARY KEY, measurement_id TEXT UNIQUE, observation_key TEXT UNIQUE)")
            w_files, a_files = [], []
            for part in range(N_PARTITIONS):
                w_fh = stack.enter_context((stage / f"worker_partition_{part:03}.csv").open("w", newline="", encoding="utf-8"))
                a_fh = stack.enter_context((stage / f"acquisition_partition_{part:03}.csv").open("w", newline="", encoding="utf-8"))
                ww = csv.DictWriter(w_fh, fieldnames=WORKER); ww.writeheader()
                aa = csv.DictWriter(a_fh, fieldnames=ACQUISITION); aa.writeheader()
                w_files.append(ww); a_files.append(aa)
            metadata_handle = stack.enter_context((stage / "private_metadata_join.csv").open("w", newline="", encoding="utf-8"))
            meta_writer = csv.DictWriter(metadata_handle, fieldnames=JOIN)
            meta_writer.writeheader()
            with sample.open("r", newline="", encoding="utf-8") as source:
                reader = csv.DictReader(source)
                if reader.fieldnames != PRIVATE_FIELDS:
                    raise ValueError("Private Mois source sample schema drift")
                for row in reader:
                    tax, photo, obs = (row[k] for k in ("inat_taxon_id","photo_id","observation_id"))
                    if (not tax.isdecimal() or not photo.isdecimal() or not obs
                        or int(tax) < 1 or int(photo) < 1):
                        raise ValueError("Invalid original source taxon/photo/observation IDs")
                    url = row["photo_url_large"]
                    mm = URL_RE.fullmatch(url)
                    if mm is None or mm.group(1) != photo or not row["photo_license"].strip():
                        raise ValueError("Original provider photo URL/licence missing or changed")
                    n, k, m = (int(row[col]) for col in
                               ("available_observations_in_stratum",
                                "sampled_observations_in_stratum",
                                "source_photos_in_observation"))
                    if not (1 <= k <= n and m >= 1):
                        raise ValueError("Impossible sampling strata or photo multiplicity")
                    pi = k / (n * m)
                    if (not math.isclose(float(row["photo_inclusion_probability"]), pi, rel_tol=1e-9)
                        or not math.isclose(float(row["photo_design_weight"]), 1/pi, rel_tol=1e-9)
                        or not math.isclose(float(row["observation_inclusion_probability"]), k/n, rel_tol=1e-9)):
                        raise ValueError("Photo/observation probabilities do not match the sealed sample")
                    ident = blind_id(photo)
                    db.execute("INSERT INTO ids VALUES (?,?,?)", (photo, ident, tax + "\x1f" + obs))
                    by_species[tax] = by_species.get(tax, 0) + 1
                    if by_species[tax] > CAP:
                        raise ValueError("The selected Mois measurement cap is exceeded")
                    p = worker_partition(ident)
                    filename = ident + ".jpg"
                    basic = {"measurement_id":ident,"image_filename":filename,
                             "photo_license":row["photo_license"]}
                    w_files[p].writerow(basic)
                    a_files[p].writerow({**basic,"photo_url_large":url})
                    meta_writer.writerow({"measurement_id":ident,**row})
                    part_counts[p] += 1
                    pilot += int(row["technical_pilot_flag"] == "1")
                    total += 1
            db.commit()
        if (total != receipt["selected_original_photo_ids"]
            or pilot != receipt["technical_pilot_selected_photo_ids"]):
            raise ValueError("Selected photo and pilot totals changed during ROI packet handoff")
        ids_sqlite.unlink()
        file_hashes = {}
        for p in range(N_PARTITIONS):
            for stem in ("worker_partition", "acquisition_partition"):
                name = f"{stem}_{p:03}.csv"
                file_hashes[name] = sha256_file(stage / name)
        file_hashes["private_metadata_join.csv"] = sha256_file(stage / "private_metadata_join.csv")
        result = {
            "schema": "mois_sample_scoped_original_roi_v4_blind_packets_v1",
            "status": "SAMPLE_ONLY_BLIND_METADATA_SEALED_PIXELS_NOT_MEASURED",
            "original_photo_sample_receipt_sha256": sha256_file(receipt_path),
            "original_selected_private_photo_csv_sha256": sha256_file(sample),
            "selected_photo_ids": total,
            "selected_original_taxa_with_photo": len(by_species),
            "technical_pilot_photo_ids": pilot,
            "source_full_population_photo_ids": receipt["source_photo_ids_authenticated"],
            "worker_partitions": N_PARTITIONS,
            "partition_photo_counts": part_counts,
            "file_sha256": file_hashes,
            "worker_only_fields": WORKER,
            "metadata_site_and_climate_hidden_from_worker": True,
            "byte_identical_to_fcp_original_photo_id_blinding_algorithm": True,
            "only_mois_sample_not_all_fcp_source_photos": True,
            "original_pixels_measured": False,
            "historical_colour_outcomes_reused": False,
            "do_not_pass_to_fcp_full_census_handoff": True,
            "original_roi_v4_model_weights_verified_here": False,
            "note": "Frozen photo ID packets are ready for a SEPARATELY authorized subset-scoped ROI-v4 executor, not FCP's 93M-photo production worker.",
        }
        (stage / "mois_roi_packet_receipt.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        stage.replace(outdir)
        return result
    except Exception:
        # Preserve a failed staged output for explicit quarantine; never
        # return a fake successful packet receipt.
        raise


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--private-selection", type=Path, required=True)
    p.add_argument("--selection-receipt", type=Path, required=True)
    p.add_argument("--out-private-packets", type=Path, required=True)
    a = p.parse_args()
    result = assemble(a.private_selection, a.selection_receipt,
                      a.out_private_packets)
    print(json.dumps({k: result[k] for k in
                      ("status", "selected_photo_ids", "selected_original_taxa_with_photo",
                       "technical_pilot_photo_ids")}, sort_keys=True))


if __name__ == "__main__":
    main()
