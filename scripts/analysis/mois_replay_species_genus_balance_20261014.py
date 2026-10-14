#!/usr/bin/env python3
"""Recalculate photo/species/genus-balanced photographed hue from original FCP ZIP.

The raw source contains copyrighted/attributed photo and location data. This
script reads the raw source locally and publishes aggregate-only results. It
cannot estimate true species frequencies, inherited morphs or fitness payoffs.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import zipfile

import numpy as np
import pandas as pd

SOURCE_SHA256 = "2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0"
BRE = "fcp42111_real_site_abiotic/breadth_flower_colour_climate_soil.csv.gz"
CEL = "fcp42111_real_site_abiotic/taxon_cell_flower_colour_climate_soil.csv.gz"
FIELDS = ["inat_taxon_id", "photo_id", "species", "morph", "measurement_status"]
CLASSES = {"white", "yellow_orange", "red_pink", "blue_purple"}
SEED = 2026101450
N_BOOT = 1999


def sha256(path: Path) -> str:
    d = hashlib.sha256()
    with path.open("rb") as stream:
        for part in iter(lambda: stream.read(1024 * 1024), b""):
            d.update(part)
    return d.hexdigest()


def load_verified_photos(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    if sha256(path) != SOURCE_SHA256:
        raise ValueError("Incorrect original archive SHA256; source scope not reproducible")
    with zipfile.ZipFile(path) as archive:
        if archive.testzip() is not None:
            raise ValueError("Corrupted original source ZIP")
        def read(member):
            with archive.open(member) as compressed, gzip.GzipFile(fileobj=compressed) as uncompressed:
                return pd.read_csv(uncompressed, usecols=FIELDS, low_memory=False)
        breadth = read(BRE)
        taxon_cell = read(CEL)
    if (len(breadth), len(taxon_cell)) != (42111, 85337):
        raise ValueError("Original photo source row count drift")
    joined = pd.concat((breadth, taxon_cell), ignore_index=True)
    if int(joined.photo_id.duplicated().sum()) != 26905:
        raise ValueError("Original photo overlap drift")
    if joined.groupby("photo_id").inat_taxon_id.nunique().max() != 1:
        raise ValueError("Original photo re-used for distinct taxa")
    joined = joined.drop_duplicates("photo_id").copy()
    if len(joined) != 100543 or joined.inat_taxon_id.nunique() != 42111:
        raise ValueError("Original unique-photo and taxon counts drift")
    if joined.groupby("inat_taxon_id").species.nunique().max() != 1:
        raise ValueError("Original single taxon ID with multiple names")
    joined["genus"] = joined.species.str.split().str[0]
    joined["classified"] = joined.measurement_status.eq("classified_four_state_morph") & joined.morph.isin(CLASSES)
    if joined.classified.sum() != 45953 or joined.genus.isna().any():
        raise ValueError("Original four-hue classifications or genus derivation drift")
    return joined, breadth


def cluster_bootstrap(genus: pd.DataFrame, indices: np.ndarray) -> np.ndarray:
    """Resample original nominal genera (not photo rows) as statistical clusters."""
    values = genus[["white_photos", "n_photos", "sum_species_white", "n_species", "mean_species_white"]].to_numpy(float)
    sums = values[indices].sum(axis=1)
    return np.column_stack((sums[:, 0] / sums[:, 1], sums[:, 2] / sums[:, 3], sums[:, 4] / indices.shape[1]))


def summarized_weights(photo: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    classified = photo.loc[photo.classified].copy()
    classified["white"] = classified.morph.eq("white").astype(np.int8)
    species = (classified.groupby(["inat_taxon_id", "genus"], sort=False)
               .agg(white=("white", "mean"), n_photos=("white", "size"), white_count=("white", "sum"))
               .reset_index())
    genus = (species.groupby("genus", sort=False)
             .agg(n_species=("white", "size"), sum_species_white=("white", "sum"),
                  mean_species_white=("white", "mean"), n_photos=("n_photos", "sum"),
                  white_photos=("white_count", "sum"))
             .reset_index())
    return species, genus


def calculate(photo: pd.DataFrame, breadth: pd.DataFrame) -> dict:
    species, genus = summarized_weights(photo)
    classifiable = photo.loc[photo.classified]
    if (len(species), len(genus), int(classifiable.morph.eq("white").sum())) != (25064, 5369, 21574):
        raise ValueError("Original unique classified taxon/genus/hue denominator drift")
    w = [classifiable.morph.eq("white").mean(), species.white.mean(), genus.mean_species_white.mean()]
    rng = np.random.default_rng(SEED)
    indices = rng.integers(0, len(genus), size=(N_BOOT, len(genus)), dtype=np.int32)
    boot = cluster_bootstrap(genus, indices)
    groups = []
    for low, high in [(1, 1), (2, 2), (3, 4), (5, 9), (10, 19), (20, 999999)]:
        grp = genus.loc[genus.n_species.between(low, high)]
        groups.append({"sampled_species_range": [low, high], "n_genera": len(grp),
                       "n_species": int(grp.n_species.sum()),
                       "genus_equal_white": float(grp.mean_species_white.mean()),
                       "species_equal_white": float(grp.sum_species_white.sum()/grp.n_species.sum())})
    thresholds = []
    min3_boot = None
    for minimum in (1, 2, 3, 5, 10, 20):
        grp = genus.loc[genus.n_species >= minimum]
        if minimum == 3:
            min3_indices = rng.integers(0, len(grp), size=(N_BOOT, len(grp)), dtype=np.int32)
            min3_boot = np.quantile(cluster_bootstrap(grp, min3_indices)[:, 2], [.025, .975]).tolist()
        thresholds.append({"min_classified_species_per_genus": minimum,
                           "n_genera": len(grp), "n_species": int(grp.n_species.sum()),
                           "genus_equal_white": float(grp.mean_species_white.mean()),
                           "species_equal_white": float(grp.sum_species_white.sum()/grp.n_species.sum())})
    breadth_classified = breadth.loc[breadth.measurement_status.eq("classified_four_state_morph") & breadth.morph.isin(CLASSES)].copy()
    breadth_classified["genus"] = breadth_classified.species.str.split().str[0]
    breadth_classified["white"] = breadth_classified.morph.eq("white").astype(np.int8)
    n_missing_species = 42111 - len(species)
    total_known_white_expectation = species.white.sum()
    return {
        "schema": "mois_original_photo_species_genus_balance_v1",
        "source_photo_archive_sha256": SOURCE_SHA256,
        "sample_scope": "One original photo-ID deduplicated breadth and taxon-cell union; four-class visible hue only; no independent flower-colour reassessment",
        "n_unique_photo_opportunities": len(photo),
        "n_classified_photos": len(classifiable),
        "n_classified_white_photos": int(classifiable.morph.eq("white").sum()),
        "n_classified_coloured_photos": int(classifiable.morph.ne("white").sum()),
        "n_nominal_species_total": int(photo.inat_taxon_id.nunique()),
        "n_species_with_classified_photo": len(species),
        "n_species_with_no_classified_photo": n_missing_species,
        "n_named_genera_total": int(photo.genus.nunique()),
        "n_genera_with_classified_photo": len(genus),
        "n_genera_one_classified_species": int(genus.n_species.eq(1).sum()),
        "point_white_proportion": dict(zip(("photo_equal", "species_equal", "genus_equal"), map(float, w))),
        "paired_whole_genus_bootstrap_95": {
            k: np.quantile(boot[:, i], [.025, .975]).tolist()
            for i, k in enumerate(("photo_equal", "species_equal", "genus_equal"))},
        "genus_equal_minus_species_equal_white_proportion": float(w[2] - w[1]),
        "paired_genus_minus_species_bootstrap_95": np.quantile(boot[:, 2] - boot[:, 1], [.025, .975]).tolist(),
        "bootstrap_seed": SEED,
        "bootstrap_draws": N_BOOT,
        "grouped_by_sampled_species_per_genus": groups,
        "min_classified_species_thresholds": thresholds,
        "min3_genus_equal_bootstrap_95": min3_boot,
        "breadth_one_photo_per_nominal_taxon": {
            "n_classified": len(breadth_classified),
            "n_genera": int(breadth_classified.genus.nunique()),
            "white_species_equal": float(breadth_classified.white.mean()),
            "white_genus_equal": float(breadth_classified.groupby("genus").white.mean().mean()),
        },
        "n_species_one_classified_photo": int(species.n_photos.eq(1).sum()),
        "n_species_multiclass_labels": int(classifiable.groupby("inat_taxon_id").morph.nunique().gt(1).sum()),
        "hypothetical_all_42111_taxon_white_mean_bounds": [
            float(total_known_white_expectation / 42111),
            float((total_known_white_expectation+n_missing_species)/42111)],
        "hypothetical_missing_species_white_mean_for_full_half": float((42111*.5-total_known_white_expectation)/n_missing_species),
        "interpretation": "Photo, species and genus weights refer to classifiable opportunistic photo hues, NOT a census of all world plant species, stable within-population phenotypes, cost-benefit fitness trade-offs or a frequency-dependent selection equilibrium.",
        "inference_limit": "Genus bootstrap resamples the observed genera as sampling clusters and does not correct missing photo hues, observation selection, genus/species phylogeny or nonrandom species coverage. No source positions, photo IDs or named genera are included in this output."
    }


def write_or_check(report: dict, output: Path, verify_only: bool) -> None:
    if verify_only:
        current = json.loads(output.read_text(encoding="utf-8"))
        def compare(a, b, path="root"):
            if isinstance(a, dict):
                assert isinstance(b, dict) and set(a) == set(b), path
                for k, v in a.items(): compare(v, b[k], path+"."+k)
            elif isinstance(a, list):
                assert isinstance(b, list) and len(a)==len(b), path
                for i, v in enumerate(a): compare(v, b[i], path+"."+str(i))
            elif isinstance(a, float):
                assert isinstance(b, (int, float)) and abs(a-b)<1e-10, (path,a,b)
            else:
                assert a==b,(path,a,b)
        compare(report, current)
        print("PASS original source species/genus balance receipt and bootstrap")
    else:
        output.write_text(json.dumps(report, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
        print("WROTE", output)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--photo-zip", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    command = parser.add_mutually_exclusive_group(required=True)
    command.add_argument("--write", action="store_true")
    command.add_argument("--check", action="store_true")
    args = parser.parse_args()
    photo, breadth = load_verified_photos(args.photo_zip)
    report = calculate(photo, breadth)
    write_or_check(report, args.output, args.check)

if __name__ == "__main__":
    main()
