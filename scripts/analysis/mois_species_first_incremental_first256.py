#!/usr/bin/env python3
"""Join real zero-ROI9 first-wave original-photo terminals into MOIS source atlas.

No fake labels and no repeat photo/observation requests. Use a blinded source
hash for never-tried species and first genuine independent fallback for failures.
A colour-bearing photograph is evidence of photographed hue, not a proven
fixed petal phenotype or GBIF distribution colour.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis.mois_species_first_palette9 import (
    PALETTE, COLUMNS, build_catalogue, class_labels, sha256,
)

SCHEMA = "mois_species_first_256_original_pixel_incremental_join_v1"
HUES = {"white", "yellow_orange", "red_pink", "blue_purple"}
BASE_PHOTOS = 1181146
BASE_STRICT = 451643
BASE_TAXA = 30048
FIRST_ATTEMPTS = 256
SOURCE_PHOTO_SHA = "3fbd0f0ec7eb347b3bbc7ab509ec9312ecd797ed932dcf7bb425959898c3835"
SOURCE_COVERAGE_SHA = "ea602238f3cd152fdd2e36e9ce4bcbc77b795a35dce802e9305f6f6cddd77fb8"
SOURCE_CENSUS_SHA = "6adab400c72932ab1079cdb66632018e90b14f20af2a53d6a8521dd8eb8441"
FIRST_SHARD_SHA = "339b6ed9b69e076b2cdb752b027411b1147fbf7bddf2db3a8619734dc9e02717"
FIRST_TERMINAL_SHA = "ff92dc3d3e94fa1fa42c3fe358fef7d9cdc8264d5b917dfd2f264d56044e8c34"

def require(cond, msg):
    if not bool(cond):
        raise ValueError("MOIS_INCREMENTAL_ORIGINAL_RGB9: " + msg)

def parse_boolean(x):
    v = x.astype(str).str.strip().str.lower()
    require(v.isin(["true", "false"]).all(), "Original strict9 Boolean status invalid")
    return v.eq("true")

def read_first_wave(source, result, receipt, *, sha_check=True):
    if sha_check:
        require(sha256(source) == FIRST_SHARD_SHA, "First real256 source shard changed")
        require(sha256(receipt) == FIRST_TERMINAL_SHA, "First actual ROI-v4 unit receipt changed")
    src = pd.read_csv(source, dtype=str, keep_default_na=False)
    measured = pd.read_csv(result, dtype=str, keep_default_na=False)
    require(len(src) == len(measured) == FIRST_ATTEMPTS and
            src.photo_id.is_unique and src.observation_id.is_unique and
            src.inat_taxon_id.is_unique and measured.photo_id.is_unique,
            "256 source taxa/photos/independent source observations lost or duplicated")
    require(set(src.photo_id) == set(measured.photo_id),
            "Source-measured original photo identities mismatch")
    report = json.loads(Path(receipt).read_text())
    require(report.get("attempted") == FIRST_ATTEMPTS and
            report.get("downloaded") >= 0 and
            report.get("photos_sha256") == sha256(result) and
            report.get("source_sha256") == sha256(source) and
            report.get("photo_colours_measured") is True,
            "Signed ROI actual terminal is not the frozen source shard")
    d = src[["photo_id", "inat_taxon_id", "observation_id"]].merge(
        measured, on="photo_id", suffixes=("_source", "_measured"),
        how="inner", validate="one_to_one")
    require((d.inat_taxon_id_source == d.inat_taxon_id_measured).all() and
            (d.observation_id_source == d.observation_id_measured).all(),
            "Measured taxon/observation conflicted with selected genuine original photo")
    numeric = d[["flower_fraction_" + x for x in PALETTE]].apply(
        pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    valid = (np.isfinite(numeric).all(axis=1) &
             (numeric >= 0).all(axis=1) &
             (numeric <= 1.00001).all(axis=1) &
             (np.abs(np.nan_to_num(numeric, nan=0).sum(axis=1) - 1) < 0.001))
    strict = parse_boolean(d.strict9_original_photo_numeric)
    expected = d.morph.isin(HUES).to_numpy(bool) & valid
    require(np.array_equal(strict.to_numpy(bool), expected),
            "Original numeric nine-palette success/unknown changed")
    frame = pd.DataFrame({"photo_id": d.photo_id,
                          "inat_taxon_id": d.inat_taxon_id_source,
                          "observation_id": d.observation_id_source,
                          "strict9": strict.to_numpy(bool)})
    for i, col in enumerate(COLUMNS):
        frame[col] = np.where(strict, numeric[:, i], np.nan)
    labels = np.full(len(frame), "", dtype=object)
    if strict.any():
        labels[strict] = class_labels(numeric[strict])
    frame["measured_palette_dominant_label"] = labels
    require(int(strict.sum()) == 91 and len(frame) == 256,
            "Expected actual signed first256 pixel outcomes drifted")
    return frame

def join_first_wave(baseline, census, coverage, first, *, expected_totals=True):
    required = {"inat_taxon_id", "photo_id", "observation_id", "strict9", *COLUMNS}
    require(required <= set(baseline.columns), "Missing baseline photo nine components")
    require(len(baseline) == BASE_PHOTOS and baseline.photo_id.is_unique and
            int(baseline.strict9.sum()) == BASE_STRICT,
            "Historical authentic original photo 1,181,146-source measurement changed")
    require(not set(baseline.photo_id.astype(str)) & set(first.photo_id.astype(str)),
            "New real256 source contains previously attempted image")
    oldobs = set(baseline.observation_id.fillna("").astype(str)) - {""}
    require(not oldobs.intersection(set(first.observation_id.astype(str))),
            "Source new photograph reuses old observation UUID")
    coverage = coverage.copy()
    coverage.inat_taxon_id = coverage.inat_taxon_id.astype(str)
    require(coverage.inat_taxon_id.is_unique and len(coverage) == 42111,
            "Source coverage 42111 nominal taxa drifted")
    source = coverage.set_index("inat_taxon_id").loc[first.inat_taxon_id]
    require(source.after_strict9.eq(0).all(),
            "A newly measured source taxon already had strict9 original image evidence")
    delta = first.groupby("inat_taxon_id", observed=True).agg(
        n_photo_attempts=("photo_id", "size"), n_strict9=("strict9", "sum"))
    coverage = coverage.set_index("inat_taxon_id")
    coverage.loc[delta.index, "after_photo_attempts"] += delta.n_photo_attempts
    coverage.loc[delta.index, "after_strict9"] += delta.n_strict9
    coverage = coverage.reset_index()
    combined = pd.concat([baseline[list(required)].copy(),
                          first[list(required)].copy()], ignore_index=True)
    catalogue = build_catalogue(combined, census, coverage)
    if expected_totals:
        require(int(combined.strict9.sum()) == BASE_STRICT + 91 and
                len(combined) == BASE_PHOTOS + FIRST_ATTEMPTS and
                int(catalogue.after_strict9.gt(0).sum()) == BASE_TAXA + 91,
                "Actual new ROI9 species-breadth gain did not reconcile")
    return catalogue, combined, coverage

def next_wave(catalogue, combined, first, candidates, *, quota=1024):
    """One photo each; retries unresolved first, then maximize NEW species."""
    require(1 <= quota <= 4096, "Unsupported measurement batch size")
    candidates = candidates.copy()
    for field in ("inat_taxon_id", "photo_id", "observation_id"):
        candidates[field] = candidates[field].fillna("").astype(str).str.strip()
    require(candidates.photo_id.is_unique and
            candidates.photo_id.str.fullmatch(r"\d+").all(),
            "Source fallback photo IDs invalid/duplicated")
    require(candidates[["inat_taxon_id", "observation_id"]].duplicated().sum() == 0,
            "Source fallback observation reused for same taxon")
    prior_photo = set(combined.photo_id.astype(str))
    prior_obs = set(combined.observation_id.fillna("").astype(str)) - {""}
    candidates = candidates.loc[
        (~candidates.photo_id.isin(prior_photo)) &
        (~candidates.observation_id.isin(prior_obs))].copy()
    catalogue = catalogue[["inat_taxon_id", "hue_status", "after_strict9"]].copy()
    candidates = candidates.merge(catalogue, on="inat_taxon_id",
                                  how="inner", validate="many_to_one")
    # User's breadth rule: after one CLEAR 9-palette category, temporarily
    # stop new photos; retain provisional_single status, not fixed species trait.
    candidates = candidates.loc[candidates.hue_status.isin(
        ["no_strict9", "mixed_image_unresolved"])].copy()
    first_taxa = set(first.inat_taxon_id)
    candidates["priority"] = np.where(candidates.inat_taxon_id.isin(first_taxa),
                                      "retry_unresolved", "new_zero9_species")
    candidates["fallback"] = pd.to_numeric(
        candidates.source_observation_rank, errors="raise").astype(int)
    candidates["stable_species_order"] = candidates.inat_taxon_id.map(
        lambda t: hashlib.sha256(("MOIS-WAVE2-FIRST|"+t).encode()).hexdigest())
    candidates = candidates.sort_values(
        ["priority", "stable_species_order", "fallback", "photo_id"],
        ascending=[False, True, True, True], kind="stable")
    # Explicit retry first, then never sampled taxa. Alphabetical priority
    # comparison is safe: 'retry...' > 'new...'.
    candidates = candidates.drop_duplicates("inat_taxon_id").head(quota)
    require(candidates.photo_id.is_unique and candidates.observation_id.is_unique and
            candidates.inat_taxon_id.is_unique, "New original photo/observation/taxon collision")
    candidates["measurement_id"] = "MOIS-ZERO9-W2-" + candidates.photo_id
    candidates["image_filename"] = candidates.measurement_id + ".jpg"
    candidates["shard"] = (np.arange(len(candidates)) // 256).astype(str)
    return candidates

def run(base_photos, source_taxa, coverage, first_source, first_pixels,
        first_receipt, global_options, output, *, quota=1024, sha_check=True):
    output=Path(output)
    require(not output.exists(), "Do not overwrite completed original data-only join")
    source_paths={"original_photo9":base_photos,"taxon_census":source_taxa,
                  "coverage":coverage}
    expected={"original_photo9":SOURCE_PHOTO_SHA, "taxon_census":SOURCE_CENSUS_SHA,
              "coverage":SOURCE_COVERAGE_SHA}
    if sha_check:
        for key, path in source_paths.items():
            require(sha256(path) == expected[key], "Frozen signed source SHA changed: "+key)
    base = pd.read_parquet(base_photos)
    census = pd.read_csv(source_taxa,dtype={"inat_taxon_id":str})
    old_cov = pd.read_csv(coverage,compression="infer",dtype={"inat_taxon_id":str})
    first = read_first_wave(first_source, first_pixels, first_receipt,
                            sha_check=sha_check)
    atlas, joined, counts = join_first_wave(base,census,old_cov,first)
    options = pd.read_csv(global_options, compression="infer", dtype=str)
    queue = next_wave(atlas,joined,first,options,quota=quota)
    output.mkdir(parents=True)
    atlas_file = output/"mois_42111_species_first_with_genuine_first256.csv.gz"
    delta_file = output/"mois_genuine_first256_strict9_append_only.parquet"
    queue_file = output/"mois_wave2_unresolved_and_new_species_photo_queue.csv"
    atlas.to_csv(atlas_file,index=False,compression="gzip",float_format="%.6f")
    first.to_parquet(delta_file,index=False)
    queue.to_csv(queue_file,index=False)
    gbif=atlas[["inat_taxon_id","species","hue_status","dominant_hue9",
                "observed_hue9_set"]].copy()
    gbif["gbif_accepted_species_key"]=""
    gbif["gbif_name_match_status"]="not_yet_matched"
    gbif.to_csv(output/"mois_updated_species_hue_for_gbif_match.csv",index=False)
    report={"schema":SCHEMA,"original_photo_attempts_after":int(len(joined)),
            "original_strict9_after":int(joined.strict9.sum()),
            "taxa_with_any_measured_strict9_after":int(atlas.after_strict9.gt(0).sum()),
            "new_real_first256_photo_attempts":FIRST_ATTEMPTS,
            "new_strict9_original_measurements":int(first.strict9.sum()),
            "new_dominant_palette9_labels":int(first.measured_palette_dominant_label.ne("").sum()),
            "new_roi9_without_confident_dominant_palette":int(
                (first.strict9 & first.measured_palette_dominant_label.eq("")).sum()),
            "new_strict9_technical_unknown":int((~first.strict9).sum()),
            "taxa_with_no_strict9_after":int(atlas.after_strict9.eq(0).sum()),
            "provisional_single_after":int(atlas.hue_status.eq("provisional_single").sum()),
            "next_wave_selected":int(len(queue)),
            "next_wave_retry":int(queue.priority.eq("retry_unresolved").sum()),
            "next_wave_new_taxa":int(queue.priority.eq("new_zero9_species").sum()),
            "no_repeated_original_photos_or_observation_UUIDs":True,
            "GBIF_occurrences_downloaded":0,
            "species_wide_fixed_colour_proven":False,
            "source_sha256":{k:sha256(v) for k,v in source_paths.items()},
            "new_first256_photo_source_sha256":sha256(first_source),
            "new_first256_terminal_sha256":sha256(first_pixels),
            "global_fallback_source_sha256":sha256(global_options)}
    report["outputs_sha256"]={p.name:sha256(p) for p in
        (atlas_file,delta_file,queue_file,output/"mois_updated_species_hue_for_gbif_match.csv")}
    (output/"receipt.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("MOIS_FIRST256_REAL_STRICT9_CATALOGUE_AND_WAVE2",
          json.dumps(report,sort_keys=True),flush=True)
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("base-photos","source-taxa","coverage","first-source",
                "first-pixels","first-receipt","global-options","output"):
        p.add_argument("--"+key,required=True)
    p.add_argument("--quota",type=int,default=1024)
    a=p.parse_args()
    run(a.base_photos,a.source_taxa,a.coverage,a.first_source,
        a.first_pixels,a.first_receipt,a.global_options,a.output,quota=a.quota)

if __name__=="__main__":
    main()
