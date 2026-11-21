#!/usr/bin/env python3
"""Species-first nine-palette MOIS catalogue from frozen REAL original-photo ROI9.

No photograph download or fresh inference. Sampling consensus is an operational
stopping rule, NOT proof of fixed petal colour, wild nativeness or evolution.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

PALETTE = ("white", "yellow", "orange", "red", "pink", "magenta", "purple", "blue", "bronze")
COLUMNS = tuple("strict9_" + c for c in PALETTE)
SCHEMA = "mois_species_first_nine_palette_taxon_screen_v1"


def require(ok, reason):
    if not bool(ok):
        raise ValueError("MOIS_SPECIES_FIRST: " + reason)


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(4 << 20), b""):
            h.update(block)
    return h.hexdigest()


def class_labels(matrix, min_fraction=0.60, min_margin=0.20):
    """Return highest-photo-palette component or empty for a blended ROI."""
    a = np.asarray(matrix, dtype=float)
    require(a.ndim == 2 and a.shape[1] == 9, "Expected nine measured components")
    ix = a.argmax(axis=1)
    top = a[np.arange(len(ix)), ix]
    second = np.partition(a, -2, axis=1)[:, -2]
    return np.where((top >= min_fraction) & (top - second >= min_margin),
                    np.asarray(PALETTE)[ix], "")


def build_catalogue(photos, census, coverage, *, min_observations=3,
                    min_fraction=0.60, min_margin=0.20):
    """All taxa retained, including genuine unknowns and zero-photo species."""
    require(min_observations >= 2, "Require >=2 independent observations")
    required = {"inat_taxon_id", "photo_id", "observation_id", "strict9", *COLUMNS}
    require(required <= set(photos.columns), "Frozen photo ROI9 schema incomplete")
    photos = photos.copy()
    for col in ("inat_taxon_id", "photo_id", "observation_id"):
        photos[col] = photos[col].fillna("").astype(str).str.strip()
    require(photos.photo_id.is_unique and photos.photo_id.str.fullmatch(r"\d+").all(),
            "Duplicate or invalid authentic original photo IDs")
    require(str(photos.strict9.dtype) in ("bool", "boolean") and not photos.strict9.isna().any(),
            "strict9 must be measured Boolean with technical unknowns retained")

    names = census[["inat_taxon_id", "species", "full_93m_original_observation_count",
                    "full_93m_original_photo_count"]].copy()
    names.inat_taxon_id = names.inat_taxon_id.astype(str)
    require(names.inat_taxon_id.is_unique and names.species.notna().all(),
            "Nominal source taxon names or identities ambiguous")
    counts = coverage[["inat_taxon_id", "after_photo_attempts", "after_strict9"]].copy()
    counts.inat_taxon_id = counts.inat_taxon_id.astype(str)
    require(counts.inat_taxon_id.is_unique, "Original source taxon coverage repeated")
    catalogue = names.merge(counts, on="inat_taxon_id", how="inner", validate="one_to_one")
    require(len(catalogue) == len(names) == len(counts), "42111 census/colour register incomplete")

    yes = photos.loc[photos.strict9].copy()
    matrix = yes[list(COLUMNS)].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    require(np.isfinite(matrix).all() and (matrix >= 0).all()
            and (matrix <= 1.00001).all()
            and np.allclose(matrix.sum(axis=1), 1, atol=0.001),
            "Measured strict nine-palette numeric vectors invalid; cannot impute missing")
    observed = yes.observation_id.ne("")
    # One observation with many photographs does NOT count as multiple trials.
    yes["observation_key"] = np.where(observed, "obs:" + yes.observation_id,
                                       "photo:" + yes.photo_id)
    yes["independent_observation_known"] = observed
    by_obs = yes.groupby(["inat_taxon_id", "observation_key"], observed=True, sort=False)
    obs = by_obs[list(COLUMNS)].mean().reset_index()
    obs["known"] = obs.observation_key.str.startswith("obs:")
    obs["hue"] = class_labels(obs[list(COLUMNS)].to_numpy(float), min_fraction, min_margin)

    rows = []
    for taxon, frame in obs.groupby("inat_taxon_id", observed=True, sort=False):
        valid_obs = frame.loc[frame.known]
        voting = valid_obs.loc[valid_obs.hue.ne("")]
        votes = voting.hue.value_counts()
        if len(votes):
            top_hue = sorted(votes[votes == votes.max()].index)[0]
            # Distinct hue sightings are real measured evidence, not a fixed genetic trait.
            classes = sorted(votes.index.tolist())
        else:
            fallback = frame.loc[frame.hue.ne("")].hue.value_counts()
            top_hue = sorted(fallback[fallback == fallback.max()].index)[0] if len(fallback) else ""
            classes = sorted(fallback.index.tolist())
        if len(votes) >= 2 and int(votes.iloc[1]) >= 2:
            status = "sampled_multicolour"
        elif len(votes) == 1 and len(valid_obs) >= min_observations and len(voting) == len(valid_obs):
            status = "sampling_consensus_single"
        elif len(votes) >= 2:
            status = "multicolour_candidate"
        elif top_hue:
            status = "provisional_single"
        else:
            status = "mixed_image_unresolved"
        vec = frame[list(COLUMNS)].mean().to_numpy(float)
        rows.append({"inat_taxon_id": taxon, "hue_status": status,
                     "dominant_hue9": top_hue, "observed_hue9_set": "|".join(classes),
                     "independent_roi9_observations": int(frame.known.sum()),
                     "distinct_palette_votes": int(len(votes)),
                     "ambiguous_observation_rois": int((frame.hue == "").sum()),
                     **{"observed_mean_" + p: float(v) for p, v in zip(PALETTE, vec)}})
    stat = pd.DataFrame(rows, columns=["inat_taxon_id", "hue_status", "dominant_hue9",
                    "observed_hue9_set", "independent_roi9_observations",
                    "distinct_palette_votes", "ambiguous_observation_rois",
                    *("observed_mean_" + x for x in PALETTE)])
    catalogue = catalogue.merge(stat, on="inat_taxon_id", how="left", validate="one_to_one")
    mask = catalogue.hue_status.isna()
    catalogue.loc[mask, "hue_status"] = "no_strict9"
    for field in ("dominant_hue9", "observed_hue9_set"):
        catalogue[field] = catalogue[field].fillna("")
    for field in ("independent_roi9_observations", "distinct_palette_votes",
                  "ambiguous_observation_rois"):
        catalogue[field] = catalogue[field].fillna(0).astype(int)
    catalogue["stop_standard_measurement"] = catalogue.hue_status.isin(
        ["sampling_consensus_single", "sampled_multicolour"])
    priority = {"no_strict9": 0, "mixed_image_unresolved": 1,
                "provisional_single": 2, "multicolour_candidate": 3,
                "sampled_multicolour": 9, "sampling_consensus_single": 9}
    catalogue["next_measure_priority"] = catalogue.hue_status.map(priority).astype(int)
    require((catalogue.after_strict9.gt(0) == catalogue.hue_status.ne("no_strict9")).all(),
            "Frozen taxon-level strict9 counts disagree with measured photo evidence")
    require(int(catalogue.after_strict9.sum()) <= int(photos.strict9.sum()),
            "Nominal source strict9 count exceeds original image evidence")
    return catalogue.sort_values("inat_taxon_id", kind="stable").reset_index(drop=True)


def plan_next_photos(catalogue, photos, candidates, *, quota=4096):
    """Source-only candidate queue: max one NEW original observation/species/wave.

    All candidate ranks are decided WITHOUT looking at the candidate's colour.
    Worker is responsible for licence and exact-source SHA validation.
    """
    require(0 < quota <= 100000, "Invalid queue capacity")
    needed = {"photo_id", "observation_id", "inat_taxon_id"}
    require(needed <= set(candidates.columns), "Candidate photo source IDs missing")
    d = candidates.copy()
    for key in needed:
        d[key] = d[key].fillna("").astype(str).str.strip()
    require(d.photo_id.str.fullmatch(r"\d+").all() and d.photo_id.is_unique,
            "Candidate source has malformed/duplicate photo IDs")
    d = d.loc[~d.photo_id.isin(set(photos.photo_id.astype(str)))].copy()
    # Do not treat 2 photographs of one recorded observation as independent sampling.
    prior_obs = set(photos.observation_id.dropna().astype(str)) - {""}
    d = d.loc[d.observation_id.ne("") & ~d.observation_id.isin(prior_obs)].copy()
    d = d.merge(catalogue[["inat_taxon_id", "hue_status", "next_measure_priority",
                             "stop_standard_measurement"]],
                on="inat_taxon_id", how="inner", validate="many_to_one")
    d = d.loc[~d.stop_standard_measurement].copy()
    d["stable_rank"] = d.photo_id.map(lambda x: hashlib.sha256(
        ("MOIS_SPECIES_FIRST_20261111|" + x).encode()).hexdigest())
    d = d.sort_values(["next_measure_priority", "stable_rank"], kind="stable")
    d = d.drop_duplicates("inat_taxon_id", keep="first").head(quota)
    return d.drop(columns=["stable_rank", "stop_standard_measurement"])


def run(photo_path, census_path, coverage_path, output, candidates=None,
        *, min_observations=3, quota=4096):
    output = Path(output)
    require(not output.exists(), "Do not overwrite an existing frozen atlas output")
    photos = pd.read_parquet(photo_path, columns=["photo_id", "observation_id",
                    "inat_taxon_id", "strict9", *COLUMNS])
    census = pd.read_csv(census_path, dtype={"inat_taxon_id": str})
    coverage = pd.read_csv(coverage_path, compression="infer", dtype={"inat_taxon_id": str})
    atlas = build_catalogue(photos, census, coverage, min_observations=min_observations)
    output.mkdir(parents=True)
    atlas_path = output / "mois_species_first_42111_nine_palette_catalogue.csv.gz"
    atlas.to_csv(atlas_path, index=False, compression="gzip", float_format="%.6f")
    match = atlas[["inat_taxon_id", "species", "hue_status", "dominant_hue9",
                   "observed_hue9_set"]].copy()
    match["gbif_accepted_species_key"] = ""
    match["gbif_name_match_status"] = "not_yet_matched"
    match.to_csv(output / "mois_species_first_for_gbif_name_match.csv", index=False)
    if candidates:
        pending = pd.read_csv(candidates, dtype=str)
        batch = plan_next_photos(atlas, photos, pending, quota=quota)
        batch.to_csv(output / "mois_species_first_next_original_photo_queue.csv", index=False)
    counts = atlas.hue_status.value_counts().to_dict()
    report = {"schema": SCHEMA, "source_sha256": {
        "photos": sha256(photo_path), "source_taxa": sha256(census_path),
        "taxon_roi9_coverage": sha256(coverage_path)},
        "nominal_source_taxa": int(len(atlas)),
        "source_photo_attempts": int(len(photos)),
        "strict_nine_palette_original_photos": int(photos.strict9.sum()),
        "taxa_with_any_strict9": int(atlas.after_strict9.gt(0).sum()),
        "taxa_no_strict9": int(atlas.after_strict9.eq(0).sum()),
        "hue_status_counts": {k: int(v) for k, v in counts.items()},
        "taxa_operationally_retired": int(atlas.stop_standard_measurement.sum()),
        "GBIF_species_occurrences_downloaded": 0,
        "RGB_photo_petal_truth_verified": False,
        "constant_hue_across_full_species_range_proven": False,
        "new_original_photos_downloaded": 0,
        "new_nine_palette_inferences": 0,
        "files_sha256": {atlas_path.name: sha256(atlas_path)}}
    (output / "mois_species_first_receipt.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n")
    print("MOIS_SPECIES_FIRST_RESULT", json.dumps(report, sort_keys=True), flush=True)
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for option in ("photos", "source-taxa", "coverage", "output"):
        p.add_argument("--" + option, required=True)
    p.add_argument("--candidates", default=None)
    p.add_argument("--quota", type=int, default=4096)
    p.add_argument("--min-observations", type=int, default=3)
    x = p.parse_args()
    run(x.photos, x.source_taxa, x.coverage, x.output, x.candidates,
        min_observations=x.min_observations, quota=x.quota)


if __name__ == "__main__":
    main()
