#!/usr/bin/env python3
"""Observation-level probability audit of existing MOIS original image attempts.

The target is one prespecified, outcome-blind canonical photo per *represented
observation* in the frozen measured-photo frame. No claim about the 93m upstream
plant population, wild/native status, or missing-photo flower colour is made.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path

import pandas as pd

from scripts.analysis.mois_1181146_original_roi9_missingness_gap_decision import (
    latitude_band, recorded_local_season,
)

SCHEMA = "mois_original1181146_observation_two_phase_botanical_audit_v1"
SEED = "MOIS_20261029_PREOUTCOME_OBSERVATION_FRAME_BIAS_AUDIT_V1"
ARMS = (
    "historical_K1_FCP_source_heterogeneous",
    "bias_refill_2048", "p2_first_4096", "bulk128_131072",
)
EXPECTED_PHOTOS = 1181146
EXPECTED_SUCCESS = 451643
DEFAULT_SUCCESS = 400
DEFAULT_UNKNOWN = 800
COLS = (
    "photo_id", "observation_id", "inat_taxon_id", "source_design",
    "observed_on", "latitude", "longitude", "strict9",
)

def require(condition, message):
    if not bool(condition):
        raise ValueError("MOIS_OBSERVATION_AUDIT: " + message)

def sha_file(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for block in iter(lambda: fh.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def hash_rank(stage, key):
    return hashlib.sha256((SEED + "|" + stage + "|" + str(key)).encode()).hexdigest()

def coarse_zone(zone):
    if zone == "N_tropical":
        return "N_tropical"
    if zone == "S_tropical":
        return "S_tropical"
    if zone == "N_temperate":
        return "N_extratropical"
    if zone == "S_temperate":
        return "S_extratropical"
    return "polar_or_unknown"

def allocate(counts, target):
    """Every nonempty stratum has positive q; capped proportional top-up."""
    counts = dict(sorted((str(k), int(v)) for k, v in counts.items() if v > 0))
    require(counts and target >= len(counts) and target <= sum(counts.values()),
            "Sample budget cannot cover every nonempty stratum")
    baseline = 2 if sum(min(2, v) for v in counts.values()) <= target else 1
    chosen = {k: min(baseline, v) for k, v in counts.items()}
    rest = target - sum(chosen.values())
    capacity = {k: counts[k] - chosen[k] for k in counts}
    total_cap = sum(capacity.values())
    if rest:
        require(total_cap >= rest, "Insufficient population after baseline")
        desired = {k: rest * capacity[k] / total_cap for k in counts}
        for k in counts:
            chosen[k] += min(capacity[k], int(math.floor(desired[k])))
        leftover = target - sum(chosen.values())
        for k in sorted(counts, key=lambda k: (-(desired[k] - math.floor(desired[k])), k)):
            if not leftover:
                break
            if chosen[k] < counts[k]:
                chosen[k] += 1
                leftover -= 1
        require(leftover == 0, "Remainder allocation was incomplete")
    require(sum(chosen.values()) == target and
            all(1 <= chosen[k] <= counts[k] for k in counts),
            "Invalid per-stratum quotas")
    return chosen

def observation_frame(photo):
    """Canonically select one original ID per observation without peeking at hue."""
    require(set(COLS).issubset(photo.columns), "Frozen source columns missing")
    d = photo.loc[:, COLS].copy()
    for c in ("photo_id", "observation_id", "inat_taxon_id", "source_design"):
        d[c] = d[c].fillna("").astype(str).str.strip()
    require(d.photo_id.is_unique and d.photo_id.str.fullmatch(r"[0-9]+").all(),
            "Original photo IDs must be nonempty unique numeric strings")
    require(d.source_design.isin(ARMS).all(), "Unknown source arm")
    require(d.strict9.notna().all() and
            d.strict9.map(lambda x: isinstance(x, (bool,))).all(),
            "strict9 must be explicit boolean; unknown is not False")
    # Accept numpy.bool_ from Arrow/Pandas without accepting 'False' strings.
    # The dtype check above also catches silent strings in production sources.
    missing_obs = int(d.observation_id.eq("").sum())
    eligible = d.loc[d.observation_id.ne("")].copy()
    require(not eligible.empty, "No source observations with identity")
    eligible["_canonical_rank"] = [
        hash_rank("canonical", o + "|" + p)
        for o, p in zip(eligible.observation_id, eligible.photo_id)
    ]
    eligible = eligible.sort_values(
        ["observation_id", "_canonical_rank", "photo_id"], kind="mergesort"
    )
    chosen = eligible.drop_duplicates("observation_id", keep="first").copy()
    require(chosen.observation_id.is_unique, "Observation frame is not independent")
    report = {
        "source_original_photo_attempts": int(len(d)),
        "photos_without_observation_uuid": missing_obs,
        "source_photos_with_uuid": int(len(eligible)),
        "distinct_represented_observations": int(len(chosen)),
        "additional_photos_of_already_represented_observations":
            int(len(eligible) - len(chosen)),
        "canonical_success_count": int(chosen.strict9.sum()),
        "canonical_unknown_count": int((~chosen.strict9).sum()),
        "canonical_photo_chosen_outcome_blind": True,
    }
    return chosen, report

def draw_audit(frame, n_success=DEFAULT_SUCCESS, n_unknown=DEFAULT_UNKNOWN):
    f = frame.copy()
    f["latitude_zone"] = latitude_band(f.latitude)
    f["local_season"] = recorded_local_season(f.observed_on, f.latitude)
    f["coarse_zone"] = f.latitude_zone.map(coarse_zone)
    f["roi_status"] = f.strict9.map({True: "success", False: "unknown"})
    f["_stratum"] = (
        f.source_design + "|" + f.coarse_zone + "|" + f.local_season +
        "|" + f.roi_status
    )
    sizes = f.groupby("_stratum", sort=True).size().to_dict()
    quotas = {}
    for status, target in (("success", n_success), ("unknown", n_unknown)):
        population = {k: v for k, v in sizes.items() if k.endswith("|" + status)}
        quotas.update(allocate(population, target))
    f["_rank"] = f.observation_id.map(lambda s: hash_rank("audit", s))
    f = f.sort_values(["_stratum", "_rank", "observation_id"], kind="mergesort")
    f["_within_stratum_index"] = f.groupby("_stratum", sort=False).cumcount()
    sampled = f.loc[f["_within_stratum_index"] <
                    f["_stratum"].map(quotas)].copy()
    sampled["stratum_N"] = sampled["_stratum"].map(sizes).astype(int)
    sampled["stratum_n"] = sampled["_stratum"].map(quotas).astype(int)
    sampled["audit_inclusion_q"] = sampled.stratum_n / sampled.stratum_N
    sampled["review_inverse_q"] = sampled.stratum_N / sampled.stratum_n
    sampled = sampled.sort_values("_rank", kind="mergesort").reset_index(drop=True)
    require(len(sampled) == n_success + n_unknown and
            sampled.photo_id.is_unique and sampled.observation_id.is_unique,
            "Audited canonical photos must be unique and equal budget")
    require(sampled.roi_status.value_counts().to_dict() ==
            {"success": n_success, "unknown": n_unknown},
            "ROI success vs unknown audit quotas wrong")
    require(sampled.audit_inclusion_q.gt(0).all() and
            sampled.audit_inclusion_q.le(1).all(), "Invalid audit probability")
    strata = pd.DataFrame([
        {
            "stratum": k, "frame_observations": int(sizes[k]),
            "sampled_observations": int(quotas[k]),
            "review_probability": quotas[k] / sizes[k],
            "review_inverse_probability": sizes[k] / quotas[k],
        } for k in sorted(sizes)
    ])
    require(strata.sampled_observations.gt(0).all(),
            "Nonempty source strata must have positive sample probability")
    return sampled, strata

def freeze(source_photo, source_receipt, out, *,
           verify_frozen=True, n_success=400, n_unknown=800):
    source_photo = Path(source_photo)
    source_receipt = Path(source_receipt)
    out = Path(out)
    require(not out.exists(), "Do not overwrite a frozen audit receipt")
    r = json.loads(source_receipt.read_text())
    source_sha = sha_file(source_photo)
    require(r.get("source_georestored_original_parquet_SHA256") == source_sha,
            "Actual photo parquet does not match signed georestoration receipt")
    if verify_frozen:
        require(r.get("schema") ==
                "mois_signed1181146_trueROI9_bulk128_original_metadata_georestore_v2" and
                r.get("same_original_image_ID_and_original_9palette_unmodified") is True,
                "Missing official georestoration provenance")
    d = pd.read_parquet(source_photo, columns=list(COLS))
    if verify_frozen:
        require(len(d) == EXPECTED_PHOTOS and int(d.strict9.sum()) == EXPECTED_SUCCESS,
                "Official four-arm 1.18m frozen census changed")
    frame, summary = observation_frame(d)
    sample, strata = draw_audit(frame, n_success, n_unknown)
    out.mkdir(parents=True)
    # Both the private audit layout and the annotation layout intentionally
    # exclude model-predicted white/colour (and climate predictor values).
    roster_cols = [
        "photo_id", "observation_id", "inat_taxon_id", "source_design",
        "latitude_zone", "coarse_zone", "local_season", "roi_status",
        "_stratum", "stratum_N", "stratum_n", "audit_inclusion_q",
        "review_inverse_q",
    ]
    sample.loc[:, roster_cols].rename(columns={
        "_stratum": "stratum"
    }).to_csv(out / "audit1200_source_observation_roster.csv", index=False)
    sample.loc[:, ["photo_id", "observation_id"]].assign(
        reviewer_id="", human_review_completed=False,
        human_open_corolla_visible="", human_corolla_colour="",
        human_corolla_colour_uncertain="", failure_cause="",
        native_or_cultivated_status="", notes=""
    ).to_csv(out / "audit1200_human_annotation_BLANK.csv", index=False)
    strata.to_csv(out / "audit1200_probability_strata.csv", index=False)
    manifest = {
        "schema": SCHEMA, "status": "REAL_SOURCE_AUDIT_IDS_ONLY_NOT_HUMAN_LABELED",
        "source_georestored_photo_parquet_sha256": source_sha,
        "source_georestoration_receipt_sha256": sha_file(source_receipt),
        "source_scope": "existing_four_arm_measured_original_photo_observations_ONLY",
        "upstream_93m_photo_population_probability_identified": False,
        "native_wild_flora_population_probability_identified": False,
        "random_review_sample_is_NOT_prior_selected_512_or_128": True,
        "audit_seed": SEED, "canonical_source_photo_rule": "SHA256 observation_uuid + photo_id minimum, independent of ROI/hue",
        "primary_inferential_unit": "represented observation with single canonical photo",
        "stratification": "source_design x coarse_zone x local_season x strict9_status",
        "sample_success_count": int((sample.roi_status == "success").sum()),
        "sample_unknown_count": int((sample.roi_status == "unknown").sum()),
        "sample_total": int(len(sample)),
        "nonempty_strata": int(len(strata)),
        "audit_status_all_human_labels_blank": True,
        "photo_pixels_downloaded": 0, "new_ROI_inferences": 0,
        "frozen_photo_summary": summary,
    }
    manifest["output_sha256"] = {
        p.name: sha_file(p) for p in sorted(out.glob("*.csv"))
    }
    (out / "audit1200_frozen_design_receipt.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
    )
    print("MOIS_ORIGINAL_OBSERVATION_AUDIT_1200_SOURCE_FROZEN " +
          json.dumps({
              "observations": summary["distinct_represented_observations"],
              "extra_photos": summary["additional_photos_of_already_represented_observations"],
              "success": manifest["sample_success_count"],
              "unknown": manifest["sample_unknown_count"],
              "strata": manifest["nonempty_strata"],
              "sha_roster": manifest["output_sha256"]["audit1200_source_observation_roster.csv"],
          }, sort_keys=True), flush=True)
    return manifest

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--photo-parquet", type=Path, required=True)
    parser.add_argument("--georestored-receipt", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    freeze(args.photo_parquet, args.georestored_receipt, args.output)

if __name__ == "__main__":
    main()
