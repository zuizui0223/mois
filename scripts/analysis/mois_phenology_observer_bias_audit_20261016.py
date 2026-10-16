#!/usr/bin/env python3
"""Audit iNaturalist Flower/Fruits annotation and observer effort in Mois.

Annotations are per OBSERVATION, not per photo. They are absent from the
2026-09-27 four-table iNaturalist Open Data snapshot. This consumes an
independently captured, checksum-pinned, *private* observation-UUID annotation
overlay after sampling; it does NOT silently call the API or treat no annotation
as "not flowering".

Term 12 currently "Flowers and Fruits" (formerly Plant Phenology):
13 Flowers, 15 Flower Buds, 14 Fruits/Seeds, 21 No Flowers/Fruits.
Multiple terms may co-occur. Contradictory 13/21 flags are quarantined.

The output is descriptive evidence about selection/annotation/classifiability
conditional on existing iNaturalist opportunities, not a correction for
unobserved plants, platform upload intensity, or evolutionary morph frequencies.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from datetime import date
from pathlib import Path

ANNOTATION_TERM_ID = "12"
ANNOTATION_VALUE_IDS = {13, 14, 15, 21}
OVERLAY_COLUMNS = [
    "observation_id", "annotation_term_id", "annotation_value_ids",
    "lookup_status", "annotation_source", "lookup_snapshot_date"
]
LOOKUPS = {"LOOKED_UP", "NOT_LOOKED_UP", "UNAVAILABLE"}
SOURCES = {"inat_api_v1", "inat_export", "authorized_manual_review"}
STATUS = {
    "FLOWERS", "BUDS_ONLY", "FRUITS_ONLY", "NO_FLOWERS_FRUITS",
    "NO_ANNOTATION", "NOT_VERIFIED"
}
JOINED_COLUMNS = {
    "observation_id", "observer_id", "inat_taxon_id", "genus",
    "observed_on", "region_stratum", "precip_stratum",
    "photo_design_weight", "observation_inclusion_probability",
    "measurement_status", "morph"
}
MIN_PUBLIC_CELL = 5


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def annotation_state(row: dict[str, str]) -> str:
    if row["annotation_term_id"] != ANNOTATION_TERM_ID:
        raise ValueError("Unexpected iNaturalist phenology annotation term: expected 12")
    lookup = row["lookup_status"]
    source = row["annotation_source"]
    if lookup not in LOOKUPS or source not in SOURCES:
        raise ValueError("Unknown annotation lookup source or status")
    if not row["lookup_snapshot_date"]:
        raise ValueError("Missing annotation lookup snapshot date")
    date.fromisoformat(row["lookup_snapshot_date"])
    value = row["annotation_value_ids"].strip()
    if lookup != "LOOKED_UP":
        if value:
            raise ValueError("Unverified annotation is not biological non-flowering evidence")
        return "NOT_VERIFIED"
    if not value:
        return "NO_ANNOTATION"
    ids = value.split("|")
    if (not all(x.isdecimal() for x in ids)
        or len(ids) != len(set(ids))
        or {int(i) for i in ids} - ANNOTATION_VALUE_IDS):
        raise ValueError("Unknown or duplicated annotation value ID")
    found = {int(i) for i in ids}
    if 21 in found and len(found) > 1:
        raise ValueError("Contradictory flowering and no-flower observation annotations")
    if 13 in found:
        return "FLOWERS"
    if 15 in found:
        return "BUDS_ONLY"
    if 14 in found:
        return "FRUITS_ONLY"
    return "NO_FLOWERS_FRUITS"


def load_overlay(path: Path, pin: str) -> tuple[dict[str,str],dict]:
    if len(pin) != 64 or any(c not in "0123456789abcdef" for c in pin):
        raise ValueError("A real byte-pinned phenology overlay SHA256 is mandatory")
    if digest(path) != pin:
        raise ValueError("Private phenology overlay SHA256 does not match frozen source")
    rows = {}
    sources = set()
    snapshots = set()
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames != OVERLAY_COLUMNS:
            raise ValueError("Expected exact observation-level phenology overlay columns")
        for r in reader:
            ident = r["observation_id"].strip()
            if not ident or ident in rows:
                raise ValueError("Missing/duplicated original iNaturalist observation UUID")
            rows[ident] = annotation_state(r)
            sources.add(r["annotation_source"])
            snapshots.add(r["lookup_snapshot_date"])
    if len(snapshots) != 1:
        raise ValueError("Mixed annotation snapshot dates; not one reproducible source")
    return rows, {
        "annotation_source_kinds": sorted(sources),
        "annotation_snapshot_date": next(iter(snapshots)),
        "annotation_overlay_sha256": pin,
        "coverage_was_fully_looked_up": all(x != "NOT_VERIFIED" for x in rows.values())
    }


def month_token(raw: str) -> str:
    if not raw:
        return "UNKNOWN_MONTH"
    try:
        return f"{date.fromisoformat(raw).month:02}"
    except ValueError:
        return "UNKNOWN_MONTH"


def mean(values, weights) -> float | None:
    z = sum(weights)
    return sum(v*w for v,w in zip(values,weights))/z if z else None


def audit(joined_csv: Path, sample_receipt: Path, overlay_csv: Path,
          overlay_sha256: str, out_private: Path, out_public: Path,
          *, synthetic: bool = False) -> dict:
    if out_private.exists() or out_public.exists():
        raise ValueError("Refuse to overwrite observation-bias audit")
    evidence = json.loads(sample_receipt.read_text(encoding="utf-8"))
    if (evidence.get("schema") != "mois_full_selected_photo_roi_terminal_collection_v1"
        or evidence.get("status") !=
          "COMPLETE_MOIS_SAMPLED_PHOTO_TRAITS_CONDITIONAL_CLASSIFIABILITY_NO_13ENV_FIT"
        or evidence.get("private_joined_photo_csv_sha256") != digest(joined_csv)
        or evidence.get("synthetic_only_test") is not synthetic
        or evidence.get("full_FCP_93m_photos_measured") is not False
        or evidence.get("all13_photo_site_predictors_joined") is not False):
        raise ValueError("Only exactly sealed Mois sampled-photo terminal evidence is allowed")
    overlay, origin = load_overlay(overlay_csv, overlay_sha256)
    by_status = defaultdict(lambda: {
        "n":0, "classified":0, "white":0, "weight":0., "weighted_class":0.,
        "weighted_white":0., "sumw2":0., "obs_weight":0.,
    })
    by_month, by_region, by_precip, by_observer = [defaultdict(int) for _ in range(4)]
    per_observer_weight = defaultdict(float)
    seen_observations = set()
    checked = 0
    unknown_month = 0
    contradicted = 0
    total = 0
    with joined_csv.open(newline="",encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if not JOINED_COLUMNS.issubset(reader.fieldnames or []):
            raise ValueError("Missing private joined photo site, observer or colour evidence")
        for row in reader:
            key = row["observation_id"]
            if key in seen_observations or key not in overlay:
                raise ValueError("Every selected observation must have EXACTLY one annotation lookup row")
            seen_observations.add(key)
            stage = overlay[key]
            if not row["inat_taxon_id"].isdigit() or not row["genus"]:
                raise ValueError("Missing immutable observation taxon/genus identity")
            w = float(row["photo_design_weight"])
            pi = float(row["observation_inclusion_probability"])
            if not all(math.isfinite(v) for v in (w,pi)) or w < 1 or not 0 < pi <= 1:
                raise ValueError("Invalid source sample probability and design weight")
            obs_weight = 1/pi  # one original photo for each chosen observation
            if row["measurement_status"]=="classified_four_state_morph":
                if row["morph"] not in ("white","yellow_orange","red_pink","blue_purple"):
                    raise ValueError("Classifiable original flower has invalid 4-colour state")
                valid = True
            elif row["morph"]=="mixed_uncertain":
                valid = False
            else:
                raise ValueError("No valid ROI terminal for original source observation")
            white = valid and row["morph"]=="white"
            d = by_status[stage]
            d["n"] += 1
            d["classified"] += int(valid)
            d["white"] += int(white)
            d["weight"] += obs_weight
            d["weighted_class"] += obs_weight*valid
            d["weighted_white"] += obs_weight*white
            d["sumw2"] += obs_weight*obs_weight
            d["obs_weight"] += obs_weight
            # No private observation/site identifiers are written to the
            # public receipt; the private report stores only aggregates.
            month = month_token(row["observed_on"])
            by_month[month] += 1
            by_region[(row["region_stratum"],stage)] += 1
            by_precip[(row["precip_stratum"],stage)] += 1
            user = row["observer_id"].strip() or "UNKNOWN_OBSERVER"
            by_observer[user] += 1
            per_observer_weight[user] += obs_weight
            unknown_month += (month=="UNKNOWN_MONTH")
            checked += stage!="NOT_VERIFIED"
            contradicted += stage=="NO_FLOWERS_FRUITS" and valid
            total += 1
    if set(overlay) != seen_observations or total!=evidence["original_photo_terminals_verified"]:
        raise ValueError("Annotation lookup universe differs from the sealed selected observations")
    n_obs = len(by_observer)
    counts = sorted(by_observer.values(),reverse=True)
    top = max(1,math.ceil(n_obs*.01))
    effective = (sum(per_observer_weight.values())**2
                 /sum(v*v for v in per_observer_weight.values()))
    summary_status = {}
    for group in sorted(STATUS):
        d=by_status[group]
        pw=d["weighted_class"]/d["weight"] if d["weight"] else None
        white=d["weighted_white"]/d["weighted_class"] if d["weighted_class"] else None
        summary_status[group] = {
            "selected_observations": d["n"],"ROI_four_colour_classifiable":d["classified"],
            "observed_white_among_classifiable":d["white"],
            "design_weighted_ROI_classifiability":pw,
            "design_weighted_white_given_ROI_classifiable":white,
            "effective_sampled_observations_kish":d["weight"]**2/d["sumw2"] if d["sumw2"] else 0,
        }
    assert sum(v["selected_observations"] for v in summary_status.values())==total
    # Anonymous aggregates only, with no observer IDs or raw coordinates.
    pub = {
        "schema":"mois_original_sample_phenology_observer_bias_audit_v1",
        "status":"PHOTO_SELECTED_OBSERVATION_BIAS_AUDIT_NOT_POPULATION_FLOWERING_FREQUENCY",
        "sample_roi_collection_receipt_sha256":digest(sample_receipt),
        "original_observation_annotation_source":origin,
        "selected_observations":total,
        "annotations_looked_up":checked,
        "annotations_NOT_verified":total-checked,
        "phenology_by_observation_status":summary_status,
        "distinct_observers_in_selected_source":n_obs,
        "fraction_of_selected_observations_from_largest_1pct_observers":sum(counts[:top])/total if total else None,
        "effective_observers_after_observation_design_weighting":effective,
        "n_observations_missing_observation_month":unknown_month,
        "selected_month_counts":dict(sorted(by_month.items())),
        "contradictions_explicit_no_flowers_but_ROI_classified":contradicted,
        "contradictions_require_blind_manual_review":contradicted>0,
        "true_floral_state_for_unannotated_observations_imputed":False,
        "annotation_is_observation_not_photo_label":True,
        "unsubmitted_plants_observation_bias_identified":False,
        "platform_level_phenology_annotation_complete":False,
        "white_environment_causal_fitness_result":False,
        "notes":(
          "Flowers (13) is positive phenology evidence, not a guarantee that the "
          "specific sampled photo passes ROI-v4. Absent annotation is not negative. "
          "No values here estimate true wild-flora colour/morph frequency."
        ),
    }
    private = {
        "schema":"mois_private_selected_observer_phenology_bias_ledger_v1",
        "source_sample_sha256":digest(joined_csv),
        "observer_n_selected":dict(sorted(by_observer.items())),
        "region_x_annotation_n_selected":[
            {"region":k[0],"annotation":k[1],"n":v}
            for k,v in sorted(by_region.items())
        ],
        "precip_x_annotation_n_selected":[
            {"precip":k[0],"annotation":k[1],"n":v}
            for k,v in sorted(by_precip.items())
        ],
        "public_aggregate_sha256":None,
        "contains_private_observer_identifiers":True,
        "contains_original_observation_identifiers":False,
    }
    out_private.parent.mkdir(parents=True,exist_ok=True)
    out_public.parent.mkdir(parents=True,exist_ok=True)
    pub_bytes=(json.dumps(pub,sort_keys=True,indent=2)+"\n").encode()
    private["public_aggregate_sha256"]=hashlib.sha256(pub_bytes).hexdigest()
    out_private.write_text(json.dumps(private,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    out_public.write_bytes(pub_bytes)
    return pub


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--private-joined-sample",type=Path,required=True)
    p.add_argument("--sample-collection-receipt",type=Path,required=True)
    p.add_argument("--observation-annotation-overlay",type=Path,required=True)
    p.add_argument("--annotation-overlay-sha256",required=True)
    p.add_argument("--out-private-aggregates",type=Path,required=True)
    p.add_argument("--out-public-aggregate",type=Path,required=True)
    a=p.parse_args()
    result=audit(a.private_joined_sample,a.sample_collection_receipt,
                 a.observation_annotation_overlay,a.annotation_overlay_sha256,
                 a.out_private_aggregates,a.out_public_aggregate)
    print(json.dumps({k:result[k] for k in (
        "status","selected_observations","annotations_looked_up",
        "annotations_NOT_verified","distinct_observers_in_selected_source",
        "contradictions_explicit_no_flowers_but_ROI_classified"
    )},sort_keys=True))


if __name__=="__main__":
    main()
