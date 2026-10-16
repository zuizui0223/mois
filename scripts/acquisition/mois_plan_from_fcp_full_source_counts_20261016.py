#!/usr/bin/env python3
"""Plan Mois flower-colour remeasurement from FCP's REAL 42,111-species count ledger.

This is a metadata-only *workload* plan. It does not download photographs,
select actual photo IDs, assess colour, or change the upstream FCP PR200
UNCAPPED biological inclusion contract. A later separate Mois measurement
sample must use a frozen photo-ID manifest and colour-blind photo selection.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import tempfile
from pathlib import Path

FRAME_SHA = "5d871bd1f1190c68fd513bb527e6c93de1606dc35ff8cfff5f6187850dc382cc"
COUNT_SHA = "ef26c16252dd785f2d08bde3b18938c7624d73a54386fa3d3974dceb09ed97a9"
SOURCE_DATE = "2026-09-27"
EXPECTED_SPECIES = 42111
EXPECTED_OPPORTUNITIES = 93209351
FIELDS = ["inat_taxon_id", "species", "n_unique_source_photo_ids", "snapshot_date"]
TARGETS = (20, 40, 100, 200)


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for b in iter(lambda: handle.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def plan(counts_path: Path, audit_path: Path, *, production=True) -> dict:
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    if (audit.get("schema") != "fcp_42111_original_photo_opportunity_audit_v1"
        or audit.get("status") != "VERIFIED_COMPLETE_COUNTS_ONLY_NOT_COLOUR_MEASUREMENT"
        or audit.get("source_frame_sha256") != FRAME_SHA
        or audit.get("snapshot_date") != SOURCE_DATE
        or audit.get("count_ledger_sha256") != sha(counts_path)
        or audit.get("original_source_taxon_IDs_and_names_independently_matched") is not True
        or audit.get("all_original_image_pixels_measured") is not False
        or audit.get("all_source_photo_id_manifests_frozen") is not False):
        raise ValueError("Unverified original full-source counts; refuse budget estimates")
    if production and (sha(counts_path) != COUNT_SHA
                       or audit.get("source_frame_species") != EXPECTED_SPECIES
                       or audit.get("all_source_photo_opportunity_ids") != EXPECTED_OPPORTUNITIES):
        raise ValueError("Original FCP ledger hash/denominator drift")
    totals = {k: 0 for k in TARGETS}
    eligible = {k: 0 for k in TARGETS}
    seen = set()
    photos = 0
    zeros = 0
    with counts_path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames != FIELDS:
            raise ValueError("Unexpected source count columns")
        for row in reader:
            sid = row["inat_taxon_id"]
            raw = row["n_unique_source_photo_ids"]
            if (not sid.isdecimal() or int(sid) <= 0 or int(sid) in seen
                or not raw.isdecimal() or not row["species"].strip()
                or row["snapshot_date"] != SOURCE_DATE):
                raise ValueError("Duplicate/malformed original source taxon or count")
            seen.add(int(sid))
            n = int(raw)
            photos += n
            zeros += (n == 0)
            for cap in TARGETS:
                totals[cap] += min(n, cap)
                eligible[cap] += (n >= cap)
    if (len(seen) != audit.get("source_frame_species")
        or photos != audit.get("all_source_photo_opportunity_ids")
        or zeros != audit.get("taxa_with_zero_photos")):
        raise ValueError("Count totals fail independent audit reconciliation")
    return {
        "schema": "mois_fcp_42111_metadata_only_measurement_workload_scenarios_v1",
        "status": "BUDGET_ONLY_NOT_PHOTO_ID_SELECTION_OR_COLOUR",
        "source_date": SOURCE_DATE,
        "source_frame_sha256": FRAME_SHA,
        "source_count_ledger_sha256": sha(counts_path),
        "source_taxa": len(seen),
        "source_photo_opportunity_records": photos,
        "source_taxa_with_no_photos": zeros,
        "scenario_targets_per_taxon": {
            str(k): {
                "maximum_candidate_photo_processing_units": totals[k],
                "species_with_at_least_target_source_opportunities": eligible[k],
                "relative_to_100543_historical_measured_photo_ids": totals[k] / 100543,
                "fraction_of_93209351_source_opportunities": totals[k] / photos,
            } for k in TARGETS
        },
        "not_a_floral_colour_result": True,
        "not_unique_photo_ids_frozen": True,
        "not_original_pixel_download": True,
        "not_genus_or_geography_balanced": True,
        "upstream_fcp_pr200_full_inclusion_remains_uncapped": True,
        "notes": [
            "This separate Mois planning cap is a computational budget, not a biological claim that FCP has capped source photos.",
            "No photo ID is selected; colour-blind stable-ID hashing requires a complete independently validated source photo-ID manifest.",
            "Sampling must retain source inclusion probabilities, taxon, observation, observer, location and missing-flower outcomes.",
            "Counts are upper bounds on processable/evaluable floral images and are not independent plant observations.",
        ],
    }


def self_test():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        p = root / "counts.csv"
        with p.open("w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(FIELDS)
            for i, n in enumerate((0, 5, 60, 300), 1):
                w.writerow((i, f"Taxon_{i}", n, SOURCE_DATE))
        a = root / "audit.json"
        a.write_text(json.dumps({
            "schema": "fcp_42111_original_photo_opportunity_audit_v1",
            "status": "VERIFIED_COMPLETE_COUNTS_ONLY_NOT_COLOUR_MEASUREMENT",
            "source_frame_sha256": FRAME_SHA, "snapshot_date": SOURCE_DATE,
            "count_ledger_sha256": sha(p),
            "original_source_taxon_IDs_and_names_independently_matched": True,
            "all_original_image_pixels_measured": False,
            "all_source_photo_id_manifests_frozen": False,
            "source_frame_species": 4, "all_source_photo_opportunity_ids": 365,
            "taxa_with_zero_photos": 1,
        }), encoding="utf-8")
        x = plan(p, a, production=False)
        assert x["source_photo_opportunity_records"] == 365
        assert x["scenario_targets_per_taxon"]["20"]["maximum_candidate_photo_processing_units"] == 45
        assert x["scenario_targets_per_taxon"]["40"]["maximum_candidate_photo_processing_units"] == 85
        assert x["scenario_targets_per_taxon"]["100"]["maximum_candidate_photo_processing_units"] == 165
        assert x["scenario_targets_per_taxon"]["200"]["maximum_candidate_photo_processing_units"] == 265
        assert x["scenario_targets_per_taxon"]["40"]["species_with_at_least_target_source_opportunities"] == 2
        try:
            plan(p, a, production=True)
            raise AssertionError("Unpinned tiny source passed production guard")
        except ValueError:
            pass
    print("PASS metadata-only budget and source gate; NO photo IDs, pixels or colours acquired")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--counts", type=Path)
    p.add_argument("--independent-audit", type=Path)
    p.add_argument("--output", type=Path)
    p.add_argument("--self-test", action="store_true")
    args = p.parse_args()
    if args.self_test:
        self_test()
        return
    if not args.counts or not args.independent_audit or not args.output:
        p.error("Provide verified --counts, --independent-audit and --output")
    r = plan(args.counts, args.independent_audit)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(r, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(r, sort_keys=True))


if __name__ == "__main__":
    main()
