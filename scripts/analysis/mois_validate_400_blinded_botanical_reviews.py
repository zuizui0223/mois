#!/usr/bin/env python3
"""Validate independent, photo-specific MOIS 400 botanical reviews.

This module NEVER predicts open flowers, attributes colour, or consults
iNaturalist observation annotations. It checks original source membership,
independent reviewer responses, and a blind third-adjudication handoff.
The coordinator token/photo/date/site map is not to be published.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

BLANK_SHA = {
    "A": "576a55e7e613faf8c93652b4ac2d72b79f8a38bec4bfe8267c8d67bbc3fed322",
    "B": "bfbbfc92106585ca260fbd5a33573ad2598209d19907a52bbce2b941b1164203",
}
MAP_SHA = "6119a9606b4d90eb80a8a90e7fd2576631d14063de320e788e7c274fd050a3c4"
COLUMNS = (
    "review_token", "species", "original_photo_url", "photo_license",
    "attribution", "reviewer_id", "open_flower_stage", "roi_organ", "notes",
)
IMMUTABLE = COLUMNS[:5]
MAP_COLUMNS = (
    "inat_taxon_id", "species", "photo_id", "observation_id", "observer_id",
    "observed_on", "photo_license", "photo_url_large", "source_rank",
    "site_latcell", "site_loncell", "season_year", "selected_season", "review_token",
)
THIRD_COLUMNS = (
    "review_token", "species", "original_photo_url", "photo_license",
    "attribution", "adjudicator_id", "adjudicated_stage",
    "adjudicated_roi_tissue", "notes",
)
STAGES = {
    "open_flower", "buds_only", "fruit_only", "no_flower_visible", "uncertain",
}
TISSUES = {
    "petal_tepal", "other_flower_organ", "spikelet_glume_bract", "fruit_seed",
    "leaf_stem", "background", "uncertain",
}
DECISION_COLUMNS = (
    "review_token", "species", "photo_id", "observation_id", "selected_season",
    "season_year", "site_latcell", "site_loncell", "r1_vote_complete",
    "r2_vote_complete", "botanical_decision", "requires_adjudication",
    "eligible_for_future_original_photo_ROI_measurement",
)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(value: bool, message: str) -> None:
    if not value:
        raise ValueError(message)


def load(path: Path, expected: tuple[str, ...],
         *, n: int | None = None) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        require(tuple(reader.fieldnames or ()) == expected,
                f"Unexpected source/reviewer columns in {path.name}")
        rows = list(reader)
    if n is not None:
        require(len(rows) == n, f"Changed frozen row denominator in {path.name}")
    require(all(row["review_token"] and all(v is not None for v in row.values())
                for row in rows), "Empty/invalid reviewer token or malformed CSV")
    tokens = [r["review_token"] for r in rows]
    require(len(tokens) == len(set(tokens)), "Duplicate blinded review token")
    return rows


def index(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    return {row["review_token"]: row for row in rows}


def vote(row: dict[str, str], *, third: bool = False) -> tuple[str, str] | None:
    a, b, user = (
        ("adjudicated_stage", "adjudicated_roi_tissue", "adjudicator_id")
        if third else ("open_flower_stage", "roi_organ", "reviewer_id")
    )
    stage, tissue, who = (row[k].strip() for k in (a, b, user))
    if not stage and not tissue:
        # A filled reviewer ID alone is a logistical assignment, not a verdict.
        return None
    require(bool(stage) and bool(tissue) and bool(who),
            "Partial stage/tissue response or anonymous botanical review")
    require(stage in STAGES and tissue in TISSUES,
            "Unrecognised flower stage or organ code")
    return stage, tissue


def phenotype_gate(v: tuple[str, str]) -> str:
    stage, organ = v
    if "uncertain" in v:
        return "unresolved_flower_or_organ"
    if stage != "open_flower":
        return "verified_not_open_flower"
    if organ == "petal_tepal":
        return "verified_open_petaloid_eligible_for_future_ROI"
    return "verified_open_nonpetaloid_not_primary"


def verify_membership(a: list[dict], b: list[dict],
                      amap: list[dict]) -> None:
    aa, bb, mm = index(a), index(b), index(amap)
    require(aa.keys() == bb.keys() == mm.keys(),
            "Both independent reviews must contain exactly the same 400 source tokens")
    for token, ar in aa.items():
        br, m = bb[token], mm[token]
        require(all(ar[c] == br[c] for c in IMMUTABLE),
                "Review A and B have different original photo/species/licence identities")
        require(ar["species"] == m["species"] and
                ar["original_photo_url"] == m["photo_url_large"] and
                ar["photo_license"] == m["photo_license"],
                "Review form does not belong to immutable coordinator photo source")


def audited_rows(template: list[dict], reviewed: list[dict], who: str) -> tuple[dict, str | None]:
    basis, measured = index(template), index(reviewed)
    require(basis.keys() == measured.keys(), f"Reviewer {who} changes original photo membership")
    for token, row in measured.items():
        require(all(row[k] == basis[token][k] for k in IMMUTABLE),
                f"Reviewer {who} altered an immutable photo/source field")
    identities = {r["reviewer_id"].strip() for r in measured.values()
                  if r["reviewer_id"].strip()}
    require(len(identities) <= 1,
            f"Reviewer {who} packet mixes more than one reviewer identity")
    for row in measured.values():
        vote(row)  # validates complete votes, but never imputes missing responses
    return measured, next(iter(identities)) if identities else None


def write_csv(path: Path, fields: tuple[str, ...], rows: list[dict]) -> None:
    with path.open("x", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows([{k:row.get(k, "") for k in fields} for row in rows])


def audit(blank_a: Path, blank_b: Path, mapping: Path, filled_a: Path,
          filled_b: Path, outdir: Path, *, adjudication: Path | None = None,
          fixture: bool = False) -> dict:
    require(not (outdir.exists() and any(outdir.iterdir())),
            "Never overwrite prior reviewer audit outputs")
    if not fixture:
        require(sha(blank_a) == BLANK_SHA["A"] and sha(blank_b) == BLANK_SHA["B"],
                "Original separately frozen blank reviewer forms changed")
        require(sha(mapping) == MAP_SHA,
                "Coordinator original photo/date/season-to-token mapping changed")
    n = None if fixture else 400
    ba, bb = load(blank_a, COLUMNS, n=n), load(blank_b, COLUMNS, n=n)
    original = load(mapping, MAP_COLUMNS, n=n)
    require(len(ba) == len(bb) == len(original),
            "Independent review source denominators differ")
    verify_membership(ba, bb, original)
    ra, identity_a = audited_rows(ba, load(filled_a, COLUMNS, n=n), "A")
    rb, identity_b = audited_rows(bb, load(filled_b, COLUMNS, n=n), "B")
    require(identity_a is None or identity_b is None or identity_a != identity_b,
            "The same botanist cannot review both independent packets")

    original_by_token = index(original)
    pending_adjudication = []
    statuses = Counter()
    decisions = []
    pairs = {}
    n_a = n_b = n_both = 0
    for token in sorted(original_by_token):
        va, vb = vote(ra[token]), vote(rb[token])
        n_a += int(va is not None)
        n_b += int(vb is not None)
        n_both += int(va is not None and vb is not None)
        if va is None and vb is None:
            state = "both_unreviewed"
        elif va is None or vb is None:
            state = "awaiting_other_independent_reviewer"
        elif va != vb:
            state = "needs_third_independent_adjudication"
            source = ra[token]
            pending_adjudication.append({
                "review_token": token, "species": source["species"],
                "original_photo_url": source["original_photo_url"],
                "photo_license": source["photo_license"],
                "attribution": source["attribution"], "adjudicator_id": "",
                "adjudicated_stage": "", "adjudicated_roi_tissue": "",
                "notes": "",
            })
        else:
            state = phenotype_gate(va)
        statuses[state] += 1
        m = original_by_token[token]
        pairs[token] = {"state": state, "va": va, "vb": vb}
        decisions.append({
            "review_token": token, "species": m["species"],
            "photo_id": m["photo_id"], "observation_id": m["observation_id"],
            "selected_season": m["selected_season"], "season_year": m["season_year"],
            "site_latcell": m["site_latcell"], "site_loncell": m["site_loncell"],
            "r1_vote_complete": va is not None, "r2_vote_complete": vb is not None,
            "botanical_decision": state,
            "requires_adjudication": state == "needs_third_independent_adjudication",
            "eligible_for_future_original_photo_ROI_measurement":
                state == "verified_open_petaloid_eligible_for_future_ROI",
        })

    if adjudication is not None:
        third = load(adjudication, THIRD_COLUMNS, n=len(pending_adjudication))
        wanted = index(pending_adjudication)
        require(set(index(third)) == set(wanted),
                "Adjudication packet changed the blinded disagreement-only membership")
        for row in third:
            token = row["review_token"]
            require(all(row[k] == wanted[token][k] for k in IMMUTABLE),
                    "Adjudication photo/source identity was modified")
            who = row["adjudicator_id"].strip()
            require(not who or who not in (identity_a, identity_b),
                    "Third adjudicator must be independent of both first reviewers")
            v = vote(row, third=True)
            if v is None:
                continue
            state = "adjudicated_" + phenotype_gate(v)
            statuses["needs_third_independent_adjudication"] -= 1
            statuses[state] += 1
            entry = next(r for r in decisions if r["review_token"] == token)
            entry["botanical_decision"] = state
            entry["requires_adjudication"] = False
            entry["eligible_for_future_original_photo_ROI_measurement"] = (
                state == "adjudicated_verified_open_petaloid_eligible_for_future_ROI"
            )

    require(sum(statuses.values()) == len(ba), "Each photo needs one exclusive status")
    n_eligible = sum(r["eligible_for_future_original_photo_ROI_measurement"]
                     for r in decisions)
    n_unresolved = sum(v for k,v in statuses.items() if (
        k in ("both_unreviewed", "awaiting_other_independent_reviewer",
              "needs_third_independent_adjudication", "unresolved_flower_or_organ",
              "adjudicated_unresolved_flower_or_organ")
    ))
    counts_season = Counter((r["species"], r["season_year"], r["selected_season"])
                            for r in decisions if
                            r["eligible_for_future_original_photo_ROI_measurement"])
    outdir.mkdir(parents=True, exist_ok=True)
    # Mapping of true source IDs and season is private, not a Git-tracked result.
    ledger = outdir / "COORDINATOR_ONLY_400_original_photo_review_decisions.csv"
    write_csv(ledger, DECISION_COLUMNS, decisions)
    disagreements = outdir / "THIRD_REVIEW_PENDING_BLINDED_ORIGINAL_PHOTOS.csv"
    write_csv(disagreements, THIRD_COLUMNS, pending_adjudication)
    redacted = {
        "schema": "mois_400_original_photo_independent_botanical_review_QC_v1",
        "status": ("ALL_SOURCE_PHOTOS_BOTANICAL_REVIEW_RESOLVED_NO_COLOUR_MEASURED"
                   if n_unresolved == 0
                   else "BOTANICAL_REVIEW_INCOMPLETE_NO_WHITE_FLOWER_RATE"),
        "original_source_400_photo_ids": len(ba),
        "original_source_named_species": len({r["species"] for r in original}),
        "photo_membership_and_licence_unchanged": True,
        "reviewer_A_complete_votes": n_a,
        "reviewer_B_complete_votes": n_b,
        "photos_with_two_independent_completed_votes": n_both,
        "disagreement_only_original_photo_count": len(pending_adjudication),
        "status_counts": dict(sorted((k,v) for k,v in statuses.items() if v)),
        "unresolved_photo_count": n_unresolved,
        "confirmed_open_petaloid_photos_for_FUTURE_original_photo_ROI": n_eligible,
        "qualified_counts_by_species_and_season_COORDINATOR_ONLY_NOT_IN_THIS_RECEIPT": (
            len(counts_season)),
        "original_photo_pixel_ROI_colour_not_measured_in_this_audit": True,
        "white_frequency_or_monthly_climate_effect_NOT_estimated": True,
        "observation_phenology_annotation_NOT_used_as_independent_picture_verification": True,
        "unreviewed_photos_NOT_recoded_as_nonflowering": True,
        "original_review_A_SHA256": sha(filled_a),
        "original_review_B_SHA256": sha(filled_b),
        "coordinator_photo_mapping_SHA256": sha(mapping),
        "coordinator_original_photo_decisions_SHA256": sha(ledger),
        "unfilled_third_independent_review_queue_SHA256": sha(disagreements),
        "adjudication_sheet_supplied": adjudication is not None,
        "adjudication_sheet_SHA256": sha(adjudication) if adjudication else None,
        "independent_biological_observations_NOT_equal_automatic_ROI_classifiability": True,
    }
    (outdir / "MOIS_400_BOTANICAL_REVIEW_STATUS_REDACTED.json").write_text(
        json.dumps(redacted, ensure_ascii=False, indent=2, sort_keys=True,
                   allow_nan=False)+"\n", encoding="utf-8"
    )
    return redacted


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    for flag in ("blank-a", "blank-b", "coordinator-map", "filled-a", "filled-b", "outdir"):
        p.add_argument("--"+flag, type=Path, required=True)
    p.add_argument("--adjudication", type=Path)
    args = p.parse_args()
    result = audit(args.blank_a, args.blank_b, args.coordinator_map,
                   args.filled_a, args.filled_b, args.outdir,
                   adjudication=args.adjudication)
    print("MOIS_400_INDEPENDENT_BOTANICAL_REVIEW_QC",
          json.dumps({"complete_A":result["reviewer_A_complete_votes"],
                      "complete_B":result["reviewer_B_complete_votes"],
                      "eligible_future_ROI":result[
                        "confirmed_open_petaloid_photos_for_FUTURE_original_photo_ROI"],
                      "unresolved":result["unresolved_photo_count"]},
                     sort_keys=True))


if __name__ == "__main__":
    main()
