#!/usr/bin/env python3
"""Stream source-authenticated stage-1 photo screens into auditable ROI queues.

Does not fetch images, run ROI inference or declare any photograph flowering.
Default SHADOW mode runs every image through the original ROI. Staged mode is
locked until a separate independent botanical validation has been supplied.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import tempfile

SCHEMA = "mois_93m_stage1_to_stage2_probability_assignment_v1"
INPUT_COLUMNS = (
    "photo_id", "observation_id", "source_arm", "selection_probability",
    "audit_stratum", "screen_status",
)
SCREEN_STATES = frozenset(("candidate", "negative", "uncertain", "error"))
PHOTO_COLUMNS = (
    "photo_id", "observation_id", "source_arm", "selection_probability",
    "audit_stratum", "screen_status", "full_roi_selected",
    "stage2_selection_probability", "joint_design_probability",
    "joint_design_inverse_weight", "assignment_reason",
)
SEED = "mois-observation-stage2-screen-negative-audit-v1"


def require(test, message):
    if not test:
        raise ValueError(message)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def draw(observation_id):
    require(bool(observation_id), "Missing original observation ID")
    payload = (SEED + "|" + observation_id).encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "big") / 2**64


def valid_probability(value):
    try:
        p = float(value)
    except (ValueError, TypeError) as exc:
        raise ValueError("Invalid source or negative-audit probability") from exc
    require(math.isfinite(p) and 0 < p <= 1, "Probability must be in (0, 1]")
    return p


def load_policy(path):
    d = json.loads(Path(path).read_text(encoding="utf-8"))
    require(d.get("schema") == "mois_93m_prespecified_negative_audit_policy_v1",
            "Unknown or unfrozen audit policy")
    rates = d.get("negative_audit_rates")
    require(isinstance(rates, dict) and bool(rates), "Negative-audit strata missing")
    return {str(k): valid_probability(v) for k, v in rates.items()}


def verify_staged_review(review_path, input_sha, policy_sha, strata):
    require(review_path is not None, "Staged exclusions require independent botanical validation")
    receipt = json.loads(Path(review_path).read_text(encoding="utf-8"))
    require(receipt.get("schema") == "mois_independent_botanical_screen_review_v1",
            "Independent botanical reference receipt missing")
    require(receipt.get("screen_input_sha256") == input_sha
            and receipt.get("audit_policy_sha256") == policy_sha,
            "Stage-1 and audit policy SHA must match reviewed versions")
    require(receipt.get("review_basis") == "independent_photo_specific_open_flower_review"
            and receipt.get("independent_review_complete") is True,
            "Screen negatives must be independently reviewed as photographs")
    reviewed = receipt.get("review_by_stratum", {})
    require(set(strata) <= set(reviewed), "Some audit strata lack independent review")
    for stratum in strata:
        r = reviewed[stratum]
        require(isinstance(r.get("reviewed_negative_count"), int)
                and r["reviewed_negative_count"] >= 100,
                "Insufficient negative-screen independent audit for " + stratum)
        # This is a deployment safeguard, not a claim about an unobserved cohort.
        lower = r.get("white_flower_recall_lower95")
        require(isinstance(lower, (int, float))
                and math.isfinite(lower) and 0.95 <= lower <= 1,
                "White-flower screen-recall lower bound is too weak")
    return sha256(review_path)


def plan(input_csv, policy_json, output_dir, *,
         mode="shadow", review_receipt=None, expected_input_sha256=None):
    """Two streamed passes, constant row memory; no full 93m-ID Python set.

    Pass 1 authenticates the scored source and counts its audit strata. Pass 2
    streams deterministic assignments from strictly increasing photo IDs.
    All screen-negatives retain a nonzero, known second-stage sampling chance.
    """
    input_csv, policy_json, output_dir = map(Path, (input_csv, policy_json, output_dir))
    require(mode in ("shadow", "staged"), "Mode must be shadow or staged")
    require(not output_dir.exists(), "Refuse to overwrite any stage assignment")
    rates = load_policy(policy_json)
    policy_sha = sha256(policy_json)
    require(expected_input_sha256 is not None
            and len(expected_input_sha256) == 64
            and all(c in "0123456789abcdef" for c in expected_input_sha256),
            "Expected signed stage-1 CSV SHA256 is mandatory")
    source_sha = sha256(input_csv)
    require(source_sha == expected_input_sha256, "Stage-1 CSV checksum differs from signed source")
    with input_csv.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        require(set(INPUT_COLUMNS) <= set(reader.fieldnames or []),
                "Stage-1 screen is missing source identity, stratum or state")
        seen_strata = set()
        for row in reader:
            require(row["audit_stratum"] in rates, "Unspecified probability stratum")
            require(row["screen_status"] in SCREEN_STATES, "Unknown screen state")
            seen_strata.add(row["audit_stratum"])
    validation_sha = None
    if mode == "staged":
        validation_sha = verify_staged_review(
            review_receipt, source_sha, policy_sha, seen_strata)
    stats = Counter()
    root = Path(tempfile.mkdtemp(prefix="mois93m-stage-", dir=output_dir.parent))
    try:
        full_path = root / "all_photo_screen_decisions.csv"
        selected_path = root / "source_ids_for_full_roi.csv"
        with (input_csv.open(newline="", encoding="utf-8") as inp,
              full_path.open("w", newline="", encoding="utf-8") as all_out,
              selected_path.open("w", newline="", encoding="utf-8") as chosen_out):
            reader = csv.DictReader(inp)
            all_writer = csv.DictWriter(all_out, fieldnames=PHOTO_COLUMNS, lineterminator="\n")
            chosen_writer = csv.DictWriter(chosen_out, fieldnames=PHOTO_COLUMNS, lineterminator="\n")
            all_writer.writeheader()
            chosen_writer.writeheader()
            previous_photo_id = -1
            for row in reader:
                pid = row["photo_id"]
                require(pid.isdecimal(), "Malformed original photo ID")
                number = int(pid)
                require(number > previous_photo_id, "Original photo IDs must be globally sorted and unique")
                previous_photo_id = number
                obs = row["observation_id"]
                require(bool(obs) and bool(row["source_arm"]), "Missing original observation/arm")
                source_p = valid_probability(row["selection_probability"])
                stratum = row["audit_stratum"]
                status = row["screen_status"]
                require(stratum in rates and status in SCREEN_STATES, "Stage-1 source changed")
                roi_p = rates[stratum] if status == "negative" and mode == "staged" else 1.0
                keep = status != "negative" or mode == "shadow" or draw(obs) < roi_p
                joint = source_p * roi_p
                out = {
                    "photo_id": pid, "observation_id": obs, "source_arm": row["source_arm"],
                    "selection_probability": row["selection_probability"],
                    "audit_stratum": stratum, "screen_status": status,
                    "full_roi_selected": str(keep).lower(),
                    "stage2_selection_probability": format(roi_p, ".17g"),
                    "joint_design_probability": format(joint, ".17g"),
                    "joint_design_inverse_weight": format(1 / joint, ".17g"),
                    "assignment_reason": (
                        "shadow_every_image" if mode == "shadow" else
                        "candidate_or_uncertain_keep" if status != "negative" else
                        "random_negative_audit" if keep else "probability_sample_not_selected"
                    ),
                }
                all_writer.writerow(out)
                if keep:
                    chosen_writer.writerow(out)
                stats["all_screened"] += 1
                stats["queued_full_roi"] += int(keep)
                stats["source_screen_" + status] += 1
                stats["negative_sampled"] += int(keep and status == "negative")
                stats["negative_not_sampled"] += int(not keep)
        require(stats["all_screened"] > 0, "No genuine screen rows")
        report = {
            "schema": SCHEMA, "mode": mode,
            "input_stage1_csv_sha256": source_sha,
            "audit_policy_sha256": policy_sha,
            "independent_review_sha256": validation_sha,
            "full_roi_assignments_without_any_pixels_downloaded": True,
            "screen_negative_not_equated_with_nonflowering_or_nonwhite": True,
            "unselected_screen_negatives_have_nonzero_randomized_roi_inclusion": True,
            "sampling_unit": "original_observation_id",
            "source_photo_ids_are_unique_sorted": True,
            "negative_audit_rates_by_stratum": rates,
            "stage2_not_an_unweighted_census": mode == "staged",
            "counts": dict(stats),
            "output_sha256": {
                full_path.name: sha256(full_path),
                selected_path.name: sha256(selected_path),
            },
        }
        (root / "stage2_assignment_receipt.json").write_text(
            json.dumps(report, sort_keys=True, indent=2) + "\n", encoding="utf-8")
        root.replace(output_dir)
        return report
    except BaseException:
        shutil.rmtree(root, ignore_errors=True)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage1", type=Path, required=True)
    parser.add_argument("--policy", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-stage1-sha256", required=True)
    parser.add_argument("--mode", choices=("shadow", "staged"), default="shadow")
    parser.add_argument("--independent-review", type=Path)
    args = parser.parse_args()
    receipt = plan(args.stage1, args.policy, args.output, mode=args.mode,
                   review_receipt=args.independent_review,
                   expected_input_sha256=args.expected_stage1_sha256)
    print("MOIS_STAGE2_ORIGINAL_PHOTO_ASSIGNMENT", json.dumps(receipt, sort_keys=True))


if __name__ == "__main__":
    main()
