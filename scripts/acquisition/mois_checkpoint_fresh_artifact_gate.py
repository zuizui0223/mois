#!/usr/bin/env python3
"""Only request a new *data-only* photo checkpoint when sealed worker IDs change.

Never downloads photo pixels, recolours a flower, deletes unknowns, mutates
previous releases, or treats unfinished acquisition targets as measurements.
The checkpoint collector independently rechecks original source and ROI SHAs.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

WAVES = {
    "p0004": 39353343403,
    "p0008": 39371199127,
    "shadow": 39369517433,
    "p0016": 39374811084,
    "p0032": 39401484641,
}
PREFIXES = {
    "p0004": "mois-93m-p0004-RGB-group-",
    "p0008": "mois-93m-p0008-RGB-group-",
    "shadow": "mois-tier03-shadow-b0-RGB-group-",
    "p0016": "mois-93m-p0016-RGB-group-",
    "p0032": "mois-93m-p0032-RGB-group-",
}
SHA = re.compile(r"sha256:[0-9a-f]{64}\Z")
CHECKPOINT = re.compile(r"mois-real-pixel-checkpoint-([0-9]+)-v1\Z")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def latest_checkpoint_tag(releases):
    """Pick the latest *published* exact checkpoint, ignoring unrelated drafts."""
    relevant = []
    for release in releases:
        tag = str(release.get("tag_name") or "")
        match = CHECKPOINT.fullmatch(tag)
        if match and release.get("draft") is False and release.get("published_at"):
            relevant.append((release["published_at"], int(match.group(1)), tag))
    return max(relevant)[2] if relevant else None


def validate_snapshots(snapshots, runs=WAVES):
    """Fail closed on altered artifact identities and unknown acquisition waves."""
    require(isinstance(snapshots, dict) and set(snapshots) == set(runs),
            "Missing or unexpected waves in source-authenticated artifact snapshot")
    for wave, run in runs.items():
        entry = snapshots[wave]
        require(isinstance(entry, dict) and entry.get("source_run_id") == run,
                "Source run ID changed for wave " + wave)
        artifacts = entry.get("artifacts")
        require(isinstance(artifacts, list), "No artifact listing for " + wave)
        seen_ids, seen_names = set(), set()
        for artifact in artifacts:
            require(isinstance(artifact, dict) and
                    type(artifact.get("id")) is int and artifact["id"] > 0 and
                    isinstance(artifact.get("name"), str) and artifact["name"] and
                    isinstance(artifact.get("digest"), str) and
                    SHA.fullmatch(artifact["digest"]) is not None and
                    type(artifact.get("size_in_bytes")) is int and
                    artifact["size_in_bytes"] > 0,
                    "Invalid sealed artifact identity/digest/size for " + wave)
            require(artifact["id"] not in seen_ids and
                    artifact["name"] not in seen_names,
                    "Artifact duplicate within original wave " + wave)
            seen_ids.add(artifact["id"])
            seen_names.add(artifact["name"])
    return snapshots


def decide_checkpoint(previous, current, runs=WAVES):
    """Compare published seal IDs to live seals, without classifying images."""
    validate_snapshots(current, runs)
    if previous is not None:
        validate_snapshots(previous, runs)
    deltas = {}
    for wave in runs:
        old_items = (previous or {}).get(wave, {}).get("artifacts", [])
        now_items = current[wave]["artifacts"]
        old = {x["id"]: x for x in old_items}
        now = {x["id"]: x for x in now_items}
        for ident, item in old.items():
            require(ident in now and now[ident] == item,
                    "Published sealed photo artifact changed or missing: " + wave)
        deltas[wave] = len(now) - len(old)
    total = sum(deltas.values())
    return {
        "schema": "mois_checkpoint_only_new_sealed_RGB_artifact_gate_v1",
        "ready": total > 0,
        "new_worker_artifacts": total,
        "new_artifacts_by_wave": deltas,
        "previous_published_snapshot_used": previous is not None,
        "new_worker_completion_is_not_new_photo_hue_count": True,
        "image_downloads_or_colour_reinference_performed": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--previous", type=Path,
                        help="Last published checkpoint receipt JSON; omit if none")
    parser.add_argument("--current", type=Path, required=True,
                        help="Live source_artifact_snapshots JSON")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    previous = (json.loads(args.previous.read_text())["source_artifact_snapshots"]
                if args.previous else None)
    current = json.loads(args.current.read_text())
    require(not args.output.exists(), "Refuse to overwrite prior checkpoint decision")
    result = decide_checkpoint(previous, current)
    args.output.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
    print("MOIS_SEALED_ORIGINAL_RGB_CHECKPOINT_GATE", json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
