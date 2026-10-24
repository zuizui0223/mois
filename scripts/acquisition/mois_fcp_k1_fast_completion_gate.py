#!/usr/bin/env python3
"""Fast fail-closed metadata gate for SHA-sealed genuine K1 original-photo jobs.

Previous K1 reader scanned *all* ~2800 GitHub release tags every 10 minutes.
Prefer the latest independently SHA-verified immutable K1 cumulative receipt.
This narrows live GitHub release lookups to newly succeeded but unsealed worker
groups, without touching photo bytes, ROI model, or inclusion probabilities.
On missing cumulative evidence, fall back to complete release pagination.
"""
from __future__ import annotations

import hashlib
import json
import re

SCHEMA="mois_fcp_K1_93m_balanced_true_original_ROI_append_group_v1"
CUM_SCHEMA="mois_fcp_K1_original_photo_cumulative_signed_group_ledger_v1"
SOURCE_RUN=40066261878
SOURCE_ARM="mois_global93m_tier06_taxon_month_lat_probability_onephoto_K1_v1"
SOURCE_SHA="4da66dc0d34f440269c2202f95b0f39ad369531fd625676e8ec69df4ae939f1c"
REUSE=3280
SOURCE_TARGET=76753
MAX_GROUPS=75

CUM_PATTERN=re.compile(
    r"mois-K1-genuine-completed-original-ROI-(\d+)groups-([0-9a-f]{16})-v1\Z")
WORKER_PATTERN=re.compile(r"warm-original-ROI \((\d+)\)\Z")
ARTIFACT_PATTERN=re.compile(r"mois-balanced-onephoto-K1-real-RGB-group-(\d+)\Z")
SEALED_TAG="mois-full93m-balanced-onephoto-K1-original-RGB-append-group-{}-v1"
SHA64=re.compile(r"sha256:[0-9a-f]{64}\Z")


def require(value,reason):
    if not value:
        raise ValueError(reason)


def newest_snapshot_tag(release_tags):
    """Newest signed group-count tag from one recent releases page, if present."""
    groups={}
    for tag in release_tags:
        m=CUM_PATTERN.fullmatch(str(tag))
        if not m:
            continue
        n=int(m.group(1))
        require(1<=n<=MAX_GROUPS,"Out-of-range K1 cumulative source release")
        if n in groups and groups[n]!=tag:
            raise ValueError("Conflicting equally sized immutable K1 cumulative evidence")
        groups[n]=tag
    return groups[max(groups)] if groups else None


def verified_archived_groups(snapshot, tag):
    """Never treat a release title alone as proof that images were processed."""
    match=CUM_PATTERN.fullmatch(str(tag))
    require(bool(match),"Unrecognized original K1 cumulative release")
    n=int(match.group(1))
    require(
        snapshot.get("schema")==CUM_SCHEMA
        and snapshot.get("source_run")==SOURCE_RUN
        and snapshot.get("total_scheduled_original_worker_groups")==MAX_GROUPS
        and snapshot.get("new_original_photo_target_total")==SOURCE_TARGET
        and snapshot.get("previously_measured_K1_photo_IDs_reused_elsewhere_NOT_NEW")==REUSE
        and snapshot.get("source_precolour_probability_photo_sampling_preserved") is True
        and snapshot.get("within_K1_original_source_group_units_disjoint") is True
        and snapshot.get("new_original_images_downloaded_or_colour_inferred_in_THIS_ledger")==0
        and snapshot.get("this_is_NOT_a_crossarm_photo_ID_deduplicated_global_total") is True
        and snapshot.get("global_wild_flora_or_blooming_petal_status_not_independently_verified") is True
        and snapshot.get("completed_and_SHA_sealed_original_worker_groups")==n,
        "Previously published K1 ROI source receipts/probability contract changed")
    groups=snapshot.get("completed_original_worker_group_IDs")
    receipts=snapshot.get("per_group")
    require(isinstance(groups,list) and isinstance(receipts,list)
            and len(groups)==len(receipts)==n
            and all(type(g) is int and 0<=g<MAX_GROUPS for g in groups)
            and groups==sorted(set(groups)),
            "Cumulative K1 original worker group source identity invalid")
    per={}
    for row in receipts:
        g=row.get("original_worker_group")
        signed=row.get("original_group_receipt_SHA256")
        require(type(g) is int and g in groups and g not in per
                and isinstance(signed,str) and re.fullmatch(r"[0-9a-f]{64}",signed)
                and type(row.get("new_original_photo_attempts")) is int
                and row["new_original_photo_attempts"]==min(1024,SOURCE_TARGET-1024*g)
                and type(row.get("four_hue_classified")) is int
                and 0<=row["four_hue_classified"]<=row["new_original_photo_attempts"],
                "Invalid recorded K1 original photo worker terminal/group identity")
        per[g]=signed
    require(sorted(per)==groups,"Signed K1 terminal receipt list does not match published groups")
    digest=hashlib.sha256("\n".join(
        f"{g}:{per[g]}" for g in groups).encode()).hexdigest()
    require(
        snapshot.get("source_original_group_receipt_aggregate_SHA256")==digest
        and digest.startswith(match.group(2)),
        "Cumulative K1 original pixel-group receipt SHA/tag mismatch")
    return set(groups)


def successful_finished_groups(jobs, artifacts):
    """Only full genuine K1 job success + immutable SHA ZIP makes a group eligible."""
    require(isinstance(jobs,list) and isinstance(artifacts,list),
            "Missing K1 original worker jobs/artifacts inventory")
    names=set()
    successes=set()
    for job in jobs:
        m=WORKER_PATTERN.fullmatch(str(job.get("name","")))
        if not m:
            continue
        group=int(m.group(1))
        require(group in range(MAX_GROUPS) and group not in names,
                "Malformed or duplicate K1 original work group")
        names.add(group)
        if job.get("conclusion")=="success":
            successes.add(group)
    require(names==set(range(MAX_GROUPS)),
            "Frozen 75-group K1 original image job roster changed")
    indexed={}
    for artifact in artifacts:
        m=ARTIFACT_PATTERN.fullmatch(str(artifact.get("name","")))
        if not m:
            continue
        g=int(m.group(1))
        require(g in range(MAX_GROUPS) and g not in indexed,
                "Duplicate/out-of-range K1 original photo worker ZIP")
        require(artifact.get("expired") is False
                and SHA64.fullmatch(str(artifact.get("digest","")))
                and type(artifact.get("size_in_bytes")) is int
                and artifact["size_in_bytes"]>0,
                "Original completed photo ZIP expired or lacked authenticated SHA")
        indexed[g]=artifact
    require(set(indexed)<=successes,
            "Sealed original photo result without matching successful K1 image worker")
    missing=successes-set(indexed)
    require(not missing,"Successful original photo work groups lack their saved ZIPs")
    return successes


def groups_needing_archive(successful, archived, has_individual_release):
    """Query ONLY new genuine groups absent from the last source-verified snapshot."""
    require(isinstance(archived,set) and archived<=successful,
            "Previously archived original-photo group not among successful K1 workers")
    ready=[]
    for group in sorted(successful-archived):
        if not has_individual_release(group):
            ready.append(group)
    return ready


def missing_from_all_releases(successful, all_tags):
    """Fallback when no signed cumulative evidence is found, never guess complete."""
    existing=set(all_tags)
    return [g for g in sorted(successful) if SEALED_TAG.format(g) not in existing]
