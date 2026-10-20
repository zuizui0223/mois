#!/usr/bin/env python3
"""Resume publication of already SHA-verified 93m source target archives.

The actual 16 worker artifacts are immutable; this script only uploads their
existing 224 tier archives plus signed receipts into an already-created DRAFT
GitHub release, by exact release ID. Never redownload or remeasure source
photos. Never mutate an existing release asset, even if corrupted.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import urllib.error
import urllib.parse
import urllib.request

PHOTO_TOTAL = 93209351
RELEASE_TAG = "mois-complete-93209351-original-photo-target-tiers-v1"
RELEASE_ID = 417967155
REPO = "zuizui0223/mois"
ARCHIVE_PATTERN = re.compile(r"tier-(?:0[0-9]|1[0-3])-group-(?:0[0-9]|1[0-5])\.tar")
RECEIPT_PATTERN = re.compile(r"group-(?:0[0-9]|1[0-5])-receipt\.json")
MANIFEST_NAME = "mois_all_original_photo_target_tiers_receipt.json"
API = "https://api.github.com"


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha256_file(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(2 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def expected_files(root):
    root = Path(root)
    m = root / MANIFEST_NAME
    report = json.loads(m.read_text(encoding="utf-8"))
    require(report.get("schema") == "mois_complete_full93m_original_photo_target_registry_v1"
            and report.get("exact_original_photo_IDs_in_all_disjoint_tiers") == PHOTO_TOTAL
            and report.get("native_vs_introduced_not_hard_filtered") is True
            and report.get("source_only_NO_photo_pixels_or_hue_measured") is True
            and sum(report["new_target_counts_by_tier"].values()) == PHOTO_TOTAL,
            "Only the source-authenticated 93209351 original photo target registry may publish")
    expected_tar = report["published_tier_archive_SHA256"]
    require(len(expected_tar) == 224, "All 16 x 14 tier archives must be source-signed")
    names = {f"tier-{tier:02d}-group-{g:02d}.tar"
             for tier in range(14) for g in range(16)}
    require(set(expected_tar) == names, "Tier archive names are incomplete or duplicated")
    tar_files = sorted(p for p in root.glob("tier-*.tar") if ARCHIVE_PATTERN.fullmatch(p.name))
    group_files = sorted(p for p in root.glob("group-*-receipt.json")
                         if RECEIPT_PATTERN.fullmatch(p.name))
    require({p.name for p in tar_files} == names and len(tar_files) == 224
            and {p.name for p in group_files} ==
            {f"group-{g:02d}-receipt.json" for g in range(16)},
            "Missing original signed worker archives or all 16 source receipts")
    files = tar_files + group_files + [m]
    hashes = {}
    for file in files:
        digest = sha256_file(file)
        if file.name in expected_tar:
            require(digest == expected_tar[file.name],
                    "Immutable real source tier archive SHA changed: " + file.name)
        hashes[file.name] = "sha256:" + digest
    return files, hashes, report


def api_request(token, method, url, body=None):
    request = urllib.request.Request(
        url, data=(json.dumps(body).encode("utf-8") if body is not None else None),
        method=method,
        headers={
            "Authorization": "Bearer " + token,
            "Accept": "application/vnd.github+json",
            "X-GitHub-Api-Version": "2022-11-28",
            **({"Content-Type": "application/json"} if body is not None else {}),
        },
    )
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            return json.loads(response.read())
    except urllib.error.HTTPError as exc:
        raise RuntimeError(f"GitHub release {method} HTTP {exc.code}: "
                           f"{exc.read(300).decode('utf-8', 'replace')}") from exc


def release_info(token, repo, release_id):
    release = api_request(token, "GET", f"{API}/repos/{repo}/releases/{release_id}")
    require(release["tag_name"] == RELEASE_TAG and release["id"] == release_id,
            "Refuse to publish another GitHub Release/tag")
    return release


def published_assets(token, repo, release_id):
    all_assets = {}
    page = 1
    while True:
        batch = api_request(token, "GET",
                            f"{API}/repos/{repo}/releases/{release_id}/assets"
                            f"?per_page=100&page={page}")
        for asset in batch:
            name = asset["name"]
            require(name not in all_assets, "Duplicate asset names on original release")
            all_assets[name] = asset
        if len(batch) < 100:
            break
        page += 1
        require(page <= 20, "Too many unknown assets on the target Release")
    return all_assets


def check_existing_assets(assets, hashes):
    for name, asset in assets.items():
        require(name in hashes, "Unexpected old asset in draft Release: " + name)
        require(asset.get("state") == "uploaded"
                and asset.get("digest") == hashes[name],
                "Existing Release asset missing or mismatching exact original byte SHA: " + name)


def upload_asset(token, repo, release_id, file, expected_hash):
    # Upload GitHub draft assets by numeric ID. Tag-name lookup returns 404 for
    # unpublished drafts; gh release upload "$tag" is therefore not valid here.
    name = urllib.parse.quote(file.name, safe="")
    url = (f"https://uploads.github.com/repos/{repo}/releases/"
           f"{release_id}/assets?name={name}")
    completed = subprocess.run(
        [
            "curl", "--fail-with-body", "--location", "--silent", "--show-error",
            "--retry", "3", "--retry-all-errors", "--connect-timeout", "30",
            "--max-time", "1200",
            "-H", "Authorization: Bearer " + token,
            "-H", "Accept: application/vnd.github+json",
            "-H", "Content-Type: application/octet-stream",
            "--data-binary", "@" + str(file), url,
        ],
        text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        timeout=1250, check=False,
    )
    if completed.returncode:
        # Authentication is not printed. A failed upload must not automatically
        # delete any incomplete GitHub asset; subsequent inspection is required.
        raise RuntimeError("Release draft upload failed for " + file.name +
                           " (curl " + str(completed.returncode) + "): " +
                           completed.stderr[-500:])
    try:
        asset = json.loads(completed.stdout)
    except json.JSONDecodeError as exc:
        raise RuntimeError("Release upload API did not return JSON: " + file.name) from exc
    require(asset.get("name") == file.name
            and asset.get("state") == "uploaded"
            and asset.get("digest") == expected_hash,
            "Upload receipt SHA not equal to frozen source bytes: " + file.name)
    return asset


def publish(root, *, repo=REPO, release_id=RELEASE_ID, token=None):
    token = token or os.environ.get("GH_TOKEN")
    require(bool(token), "No GitHub token for authenticated Release publication")
    files, hashes, report = expected_files(root)
    require(repo == REPO and release_id == RELEASE_ID,
            "Wrong repository or existing draft release identity")
    release = release_info(token, repo, release_id)
    assets = published_assets(token, repo, release_id)
    check_existing_assets(assets, hashes)
    if not release["draft"]:
        require(set(assets) == set(hashes),
                "Published source target Release is missing original assets")
        print("MOIS_93209351_FULL_TARGET_RELEASE_ALREADY_PUBLIC_AND_SHA_VERIFIED",
              len(assets), flush=True)
        return
    print("MOIS_SOURCE_ONLY_RELEASE_RESUME_FROM_PINNED_SUCCESSFUL_16_WORKER_ARTIFACTS",
          json.dumps({"verified_photo_targets": PHOTO_TOTAL,
                      "assets_expected": len(files), "assets_present": len(assets)}),
          flush=True)
    for i, file in enumerate(files):
        if file.name not in assets:
            upload_asset(token, repo, release_id, file, hashes[file.name])
            # Keep an explicit completion receipt for recoverability.
            if i % 16 == 0 or file.name == MANIFEST_NAME:
                print("MOIS_SIGNED_SOURCE_ONLY_RELEASE_ASSET_COMPLETED",
                      i + 1, "/", len(files), file.name, flush=True)
    final_assets = published_assets(token, repo, release_id)
    check_existing_assets(final_assets, hashes)
    require(set(final_assets) == set(hashes),
            "Full 93m source target Release is not complete; remain in draft")
    updated = api_request(token, "PATCH",
                          f"{API}/repos/{repo}/releases/{release_id}",
                          {"draft": False})
    require(updated["tag_name"] == RELEASE_TAG and updated["draft"] is False,
            "Signed source-only target Release was not published")
    print("MOIS_ALL_93209351_SOURCE_TARGET_TIERS_PUBLISHED_SHA_VERIFIED_NOT_RGB",
          json.dumps({"photos": PHOTO_TOTAL, "assets": len(final_assets),
                      "tiers": report["new_target_counts_by_tier"]},
                     sort_keys=True), flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--release-id", type=int, default=RELEASE_ID)
    args = p.parse_args()
    publish(args.source, release_id=args.release_id)


if __name__ == "__main__":
    main()
