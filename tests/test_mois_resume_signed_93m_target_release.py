"""Offline identity and resumability tests; never upload GitHub assets here."""
import json
from pathlib import Path

import pytest

from scripts.acquisition import mois_resume_signed_93m_target_release as r


def mock_source(tmp_path):
    hashes = {}
    for tier in range(14):
        for group in range(16):
            path = tmp_path / f"tier-{tier:02d}-group-{group:02d}.tar"
            path.write_bytes(f"{tier}-{group}".encode())
            hashes[path.name] = r.sha256_file(path)
    for group in range(16):
        (tmp_path / f"group-{group:02d}-receipt.json").write_text(
            json.dumps({"group": group}))
    manifest = {
        "schema": "mois_complete_full93m_original_photo_target_registry_v1",
        "exact_original_photo_IDs_in_all_disjoint_tiers": r.PHOTO_TOTAL,
        "native_vs_introduced_not_hard_filtered": True,
        "source_only_NO_photo_pixels_or_hue_measured": True,
        "new_target_counts_by_tier": {"0": r.PHOTO_TOTAL},
        "published_tier_archive_SHA256": hashes,
    }
    (tmp_path / r.MANIFEST_NAME).write_text(json.dumps(manifest))
    return manifest


def test_exact_241_signed_assets_are_required(tmp_path):
    mock_source(tmp_path)
    files, hashes, manifest = r.expected_files(tmp_path)
    assert len(files) == len(hashes) == 241
    assert files[-1].name == r.MANIFEST_NAME
    assert all(x.startswith("sha256:") for x in hashes.values())
    assert sum(manifest["new_target_counts_by_tier"].values()) == r.PHOTO_TOTAL


def test_modified_tier_archive_is_rejected_before_any_upload(tmp_path):
    mock_source(tmp_path)
    (tmp_path / "tier-02-group-10.tar").write_bytes(b"corrupt")
    with pytest.raises(ValueError, match="SHA changed"):
        r.expected_files(tmp_path)


def test_incomplete_worker_receipts_are_not_published(tmp_path):
    mock_source(tmp_path)
    (tmp_path / "group-15-receipt.json").unlink()
    with pytest.raises(ValueError, match="Missing original"):
        r.expected_files(tmp_path)


def test_source_photo_count_or_native_filter_mismatch_rejected(tmp_path):
    m = mock_source(tmp_path)
    m["native_vs_introduced_not_hard_filtered"] = False
    (tmp_path / r.MANIFEST_NAME).write_text(json.dumps(m))
    with pytest.raises(ValueError, match="source-authenticated"):
        r.expected_files(tmp_path)


def test_existing_assets_must_match_exact_source_sha():
    sha = "sha256:" + "a"*64
    r.check_existing_assets({"good.tar": {
        "name": "good.tar", "digest": sha, "state": "uploaded",
    }}, {"good.tar": sha})
    with pytest.raises(ValueError, match="mismatching"):
        r.check_existing_assets({"good.tar": {
            "name": "good.tar", "digest": "sha256:" + "b"*64, "state": "uploaded",
        }}, {"good.tar": sha})
    with pytest.raises(ValueError, match="Unexpected"):
        r.check_existing_assets({"unknown.tar": {
            "name": "unknown.tar", "digest": sha, "state": "uploaded",
        }}, {"good.tar": sha})
