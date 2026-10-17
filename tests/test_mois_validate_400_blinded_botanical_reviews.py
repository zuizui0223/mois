"""Synthetic, no-actual-botanical-votes tests for frozen MOIS review adjudication."""
from __future__ import annotations

import csv
from pathlib import Path

import pytest

from scripts.analysis.mois_validate_400_blinded_botanical_reviews import (
    COLUMNS, MAP_COLUMNS, THIRD_COLUMNS, audit,
)


def csvwrite(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def csvread(path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def sources(tmp: Path):
    x = []
    for idx, photo in enumerate(["101", "202"]):
        x.append({
            "review_token": f"token{idx}", "species": "Testus flowers",
            "original_photo_url": f"https://example.invalid/photo/{photo}",
            "photo_license": "cc-by", "attribution": "Synthetic test photographer",
            "reviewer_id": "", "open_flower_stage": "", "roi_organ": "",
            "notes": "",
        })
    maps = []
    for idx, v in enumerate(x):
        maps.append({
            "inat_taxon_id": "999", "species": v["species"],
            "photo_id": str(101 * (idx + 1)), "observation_id": f"uuid-{idx}",
            "observer_id": f"photo-source-{idx}",
            "observed_on": "2025-06-01", "photo_license": "cc-by",
            "photo_url_large": v["original_photo_url"],
            "source_rank": f"hash-{idx}",
            "site_latcell": "542", "site_loncell": "425",
            "season_year": "2025", "selected_season": "summer",
            "review_token": v["review_token"],
        })
    ba, bb, mapping, ra, rb = (tmp / n for n in
            ["blank_a.csv", "blank_b.csv", "mapping.csv", "a.csv", "b.csv"])
    for file, fields, rows in (
        (ba, COLUMNS, x), (bb, COLUMNS, list(reversed(x))),
        (mapping, MAP_COLUMNS, maps), (ra, COLUMNS, x),
        (rb, COLUMNS, list(reversed(x))),
    ):
        csvwrite(file, fields, rows)
    return ba, bb, mapping, ra, rb


def vote(rows, token, reviewer, stage="open_flower", tissue="petal_tepal"):
    target = next(r for r in rows if r["review_token"] == token)
    target.update(reviewer_id=reviewer, open_flower_stage=stage,
                  roi_organ=tissue)


def test_all_blank_forms_do_not_produce_biological_claim(tmp_path: Path):
    a, b, m, ar, br = sources(tmp_path)
    r = audit(a, b, m, ar, br, tmp_path / "out", fixture=True)
    assert r["original_source_400_photo_ids"] == 2
    assert r["status"] == "BOTANICAL_REVIEW_INCOMPLETE_NO_WHITE_FLOWER_RATE"
    assert r["status_counts"] == {"both_unreviewed": 2}
    assert r["reviewer_A_complete_votes"] == r["reviewer_B_complete_votes"] == 0
    assert r["confirmed_open_petaloid_photos_for_FUTURE_original_photo_ROI"] == 0
    assert r["white_frequency_or_monthly_climate_effect_NOT_estimated"] is True
    assert len(csvread(tmp_path / "out" / "THIRD_REVIEW_PENDING_BLINDED_ORIGINAL_PHOTOS.csv")) == 0


def test_votes_align_and_conflicts_reach_independent_third_queue(tmp_path: Path):
    a, b, m, ar, br = sources(tmp_path)
    aa, bb = csvread(ar), csvread(br)
    vote(aa, "token0", "reviewer_A")
    vote(bb, "token0", "reviewer_B")
    vote(aa, "token1", "reviewer_A", stage="fruit_only", tissue="fruit_seed")
    vote(bb, "token1", "reviewer_B", stage="open_flower", tissue="petal_tepal")
    csvwrite(ar, COLUMNS, aa)
    csvwrite(br, COLUMNS, bb)
    r = audit(a, b, m, ar, br, tmp_path / "open", fixture=True)
    assert r["photos_with_two_independent_completed_votes"] == 2
    assert r["confirmed_open_petaloid_photos_for_FUTURE_original_photo_ROI"] == 1
    assert r["disagreement_only_original_photo_count"] == 1
    third_queue = tmp_path / "open" / "THIRD_REVIEW_PENDING_BLINDED_ORIGINAL_PHOTOS.csv"
    third = csvread(third_queue)
    assert len(third) == 1 and third[0]["review_token"] == "token1"
    assert "observed_on" not in third[0] and "selected_season" not in third[0]
    third[0].update(adjudicator_id="reviewer_C", adjudicated_stage="fruit_only",
                    adjudicated_roi_tissue="fruit_seed")
    filled_third = tmp_path / "adjudication.csv"
    csvwrite(filled_third, THIRD_COLUMNS, third)
    resolved = audit(a, b, m, ar, br, tmp_path / "resolved",
                     adjudication=filled_third, fixture=True)
    assert resolved["unresolved_photo_count"] == 0
    assert resolved["confirmed_open_petaloid_photos_for_FUTURE_original_photo_ROI"] == 1
    assert resolved["status"] == "ALL_SOURCE_PHOTOS_BOTANICAL_REVIEW_RESOLVED_NO_COLOUR_MEASURED"
    third[0]["adjudicator_id"] = "reviewer_B"
    csvwrite(filled_third, THIRD_COLUMNS, third)
    with pytest.raises(ValueError, match="Third adjudicator"):
        audit(a, b, m, ar, br, tmp_path / "invalidthird",
              adjudication=filled_third, fixture=True)


def test_identity_tamper_and_reuse_are_rejected(tmp_path: Path):
    a, b, m, ar, br = sources(tmp_path)
    arr, brr = csvread(ar), csvread(br)
    vote(arr, "token0", "same_reviewer")
    vote(brr, "token0", "same_reviewer")
    csvwrite(ar, COLUMNS, arr)
    csvwrite(br, COLUMNS, brr)
    with pytest.raises(ValueError, match="same botanist"):
        audit(a, b, m, ar, br, tmp_path / "reuse", fixture=True)
    vote(brr, "token0", "independent_B")
    brr[0]["original_photo_url"] = "https://attacker.example.org/not-original"
    csvwrite(br, COLUMNS, brr)
    with pytest.raises(ValueError, match="immutable"):
        audit(a, b, m, ar, br, tmp_path / "tamper", fixture=True)


def test_partial_or_invalid_vote_not_inferred_as_nonflower(tmp_path: Path):
    a, b, m, ar, br = sources(tmp_path)
    aa = csvread(ar)
    target = next(r for r in aa if r["review_token"] == "token0")
    target.update(reviewer_id="A", open_flower_stage="open_flower", roi_organ="")
    csvwrite(ar, COLUMNS, aa)
    with pytest.raises(ValueError, match="Partial stage"):
        audit(a, b, m, ar, br, tmp_path / "partial", fixture=True)
    target["roi_organ"] = "not_a_petaloid_code"
    csvwrite(ar, COLUMNS, aa)
    with pytest.raises(ValueError, match="Unrecognised"):
        audit(a, b, m, ar, br, tmp_path / "illegal", fixture=True)


def test_reviewers_with_identical_uncertain_codes_are_not_eligible(tmp_path: Path):
    a, b, m, ar, br = sources(tmp_path)
    aa, bb = csvread(ar), csvread(br)
    vote(aa, "token0", "A", stage="uncertain", tissue="uncertain")
    vote(bb, "token0", "B", stage="uncertain", tissue="uncertain")
    csvwrite(ar, COLUMNS, aa)
    csvwrite(br, COLUMNS, bb)
    r = audit(a, b, m, ar, br, tmp_path / "unresolved", fixture=True)
    assert r["status_counts"]["unresolved_flower_or_organ"] == 1
    assert r["unresolved_photo_count"] == 2
    assert r["confirmed_open_petaloid_photos_for_FUTURE_original_photo_ROI"] == 0
