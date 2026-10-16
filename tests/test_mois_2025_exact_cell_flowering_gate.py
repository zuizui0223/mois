"""Synthetic constraints for MOIS original month-specific site-climate design gate.

These examples are not measured flowering phenotypes or ecological results.
"""
import csv
import json
from pathlib import Path

import pytest

from scripts.analysis.mois_2025_exact_cell_climate_flower_stage_gate import (
    MONTH_TO_SEASON, VARS, average, run, seasonal_support, sha,
)


def test_calendar_is_original_2025_northern_temperate_months_only():
    assert [k for k, v in MONTH_TO_SEASON.items() if v == "summer"] == [
        "2025-06", "2025-07", "2025-08"]
    assert [k for k, v in MONTH_TO_SEASON.items() if v == "autumn"] == [
        "2025-09", "2025-10", "2025-11"]
    assert len(VARS) == 10


def test_photo_level_average_is_not_equal_weighted_by_occupied_cells():
    d = [{"value": "20"}, {"value": "20"}, {"value": "5"}]
    assert average(d, "value") == 15.0


def test_reject_fake_empty_or_bad_original_denominator():
    with pytest.raises(ValueError, match="Original 64/40"):
        seasonal_support([], [])
    with pytest.raises(ValueError, match="Original 64/40"):
        seasonal_support([{} for _ in range(64)], [{} for _ in range(39)])


def test_no_unknown_climate_receipt_is_accepted(tmp_path: Path):
    p = tmp_path / "photos.csv"
    k = tmp_path / "keys.csv"
    s = tmp_path / "source.json"
    c = tmp_path / "climate.csv"
    r = tmp_path / "receipt.json"
    p.write_text("photo_id\n1\n", encoding="utf-8")
    k.write_text("year_month\n2025-07\n", encoding="utf-8")
    c.write_text("year_month\n2025-07\n", encoding="utf-8")
    s.write_text(json.dumps({
        "status": "REAL_DATE_AND_ROI_ATTEMPTS_NO_BOTANICAL_WHITE_CLIMATE_EFFECT",
        "original_photo_attempts_sha256": sha(p),
        "original_month_climate_grid_csv_sha256": sha(k),
    }), encoding="utf-8")
    r.write_text(json.dumps({
        "status": "AUTHENTIC_2025_MONTHLY_CLIMATE_VALUES_NOT_WHITE_FLOWER_ASSOCIATION",
        "input_original_2025_photo_climate_request_sha256": sha(k),
        "input_original_2025_roi_source_receipt_sha256": "tampered",
        "decoded_month_grid_climate_sha256": sha(c),
        "genuine_original_monthly_photo_observations": 64,
        "provider_request_count": 126,
        "no_botanically_confirmed_open_flower_colour_data": True,
    }), encoding="utf-8")
    with pytest.raises(ValueError, match="TerraClimate receipt mismatch"):
        run(p, k, s, c, r, tmp_path / "out.json")
