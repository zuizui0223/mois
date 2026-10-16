"""Recompute original-photo wet/dry missingness bounds from the frozen receipt."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/analysis/mois_photo_nonresponse_bounds_20261014.py"
RESULT = ROOT / "results/mois_photo_nonresponse_bounds_20261014.json"


def test_sharp_bounds_reproduce_source_aggregates_and_scientific_limits():
    spec = importlib.util.spec_from_file_location("photo_bounds", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    receipt = module.build_receipt()
    assert receipt == json.loads(RESULT.read_text(encoding="utf-8"))
    assert receipt["n_georeferenced_rain_photo_opportunities"] == 100255
    assert receipt["dry_decile"]["n_classified_white"] == 2082
    assert receipt["wet_decile"]["n_classified_white"] == 2079
    assert receipt["dry_decile"]["n_unclassified"] == 5103
    assert receipt["wet_decile"]["n_unclassified"] == 6086
    assert receipt["hypothetical_all_photo_wet_minus_dry_sharp_bounds_pp"] == [-50.440677, 61.085163]
    assert receipt["observed_classified_only_wet_minus_dry_pp"] == 11.454238
    assert receipt["original_mnar_tipping_wet_unknown_white_if_dry_unknown_white_is_0_5"] == 0.412637
    assert "not bootstrap confidence intervals" in receipt["scope_restriction"]


def test_manuscript_bounds_are_explicitly_raw_not_adjusted():
    manuscript = (ROOT / "docs/MOIS_GEB_RESEARCH_ARTICLE_DRAFT_20261014.md").read_text(encoding="utf-8")
    si = (ROOT / "docs/MOIS_GEB_SUPPORTING_INFORMATION_20261014.md").read_text(encoding="utf-8")
    assert "−50.44 to +61.09 percentage points" in manuscript
    assert "do not apply to the adjusted congeneric slope" in manuscript
    assert "−50.4407 to +61.0852 percentage points" in si
    assert "not the congeneric genus-adjusted BIO12 slope" in si
    assert "biologically meaningful latent binary hue" in si
