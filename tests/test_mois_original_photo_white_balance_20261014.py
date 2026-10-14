"""Globally photographed hue near-parity is not evidence of evolutionary balancing selection."""
from __future__ import annotations
import importlib.util
import json
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]

def test_source_derived_white_balance_parity_threshold():
    p = ROOT / "scripts/analysis/mois_original_photo_white_balance_20261014.py"
    spec = importlib.util.spec_from_file_location("mois_original_white_balance", p)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    result = module.build()
    assert result == json.loads((ROOT / "results/mois_original_photo_white_balance_20261014.json").read_text())
    assert result["photos"]["white"] == 21574
    assert result["photos"]["coloured"] == 24379
    assert result["photos"]["unclassified"] == 54590
    assert result["photos"]["classified_coloured_to_white_relabels_to_exceed_half"] == 1403
    assert abs(result["hypothetical_all_photo_balance"]["fraction_white_among_unclassified_needed_for_exact_50_percent"] - 0.52569) < 0.0001
    assert result["latitudinal"]["within_species_photo_label_shuffle_one_sided_p_unusually_even"] == 0.214
    assert result["spatial_components"]["null_one_sided_p_negative_compensation"] == 0.997

def test_claim_gate_near_parity_is_not_fitness_or_species_balance():
    report = (ROOT / "results/mois_original_photo_white_balance_20261014.json").read_text(encoding="utf-8")
    result = json.loads(report)
    assert not result["conclusions"]["observed_classified_photos_exactly_half"]
    assert not result["conclusions"]["evidence_for_frequency_dependent_balancing_selection"]
    assert not result["conclusions"]["evidence_for_pigment_cost_or_fitness_crossing"]
    assert not result["conclusions"]["biological_species_half_ratio_estimable_from_aggregate"]
    author = (ROOT / "docs/MOIS_GEB_RESEARCH_ARTICLE_DRAFT_20261014.md").read_text(encoding="utf-8")
    si = (ROOT / "docs/MOIS_GEB_SUPPORTING_INFORMATION_20261014.md").read_text(encoding="utf-8")
    memo = (ROOT / "docs/MOIS_WHITE_CHROMATIC_BALANCE_AND_FITNESS_HYPOTHESES_20261014.md").read_text(encoding="utf-8")
    assert "52.57%" in author and "50:50" in author
    assert "52.57%" in si and "0.997" in si
    assert "relative fitness" in memo.lower() and "negative frequency-dependent" in memo.lower()
