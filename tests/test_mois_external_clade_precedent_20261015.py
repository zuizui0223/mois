"""Guard Mois against claiming unperformed independent replication or first global aridity trend.

All evidence on the external 2026 precedent is publisher literature, NOT an
independent downloaded Dryad/Phaidra input. The separate script --self-test
exercises artificial records only.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MS = ROOT / "docs/MOIS_GEB_RESEARCH_ARTICLE_DRAFT_20261014.md"
NOTE = ROOT / "docs/MOIS_EXTERNAL_CLADE_HUE_AND_PRECEDENT_AUDIT_20261015.md"
REFERENCE = ROOT / "docs/MOIS_GEB_REFERENCE_AND_GUIDELINE_AUDIT_20261014.md"
SCRIPT = ROOT / "scripts/analysis/mois_external_clade_flower_colour_benchmark_20261015.py"


def test_global_precedent_citation_and_distinct_estimands():
    article = MS.read_text(encoding="utf-8")
    assert "10.1002/ajb2.70044" in article
    assert "Dellinger et al. (2026)" in article
    assert "2,815 species across 51 clades" in article
    assert "approximately 82% to 50–55%" in article
    assert "not a direct replication" in article
    assert "not a universal genus response" in article
    assert "0.21%" in article


def test_external_color_source_is_protocol_not_fictitious_result():
    note = NOTE.read_text(encoding="utf-8")
    refs = REFERENCE.read_text(encoding="utf-8")
    src = SCRIPT.read_text(encoding="utf-8")
    for token in ("10.1002/ajb2.70044", "10.5061/dryad.r4xgxd2sc", "2,960", "2,815", "51"):
        assert token in note
    assert "NOT YET RUN" in note
    assert "no external CSV bytes" in note
    assert "NO downloaded dataset" not in note  # this is a lab note, not synthetic CI output
    assert "10.1002/ajb2.70044" in refs
    assert "Source download redirects" in refs
    assert "flower_color" in src
    assert "--self-test" in src
    assert "source_row_count" in src
    assert "clade_mean_binary_white_variance_fraction" in src
    assert "Species appears in more than one source row" in src
