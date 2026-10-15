"""Keep mois focused on FCP-derived flower-colour × environment inference.

Cited publication titles and original source schemas can legitimately contain
'fruit'; that does not license new fruit phenotype models in this article.
"""
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
ARTICLE=ROOT/"docs/MOIS_GEB_RESEARCH_ARTICLE_DRAFT_20261014.md"
SUPPLEMENT=ROOT/"docs/MOIS_GEB_SUPPORTING_INFORMATION_20261014.md"
READINESS=ROOT/"docs/MOIS_GEB_READINESS_AND_COVER_LETTER_20261014.md"
README=ROOT/"README.md"


def test_fcp_derived_article_keeps_floral_colour_only():
    raw=ARTICLE.read_text(encoding="utf8")
    body=raw[raw.index("## 1 | Introduction"):raw.index("## Data and code availability")]
    assert "rainfall" in body.lower() and "flower" in body.lower()
    assert "2,960" in body and "27.48%" in body
    for forbidden in ("fruit", "果実", "S2.19", "cross-organ", "red/black"):
        assert forbidden not in body.lower(), forbidden
    assert len(body.split())<=5000
    abstract=raw[raw.index("## Abstract"):raw.index("## 1 | Introduction")]
    assert len(abstract.split())<=300


def test_supplement_and_working_scope_exclude_cross_organ_claims():
    si=SUPPLEMENT.read_text(encoding="utf8")
    readiness=READINESS.read_text(encoding="utf8")
    readme=README.read_text(encoding="utf8")
    assert "### S2.18 |" in si
    assert "### S2.19 |" not in si
    assert "fruit →" not in si and "red/black fruit" not in si
    assert "mois_external_paired_flower_fruit_20261016.json" not in si+readiness+readme
    assert "mois_external_cross_organ_clade_20261016.py" not in si+readiness+readme
    assert "S2.19" not in readiness
    assert "flower-colour-only" in readiness
    assert "independent" in readme.lower() and "BIO12" in readme


def test_unrelated_fruit_work_is_not_shipped_in_mois_review_branch():
    paths=[
        ROOT/"scripts/analysis/mois_external_cross_organ_clade_20261016.py",
        ROOT/"results/mois_external_paired_flower_fruit_20261016.json",
        ROOT/"tests/test_mois_external_paired_organs_20261016.py",
    ]
    assert all(not path.exists() for path in paths)
