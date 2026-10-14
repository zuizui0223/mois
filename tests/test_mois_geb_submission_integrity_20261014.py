"""Journal and source-science guards for a still-unsubmitted GEB manuscript.

No new photographic hue labels, environmental models, original bootstrap draws
or statistical significance tests are performed by this manuscript audit.
"""
from __future__ import annotations
import csv
import math
import re
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
MS=ROOT/"docs/MOIS_GEB_RESEARCH_ARTICLE_DRAFT_20261014.md"
SI=ROOT/"docs/MOIS_GEB_SUPPORTING_INFORMATION_20261014.md"
BLIND=ROOT/"docs/MOIS_GEB_BLINDED_REVIEW_COPY_20261014.md"
NUMBERED=ROOT/"docs/MOIS_GEB_LINE_NUMBERED_REVIEW_COPY_20261014.txt"
REFERENCE=ROOT/"docs/MOIS_GEB_REFERENCE_AND_GUIDELINE_AUDIT_20261014.md"
TITLE=ROOT/"docs/MOIS_GEB_IDENTIFYING_TITLE_PAGE_TEMPLATE_20261014.md"
COVER=ROOT/"docs/MOIS_GEB_COVER_LETTER_FOR_PDF_20261014.md"
CSV=ROOT/"results/fcp_macro_direction_20261013/direction_ledger.csv"


def text(path):
    return path.read_text(encoding="utf8")


def test_four_hue_coefficients_and_intervals_retain_orig_source_family():
    with CSV.open(newline="",encoding="utf8") as f:
        rows=list(csv.DictReader(f))
    assert len(rows)==1280
    for axis in ("wc_bio12","wc_vapr_annual_mean"):
        matching=[r for r in rows if all([
            r["mode"]=="climate",r["quality"]=="all",r["region"]=="world",
            r["adjustment"]=="genus_geo",r["weight"]=="genus_equal",r["anchor"]==axis])]
        assert len(matching)==4
        assert set(r["hue"] for r in matching)=={"white","yellow_orange","red_pink","blue_purple"}
        assert math.isclose(sum(float(r["direction_pp_per_env_global_sd"]) for r in matching),0,abs_tol=1e-11)
        for r in matching:
            assert (int(r["n_species"]),int(r["n_genera"]))==(20694,1934)
            lo,hi=(float(r[z]) for z in ("maxT_ci_lo_pp","maxT_ci_hi_pp"))
            pt=float(r["direction_pp_per_env_global_sd"])
            assert lo<=pt<=hi
            assert (lo>0 or hi<0) == (r["maxT_excludes_zero"]=="True")
            assert (r["maxT_excludes_zero"]=="True") == (r["hue"]=="white")
    s=text(MS);supp=text(SI)
    assert s.count("**Figure 4.")==1
    assert "no single nonwhite hue" not in s or "each" in s
    for token in ("−1.221","−0.446","−0.588","−1.385","−1.509","−1.057"):
        assert token in s,token
    assert "four-colour" in supp.lower()
    assert (ROOT/"figures/mois_original_four_hue_moisture.svg").exists()


def test_research_article_journal_submission_envelope_stays_bounded():
    s=text(MS)
    assert len(s.split("**Running title:** ",1)[1].splitlines()[0])<40
    abst=s[s.index("## Abstract"):s.index("## 1 | Introduction")]
    assert 1<=len(abst.split())<=300
    assert [abst.index("**"+x+":**") for x in ("Aim","Location","Time period",
        "Major taxa studied","Methods","Results","Main conclusions")]==sorted(
        abst.index("**"+x+":**") for x in ("Aim","Location","Time period",
        "Major taxa studied","Methods","Results","Main conclusions"))
    main=s[s.index("## 1 | Introduction"):s.index("## Data and code availability")]
    assert 3000<=len(main.split())<=5000
    wordline=re.search(r"^\*\*Keywords:\*\* (.+)$",s,re.M).group(1)
    words=[x.strip() for x in wordline.split(",")]
    assert 6<=len(words)<=10 and words==sorted(words,key=str.casefold)
    assert [str(i) for i in range(1,5)]==[m.group(1) for m in re.finditer(r"^\*\*Figure (\d+)\.",s,re.M)]
    assert [str(i) for i in range(1,3)]==[m.group(1) for m in re.finditer(r"^\*\*Table (\d+)\.",s,re.M)]
    assert "not submitted by this workflow" in s


def test_review_copy_anonymizes_public_repository_and_keeps_blocker_visible():
    a=text(MS);b=text(BLIND);line=text(NUMBERED)
    assert a!=b
    assert "# The geographic and measurement dependence" in b
    assert "## 1 | Introduction" in b and "## Main tables" in b
    assert b.count("## Data and code availability")==1
    assert "NOT YET PROVIDED" in b and "reviewer" in b.lower()
    assert not re.search(r"zuizui0223|github\.com|(?:^|/)fcp/|mois_geb|fcp_macro",b,re.I)
    assert not re.search(r"\b[0-9a-f]{64}\b",b,re.I)
    assert "## Acknowledgements, author contributions and declarations" not in b
    assert "## References" in b and b.index("## References")<b.index("## Data and code availability")
    rows=line.splitlines()
    assert len(rows)>180
    assert all(row.startswith(str(i+1).zfill(4)+"  ") for i,row in enumerate(rows))
    assert "AUTHOR(S) TO CONFIRM" in text(TITLE)
    assert "ONE CORRESPONDING AUTHOR" in text(TITLE).upper()


def test_reference_DOIs_are_verified_and_nothing_claims_a_source_archive_DOI():
    s=text(MS);g=text(REFERENCE)
    references=["10.1111/nph.16737","10.1093/biosci/biab093","10.1002/joc.5086",
        "10.3389/fpls.2020.00991","10.5194/soil-7-217-2021",
        "10.1111/jbi.14751","10.1109/CVPR52733.2024.01525"]
    for doi in references:
        assert doi in s and doi in g,doi
    assert "NOT YET PROVIDED" in text(BLIND)
    assert "GitHub" in g and "stable" in g.lower()
    assert "not present in the new mois repository" in text(SI)


def test_cover_letter_interest_and_submission_require_human_approval():
    cover=text(COVER)
    interest=cover[cover.index("Our study addresses"):cover.index("\n\nThe closely related")]
    assert len(interest.split())<250
    assert "pigment" in interest.lower()
    assert "not submitted" in cover.lower()
    assert "HUMAN" in cover and "TO CONFIRM" in cover
    audit=text(REFERENCE)
    assert "HOLD" in audit
    assert "PRIVATE" in audit.upper()
    assert "source archive" in audit.lower()
