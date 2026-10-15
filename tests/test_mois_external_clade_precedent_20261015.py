"""Guard Mois against claiming unperformed independent replication or first global aridity trend.

The observed external CSV was retrieved from a previously saved Library copy;
this CI has only its independent aggregate receipt, not the original 2,960 rows.
Synthetic code tests and numerical claim checks do not recreate the raw data.
"""
from pathlib import Path
import json
import math

ROOT = Path(__file__).resolve().parents[1]
MS = ROOT / "docs/MOIS_GEB_RESEARCH_ARTICLE_DRAFT_20261014.md"
NOTE = ROOT / "docs/MOIS_EXTERNAL_CLADE_HUE_AND_PRECEDENT_AUDIT_20261015.md"
REFERENCE = ROOT / "docs/MOIS_GEB_REFERENCE_AND_GUIDELINE_AUDIT_20261014.md"
SCRIPT = ROOT / "scripts/analysis/mois_external_clade_flower_colour_benchmark_20261015.py"
EXTERNAL = ROOT / "results/mois_external_clade_white_source_20261015.json"


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
    assert "Independent external white-colour clustering — executed" in note
    assert "27.4843%" in note and "1.6845%" in note
    assert "no independent climate–colour coefficient" in note
    assert "publisher-side manifest SHA" in note
    assert "10.1002/ajb2.70044" in refs
    assert "Source download redirects" in refs
    assert "flower_color" in src
    assert "--self-test" in src
    assert "source_row_count" in src
    assert "clade_mean_binary_white_variance_fraction" in src
    assert "Species appears in more than one source row" in src


def test_real_external_aggregates_are_consistent_but_not_a_climate_replication():
    q = json.loads(EXTERNAL.read_text(encoding="utf-8"))
    assert q["external_data_doi"] == "10.5061/dryad.r4xgxd2sc"
    assert q["external_csv_sha256"] == "a253308785e4cbd0e361b3ca04cdfdae29c523eefa843375cebf8030ee0874af"
    assert q["external_csv_size_bytes"] == 118867
    assert (q["source_row_count"], q["qualified_classifiable_species"], q["qualified_clades"]) == (2960, 2960, 51)
    assert q["source_rows_without_classifiable_white_category"] == 0
    assert q["n_white_species_in_all_classifiable"] == 1449
    assert sum(q["external_eight_hue_counts"].values()) == 2960
    assert len(q["external_eight_hue_counts"]) == 8
    assert set(q["original_external_label_values"]) == set(q["external_eight_hue_counts"])
    assert math.isclose(q["white_species_equal_all_classifiable"], 1449 / 2960, abs_tol=1e-14)
    assert math.isclose(q["white_equal_qualified_clade_means"], .4819389233555799, abs_tol=1e-12)
    theta = q["clade_mean_binary_white_variance_fraction"]
    assert math.isclose(theta["observed"], .27484321871615275, abs_tol=1e-12)
    assert math.isclose(theta["global_species_colour_shuffle_mean"], .01684536208634321, abs_tol=1e-12)
    assert theta["global_species_colour_shuffle_95"][1] < theta["observed"]
    assert theta["one_sided_upper_tail_permutation_p"] == .0001
    assert (q["permutations"], q["random_seed"]) == (9999, 202610152328)
    p3 = q["nearly_uniform_white_or_nonwhite_clade_share_ge3"]
    assert math.isclose(p3["observed"], 12 / 51, abs_tol=1e-12)
    assert math.isclose(q["nearly_uniform_white_or_nonwhite_clade_share_ge5"], 12 / 49, abs_tol=1e-12)
    assert q["nearly_uniform_white_or_nonwhite_clade_share_ge5_null"]["one_sided_upper_tail_permutation_p"] == .0001
    boot = q["clade_resampling_diagnostic"]
    assert boot["draws"] == 9999 and boot["seed"] == 202610152329
    assert boot["selected_clade_bootstrap95_white_species_equal"][0] < .5 < boot["selected_clade_bootstrap95_white_species_equal"][1]
    assert "no independent moisture gradient" in q["comparison_warning"]
    assert not any(x in EXTERNAL.read_text(encoding="utf-8") for x in ("Annona_reticulata", "Solanum_lycopersicum", "gps_coordinates"))
