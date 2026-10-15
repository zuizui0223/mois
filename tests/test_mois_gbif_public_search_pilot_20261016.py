"""Immutable, aggregate-only source audit for preliminary GBIF taxonomic joins.

Search counts are public indexed occurrence counts, NOT quality-filtered,
spatially independent specimen sites or a citable GBIF bulk DOI.
"""
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/"results/mois_independent_gbif_specimen_search_pilot_20261016.json"
SCRIPT=ROOT/"scripts/analysis/mois_independent_gbif_search_coverage_pilot_20261016.py"
EXTRACTOR=ROOT/"scripts/analysis/mois_external_independent_bio12_transfer_20261016.py"
NOTE=ROOT/"docs/MOIS_INDEPENDENT_GBIF_WORLDCLIM_WHITE_REPLICATION_GATE_20261016.md"


def test_public_gbif_search_receipt_has_exact_immutable_scope():
    q=json.loads(SRC.read_text(encoding="utf8"))
    assert q["source_traits_sha256"]=="a253308785e4cbd0e361b3ca04cdfdae29c523eefa843375cebf8030ee0874af"
    assert q["selection_ordered_binomial_names_sha256"]=="d025092677a9377265d48d1484fa8ab9555788a76b3c3b08c04fdb2729f6b0be"
    assert (q["n_clades_sampled"],q["n_species_sampled"])==(15,30)
    assert len(q["anonymous_pairs"])==15
    assert q["totals"]["nonzero_named_searches"]==27
    assert q["totals"]["searches_with_at_least_5_raw_indexed_occurrences"]==23
    assert q["totals"]["searches_with_zero"]==3
    assert q["totals"]["first_result_exact_binomial_matching_source"]==19
    assert q["totals"]["first_result_accepted_species_mismatch"]==8
    assert q["totals"]["paired_clades_both_names_raw_indexed_count_at_least_5"]==8
    assert q["totals"]["paired_clades_both_first_accepted_species_exact"]==8
    assert q["totals"]["paired_clades_both_first_exact_and_raw_count_at_least_5"]==6
    assert q["totals"]["n_hits_over_all_30_search_queries"]==1634
    a=[t for pair in q["anonymous_pairs"] for t in (pair["white"],pair["nonwhite"])]
    assert len(a)==30 and sum(t["indexed_specimens"] for t in a)==1634
    assert sum(t["first_result_species_exact"] is True for t in a)==19
    assert sum(t["first_result_species_exact"] is False for t in a)==8
    assert sum(t["first_result_species_exact"] is None for t in a)==3
    assert sum(t["indexed_specimens"]>=5 for t in a)==23
    assert all("source species" not in repr(t).lower() for t in a)
    assert "NOT an authenticated" in q["search_basis"]
    assert any("WorldClim BIO12" in note for note in q["special_cautions"])


def test_strict_binomial_source_and_taxonomic_mismatch_gates():
    py=SCRIPT.read_text(encoding="utf8")
    cl=EXTRACTOR.read_text(encoding="utf8")
    note=NOTE.read_text(encoding="utf8")
    for literal in ("PILOT_NAMES_SHA256","source_clade_sizes_in_sample_order","fetch_one","def self_test","basisOfRecord","limit"):
        assert literal in py,literal
    for literal in ("is_hybrid_source_name","n_excluded_unresolved_hybrid_names",
                    "accepted_species_mismatch_requires_crosswalk","SOURCE_EXPECTED_HYBRID_NAMES",
                    "Unrecognised source flower colour"):
        assert literal in cl,literal
    assert "Rosa_X_odorata" in note and "Garcinia_X_mangostana" in note
    assert "2,958 strict nonhybrid" in note
    assert "not a white–rainfall slope" in note.lower()
    assert not (ROOT/"results/mois_external_independent_bio12_source_20261016.json").exists()
