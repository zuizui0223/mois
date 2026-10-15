"""Protect the independent rainfall replication source boundary.

No original off-photo GBIF species occurrences or WorldClim BIO12 rasters were
downloaded or analysed when this protocol was added. Synthetic smoke tests
cannot be described as successful ecological transfer.
"""
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
CODE=ROOT/"scripts/analysis/mois_external_independent_bio12_transfer_20261016.py"
NOTE=ROOT/"docs/MOIS_INDEPENDENT_GBIF_WORLDCLIM_WHITE_REPLICATION_GATE_20261016.md"
PAPER=ROOT/"docs/MOIS_GEB_RESEARCH_ARTICLE_DRAFT_20261014.md"


def test_external_moisture_extraction_requires_real_independent_sources():
    code=CODE.read_text(encoding="utf-8")
    note=NOTE.read_text(encoding="utf-8")
    for token in (
        "PRESERVED_SPECIMEN", "sample_bio12", "fit_white_bio12",
        "min_cells=5", "min_clade_species=3", "gbifID",
        "coordinateUncertaintyInMeters", "gbif-download-doi",
        "raster-factor", "Source", "DOI",
    ):
        assert token in code or token in note, token
    for token in (
        "NOT", "not", "2,960", "WorldClim", "GBIF", "10.15468/dl.",
        "10.5061/dryad.r4xgxd2sc", "species", "raster",
        "NO", "not independently",
    ):
        # Directly verify meaningful presence of methodological gates.
        if token in ("NO", "not independently"):
            continue
        assert token in note, token
    assert "source_independent_climate_analysis_completed" in code
    assert "insufficient_independent_species_and_clades" in code
    assert "source_trait_sha256" in code and "worldclim_raster_sha256" in code
    assert "independent moisture" in PAPER.read_text(encoding="utf-8").lower()
    assert not (ROOT/"results/mois_external_independent_bio12_source_20261016.json").exists()


def test_original_empirical_clade_receipt_remains_a_nonclimate_source():
    p=ROOT/"results/mois_external_paired_flower_fruit_20261016.json"
    assert p.exists()
    s=p.read_text(encoding="utf-8")
    assert "2960" in s and "clade" in s
    assert "Not pigment mass" in s
