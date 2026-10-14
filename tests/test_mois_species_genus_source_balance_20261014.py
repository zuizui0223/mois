"""Source-only species/genus hue estimator fixtures and committed aggregate scientific gates.

CI does not include original photo source ZIP; full raw-data replay was separately
executed with SHA verification and may be rechecked by authorized owners only.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/analysis/mois_replay_species_genus_balance_20261014.py"
RESULT = ROOT / "results/mois_original_photo_species_genus_balance_20261014.json"


def module():
    spec = importlib.util.spec_from_file_location("species_genus_balance", SCRIPT)
    loaded = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(loaded)
    return loaded


def test_synthetic_levels_do_not_collapse_into_photo_share():
    m = module()
    example = pd.DataFrame([
        {"inat_taxon_id":1,"genus":"A","morph":"white","classified":True},
        {"inat_taxon_id":1,"genus":"A","morph":"white","classified":True},
        {"inat_taxon_id":1,"genus":"A","morph":"white","classified":True},
        {"inat_taxon_id":2,"genus":"A","morph":"red_pink","classified":True},
        {"inat_taxon_id":3,"genus":"B","morph":"blue_purple","classified":True},
        {"inat_taxon_id":3,"genus":"B","morph":"mixed_uncertain","classified":False},
    ])
    species, genera = m.summarized_weights(example)
    assert len(species) == 3 and len(genera) == 2
    assert example.loc[example.classified].morph.eq("white").mean() == 3/5
    assert np.isclose(species.white.mean(), 1/3)
    assert np.isclose(genera.mean_species_white.mean(), 1/4)
    # Same genus redraw weights all three estimators in a paired resample.
    out = m.cluster_bootstrap(genera, np.array([[0, 1], [1, 0]], dtype=np.int32))
    assert out.shape == (2, 3)
    np.testing.assert_allclose(out[:,0], [3/5, 3/5])
    np.testing.assert_allclose(out[:,1], [1/3, 1/3])
    np.testing.assert_allclose(out[:,2], [1/4, 1/4])


def test_source_aggregate_denominators_and_weighting_support():
    report = json.loads(RESULT.read_text(encoding="utf-8"))
    assert report["schema"] == "mois_original_photo_species_genus_balance_v1"
    assert report["n_classified_white_photos"] == 21574
    assert report["n_classified_coloured_photos"] == 24379
    assert report["n_classified_photos"] == 45953
    assert report["n_nominal_species_total"] == 42111
    assert report["n_species_with_classified_photo"] == 25064
    assert report["n_species_with_no_classified_photo"] == 17047
    assert report["n_genera_with_classified_photo"] == 5369
    assert report["n_genera_one_classified_species"] == 2529
    assert sum(x["n_genera"] for x in report["grouped_by_sampled_species_per_genus"]) == 5369
    assert sum(x["n_species"] for x in report["grouped_by_sampled_species_per_genus"]) == 25064
    assert report["n_classified_white_photos"]/report["n_classified_photos"] == report["point_white_proportion"]["photo_equal"]
    values=report["point_white_proportion"]
    assert values["genus_equal"] > values["species_equal"]
    assert np.isclose(values["genus_equal"]-values["species_equal"],report["genus_equal_minus_species_equal_white_proportion"])
    assert report["paired_genus_minus_species_bootstrap_95"][0] > 0
    assert report["min_classified_species_thresholds"][3]["min_classified_species_per_genus"] == 5
    assert np.isclose(report["min_classified_species_thresholds"][3]["genus_equal_white"],0.46625925585926914)
    assert report["bootstrap_draws"] == 1999
    assert report["bootstrap_seed"] == 2026101450
    assert np.isclose(report["hypothetical_all_42111_taxon_white_mean_bounds"][0],0.2786614250380214)
    assert np.isclose(report["hypothetical_all_42111_taxon_white_mean_bounds"][1],0.6834725195263973)


def test_manuscript_keeps_new_result_descriptive_not_fitness():
    ms=(ROOT/"docs/MOIS_GEB_RESEARCH_ARTICLE_DRAFT_20261014.md").read_text()
    si=(ROOT/"docs/MOIS_GEB_SUPPORTING_INFORMATION_20261014.md").read_text()
    memo=(ROOT/"docs/MOIS_SOURCE_REPLAYED_SPECIES_GENUS_WHITE_BALANCE_20261014.md").read_text()
    assert "species-equal and genus-equal white values were **46.82%**" in ms
    assert "### S2.4 | Original source:" in si
    for token in ("25,064","5,369","49.891%","46.626%","+3.072","17,047"):
        assert token in si and token in memo
    assert "NOT" in memo and "relative fitness" in memo
    assert "not evolutionary equilibria" in ms
