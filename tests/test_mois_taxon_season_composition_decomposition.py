"""Synthetic taxonomy-mix counterfactual guards; NEVER observed ecology."""
from __future__ import annotations
import numpy as np
import pandas as pd
import pytest

from scripts.analysis.mois_taxon_season_composition_decomposition import (
    bootstrap_taxa, summarize_zone, symmetric_photo_decomposition,
)


def test_symmetric_taxon_mix_exactly_accounts_for_fake_photo_summer_difference():
    # Summer photos overwhelmingly capture the white taxon, winter photos the
    # chromatic taxon. No taxon changes its measured white fraction whatsoever.
    a=np.array([[80, 20, 80, 20], [20, 80, 0, 0]],dtype=float)
    within, composition, total, taxon_equal=symmetric_photo_decomposition(a)
    assert within == pytest.approx(0)
    assert composition == pytest.approx(60)
    assert total == pytest.approx(60)
    assert taxon_equal == pytest.approx(0)


def test_identical_source_taxon_photo_weights_cannot_create_composition_effect():
    a=np.array([[20, 20, 10, 0], [30, 30, 15, 30]],dtype=float)
    within, composition, total, equal=symmetric_photo_decomposition(a)
    assert composition==pytest.approx(0)
    assert within==pytest.approx(total)
    assert within != pytest.approx(0)


def test_whole_taxon_bootstrap_returns_finite_additive_components():
    a=np.array([[5, 4, 1, 0], [4, 8, 2, 7], [6, 6, 3, 3]],dtype=float)
    point,interval=bootstrap_taxa(a,n_draws=100,seed=33)
    assert np.all(np.isfinite(interval))
    assert point[0]+point[1]==pytest.approx(point[2])
    assert interval.shape==(4,2)
    assert np.all(interval[:,0]<=interval[:,1])


def test_reject_missing_one_season_or_duplicate_original_observation_id():
    with pytest.raises(ValueError,match="Invalid photographic"):
        symmetric_photo_decomposition(np.array([[10, 0, 3, 0]]))
    q=[]
    for tax in ("a", "b"):
        for season in ("summer", "winter"):
            for n in (1,2):
                q.append({"inat_taxon_id":tax,"observation_id":f"{tax}{season}{n}",
                          "season":season,"white":season=="summer"})
    tab=pd.DataFrame(q)
    results=summarize_zone(tab,min_per_season=2)
    assert results["n_common_taxa"]==2
    tab.loc[1,"observation_id"]=tab.loc[0,"observation_id"]
    with pytest.raises(ValueError,match="Duplicate or missing original observation"):
        summarize_zone(tab,min_per_season=2)
