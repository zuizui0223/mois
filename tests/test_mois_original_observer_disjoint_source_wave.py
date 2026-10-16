"""Synthetic-only regression checks for photo-opportunity observer audit (no petals)."""
import pytest

from scripts.analysis.mois_observer_disjoint_photo_wave_20261017 import (
    balanced_observer_capacity, rank_photo, season,
)


def test_observer_independence_cannot_be_replaced_by_observation_counts():
    # 100 observations by just 2 photographers cannot support 40 independent
    # observers per original photo-season, no matter how many photos exist.
    assert balanced_observer_capacity({"u1", "u2"}, {"u3", "u4"}) == 2
    assert balanced_observer_capacity({"u1", "u2"}, {"u1", "u2"}) == 1
    assert balanced_observer_capacity(set(), {"u2"}) == 0


def test_cross_season_observer_overlap_must_be_allocated_only_once():
    a = {"a", "b", "overlap1", "overlap2"}
    b = {"c", "d", "overlap1", "overlap2"}
    assert balanced_observer_capacity(a, b) == 3
    assert balanced_observer_capacity({"a", "b", "c"}, {"x", "y", "z"}) == 3


def test_local_seasons_swap_southern_hemisphere_and_december_is_winter_north():
    assert season(1, 40) == "winter"
    assert season(7, 40) == "summer"
    assert season(1, -40) == "summer"
    assert season(7, -40) == "winter"
    assert season(7, 0) == "tropical_no_four_season_label"


def test_sha_photo_rank_is_deterministic_and_colour_independent():
    assert rank_photo("12345") == rank_photo("12345")
    assert rank_photo("12345") != rank_photo("12346")
    assert len(rank_photo("12345")) == 64
