"""Guard source-pixel cohort reconciliation from accidental prevalence pooling."""
import pandas as pd
import pytest

from scripts.analysis import mois_reconcile_real_disjoint_fourth512_targeted256 as m


def table(ids, labels, arm="baseline"):
    d = pd.DataFrame({
        "photo_id": list(ids),
        "morph": list(labels),
        "qualified_date": ["2025-01-15"] * len(ids),
        "old_mois_or_new_fcp": [arm] * len(ids),
    })
    return d


def test_frozen_denominators_not_source_only_3point43million():
    assert (m.BASE_N, m.BASE_HUES, m.BASE_DATED) == (309439, 149662, 149653)
    assert (m.FOURTH_N, m.FOURTH_HUES) == (512, 203)
    assert (m.TARGET_N, m.TARGET_HUES) == (256, 71)
    assert (m.TOTAL_N, m.TOTAL_HUES, m.TOTAL_DATED) == (310207, 149936, 149927)


def test_disjoint_merge_preserves_original_pixel_hues_and_dates():
    old = table(["1", "2"], ["white", "mixed_uncertain"])
    fourth = table(["3", "4"], ["white", "blue_purple"], "fourth")
    targeted = table(["5"], ["red_pink"], "targeted")
    merged = m.join_disjoint(old, fourth, targeted)
    assert merged.photo_id.tolist() == ["1", "2", "3", "4", "5"]
    assert merged.iloc[:2].reset_index(drop=True).equals(old)
    assert m.counted(merged) == {
        "measured_source_photo_IDs": 5,
        "automatic_four_hue_ROIs": 4,
        "dated_and_automatic_four_hue_ROIs": 4,
        "white_hue_ROIs": 2,
        "technical_hue_unknown_ROIs": 1,
    }


def test_source_overlap_and_modified_historical_colour_fail_fast():
    old = table(["1", "2"], ["white", "red_pink"])
    good = table(["3", "4"], ["white", "blue_purple"])
    with pytest.raises(ValueError, match="overlap"):
        m.join_disjoint(old, good, table(["3"], ["white"]))
    with pytest.raises(ValueError, match="baseline"):
        m.join_disjoint(old, table(["1"], ["white"]), table(["8"], ["white"]))
    tampered = pd.concat((old.copy(), good), ignore_index=True)
    tampered.loc[0, "morph"] = "yellow_orange"
    with pytest.raises(ValueError, match="historical"):
        m.original_prefix_and_new_photos(old, tampered, 2, 2, "test")


def test_no_unknown_assigned_nonwhite_colour_and_unmeasured_no_entry():
    old = table(["1", "2"], ["white", "mixed_uncertain"])
    extra = table(["3"], ["mixed_uncertain"])
    assert m.counted(extra)["automatic_four_hue_ROIs"] == 0
    with pytest.raises(ValueError, match="classified ROIs"):
        m.original_prefix_and_new_photos(old, pd.concat(
            (old, extra), ignore_index=True), 1, 1, "unknown-only")


def test_matched_source_support_separates_one_vs_two_real_observations():
    n1 = "same_original_taxon_025degree_cell_same_seasonyear_both_summer_winter_min1_obs"
    n2 = "same_original_taxon_025degree_cell_same_seasonyear_both_summer_winter_min2_obs"
    before = {
        "north_extratropical": {n1: 40, n2: 2},
        "south_extratropical": {n1: 34, n2: 0},
    }
    after = {
        "north_extratropical": {n1: 87, n2: 2},
        "south_extratropical": {n1: 54, n2: 1},
    }
    rows = m.support_rows(before, after)
    assert sum(row["gain_min1_source_groups"] for row in rows) == 67
    assert sum(row["gain_min2_source_groups"] for row in rows) == 1
    after["south_extratropical"][n2] = -1
    with pytest.raises(ValueError, match="cannot erase"):
        m.support_rows(before, after)
