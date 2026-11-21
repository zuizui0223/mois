"""Hand-derived invariants for the native cohort permutation analysis.

Each fixture describes actual species/cell memberships, not mock source files.
The tests catch incoherent species draws, changed chromatic membership, tail
ties, and an incorrect independent-samples treatment of paired contrasts.
"""
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest


MODULE = "scripts.analysis.mois_native_hue9_whitefixed8_georegion_null_20261121"
SCRIPT = Path(__file__).resolve().parents[1] / (
    "scripts/analysis/mois_native_hue9_whitefixed8_georegion_null_20261121.py"
)


def analysis():
    # Missing implementation is an assertion failure during the initial RED
    # run, rather than an import/collection error hiding the intended tests.
    assert SCRIPT.exists(), "The native geographic null implementation is missing"
    return __import__(MODULE, fromlist=["prepare_assemblage"])


def test_same_species_draw_is_shared_by_cells_and_overlapping_cohorts():
    """An independent cell/cohort shuffle breaks these exact draw identities."""
    m = analysis()
    labels = np.array([0, 1, 1], dtype=np.int8)
    cell = np.array([0, 0, 1, 1, 2, 2], dtype=np.int32)
    species = np.array([0, 2, 1, 2, 0, 2], dtype=np.int32)
    strict = m.prepare_assemblage("strict", cell, species, labels, 4, 9)
    moderate = m.prepare_assemblage(
        "moderate", np.r_[cell, 3], np.r_[species, 2], labels, 4, 9
    )
    global_source = m.prepare_assemblage(
        "global", np.r_[cell, 3], np.r_[species, 2], labels, 4, 9
    )
    for c in (strict, moderate):
        c.indices = {"a": np.array([0]), "b": np.array([1]),
                     "same_species_elsewhere": np.array([2])}
    result = m.simulate_model(
        global_source, {"strict": strict, "moderate": moderate},
        [np.array([0, 1])], permutations=31, seed=17, model="toy"
    )
    draws = result["draws"]
    log2 = 0.6931471805599453
    np.testing.assert_allclose(draws["strict"]["a"] + draws["strict"]["b"],
                               log2, rtol=0, atol=1e-14)
    np.testing.assert_array_equal(draws["strict"]["a"],
                                  draws["strict"]["same_species_elsewhere"])
    np.testing.assert_array_equal(draws["strict"]["a"], draws["moderate"]["a"])
    assert len(np.unique(draws["strict"]["a"])) == 2
    np.testing.assert_array_equal(strict.richness, [2, 2, 2, 0])


def test_white_fixed_chromatic_counts_use_original_nonwhite_pairs_and_zero_to_seven():
    """Changing membership or forgetting the 1..8 -> 0..7 shift fails here."""
    m = analysis()
    labels = np.array([0, 1, 8, 1], dtype=np.int8)
    c = m.prepare_assemblage(
        "native8", np.array([0, 0, 0, 1, 1, 2]),
        np.array([0, 1, 2, 1, 3, 0]), labels, 4, 8
    )
    np.testing.assert_array_equal(c.counts, [
        [1, 0, 0, 0, 0, 0, 0, 1],
        [2, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
    ])
    np.testing.assert_array_equal(c.richness, [2, 2, 0, 0])
    np.testing.assert_allclose(c.observed, [0.6931471805599453, 0, 0, 0],
                               rtol=0, atol=1e-14)
    changed = m.count_assemblage(c, np.array([0, 8, 1, 1], dtype=np.int8))
    np.testing.assert_array_equal(changed[1], [1, 0, 0, 0, 0, 0, 0, 1])
    np.testing.assert_array_equal(changed.sum(axis=1), [2, 2, 0, 0])
    with pytest.raises(ValueError, match="[Ww]hite"):
        m.count_assemblage(c, np.array([1, 8, 0, 1], dtype=np.int8))


def test_all_white_assemblage_has_zero_nonwhite_entropy_and_empty_cells():
    m = analysis()
    c = m.prepare_assemblage("white", np.array([0]), np.array([0]),
                             np.array([0], dtype=np.int8), 2, 8)
    np.testing.assert_array_equal(c.counts, np.zeros((2, 8), dtype=int))
    np.testing.assert_array_equal(c.observed, [0, 0])
    np.testing.assert_array_equal(c.richness, [0, 0])


@pytest.mark.parametrize("observed,want_lower,want_upper,want_two", [
    (2.0, 0.8, 0.8, 1.0),
    (0.0, 0.2, 1.0, 0.4),
    (4.0, 1.0, 0.2, 0.4),
])
def test_rank_tails_count_near_ties_and_never_return_zero(
    observed, want_lower, want_upper, want_two
):
    """Strict inequalities, no +1, or a centered main tail all fail this."""
    m = analysis()
    result = m.permutation_statistics(observed, np.array([1, 2, 2 + 5e-13, 3]))
    assert result["p_rank_lower"] == want_lower
    assert result["p_rank_upper"] == want_upper
    assert result["p_rank_two_sided"] == want_two
    assert result["null_mean_monte_carlo_se"] == pytest.approx(
        0.408248290463863, abs=1e-12
    )


def test_matched_mobility_has_descriptive_summaries_without_inferential_p():
    m = analysis()
    result = m.permutation_statistics(
        1.0, np.array([2.0, 3.0, 4.0]), descriptive_only=True
    )
    assert result["observed_minus_null"] == -2.0
    assert result["p_rank_lower"] is None
    assert result["p_rank_two_sided"] is None
    assert result["legacy_simulated_mean_centered_p_two_sided_diagnostic_only"] is None
    assert result["inference_status"] == "descriptive_only_outcome_selected_mobility"


def test_direct_contrast_preserves_covariance_of_synchronous_draws():
    """Independent variances give positive variance for this constant contrast."""
    m = analysis()
    result, draws = m.paired_contrast(
        4.0, 1.0, np.array([2.0, 4.0, 6.0]), np.array([1.0, 3.0, 5.0])
    )
    np.testing.assert_array_equal(draws, [1.0, 1.0, 1.0])
    assert result["observed_contrast"] == 3.0
    assert result["null_mean_contrast"] == 1.0
    assert result["residual_difference"] == 2.0
    assert result["null_sd"] == 0.0
    assert result["null_mean_monte_carlo_se"] == 0.0
    assert result["p_rank_two_sided"] == 0.5


def test_holm_family_handles_eight_dependent_tests_monotonically():
    m = analysis()
    actual = m.adjust_pvalues(np.array([.001, .01, .02, .03, .04, .2, .5, .9]))
    np.testing.assert_allclose(actual["holm"], [.008, .07, .12, .15, .16, .6, 1, 1],
                               rtol=0, atol=1e-14)


@pytest.mark.parametrize("cell,species,labels", [
    ([-1], [0], [1]),
    ([2], [0], [1]),
    ([0], [-1], [1]),
    ([0], [1], [1]),
    ([0], [0], [9]),
    ([0, 0], [0, 0], [1]),
])
def test_invalid_pairs_or_labels_are_rejected(cell, species, labels):
    m = analysis()
    with pytest.raises(ValueError):
        m.prepare_assemblage("invalid", np.array(cell), np.array(species),
                             np.array(labels), 2, 9)


def test_sector_weights_and_polar_hemispheres_use_fixed_source_coordinates():
    """Missing polar latitude cannot erase a known signed-northing hemisphere."""
    m = analysis()
    richness = np.array([30, 30, 30, 30, 30, 30, 30, 100, 29])
    zone = np.array([0, 1, 2, 3, 4, 5, 0, 2, 4])
    lat = np.array([70, -10, 50, -10, np.nan, np.nan, 66.562, 66.561999, 80])
    indices = m.summary_indices(richness, zone, lat)
    np.testing.assert_array_equal(indices[m.NORTH], [0, 2, 4, 6, 7])
    np.testing.assert_array_equal(indices[m.SOUTH], [1, 3, 5])
    np.testing.assert_array_equal(indices[m.ARCTIC], [0, 6])
    np.testing.assert_array_equal(indices[m.BOREAL], [2, 7])
    np.testing.assert_array_equal(indices["richness_ge100"], [7])
    c = m.prepare_assemblage("sectors", np.arange(9), np.zeros(9, dtype=int),
                             np.array([1]), 9, 9)
    c.indices = indices
    result = m.group_means(c, np.array([0, 10, 20, 30, 40, 50, 6, 28, 0]))
    assert result[m.EQUAL_SECTORS] == pytest.approx(26.166666666666668)
    assert result["richness_ge30"] == 23.0


def test_native_mobility_counts_unique_species_and_all_their_native_pairs():
    m = analysis()
    labels = np.array([0, 1, 2, 2])
    global_source = m.prepare_assemblage(
        "global", np.array([0, 0, 1, 1, 2, 2]),
        np.array([0, 1, 1, 2, 0, 3]), labels, 3, 9
    )
    native = m.prepare_assemblage("native", np.array([0, 0, 2]),
                                  np.array([0, 1, 0]), labels, 3, 9)
    result = m.exchangeability_support([np.array([0, 2])], global_source,
                                       native, np.zeros(4, dtype=int))
    assert result["global_movable_species"] == 2
    assert result["native_unique_species"] == 2
    assert result["native_unique_movable_species"] == 1
    assert result["native_species_cell_pairs_eligible_for_label_change"] == 2
    assert result["native_species_cell_pair_fraction_eligible_for_label_change"] == 2 / 3


def test_reproducibility_gate_rejects_drift_and_shape_changes():
    m = analysis()
    assert m.compare_gate("rounding", np.array([1e-11]), np.array([0.0]))["passed"]
    with pytest.raises(ValueError, match="differs"):
        m.compare_gate("drift", np.array([1.1e-10]), np.array([0.0]))
    with pytest.raises(ValueError, match="shape"):
        m.compare_gate("shape", np.array([0.0, 0.0]), np.array([0.0]))


def test_existing_output_directory_is_rejected_before_reading_sources(tmp_path):
    m = analysis()
    protected = tmp_path / "receipt.json"
    protected.write_text("frozen result")
    with pytest.raises(ValueError, match="overwrit"):
        m.run(SimpleNamespace(out=str(tmp_path), permutations=999))
    assert protected.read_text() == "frozen result"


def test_run_rejects_too_few_permutations_before_reading_sources(tmp_path):
    m = analysis()
    with pytest.raises(ValueError, match="199"):
        m.run(SimpleNamespace(out=str(tmp_path / "unused"), permutations=198))
    assert not (tmp_path / "unused").exists()
