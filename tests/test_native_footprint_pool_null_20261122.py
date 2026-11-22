"""Checks of the new native-footprint estimand, not snapshots of real results."""
import unittest
from types import SimpleNamespace

import numpy as np

from scripts.analysis.mois_native_footprint_pool_null_20261122 import (
    build_native_pool, shuffle_groups, draw_labels, entropy_counts,
    permutation_statistics, paired_contrast, holm_adjust, simulate_native,
)


class NativeFootprintTests(unittest.TestCase):
    def test_native_counts_determine_width_and_no_footprint_is_excluded(self):
        coverage = [1, 2, 5, 6, 20, 21, 100, 101, 0]
        sp = np.repeat(np.arange(len(coverage)), coverage)
        cell = np.concatenate([np.arange(n) for n in coverage])
        pool = build_native_pool(np.zeros(9, dtype=int), np.arange(9),
                                 cell, sp, np.zeros(101, dtype=int))
        np.testing.assert_array_equal(pool['coverage'], coverage)
        np.testing.assert_array_equal(pool['breadth_bin'], [0, 1, 1, 2, 2, 3, 3, 4, -1])
        self.assertEqual(pool['affinity'][8], 7)
        self.assertEqual(pool['eligible'].sum(), 8)
        self.assertEqual(pool['breadth_bin'][0], 0)  # Even if its global range is large.

    def test_native_affinity_threshold_and_multizone(self):
        zone = np.r_[np.zeros(13, dtype=int), np.ones(7, dtype=int),
                     np.zeros(12, dtype=int), np.ones(8, dtype=int)]
        pool = build_native_pool(np.zeros(3, dtype=int), np.array([1, 2, 3]),
                                 np.arange(40), np.repeat([0, 1], 20), zone)
        np.testing.assert_array_equal(pool['affinity'], [0, 6, 7])
        np.testing.assert_allclose(pool['top_sector_share'], [.65, .6, 0])

    def test_duplicate_species_cell_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            build_native_pool(np.array([0]), np.array([1]),
                              np.array([0, 0]), np.array([0, 0]), np.array([0]))

    def test_groups_preserve_genus_native_width_and_affinity(self):
        labels = np.array([1, 2, 3, 4, 5, 6])
        genus = np.array([0, 0, 0, 0, 1, 1])
        cell = np.array([0, 0, 1, 1, 0, 0])
        pool = build_native_pool(genus, labels, cell, np.arange(6), np.array([0, 1]))
        groups = shuffle_groups(pool, labels, 9, geography=True)
        self.assertEqual([g.tolist() for g in groups], [[0, 1], [2, 3], [4, 5]])
        broad = shuffle_groups(pool, labels, 9, geography=False)
        self.assertEqual([g.tolist() for g in broad], [[0, 1, 2, 3], [4, 5]])

    def test_whites_and_absent_native_species_never_exchange(self):
        labels = np.array([0, 1, 2, 3])
        pool = build_native_pool(np.zeros(4, dtype=int), labels,
                                 np.zeros(3, dtype=int), np.arange(3), np.array([0]))
        groups = shuffle_groups(pool, labels, 8, geography=True)
        self.assertEqual([g.tolist() for g in groups], [[1, 2]])
        rng = np.random.default_rng(9)
        for _ in range(30):
            result = draw_labels(labels, groups, rng, 8)
            self.assertEqual(result[0], 0)
            self.assertEqual(result[3], 3)
            self.assertEqual(sorted(result[1:3]), [1, 2])

    def test_draw_is_one_label_per_species_across_cells(self):
        original = np.array([1, 2, 3])
        labels = draw_labels(original, [np.arange(3)], np.random.default_rng(3), 9)
        cell = np.array([0, 1, 2, 2])
        sp = np.array([0, 0, 1, 2])
        counts, entropy = entropy_counts(cell, sp, labels, 3, 9, original)
        np.testing.assert_array_equal(counts[0], counts[1])
        self.assertEqual(counts[0, labels[0]], 1)
        np.testing.assert_allclose(entropy, [0, 0, np.log(2)])

    def test_white_identity_violation_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'White'):
            entropy_counts(np.array([0]), np.array([1]), np.array([1, 0]),
                           1, 8, np.array([0, 1]))

    def test_overlapping_exchange_groups_are_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Overlapping'):
            draw_labels(np.array([1, 2, 3]), [np.array([0, 1]), np.array([1, 2])],
                        np.random.default_rng(3), 9)

    def test_rank_ties_and_finite_monte_carlo_correction(self):
        stats = permutation_statistics(1.0, np.array([0., 1., 2., 3.]))
        self.assertEqual(stats['p_rank_lower'], 3/5)
        self.assertEqual(stats['p_rank_upper'], 4/5)
        self.assertEqual(stats['p_rank_two_sided'], 1.0)
        extreme = permutation_statistics(-1.0, np.array([0., 1., 2., 3.]))
        self.assertEqual(extreme['p_rank_two_sided'], 2/5)

    def test_contrast_uses_paired_draws_not_differences_of_pvalues(self):
        a, b = np.array([1., 3., 2., 4.]), np.array([.5, 2., 1., 3.])
        stats, draws = paired_contrast(2.0, 1.0, a, b)
        np.testing.assert_allclose(draws, [.5, 1., 1., 1.])
        self.assertAlmostEqual(stats['observed_minus_null'], .125)
        self.assertEqual(stats['p_rank_lower'], 1.0)
        self.assertEqual(stats['p_rank_upper'], .8)

    def test_holm_four_controls_arbitrary_dependence(self):
        np.testing.assert_allclose(holm_adjust([.002, .012, .03, .4]), [.008, .036, .06, .4])

    def test_simulation_keeps_richness_and_cell_expectation_identity(self):
        original = np.array([1, 1, 2])
        cell, sp = np.array([0, 0, 1, 1]), np.array([0, 1, 0, 2])
        counts, observed = entropy_counts(cell, sp, original, 2, 9, original)
        c = SimpleNamespace(original_labels=original, cell=cell, species=sp,
                            n_cells=2, n_hues=9, richness=counts.sum(axis=1),
                            observed=observed, indices={'richness_ge30': np.array([0, 1])})
        result = simulate_native(c, [np.arange(3)], 31, 8, 'toy')
        self.assertAlmostEqual(result['expected'].mean(), result['draws']['richness_ge30'].mean())
        again = simulate_native(c, [np.arange(3)], 31, 8, 'toy')
        np.testing.assert_array_equal(result['expected'], again['expected'])


if __name__ == '__main__':
    unittest.main()
