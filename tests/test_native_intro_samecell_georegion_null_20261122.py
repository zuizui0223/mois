"""Substantive invariants for same-cell synchronous native/introduced nulls."""
import unittest
from types import SimpleNamespace

import numpy as np

from scripts.analysis.mois_native_intro_samecell_georegion_null_20261122 import (
    PRIMARY, counts_for_cohort, make_pair_specs, permutation_stats,
    holm_adjust, weighted_statistic, simulate_paired,
)


def cohort(labels, cell, species, n_cells=2, n_hues=9):
    c = SimpleNamespace(original_labels=np.asarray(labels, dtype=np.int8),
                        cell=np.asarray(cell, dtype=np.int32),
                        species=np.asarray(species, dtype=np.int32),
                        n_cells=n_cells, n_hues=n_hues)
    c.counts = counts_for_cohort(c, c.original_labels)
    c.richness = c.counts.sum(axis=1)
    return c


class SameCellTests(unittest.TestCase):
    def test_pair_mask_intersection_and_inclusive_threshold(self):
        n = np.array([30, 200, 29, 120, 45])
        i = np.array([40, 29, 200, 100, 35])
        specs = make_pair_specs(n, i, np.array([0, 0, 1, 2, 3]))
        np.testing.assert_array_equal(specs[PRIMARY]['indices'], [0, 3, 4])
        np.testing.assert_array_equal(
            specs['both_richness_ge100_equal_cells']['indices'], [3])

    def test_minority_richness_weight_not_pooled_weight(self):
        specs = make_pair_specs(np.array([30, 150]), np.array([90, 100]), np.array([0, 1]))
        s = specs['both_richness_ge30_min_richness_weight']
        np.testing.assert_allclose(s['weights'], [30 / 130, 100 / 130])
        self.assertAlmostEqual(weighted_statistic(np.array([1., 3.]), s), 330 / 130)

    def test_nine_on_fixed_eight_mask(self):
        specs = make_pair_specs(np.array([30, 60, 100]), np.array([40, 70, 110]),
                                np.array([0, 1, 2]), common8=np.array([False, True, False]))
        np.testing.assert_array_equal(
            specs['both_nonwhite_richness_ge30_equal_cells']['indices'], [1])

    def test_eight_mask_cannot_include_nine_ineligible_cell(self):
        with self.assertRaisesRegex(ValueError, 'subset'):
            make_pair_specs(np.array([10, 30]), np.array([100, 100]), np.array([0, 1]),
                            common8=np.array([True, True]))

    def test_no_jointly_eligible_cells_is_explicit(self):
        with self.assertRaisesRegex(ValueError, 'eligible'):
            make_pair_specs(np.array([29, 100]), np.array([100, 29]), np.array([0, 1]))

    def test_region_balanced_weights_equalize_six_sectors(self):
        zone = np.repeat(np.arange(6), np.arange(30, 36))
        specs = make_pair_specs(np.full(len(zone), 100), np.full(len(zone), 100), zone)
        s = specs['six_coordinate_sectors_equal_weight']
        for z in range(6):
            self.assertAlmostEqual(s['weights'][zone[s['indices']] == z].sum(), 1 / 6)

    def test_count_retains_common_species_hue_in_both_statuses(self):
        labels = [0, 1, 2]
        native = cohort(labels, [0, 0, 1], [0, 1, 2])
        intro = cohort(labels, [0, 1, 1], [2, 0, 1])
        permutation = np.array([2, 0, 1], dtype=np.int8)
        cn, ci = counts_for_cohort(native, permutation), counts_for_cohort(intro, permutation)
        self.assertEqual(cn[1, 1], 1)
        self.assertEqual(ci[0, 1], 1)
        np.testing.assert_array_equal(cn.sum(1), native.richness)
        np.testing.assert_array_equal(ci.sum(1), intro.richness)

    def test_fixed_white_rejects_moved_identity(self):
        c = cohort([0, 1, 2], [0, 0], [1, 2], n_hues=8)
        with self.assertRaisesRegex(ValueError, 'White'):
            counts_for_cohort(c, np.array([1, 0, 2]))

    def test_chromatic_labels_remap_correctly_and_empty_cells_zero(self):
        c = cohort([0, 1, 8], [0, 0], [1, 2], n_hues=8)
        self.assertEqual(c.counts[0, 0], 1)
        self.assertEqual(c.counts[0, 7], 1)
        self.assertEqual(c.counts[1].sum(), 0)

    def test_synchronous_identical_membership_gives_zero_contrast_every_draw(self):
        labels = [0, 1, 2, 3]
        a = cohort(labels, [0, 0, 1], [0, 1, 2])
        b = cohort(labels, [0, 0, 1], [0, 1, 2])
        specs = {'toy': {'indices': np.array([0, 1]), 'weights': np.array([.5, .5]),
                         'weighting': 'equal cells'}}
        result = simulate_paired({'strict_native': a, 'strict_introduced': b},
                                 [np.arange(4)], {'strict': specs}, 30, 100,
                                 model='toy', log_every=0)
        v = result['draws']['strict']['toy']
        np.testing.assert_array_equal(v[:, 0], v[:, 1])
        np.testing.assert_array_equal(v[:, 2], np.zeros(30))

    def test_overlapping_donor_groups_rejected(self):
        a = cohort([0, 1, 2], [0, 0], [0, 1])
        s = {'x': {'indices': np.array([0]), 'weights': np.array([1.]), 'weighting': 'equal'}}
        with self.assertRaisesRegex(ValueError, 'Overlapping'):
            simulate_paired({'strict_native': a, 'strict_introduced': a},
                            [np.array([0, 1]), np.array([1, 2])], {'strict': s}, 2, 1)

    def test_rank_p_has_identity_correction_and_ties(self):
        r = permutation_stats(0., np.array([0., 1., 2., 3.]))
        self.assertEqual(r['p_rank_lower'], 2 / 5)
        self.assertEqual(r['p_rank_upper'], 1.)
        self.assertEqual(r['p_rank_two_sided'], 4 / 5)
        r = permutation_stats(-1., np.array([0., 1., 2., 3.]))
        self.assertEqual(r['p_rank_two_sided'], 2 / 5)

    def test_identical_null_and_observation_p_is_one(self):
        r = permutation_stats(0., np.zeros(9))
        self.assertEqual(r['p_rank_two_sided'], 1.)
        self.assertEqual(r['null_mean_monte_carlo_se'], 0.)

    def test_native_minus_introduced_sign_and_holm(self):
        result = permutation_stats(.2 - .1, np.array([.05, .06, .07]))
        self.assertGreater(result['observed_minus_null'], 0)
        np.testing.assert_allclose(holm_adjust([.01, .04]), [.02, .04])
        np.testing.assert_allclose(holm_adjust([.04, .01]), [.04, .02])


if __name__ == '__main__':
    unittest.main()
