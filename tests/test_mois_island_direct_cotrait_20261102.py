"""Direct flower trait co-occurrence, rank support, identity tests."""
import unittest
import pandas as pd
import numpy as np
from scripts.analysis import mois_island_direct_cotrait_20261102 as c


class TestDirectCotrait(unittest.TestCase):
    def test_missing_colour_never_becomes_other_colour(self):
        df = pd.DataFrame({"accepted_id":["1","2","3"], "green_brown":[1.,0.,np.nan],
                           "family":["A","A","A"],"genus":["X","X","X"],
                           "small_flower":[1.,0.,0.]})
        r = c.paired_summary(df,"small_flower","test")
        self.assertEqual(r["jointly_evaluable_accepted_species"],2)
        self.assertEqual(r["other_colour_evaluable"],1)

    def test_unmatched_genus_gives_null_not_zero(self):
        df = pd.DataFrame({"accepted_id":["1","2"],"green_brown":[1.,0.],
                           "genus":["A","B"],"family":["X","Y"],"small_flower":[1.,0.]})
        r = c.paired_summary(df,"small_flower","test")
        self.assertEqual(r["both_colour_genus_count"],0)
        self.assertIsNone(r["within_genus_equal_diff_pp"])

    def test_within_genus_equal_vs_species_weighting(self):
        df = pd.DataFrame({"accepted_id":["1","2","3","4","5","6"],
                           "green_brown":[1.,0.,0.,0.,1.,0.],
                           "genus":["G1","G1","G1","G1","G2","G2"],
                           "family":["F1"]*4+["F2"]*2,
                           "small_flower":[1.,0.,0.,0.,0.,1.]})
        r = c.paired_summary(df,"small_flower","test")
        self.assertEqual(r["both_colour_genus_count"],2)
        self.assertAlmostEqual(r["within_genus_equal_diff_pp"],0)
        self.assertAlmostEqual(r["within_family_equal_diff_pp"],0)
        self.assertEqual(r["within_genus_species_support"],6)

    def test_mixed_colour_and_underspecified_form_not_negative(self):
        self.assertNotIn("green_brown_inconspicuous|white",c.COLOR_STATES)
        self.assertNotIn("other_described",c.RESPONSE_STATES["open_radial_form"][1])
        self.assertEqual(c.RESPONSE_STATES["open_radial_form"][1]["open_radial"],1)

    def test_duplicate_species_within_scope_rejected(self):
        f = pd.DataFrame({"accepted_id":["1","1"],"green_brown":[1.,0.],
                          "genus":["A","A"],"family":["A","A"],"small_flower":[1.,0.]})
        with self.assertRaisesRegex(ValueError,"twice"):
            c.paired_summary(f,"small_flower","global")

    def test_sc_not_equated_to_selfing_and_mixed_outcome_not_forced(self):
        self.assertEqual(c.RESPONSE_STATES["self_compatible"][1]["SC"],1)
        self.assertNotIn("mixed_or_variable",c.RESPONSE_STATES["self_compatible"][1])
        self.assertNotIn("mixed_mating",c.RESPONSE_STATES["predominantly_selfing"][1])

    def test_no_positives_zero_given_both_colours_known(self):
        f = pd.DataFrame({"accepted_id":["1","2"],"green_brown":[1.,0.],
                          "small_flower":[0.,0.],"genus":["G","G"],"family":["F","F"]})
        r = c.paired_summary(f,"small_flower","test")
        self.assertEqual(r["species_equal_diff_pp"],0)
        self.assertEqual(r["within_genus_equal_diff_pp"],0)

    def test_missing_positive_class_returns_null_difference(self):
        f = pd.DataFrame({"accepted_id":["1","2"],"green_brown":[0.,0.],
                          "small_flower":[1.,0.],"genus":["G","G"],"family":["F","F"]})
        r = c.paired_summary(f,"small_flower","test")
        self.assertIsNone(r["species_equal_diff_pp"])
        self.assertIsNone(r["within_family_equal_diff_pp"])


if __name__ == "__main__":
    unittest.main()
