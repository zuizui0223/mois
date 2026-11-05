"""Independent source-provenance and support tests for GloBI white ROI × ALAN."""
import unittest
import numpy as np
import pandas as pd
from scripts.analysis import mois_alan_globi_exclusive_guild_roi_2020 as g

class TestGuildSupport(unittest.TestCase):
    def fake(self):
        return pd.DataFrame([
            dict(cell_id="1",accepted_species="A a",evidence_group="moth_only",strict9_white=.8,direct_colour="white"),
            dict(cell_id="1",accepted_species="B b",evidence_group="butterfly_only",strict9_white=.2,direct_colour="white"),
            dict(cell_id="1",accepted_species="C c",evidence_group="moth_only",strict9_white=.6,direct_colour="nonwhite"),
            dict(cell_id="1",accepted_species="D d",evidence_group="butterfly_only",strict9_white=.1,direct_colour="nonwhite"),
            dict(cell_id="1",accepted_species="E e",evidence_group="both",strict9_white=.9,direct_colour="white"),
        ])
    def test_both_guild_documentation_does_not_become_nocturnal_only(self):
        src=pd.DataFrame({"plant_name":["A a","B b","C c","C c"],
            "guild":["moth_family_visit","butterfly_family_visit","moth_family_visit","butterfly_family_visit"]})
        result=g.evidence_groups(src)
        self.assertEqual(result["A a"],"moth_only")
        self.assertEqual(result["B b"],"butterfly_only")
        self.assertEqual(result["C c"],"both")
    def test_white_two_groups_and_four_class_counts(self):
        summary,checks=g.cell_colour_support(self.fake())
        self.assertEqual(len(summary),1)
        self.assertEqual(checks["cells_with_direct_species_white_both_evidence_groups"],1)
        self.assertEqual(checks["cells_with_direct_white_nonwhite_by_both_evidence_groups"],1)
        self.assertEqual(checks["cells_with_at_least_2_distinct_species_in_each_group"],1)
        self.assertAlmostEqual(summary.moth_minus_butterfly_white_pp.iloc[0],55.0)
    def test_technical_unknown_cannot_be_called_nonwhite(self):
        sample=self.fake()
        sample.loc[sample.accepted_species.eq("D d"),"direct_colour"]=np.nan
        _,result=g.cell_colour_support(sample)
        self.assertEqual(result["cells_with_direct_white_nonwhite_by_both_evidence_groups"],0)
        self.assertEqual(result["cells_with_direct_species_white_both_evidence_groups"],1)
    def test_species_duplicate_in_one_cell_fails(self):
        d=pd.concat([self.fake(),self.fake().iloc[[0]]],ignore_index=True)
        with self.assertRaisesRegex(ValueError,"Repeated same accepted species"):
            g.cell_colour_support(d)
    def test_no_same_group_assigned_twice(self):
        d=pd.DataFrame({"plant_name":["X x","X x"],"guild":["moth_family_visit"]*2})
        self.assertEqual(g.evidence_groups(d)["X x"],"moth_only")
    def test_bright_cell_without_species_both_guilds_not_selected(self):
        d=self.fake().loc[self.fake().accepted_species.eq("A a")]
        result,checks=g.cell_colour_support(d)
        self.assertEqual(len(result),0)
        self.assertEqual(checks["cells_with_direct_species_white_both_evidence_groups"],0)
    def test_bad_source_channel_fail_closed(self):
        d=pd.DataFrame({"plant_name":["A a"],"guild":["unconfirmed"]})
        with self.assertRaisesRegex(ValueError,"Unspecified"):
            g.evidence_groups(d)
    def test_all_photo_white_is_continuous_not_binarized(self):
        d=self.fake();d.loc[0,"strict9_white"]=.71
        result,_=g.cell_colour_support(d)
        self.assertNotEqual(result.moth_minus_butterfly_white_pp.iloc[0],55.0)

if __name__=="__main__":unittest.main()
