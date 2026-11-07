"""Source-positive pollinator evidence preserves multi-label overlap and missingness."""
import unittest
import pandas as pd
from scripts.analysis import mois_alan_positive_visitor_evidence_support_20261105 as m

class PositiveEvidenceTests(unittest.TestCase):
    def make(self):
        d=[
         {"longitude":140.1,"latitude":35.1,"accepted_species":"A alba",
          "color_state":"white","moth_record":1,"butterfly_record":1},
         {"longitude":140.1,"latitude":35.1,"accepted_species":"B rubra",
          "color_state":"nonwhite","moth_record":1,"butterfly_record":0},
         {"longitude":140.1,"latitude":35.1,"accepted_species":"C flava",
          "color_state":"nonwhite","moth_record":0,"butterfly_record":1},
         {"longitude":140.1,"latitude":35.1,"accepted_species":"D alba",
          "color_state":"white","moth_record":1,"butterfly_record":0},
         {"longitude":140.1,"latitude":35.1,"accepted_species":"E alba",
          "color_state":"white","moth_record":0,"butterfly_record":1},
        ]
        return pd.DataFrame(d)
    def test_four_evidence_classes_do_not_imply_guild_exclusivity(self):
        r=m.spatial_counts(self.make(),.25)
        self.assertEqual(r["four_evidence_color_classes_min1_cells"],1)
        self.assertEqual(r["four_evidence_color_classes_min2_cells"],0)
        self.assertEqual(r["four_classes_no_shared_white_species_cells"],0)
        self.assertEqual(r["white_species_shared_moth_butterfly_evidence_cell_rows"],1)
    def test_no_overlapping_white_species_is_distinct_from_supported_four_groups(self):
        d=self.make().loc[lambda x:x.accepted_species.ne("A alba")]
        r=m.spatial_counts(d,.25)
        self.assertEqual(r["four_evidence_color_classes_min1_cells"],1)
        self.assertEqual(r["four_classes_no_shared_white_species_cells"],1)
    def test_unknown_trait_is_not_a_negative_or_white(self):
        d=self.make()
        d.loc[d.accepted_species.eq("C flava"),"color_state"]="unknown"
        r=m.spatial_counts(d,.25)
        self.assertEqual(r["four_evidence_color_classes_min1_cells"],0)
    def test_duplicate_photo_same_accepted_species_in_cell_collapsed(self):
        d=self.make()
        d=pd.concat([d,d.iloc[[1]]],ignore_index=True)
        r=m.spatial_counts(d,.25)
        self.assertEqual(r["four_evidence_color_classes_min1_cells"],1)
    def test_declared_guild_labels_do_not_claim_true_nocturnal_pollination(self):
        self.assertNotIn("nocturnal",m.KNOWN_COLOUR)
        self.assertNotIn("diurnal",m.KNOWN_COLOUR)
    def test_signed_raw_source_shas_present(self):
        self.assertTrue(all(len(h)==64 for h in m.SOURCE.values()))
    def test_invalid_guild_evidence_does_not_pass(self):
        with self.assertRaisesRegex(ValueError,"Not a supported"):
            m.evidence_flags(pd.DataFrame([{"plant_name":"A alba","guild":"unknown"}]))
    def test_same_plant_with_both_visit_records_keeps_both_flags(self):
        d=pd.DataFrame([{"plant_name":"A alba","guild":"moth_family_visit"},
                        {"plant_name":"A alba","guild":"butterfly_family_visit"}])
        q=m.evidence_flags(d)
        self.assertEqual(int(q.loc["A alba","moth_record"]),1)
        self.assertEqual(int(q.loc["A alba","butterfly_record"]),1)

if __name__=="__main__":unittest.main()
