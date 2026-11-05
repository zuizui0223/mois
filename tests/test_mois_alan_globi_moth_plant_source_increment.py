"""Fail-closed GloBI plant-moth category and photo provenance test."""
import unittest
import pandas as pd
from scripts.analysis import mois_alan_globi_moth_plant_source_increment as g


class TestSourceGloBI(unittest.TestCase):
    def row(self,family="Sphingidae",relation="visitsFlowersOf"):
        return {"source_taxon_name":"Xylophanes tersa",
                "source_taxon_path":"Animalia | Arthropoda | Insecta | Lepidoptera | "+
                           family+" | Xylophanes | Xylophanes tersa",
                "target_taxon_name":"Pentas lanceolata",
                "target_taxon_path":"Eukaryota | Viridiplantae | Angiosperms | Rubiaceae | Pentas | Pentas lanceolata",
                "interaction_type":relation,
                "study_title":"GloBI original study reference"}
    def test_correct_moth_family_evidence_survives(self):
        v=g.clean_record(self.row(),"Sphingidae","visitsFlowersOf")
        self.assertEqual(v["guild"],"moth_family_visit")
        self.assertEqual(v["plant_name"],"Pentas lanceolata")
    def test_butterfly_is_different_comparison_not_negative_moth(self):
        v=g.clean_record(self.row("Papilionidae"),"Papilionidae","visitsFlowersOf")
        self.assertEqual(v["guild"],"butterfly_family_visit")
    def test_global_api_false_positive_family_is_rejected(self):
        self.assertIsNone(g.clean_record(self.row("Mecoptera"),"Noctuidae","visitsFlowersOf"))
    def test_name_must_be_strict_botanical_binomial(self):
        v=self.row();v["target_taxon_name"]="Pentas"
        self.assertIsNone(g.clean_record(v,"Sphingidae","visitsFlowersOf"))
    def test_pollinates_cannot_come_from_visit_statement(self):
        self.assertIsNone(g.clean_record(self.row(),"Sphingidae","pollinates"))
    def test_nonplant_target_rejected(self):
        v=self.row();v["target_taxon_path"]="Animalia | Insecta | Lepidoptera"
        self.assertIsNone(g.clean_record(v,"Sphingidae","visitsFlowersOf"))
    def test_missing_original_study_not_silently_verified(self):
        v=self.row();v["study_title"]=""
        self.assertIsNone(g.clean_record(v,"Sphingidae","visitsFlowersOf"))
    def test_source_taxonomy_does_not_inflate_ids(self):
        tax=pd.DataFrame({"inat_taxon_id":["1","2","3"],
            "wcvp_match_state":["unique_exact_name_to_accepted_id",
               "multiple_accepted_ids_ambiguous","unique_exact_name_to_accepted_id"],
            "wcvp_accepted_plant_name_id":["11","12","11"]})
        ac=pd.DataFrame({"wcvp_accepted_plant_name_id":["11"],"wcvp_accepted_name":["Pentas lanceolata"]})
        r=g.source_taxonomy(tax,ac)
        self.assertEqual(len(r),2)
        self.assertEqual(set(r.wcvp_accepted_name),{"Pentas lanceolata"})
    def test_missing_no_record_not_negative_guild(self):
        self.assertIn("moth_family_visit",g.clean_record(self.row(),"Sphingidae","visitsFlowersOf")["guild"])
        self.assertNotIn("diurnal_only",("moth_family_visit","butterfly_family_visit"))


if __name__=="__main__":
    unittest.main()
