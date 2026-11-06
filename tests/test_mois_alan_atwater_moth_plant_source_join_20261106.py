"""Atwater direct moth flower visits/pollen carriage: source QA, no false night-pollination labels."""
import unittest
import pandas as pd
from scripts.analysis import mois_alan_atwater_moth_plant_source_join_20261106 as a

class AtwaterStudyTests(unittest.TestCase):
    def test_table_2_exact_plant_count_and_two_interaction_tiers(self):
        d=a.public_paper_ledger()
        self.assertEqual(len(d),7)
        self.assertTrue(d.plant_name.is_unique)
        self.assertEqual(int(d.flower_visitor_count_in_paper_table2.sum()),151)
        self.assertEqual(int(d.pollen_transport_count_in_paper_table2.sum()),28)
        self.assertEqual(int((d.pollen_transport_count_in_paper_table2>0).sum()),6)
        self.assertEqual(int((d.pollen_transport_count_in_paper_table2==0).sum()),1)
    def test_Liatris_typographical_note_is_not_unknown_taxon_promotion(self):
        d=a.public_paper_ledger()
        self.assertEqual(d.loc[d.plant_name.str.startswith("Liatris"),"plant_name"].iloc[0],"Liatris tenuifolia")
        self.assertIn("OCR",d.loc[d.plant_name.str.startswith("Liatris"),"source_taxonomy_note"].iloc[0])
    def test_no_claim_of_stigma_deposition_from_pollen_carriage(self):
        d=a.public_paper_ledger()
        self.assertTrue(d.source_evidence_tier.str.contains("NOT effective").all())
    def test_mock_source_unknown_colour_stays_unknown(self):
        c=pd.DataFrame([
         dict(photo_id="1",observation_id="o1",inat_taxon_id="taxon1",
              source_subarm="historical_main",strict9=True,
              latitude=35.,longitude=-81.,observed_on="2020-05-10"),
         dict(photo_id="2",observation_id="o2",inat_taxon_id="taxon2",
              source_subarm="historical_main",strict9=True,
              latitude=35.1,longitude=-81.1,observed_on="2024-05-10")])
        names=pd.DataFrame([
          dict(inat_taxon_id="taxon1",wcvp_match_state="unique_exact_name_to_accepted_id",
               wcvp_accepted_plant_name_id="idA"),
          dict(inat_taxon_id="taxon2",wcvp_match_state="unique_exact_name_to_accepted_id",
               wcvp_accepted_plant_name_id="idB")])
        accepted=pd.DataFrame([
          {"wcvp_accepted_plant_name_id":"idA","wcvp_accepted_name":"Ageratina aromatica"},
          {"wcvp_accepted_plant_name_id":"idB","wcvp_accepted_name":"Dalea pinnata"}])
        direct=pd.DataFrame([
          {"accepted_species":"Ageratina aromatica","trait_name":"flower_primary_color",
           "normalized_value":"mixed_or_variable","quality":"high","resolution_status":"resolved",
           "classification":"direct"},
          {"accepted_species":"Dalea pinnata","trait_name":"flower_primary_color",
           "normalized_value":"white","quality":"high","resolution_status":"resolved",
           "classification":"direct"}])
        # Source static VIIRS frame expected numbers in strict gate; use 1044 synthetic IDs
        original = pd.DataFrame({"cell_id":[f"{i}" for i in range(1044)],
                 "source_complete":[True]*1044})
        v24=pd.DataFrame({"cell_id":[f"{i}" for i in range(302)]})
        matched,yrs=a.join_sources(c,names,accepted,direct,original,v24)
        self.assertEqual(len(matched),2)
        self.assertEqual(matched.known_single_colour.eq("white").sum(),1)
        self.assertEqual(matched.known_single_colour.isna().sum(),1)
    def test_no_study_moth_dependence_encoded(self):
        self.assertNotIn("confirmed_nocturnal_pollination",a.public_paper_ledger().columns)

if __name__=="__main__":unittest.main()
