"""No-network tests for cautious GBIF COL XR species/taxon reconciliation."""
import unittest
import pandas as pd

from scripts.analysis.mois_gbif_colxr_nine_palette_name_pilot import (
    inspect_response, make_pilot, COL_XR
)


class GbifColXrMatchTests(unittest.TestCase):
    def test_only_exact_accepted_same_binomial_species_passes(self):
        accepted={"usage":{"key":"ABC1","canonicalName":"Campanula punctata",
                           "rank":"SPECIES","status":"ACCEPTED"},
                  "matchType":"EXACT"}
        result=inspect_response(accepted,"Campanula punctata")
        self.assertEqual(result["status"],"colxr_exact_accepted_species")
        self.assertEqual(result["key"],"ABC1")
        synonym={"usage":{"key":"old","canonicalName":"Campanula punctata",
                           "rank":"SPECIES","status":"SYNONYM"},
                 "acceptedUsage":{"key":"NEW1","canonicalName":"Campanula changed"},
                 "matchType":"EXACT"}
        self.assertEqual(inspect_response(synonym,"Campanula punctata")["status"],
                         "synonym_or_doubtful_review")
        highrank={"usage":{"key":"GG","canonicalName":"Campanula",
                          "rank":"GENUS","status":"ACCEPTED"},"matchType":"EXACT"}
        self.assertEqual(inspect_response(highrank,"Campanula punctata")["status"],
                         "not_species_rank_review")
        fuzzy=dict(accepted,matchType="FUZZY")
        self.assertEqual(inspect_response(fuzzy,"Campanula punctata")["status"],
                         "fuzzy_or_nonexact_review")
        self.assertEqual(inspect_response(accepted,"Campanula unrelated")["status"],
                         "canonical_disagreement_review")

    def test_real_v2_diagnostics_nested_match_type(self):
        real_v2={"usage":{"key":"34RMK","canonicalName":"Dendrophthoe curvata",
                          "rank":"SPECIES","status":"ACCEPTED"},
                 "diagnostics":{"matchType":"EXACT","confidence":97}}
        r=inspect_response(real_v2,"Dendrophthoe curvata")
        self.assertEqual(r["status"],"colxr_exact_accepted_species")
        self.assertEqual(r["key"],"34RMK")

    def test_pilot_species_id_unique_and_hue_stratified(self):
        rows=[]
        for p in ("white","yellow","red"):
            for i in range(5):
                rows.append({"inat_taxon_id":str(3000+len(rows)),
                             "species":f"Plantus_{p} epithet{i}",
                             "dominant_hue9":p,"hue_status":"provisional_single"})
        for i,r in enumerate(rows):
            r["species"]="Plantus epithet"+chr(97+i)
        rows.append({"inat_taxon_id":"501","species":"Plantus",
                     "dominant_hue9":"white","hue_status":"provisional_single"})
        rows.append({"inat_taxon_id":"502","species":"Plantus epithet",
                     "dominant_hue9":"","hue_status":"mixed_image_unresolved"})
        pilot=make_pilot(pd.DataFrame(rows),per_hue=2)
        self.assertEqual(len(pilot),6)
        self.assertEqual(pilot.dominant_hue9.value_counts().to_dict(),
                         {"white":2,"yellow":2,"red":2})
        self.assertTrue(pilot.inat_taxon_id.is_unique)


if __name__=="__main__":
    unittest.main()
