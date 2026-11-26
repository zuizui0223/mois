"""No-network tests of exact GBIF occurrence-media provenance and field/specimen gate."""
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import pandas as pd

from scripts.acquisition import mois_CCGBIF_noninat_field_vs_specimen as m
from scripts.acquisition.mois_zero9_gbif_other_cc_media_pilot import source_sha


class GBIFRecordTypeTests(unittest.TestCase):
    def row(self):
        return SimpleNamespace(
            inat_taxon_id="100",source_binomial="Plantus alba",
            gbif_occurrence_key="19",gbif_dataset_key="ds",
            gbif_COL_XR_accepted_species_key="G123",
            image_url="https://d2jcv3kl45hlgi.cloudfront.net/a.jpg",
            image_license="https://creativecommons.org/licenses/by/4.0/",
            media_host="d2jcv3kl45hlgi.cloudfront.net",
        )

    def data(self,basis="HUMAN_OBSERVATION",**changes):
        result={"key":19,"datasetKey":"ds","species":"Plantus alba",
                "basisOfRecord":basis,
                "media":[{"identifier":self.row().image_url}],
                "eventDate":"2020-01-01"}
        result.update(changes)
        return result

    def test_field_is_not_equivalent_to_true_flower_visibility(self):
        a=m.inspect_row(self.row(),self.data())
        self.assertTrue(a["same_source_species_and_dataset_proven"])
        self.assertEqual(a["field_vs_specimen_gate"],
                         "potential_field_observation_not_flower_verified")
        self.assertFalse(a["flower_visible_verified"])
        self.assertFalse(a["image_pixels_downloaded"])

    def test_herbarium_photos_not_counted_as_field_flower(self):
        a=m.inspect_row(self.row(),self.data("PRESERVED_SPECIMEN"))
        self.assertEqual(a["field_vs_specimen_gate"],
                         "preserved_or_material_specimen_not_field_flower")

    def test_taxon_or_media_URL_mismatch_cannot_pass(self):
        a=m.inspect_row(self.row(),self.data(species="Plantus rubra"))
        self.assertEqual(a["field_vs_specimen_gate"],
                         "occurrence_url_or_taxonomy_mismatch_review")
        b=m.inspect_row(self.row(),self.data(media=[{"identifier":"https://other.org/a.jpg"}]))
        self.assertEqual(b["field_vs_specimen_gate"],
                         "occurrence_url_or_taxonomy_mismatch_review")

    def test_no_partially_signed_source_strata_combined(self):
        with tempfile.TemporaryDirectory() as td:
            roots=[]
            for i in range(4):
                p=Path(td)/f"part{i}";p.mkdir()
                d=pd.DataFrame([
                  {"inat_taxon_id":str(i+1),"gbif_occurrence_key":str(20+i),
                   "field_vs_specimen_gate":
                    "potential_field_observation_not_flower_verified"},
                  {"inat_taxon_id":str(i+100),"gbif_occurrence_key":str(i+120),
                   "field_vs_specimen_gate":
                    "preserved_or_material_specimen_not_field_flower"}
                ])
                a=p/f"mois_GBIF_CC_media_occurrence_basis_group{i}.csv"
                d.to_csv(a,index=False)
                (p/"receipt.json").write_text(json.dumps({
                    "schema":m.SCHEMA,"stratum_index":i,
                    "result_SHA256":source_sha(a)})+"\n")
                roots.append(p)
            result=m.combine(roots,Path(td)/"all")
            self.assertEqual(result["actual_GBIF_occurrence_metadata_calls"],8)
            self.assertEqual(result["field_image_candidates_not_botanically_verified"],4)
            self.assertEqual(result["preserved_specimen_record_count"],4)
            self.assertEqual(result["numeric9_ROI_measured"],0)


if __name__=="__main__":
    unittest.main()
