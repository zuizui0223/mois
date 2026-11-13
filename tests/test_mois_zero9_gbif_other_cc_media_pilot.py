"""Offline tests of public GBIF external CC-image metadata source pilot.

These tests explicitly do not inspect an image or assign a biological flower hue.
"""
import unittest
import pandas as pd
from scripts.acquisition import mois_zero9_gbif_other_cc_media_pilot as m


class GBIFAlternativeSourceTests(unittest.TestCase):
    def test_source_only_chooses_zero_numeric_photo_taxa_with_strict_binomial(self):
        atlas=pd.DataFrame({
            "inat_taxon_id":[str(i+1) for i in range(42111)],
            "species":["Plantus alba","Plantus rubra","Plantus","Plantus rosea"]+
                       [f"Plantus epithet{chr(97+i%26)}" for i in range(42107)],
            "dominant_hue9":["","","","white"]+[""]*42107,
            "hue_status":["no_strict9","mixed_image_unresolved","no_strict9",
                          "provisional_single"]+["no_strict9"]*42107,
            "after_strict9":[0,1,0,1]+[0]*42107,
        })
        selected=m.source_taxa(atlas,sample=10)
        self.assertEqual(len(selected),10)
        self.assertFalse(set(selected.inat_taxon_id)&{"2","3","4"})
        self.assertTrue(selected.species.str.contains(" ").all())
        self.assertTrue(selected.priority.is_unique)
        self.assertTrue(selected.inat_taxon_id.is_unique)

    def test_exact_name_only_accepts_same_accepted_binomial_species(self):
        src={"usage":{"key":"GBIF:123","canonicalName":"Plantus alba",
                      "rank":"SPECIES","status":"ACCEPTED"},
             "diagnostics":{"matchType":"EXACT"}}
        self.assertEqual(m.is_exact_accepted_name(src,"Plantus alba"),
                         ("GBIF:123","accepted_exact"))
        synonym={"usage":dict(src["usage"],status="SYNONYM"),
                 "diagnostics":{"matchType":"EXACT"}}
        self.assertEqual(m.is_exact_accepted_name(synonym,"Plantus alba")[0],"")
        self.assertEqual(m.is_exact_accepted_name(src,"Plantus rubra")[0],"")
        self.assertEqual(m.is_exact_accepted_name(dict(src,diagnostics={"matchType":"FUZZY"}),
                                                 "Plantus alba")[0],"")
        self.assertEqual(m.is_exact_accepted_name(src,"Plantus alba var. alba")[0],"")

    def test_media_photo_own_license_and_non_inaturalist_provider(self):
        self.assertTrue(m.safe_image("https://images.example.org/specimen/flower01.jpg",
                                    "https://creativecommons.org/licenses/by/4.0/"))
        self.assertTrue(m.safe_image("https://cdn.example.org/flower.jpg",
                                    "http://creativecommons.org/publicdomain/zero/1.0/"))
        self.assertFalse(m.safe_image("https://inaturalist-open-data.s3.amazonaws.com/photos/1/large.jpg",
                                    "https://creativecommons.org/licenses/by/4.0/"))
        self.assertFalse(m.safe_image("https://www.inaturalist.org/photos/1",
                                     "https://creativecommons.org/licenses/by/4.0/"))
        self.assertFalse(m.safe_image("https://images.example.org/flower.jpg",""))
        self.assertFalse(m.safe_image("file:///tmp/photo.jpg",
                                     "https://creativecommons.org/licenses/by/4.0/"))

    def test_metadata_only_keeps_provider_and_rejects_occurrence_name_conflict(self):
        good={"key":111,"datasetKey":"ds1","species":"Plantus alba",
              "media":[{"type":"StillImage","identifier":"https://herb.example.org/1.jpg",
                        "license":"https://creativecommons.org/licenses/by/4.0/",
                        "creator":"botany team"},
                       {"type":"StillImage",
                        "identifier":"https://inaturalist-open-data.s3.amazonaws.com/photos/22/large.jpg",
                        "license":"https://creativecommons.org/licenses/by/4.0/"}]}
        bad={"key":112,"species":"Plantus rubra","media":good["media"]}
        res,review=m.candidate_media([good,bad],key="123",
                                     source_name="Plantus alba",inat_taxon_id="99")
        self.assertEqual(len(res),1)
        self.assertEqual(res[0]["image_url"],"https://herb.example.org/1.jpg")
        self.assertFalse(res[0]["flower_visible_verified"])
        self.assertFalse(res[0]["colour_measured"])
        self.assertTrue("occurrence_name_conflict" in review)
        self.assertTrue("unsupported_source_or_media_license" in review)


if __name__=="__main__":
    unittest.main()
