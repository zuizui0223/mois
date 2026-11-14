"""No-network signed GBIF field-only external FCP image source tests."""
import json
import tempfile
import unittest
from pathlib import Path
import pandas as pd

from scripts.acquisition import mois_CCGBIF_human_field_photo_FCP_source as m
from scripts.acquisition.mois_zero9_gbif_other_cc_media_pilot import source_sha


class FieldGBIFOriginalPhotoSourceTests(unittest.TestCase):
    def fixture(self,root):
        d=Path(root)/"field";d.mkdir()
        r=[
           {"inat_taxon_id":"11","source_binomial":"Plantus alba",
            "gbif_occurrence_key":"201","gbif_dataset_key":"ds1",
            "gbif_COL_XR_accepted_species_key":"1001",
            "image_url":"https://d2jcv3kl45hlgi.cloudfront.net/a.jpg",
            "image_license":"https://creativecommons.org/licenses/by/4.0/",
            "gbif_basisOfRecord":"HUMAN_OBSERVATION",
            "GBIF_media_attached_to_correct_species":True},
           {"inat_taxon_id":"11","source_binomial":"Plantus alba",
            "gbif_occurrence_key":"202","gbif_dataset_key":"ds1",
            "gbif_COL_XR_accepted_species_key":"1001",
            "image_url":"https://oxalis.br.fgov.be/b.jpg",
            "image_license":"https://creativecommons.org/licenses/by/4.0/",
            "gbif_basisOfRecord":"HUMAN_OBSERVATION",
            "GBIF_media_attached_to_correct_species":True},
           {"inat_taxon_id":"12","source_binomial":"Plantus rubra",
            "gbif_occurrence_key":"203","gbif_dataset_key":"ds2",
            "gbif_COL_XR_accepted_species_key":"1002",
            "image_url":"https://iiif.rbge.org.uk/a.jpg",
            "image_license":"https://creativecommons.org/publicdomain/zero/1.0/",
            "gbif_basisOfRecord":"HUMAN_OBSERVATION",
            "GBIF_media_attached_to_correct_species":True}
        ]
        p=d/m.FIELD_FILE
        pd.DataFrame(r).to_csv(p,index=False,compression="gzip")
        (d/"receipt.json").write_text(json.dumps({
           "schema":m.FIELD_SCHEMA,
           "source_accepted_species_without_iNat_numeric_ROI":289,
           "new_image_pixels_downloaded":0,
           "new_nine_fraction_ROI9_measurements":0,
           "outputs_sha256":{"field":source_sha(p)}
        })+"\n")
        return d

    def test_authenticated_photo_once_per_species_and_no_repeated_observations(self):
        with tempfile.TemporaryDirectory() as root:
            d=self.fixture(root)
            approved={"d2jcv3kl45hlgi.cloudfront.net","oxalis.br.fgov.be",
                      "iiif.rbge.org.uk"}
            table=m.select(d,approved,quota=2)
            self.assertEqual(len(table),2)
            self.assertEqual(table.inat_taxon_id.nunique(),2)
            self.assertTrue(table.observation_id.is_unique)
            self.assertTrue(table.image_url.is_unique)
            self.assertTrue(table.photo_url_large.str.startswith("https://").all())
            info=m.run(d,Path(root)/"frozen",list(approved),quota=2)
            self.assertFalse(info["photo_colour_measured"])
            self.assertFalse(info["downstream_appending_to_iNaturalist_original_frame_allowed"])
            self.assertFalse(info["botanical_flower_visibility_checked"])
            self.assertEqual(info["original_external_image_count_planned"],2)

    def test_excluded_external_images_not_reused(self):
        with tempfile.TemporaryDirectory() as root:
            d=self.fixture(root)
            approved={"d2jcv3kl45hlgi.cloudfront.net","oxalis.br.fgov.be"}
            old={"https://d2jcv3kl45hlgi.cloudfront.net/a.jpg"}
            t=m.select(d,approved,quota=2,exclude_URLs=old)
            self.assertEqual(len(t),1)
            self.assertEqual(t.image_url.iloc[0],"https://oxalis.br.fgov.be/b.jpg")

    def test_signed_field_only_source_cannot_promote_herbarium_or_fuzzy(self):
        with tempfile.TemporaryDirectory() as root:
            d=self.fixture(root)
            file=d/m.FIELD_FILE
            dd=pd.read_csv(file)
            dd.loc[0,"gbif_basisOfRecord"]="PRESERVED_SPECIMEN"
            dd.to_csv(file,index=False,compression="gzip")
            # A signature-changing edit is rejected before record-level checks.
            with self.assertRaisesRegex(ValueError,"hash missing or changed"):
                m.select(d,{"d2jcv3kl45hlgi.cloudfront.net"})


if __name__=="__main__":
    unittest.main()
