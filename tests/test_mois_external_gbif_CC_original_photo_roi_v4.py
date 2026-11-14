"""No-network safeguards for independently sourced original GBIF CC photo ROI9."""
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from scripts.acquisition import mois_external_gbif_CC_original_photo_roi_v4 as m
from scripts.acquisition.mois_zero9_gbif_other_cc_media_pilot import (
    BASELINE, source_sha
)
from scripts.acquisition.mois_zero9_gbif_other_cc_media_union import SCHEMA


class ExternalCCGBIFPhotoROITests(unittest.TestCase):
    def source_fixture(self,root):
        s=Path(root)/"external";s.mkdir()
        d=pd.DataFrame([
            {"inat_taxon_id":"100","source_binomial":"Plantus alba",
             "gbif_COL_XR_accepted_species_key":"g1","gbif_occurrence_key":"o1",
             "gbif_dataset_key":"ds1",
             "image_url":"https://d2jcv3kl45hlgi.cloudfront.net/abc.jpg",
             "image_license":"http://creativecommons.org/licenses/by/4.0/"},
            {"inat_taxon_id":"100","source_binomial":"Plantus alba",
             "gbif_COL_XR_accepted_species_key":"g1","gbif_occurrence_key":"o2",
             "gbif_dataset_key":"ds1",
             "image_url":"https://oxalis.br.fgov.be/images/BR00001.jpg",
             "image_license":"http://creativecommons.org/publicdomain/zero/1.0/"},
            {"inat_taxon_id":"101","source_binomial":"Plantus rosea",
             "gbif_COL_XR_accepted_species_key":"g2","gbif_occurrence_key":"o3",
             "gbif_dataset_key":"ds2",
             "image_url":"https://ibss-images.calacademy.org/photo.jpg",
             "image_license":"https://creativecommons.org/licenses/by/4.0/"}
        ])
        img=s/"mois_CC_noninat_GBIF_unique_image_metadata_304taxa.csv.gz"
        d.to_csv(img,index=False,compression="gzip")
        manifest={"schema":SCHEMA,
                  "source_original_42111_taxon_photo_atlas_sha256":BASELINE,
                  "stratified_hash_draw_taxa":304,
                  "photo_pixels_downloaded":0,
                  "new_strict9_numerical_photo_measurements":0,
                  "outputs_sha256":{img.name:source_sha(img)}}
        (s/"receipt.json").write_text(json.dumps(manifest)+"\n")
        return s

    def test_exact_trusted_https_provider_only(self):
        with tempfile.TemporaryDirectory() as root:
            src=self.source_fixture(root)
            rows=m.approved_media_rows(src,{
                 "d2jcv3kl45hlgi.cloudfront.net","oxalis.br.fgov.be",
                 "ibss-images.calacademy.org"},quota=2)
            self.assertEqual(len(rows),2)
            self.assertTrue(rows.inat_taxon_id.is_unique)
            self.assertTrue(rows.observation_id.is_unique)
            self.assertTrue(rows.photo_url_large.str.startswith("https://").all())
            self.assertTrue(rows.photo_license.str.contains("creativecommons").all())
            self.assertTrue(rows.photo_id.str.startswith("GBIF-CC-media:").all())
            self.assertTrue(rows.image_filename.str.endswith(".jpg").all())
            self.assertEqual(set(rows.inat_taxon_id),{"100","101"})
            m.source(src,Path(root)/"frozen",quota=2,allowed_hosts=[
                 "d2jcv3kl45hlgi.cloudfront.net","oxalis.br.fgov.be",
                 "ibss-images.calacademy.org"])
            cohort=json.loads((Path(root)/"frozen"/"manifest.json").read_text())
            self.assertEqual(cohort["original_external_image_count_planned"],2)
            self.assertFalse(cohort["photo_colour_measured"])
            self.assertFalse(cohort["downstream_appending_to_iNaturalist_original_frame_allowed"])

    def test_unsafe_host_and_reused_original_url_not_approved(self):
        hosts={"d2jcv3kl45hlgi.cloudfront.net"}
        self.assertEqual(m.validate_https_public_host(
            "https://d2jcv3kl45hlgi.cloudfront.net/a.jpg",hosts),
            "d2jcv3kl45hlgi.cloudfront.net")
        for url in (
            "http://d2jcv3kl45hlgi.cloudfront.net/a.jpg",
            "https://127.0.0.1/image.jpg",
            "https://d2jcv3kl45hlgi.cloudfront.net:8443/a.jpg",
            "https://inaturalist-open-data.s3.amazonaws.com/photos/10/large.jpg",
            "https://evil.example.org/flower.jpg",
        ):
            with self.assertRaises(ValueError):
                m.validate_https_public_host(url,hosts)
        with tempfile.TemporaryDirectory() as directory:
            src=self.source_fixture(directory)
            with self.assertRaisesRegex(ValueError,"No GBIF licensed photo URL"):
                m.approved_media_rows(src,{"upload.wikimedia.org"},quota=2)

    def test_source_SHA_tamper_is_detected(self):
        with tempfile.TemporaryDirectory() as root:
            src=self.source_fixture(root)
            file=src/"mois_CC_noninat_GBIF_unique_image_metadata_304taxa.csv.gz"
            with file.open("ab") as f:f.write(b"forged")
            with self.assertRaisesRegex(ValueError,"digest changed"):
                m.approved_media_rows(src,{"ibss-images.calacademy.org"})

    def test_external_measurement_cannot_be_automatically_added_to_original_frame(self):
        self.assertNotEqual(m.COHORT,"mois_full_source_observation_sampling_cohort_v1")


if __name__=="__main__":
    unittest.main()
