"""Synthetic no-network source/worker SHA tests for new-1909 original ROI9."""
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from unittest.mock import patch

from scripts.acquisition import mois_postsnap1909_remaining_true_original_ROI9 as m
from scripts.analysis.mois_species_first_palette9 import PALETTE,sha256


class Remaining1909GenuineOriginals(unittest.TestCase):
    def metadata(self,root):
        src=Path(root)/"candidates";src.mkdir()
        a=pd.DataFrame([
          {"inat_taxon_id":"1","photo_id":"950000001",
           "photo_url_large":"https://inaturalist-open-data.s3.amazonaws.com/photos/950000001/large.jpg",
           "photo_license":"cc-by","observation_id":"obs-new-1",
           "source_observation_created_at":"2026-11-05T08:00:00Z"},
          {"inat_taxon_id":"2","photo_id":"950000002",
           "photo_url_large":"https://inaturalist-open-data.s3.amazonaws.com/photos/950000002/large.jpg",
           "photo_license":"cc0","observation_id":"obs-new-2",
           "source_observation_created_at":"2026-11-06T08:00:00Z"},
          {"inat_taxon_id":"3","photo_id":"950000003",
           "photo_url_large":"https://inaturalist-open-data.s3.amazonaws.com/photos/950000003/large.jpg",
           "photo_license":"cc-by","observation_id":"obs-old",
           "source_observation_created_at":"2026-11-06T08:00:00Z"},
          {"inat_taxon_id":"4","photo_id":"950000004",
           "photo_url_large":"https://inaturalist-open-data.s3.amazonaws.com/photos/950000004/large.jpg",
           "photo_license":"all-rights-reserved","observation_id":"obs-new-4",
           "source_observation_created_at":"2026-11-06T08:00:00Z"},
        ])
        photo=src/"MOIS_post_snapshot_new1909_unmeasured_taxa_CC_photo_ID_opportunities.csv"
        a.to_csv(photo,index=False)
        (src/"receipt.json").write_text(json.dumps({
          "schema":m.SOURCE_SCHEMA,
          "new_source_taxa_API_queried":1909,
          "original_source_zero_numeric_at_pre20261028":2229,
          "actual_new_photo_pixels_downloaded":0,
          "new_strict9_numerical_ROI":0,
          "results_SHA256":{photo.name:sha256(photo)}
        }))
        return src

    def test_candidates_exclude_measured_taxon_prior_UUID_and_ineligible_photo(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);src=self.metadata(root)
            atlas=pd.DataFrame({"inat_taxon_id":["1","2","3","4"],
                               "after_strict9":[0,1,0,0]})
            old=pd.DataFrame({"photo_id":["950000000"],
                              "observation_id":["obs-old"]})
            frame,count=m.pick_source(src,atlas,old,quota=4)
            self.assertEqual(count,4)
            self.assertEqual(len(frame),1)
            self.assertEqual(frame.inat_taxon_id.item(),"1")
            self.assertEqual(frame.photo_id.item(),"950000001")

    def test_source_photo_and_original_FCP_worker_terminal_QC(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);src=root/"frozen";(src/"shards").mkdir(parents=True)
            rows=pd.DataFrame([
                {"inat_taxon_id":"1","photo_id":"950000001",
                 "observation_id":"obs1","photo_url_large":"https://inaturalist-open-data.s3.amazonaws.com/photos/950000001/large.jpg",
                 "photo_license":"cc-by","measurement_id":"m1"},
                {"inat_taxon_id":"2","photo_id":"950000002",
                 "observation_id":"obs2","photo_url_large":"https://inaturalist-open-data.s3.amazonaws.com/photos/950000002/large.jpg",
                 "photo_license":"cc-by","measurement_id":"m2"}])
            f=src/"shards"/"shard-000.csv";rows.to_csv(f,index=False)
            manifest={"schema":m.WORKER_SOURCE_SCHEMA,"source_arm":m.SOURCE_ARM,
                "source_previous_attempts":m.PRIOR_ATTEMPTS,
                "new_photo_attempts_planned":2,
                "source_historical_original_photo_and_observation_IDs_checked":True,
                "photo_shard_sha256":{"0":sha256(f)}}
            (src/"manifest.json").write_text(json.dumps(manifest))
            worker=root/"worker"/"shard_000";worker.mkdir(parents=True)
            terminal=rows.copy()
            terminal["morph"]=["white","mixed_uncertain"]
            terminal["image_sha256"]=["1"*64,"2"*64]
            for p in PALETTE:
                terminal["flower_fraction_"+p]=[
                    "1" if p=="white" else "0",""]
            photo=worker/"photos.csv";terminal.to_csv(photo,index=False)
            (worker/"receipt.json").write_text(json.dumps({
                "schema":"mois_full_original93m_real_RGB_photo_unit_v1",
                "photo_colours_measured":True,"attempted":2,
                "classifiable":1,"source_sha256":sha256(f),
                "photos_sha256":sha256(photo)}))
            data,proof=m.genuine_frame(src,root/"worker")
            self.assertEqual(len(data),2)
            self.assertEqual(int(data.strict9.sum()),1)
            self.assertEqual(data.dominant_label9.tolist(),["white",""])
            self.assertTrue(np.isnan(data.loc[1,"strict9_white"]))
            self.assertEqual(len(proof),3)
            note=m.audit(src,root/"worker",root/"audited")
            self.assertEqual(note["new_technical_unknowns"],1)
            self.assertFalse(note["image_SHA_collisions_with_all_1231255_prior_originals_certified"])


if __name__=="__main__":
    unittest.main()
