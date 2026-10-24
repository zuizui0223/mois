"""Data-only original FCP ROI-v4 source/terminal/photo identity evidence tests."""
from __future__ import annotations
import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import pandas as pd

from scripts.acquisition import mois_fcp_u100_depth1_real_photo_collector as c
from scripts.acquisition.mois_queue_from_full93m_target_tiers import WORKER_COLUMNS
from scripts.analysis.mois_reassemble_original_new4096_measurements import METHOD


def file_csv(file,rows):
    with file.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


class GenuineFCPFirstPhotoCollectorTests(unittest.TestCase):
    def setUp(self):
        self.temp=TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.source=self.root/"source"
        self.workers=self.root/"workers"
        self.out=self.root/"original-photos"
        (self.source/"shards").mkdir(parents=True)
        (self.workers/"shard_000/sealed").mkdir(parents=True)
        self.inputs=[]
        for i,photo in enumerate(["901","902"]):
            row={k:"" for k in WORKER_COLUMNS}
            row.update(measurement_id=f"MOIS-SIGNED-{i}",photo_id=photo,
                inat_taxon_id=str(i+10),species="Plantus sample",
                observation_id=f"real-uuid-{i}",observer_id=f"observer{i}",
                observed_on="2025-06-10",latitude="35",longitude="140",
                photo_license="cc-by",photo_url_large=
                    f"https://inaturalist-open-data.s3.amazonaws.com/photos/{photo}/large.jpg",
                image_filename=f"MOIS-SIGNED-{i}.jpg",
                source_arm="mois_full93m_fcp_u100_prefix1_one_observation_per_eligible_taxon_v1",
                selection_probability=format(1/(100+i)/2,".17g"),
                source_snapshot_date="2026-09-27",shard="0")
            self.inputs.append(row)
        self.roster=self.source/"shards/shard-000.csv"
        file_csv(self.roster,self.inputs)
        self.manifest={
            "schema":"mois_full_source_observation_sampling_cohort_v1",
            "fcp_prefix_source_schema":
                "mois_FCP_U100_full93m_progressive_prefix1_prior_and_active_K1_photo_reuse_v1",
            "selected_FCP_U100_eligible_source_taxa":6,
            "full_frame_photo_ids":93209351,
            "new_photo_attempts_planned":2,
            "previously_measured_original_photo_ids_reused":1,
            "inflight_K1_photo_IDs_reserved_without_redownload":1,
            "technical_pilot_photo_IDs_not_refetched":1,
            "unsupported_original_photo_format_selected_as_unknown":1,
            "photo_already_reserved_by_K1_never_fetched_twice":True,
            "FCP_top100_photo_prefix_inclusion_probability_corrected_to_1_over_N_M":True,
            "colour_measured":False,
            "source_arm":self.inputs[0]["source_arm"],
            "work_units":1,"worker_ids":[0],"group_ids":[0],
            "photo_shard_sha256":{"0":c.sha(self.roster)},
        }
        (self.source/"manifest.json").write_text(json.dumps(self.manifest))
        terminal=[
            {"measurement_id":"MOIS-SIGNED-0","morph":"white",
             "measurement_status":"classified_four_state_morph",
             "image_sha256":"a"*64},
            {"measurement_id":"MOIS-SIGNED-1","morph":"mixed_uncertain",
             "measurement_status":"not_evaluable_roi_or_flip_gate",
             "image_sha256":"b"*64},
        ]
        sealed=self.workers/"shard_000/sealed/terminal_p000_b0000000.csv"
        file_csv(sealed,terminal)
        photo_rows=[
            {**row,**data,"classified_four_hue":str(data["morph"]=="white").lower()}
            for row,data in zip(self.inputs,terminal)
        ]
        photos=self.workers/"shard_000/photos.csv"
        file_csv(photos,photo_rows)
        receipt={
            "schema":"mois_full_original93m_real_RGB_photo_unit_v1",
            "shard":0,"attempted":2,"downloaded":2,"classifiable":1,
            "acquisition_failed":0,"unclassifiable_acquired":1,
            "source_sha256":c.sha(self.roster),
            "method":METHOD,"terminal_sha256":c.sha(sealed),
            "photos_sha256":c.sha(photos),
            "photo_colours_measured":True,
            "raw_licensed_image_bytes_persisted":False,
        }
        (self.workers/"shard_000/receipt.json").write_text(json.dumps(receipt))
        self.diag=self.workers/"warm_roi_group_000_diagnostic.json"
        self.diag.write_text(json.dumps({
            "group":0,
            "frozen_source_blob":c.MODEL_BLOB,
            "method_or_thresholds_modified":False,
            "artificial_flower_labels_created":False,
            "frozen_yolo_and_onnx_initializations":1,
            "group_wall_seconds":3.0,
            "new_photos_downloaded":2,
            "photo_prefetch":{
                "loader_seconds":1.5,
                "consumer_wait_seconds":0.2,
            },
        }))

    def collect(self):
        return c.seal(self.source,self.workers,self.out,
                      expected_source_sha=c.sha(self.source/"manifest.json"),
                      expected=2,expected_source_n=6,
                      prior_reused=1,k1_reserved=1)

    def test_exact_source_original_terminal_and_technical_unknowns(self):
        r=self.collect()
        self.assertEqual(r["new_unique_original_photo_attempts"],2)
        self.assertEqual(r["new_automatically_four_hue_ROI_classified"],1)
        self.assertEqual(r["new_acquired_unclassifiable"],1)
        self.assertEqual(r["new_image_acquisition_failed"],0)
        self.assertEqual(r["prior_real_original_photo_ROI_reuse_without_new_inference"],1)
        self.assertEqual(r["overlapping_inflight_K1_photos_not_refetched"],1)
        self.assertEqual(r["new_photo_fetch_or_ROI_inference_in_this_collector"],0)
        self.assertTrue(r["independent_terminal_and_photo_source_ROI_SHA_verified"])
        self.assertEqual(r["photo_group_diagnostics"][0]["fixed_model_session_initializations"],1)
        x=pd.read_parquet(self.out/"new_FCP_U100_prefix1_original_photo_ROIs.parquet")
        self.assertEqual(len(x),2)
        self.assertEqual(set(x.morph),{"white","mixed_uncertain"})
        self.assertEqual(set(x.photo_id.astype(str)),{"901","902"})
        self.assertEqual(c.sha(self.out/"new_FCP_U100_prefix1_original_photo_ROIs.parquet"),
                         r["output_sha256"]["new_FCP_U100_prefix1_original_photo_ROIs.parquet"])

    def test_corrupt_original_terminal_sha_rejected(self):
        f=self.workers/"shard_000/sealed/terminal_p000_b0000000.csv"
        f.write_text(f.read_text()+"\ncorrupt")
        with self.assertRaisesRegex(ValueError,"Original terminal SHA"):
            self.collect()
        self.assertFalse(self.out.exists())

    def test_missing_warm_model_diagnostic_is_unverified(self):
        self.diag.unlink()
        with self.assertRaisesRegex(ValueError,"diagnostic absent"):
            self.collect()

    def test_missing_worker_photo_receipt_cannot_be_called_complete(self):
        (self.workers/"shard_000/receipt.json").unlink()
        with self.assertRaisesRegex(ValueError,"not all fully sealed"):
            self.collect()
        self.assertFalse(self.out.exists())

    def test_sample_photo_source_sha_tamper_refused(self):
        self.roster.write_text(self.roster.read_text()+"\n")
        with self.assertRaisesRegex(ValueError,"Source SHA mismatch"):
            self.collect()
        self.assertFalse(self.out.exists())

    def test_immutable_original_evidence_not_rewritten(self):
        self.collect()
        with self.assertRaisesRegex(ValueError,"overwrite"):
            self.collect()


if __name__=="__main__":
    unittest.main()
