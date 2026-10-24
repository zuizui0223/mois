"""Test genuine original-photo ID collision gates across old MOIS and new K1."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import pandas as pd
from scripts.acquisition import mois_k1_historical_photo_id_union as m
from scripts.acquisition.mois_completed_pixel_checkpoint import sha


class HistoricalAndK1PhotoIDUnionTests(unittest.TestCase):
    def setUp(self):
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.group = self.root/"k1-original-group0"
        self.group.mkdir()
        self.output = self.root/"output"
        self.old = {
            "101":("main","old-obs","white","classified_four_state_morph","a"*64)
        }
        self.old_obs = {"old-obs":("main","101")}
        self.old_counts = {
            "main":(1,1),
            "prior_group00":(0,0),
            "prior_R1":(0,0),
            "prior_R2":(0,0)
        }
        self.prior_overlap = {"unique_old_photo_IDs":1,"signed_arm_photo_rows":1}
        self.k1 = [
            dict(photo_id="201",observation_id="old-obs",
                 morph="white",measurement_status="classified_four_state_morph",
                 image_sha256="b"*64,source_arm=m.SOURCE_ARM,
                 selection_probability="0.00001"),
            dict(photo_id="202",observation_id="new-obs",
                 morph="mixed_uncertain",measurement_status="image_acquisition_failed",
                 image_sha256="",source_arm=m.SOURCE_ARM,
                 selection_probability="0.00002"),
        ]
        self.make_files()

    def make_files(self):
        df=pd.DataFrame(self.k1)
        df.to_parquet(self.group/"original_photo_ROIs.parquet",index=False)
        cov={
            "schema":"mois_fcp_K1_signed_original_taxon_recorded_month_latitude_coverage_v1",
            "sampled_original_photo_inclusion_probability_varies_per_photo":True,
            "one_selected_photo_per_original_observation":True,
            "individual_photo_open_flowers_independently_verified":False,
        }
        (self.group/"original_photo_source_coverage.json").write_text(json.dumps(cov))
        with (self.group/"original_photo_hue_counts.csv").open(
                "w",encoding="utf-8",newline="") as stream:
            writer=csv.writer(stream)
            writer.writerow(("morph","original_photo_attempts"))
            for hue in (*m.HUES,"mixed_uncertain"):
                writer.writerow((hue, sum(x["morph"]==hue for x in self.k1)))
        receipt={
            "original_hue_counts":{
                hue:sum(x["morph"]==hue for x in self.k1)
                for hue in (*m.HUES,"mixed_uncertain")
            },
            "artifact_SHA256":{name:sha(self.group/name) for name in m.SOURCE_SHAS}
        }
        self.receipt=self.group/"original_K1_group_receipt.json"
        self.receipt.write_text(json.dumps(receipt))
        self.raw_sha=sha(self.receipt)
        self.audit={
            "original_worker_group":0,
            "original_group_receipt_SHA256":self.raw_sha,
            "new_original_photo_attempts":2,
            "four_hue_classified":1,
            "image_acquisition_failed":1,
            "acquired_but_unclassifiable":0,
        }
        fingerprint=hashlib.sha256(f"0:{self.raw_sha}".encode()).hexdigest()
        self.cumulative={
            "schema":m.CUMULATIVE_SCHEMA,
            "source_run":m.SOURCE_RUN,
            "new_original_photo_target_total":m.NEW_TARGET,
            "total_scheduled_original_worker_groups":m.NUM_GROUPS,
            "previously_measured_K1_photo_IDs_reused_elsewhere_NOT_NEW":m.REUSED_OLD_ORIGINALS,
            "new_original_images_downloaded_or_colour_inferred_in_THIS_ledger":0,
            "this_is_NOT_a_crossarm_photo_ID_deduplicated_global_total":True,
            "per_group":[{
                "original_worker_group":0,
                "original_group_receipt_SHA256":self.raw_sha,
                "new_original_photo_attempts":2,
                "four_hue_classified":1,
            }],
            "completed_and_SHA_sealed_original_worker_groups":1,
            "completed_original_worker_group_IDs":[0],
            "source_original_group_receipt_aggregate_SHA256":fingerprint,
            "completed_real_new_original_photo_attempts":2,
            "completed_automatic_four_hue_classifications":1,
            "completed_original_acquisition_failures":1,
            "completed_acquired_but_automatic_four_hue_unclassifiable":0,
        }

    def run_merge(self):
        with patch.object(m,"PRIOR_COUNTS",self.old_counts), \
             patch.object(m,"audit_receipt",return_value=self.audit):
            return m.merge_verified(self.old,self.old_obs,[],
                self.prior_overlap,self.cumulative,[self.group],self.output,
                expected_prior=1)

    def test_exact_source_photo_ID_union_and_observation_nonindependence(self):
        r=self.run_merge()
        self.assertEqual(r["historical_actual_unique_original_photo_IDs"],1)
        self.assertEqual(r["new_K1_sha_verified_original_photo_IDs_disjoint_from_history"],2)
        self.assertEqual(r["exact_global_historical_plus_K1_distinct_original_photo_IDs"],3)
        self.assertEqual(r["exact_global_historical_plus_K1_automatic_four_hue_photo_ROIs"],2)
        self.assertEqual(r["K1_vs_old_crossarm_literal_original_observation_ID_overlaps"],1)
        self.assertTrue(r["K1_vs_historical_observation_level_independence_NOT_assumed"])
        self.assertFalse(r["unweighted_globally_representative_white_flower_rate_estimated"])
        self.assertEqual(r["new_original_image_downloads_or_colour_reinference_in_this_union"],0)

    def test_old_photo_repeated_in_K1_is_rejected_before_publishing(self):
        self.k1[0]["photo_id"]="101"
        self.make_files()
        with self.assertRaisesRegex(ValueError,"photo IDs repeat"):
            self.run_merge()
        self.assertFalse(self.output.exists())

    def test_original_roi_parquet_byte_tampering_is_rejected(self):
        file=self.group/"original_photo_ROIs.parquet"
        file.write_bytes(file.read_bytes()+b"bad bytes")
        with self.assertRaisesRegex(ValueError,"artifact SHA changed"):
            self.run_merge()
        self.assertFalse(self.output.exists())

    def test_unknown_cannot_become_unproven_white(self):
        self.k1[1]["morph"]="white"
        self.make_files()
        with self.assertRaisesRegex(ValueError,"original terminal changed"):
            self.run_merge()
        self.assertFalse(self.output.exists())

    def test_missing_group_or_tampered_cumulative_sha_fails_closed(self):
        self.cumulative["source_original_group_receipt_aggregate_SHA256"]="0"*64
        with self.assertRaisesRegex(ValueError,"receipt set SHA"):
            self.run_merge()
        self.assertFalse(self.output.exists())
        self.make_files()
        with patch.object(m,"PRIOR_COUNTS",self.old_counts):
            with self.assertRaisesRegex(ValueError,"Not all actually finished"):
                m.merge_verified(self.old,self.old_obs,[],self.prior_overlap,
                                 self.cumulative,[],self.output,expected_prior=1)

    def test_original_classification_and_historical_count_cannot_drift(self):
        self.old["101"]=("main","old-obs","mixed_uncertain","not_evaluable_roi","")
        with self.assertRaisesRegex(ValueError,"four-hue total"):
            self.run_merge()
        self.assertFalse(self.output.exists())


if __name__=="__main__":
    unittest.main()
