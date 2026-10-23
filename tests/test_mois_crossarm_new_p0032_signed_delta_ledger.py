"""Actual original MOIS photo/observation cross-arm delta regression tests."""
from __future__ import annotations

import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
import pandas as pd

from scripts.acquisition import mois_crossarm_new_p0032_signed_delta_ledger as ledger
from scripts.acquisition.mois_completed_pixel_checkpoint import sha


class SignedRealPhotoCrossarmDeltaTests(unittest.TestCase):
    def setUp(self):
        self.tmp=TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.oldregistry=self.root/"oldregistry.json"
        self.new=self.root/"p0032.parquet"
        self.newreceipt=self.root/"newreceipt.json"
        self.first=self.root/"first.parquet"
        self.firstreceipt=self.root/"first_receipt.json"
        self.broad=self.root/"broad.parquet"
        self.broadreceipt=self.root/"broad_receipt.json"
        self.out=self.root/"crossarm.json"
        self.frames=[
            (self.first,self.firstreceipt,[("9","first-a"),("10","first-b")],1),
            (self.broad,self.broadreceipt,
             [("11","broad-a"),("12","broad-b"),("13","broad-b")],2),
        ]
        for path,receipt,rows,h in self.frames:
            pd.DataFrame(rows,columns=["photo_id","observation_id"]).to_parquet(path,index=False)
            receipt.write_text(json.dumps({
                "schema":"mois_global93m_observation_baseline_real_RGB_photo_results_v1",
                "attempted":len(rows),"classifiable":h,
                "measured_photo_parquet_sha256":sha(path),
            }))
        self.oldregistry.write_text(json.dumps({
            "schema":"mois_original_pixel_immutable_cross_arm_group_registry_v1",
            "main_checkpoint_run":ledger.BASE_RUN,
            "main_checkpoint_photo_attempts":8,
            "main_checkpoint_four_hue_classified":4,
            "previous_tier06_source_photo_attempts":2,
            "previous_tier06_four_hue_classified":1,
            "broad_completed_photo_attempts":3,
            "broad_completed_four_hue_classified":2,
            "broad_pending_photo_attempts":0,
            "original_photo_attempts_distinct_across_all_verified_arms":13,
            "four_hue_classifications_across_all_verified_arms":7,
            "original_photo_or_observation_ID_cross_arm_overlap":0,
            "all_original_photo_per_partition_SHA_receipts_verified":True,
            "broad_sealed_groups":[{"group":i} for i in range(5)],
        }))
        self.original_delta=pd.DataFrame([
            {"photo_id":"14","observation_id":"new-a","morph":"white",
             "measurement_status":"classified_four_state_morph",
             "selection_probability":0.001,"source_arm":"p0032-original","image_sha256":"f"*64},
            {"photo_id":"15","observation_id":"new-a","morph":"mixed_uncertain",
             "measurement_status":"not_evaluable_roi_or_flip_gate",
             "selection_probability":0.001,"source_arm":"p0032-original","image_sha256":"e"*64},
        ])
        self.original_delta.to_parquet(self.new,index=False)
        self.updated_delta_receipt()

    def updated_delta_receipt(self):
        self.newreceipt.write_text(json.dumps({
            "schema":ledger.DELTA_SCHEMA,
            "prior_immutable_main_checkpoint_run":ledger.BASE_RUN,
            "prior_immutable_main_photo_attempts":8,
            "prior_immutable_main_four_hue_classified":4,
            "completed_true_original_photo_attempts":2,
            "new_true_original_four_hue_classified":1,
            "new_sealed_p0032_groups":ledger.DELTA_GROUPS,
            "previous_photo_ids_and_new_pixel_ids_exactly_disjoint":True,
            "unchanged_original_roi_model_and_image_terminal_byte_SHA_verified":True,
            "new_original_photos_downloaded_from_S3":0,
            "new_roi_or_hue_inference_performed":False,
            "previous_photo_dataset_rewritten":False,
            "combined_disjoint_main_plus_delta_photo_attempts":10,
            "combined_disjoint_main_plus_delta_four_hue_classified":5,
            "source_arm":"p0032-original",
            "immutable_output_sha256":{self.new.name:sha(self.new)},
        }))

    def run_audit(self):
        return ledger.merge_cross_arm(
            self.oldregistry,self.new,self.newreceipt,
            self.first,self.firstreceipt,
            self.broad,self.broadreceipt,self.out,
            old_total=13,old_hues=7,delta_n=2,delta_hues=1,
            expected_new=15,expected_new_hues=8,
            main_run=ledger.BASE_RUN,base_n=8,base_hues=4,
            old_tier06_n=2,old_tier06_hues=1,
            broad_n=3,broad_hues=2)

    def test_sha_signed_original_source_ids_and_incremental_totals(self):
        r=self.run_audit()
        self.assertEqual(r["full_authenticated_distinct_original_photo_attempts"],15)
        self.assertEqual(r["full_authenticated_four_hue_original_image_classifications"],8)
        self.assertEqual(r["newly_completed_original_p0032_distinct_observations"],1)
        self.assertEqual(r["photo_ID_overlap_new_p0032_vs_other_original_sampling_arms"],0)
        self.assertEqual(r["observation_ID_overlap_new_p0032_vs_other_original_sampling_arms"],0)
        self.assertFalse(r["unweighted_global_flora_white_rate_inferred"])
        self.assertFalse(r["historical_six_arm_parquet_rewritten"])

    def test_cross_arm_same_original_observation_is_not_allowed(self):
        d=pd.read_parquet(self.new)
        d.loc[0,"observation_id"]="broad-b"
        d.to_parquet(self.new,index=False)
        self.updated_delta_receipt()
        with self.assertRaisesRegex(ValueError,"overlap another existing sampling arm"):
            self.run_audit()
        self.assertFalse(self.out.exists())

    def test_duplicate_original_photo_with_old_arm_rejected(self):
        d=pd.read_parquet(self.new)
        d.loc[1,"photo_id"]="10"
        d.to_parquet(self.new,index=False)
        self.updated_delta_receipt()
        with self.assertRaisesRegex(ValueError,"overlap another existing sampling arm"):
            self.run_audit()

    def test_acquisition_missingness_never_relabelled_as_white(self):
        d=pd.read_parquet(self.new)
        d.loc[1,"morph"]="white"
        d.to_parquet(self.new,index=False)
        self.updated_delta_receipt()
        with self.assertRaisesRegex(ValueError,"hue terminal status mismatch"):
            self.run_audit()

    def test_prior_crossarm_identity_or_sha_drift_is_rejected(self):
        self.first.write_bytes(self.first.read_bytes()+b"tamper")
        with self.assertRaisesRegex(ValueError,"not authenticated"):
            self.run_audit()

    def test_refuse_overwrite_signed_crossarm_snapshot(self):
        self.run_audit()
        with self.assertRaisesRegex(ValueError,"cannot be rewritten"):
            self.run_audit()


if __name__=="__main__":
    unittest.main()
