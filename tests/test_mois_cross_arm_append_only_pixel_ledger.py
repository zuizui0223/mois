"""Cross-arm real original photo-ID/observation-ID SHA-only union audit."""
from __future__ import annotations

import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import pandas as pd

from scripts.acquisition import mois_cross_arm_append_only_pixel_ledger as ledger


class SignedCrossArmOriginalLedgerTests(unittest.TestCase):
    def setUp(self):
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.main = self.root/"main.parquet"
        self.main_receipt = self.root/"main_receipt.json"
        self.old = self.root/"old.parquet"
        self.old_receipt = self.root/"old_receipt.json"
        self.source = self.root/"source.json"
        self.output = self.root/"union.json"
        self.groups = [self.root/"group0",self.root/"group1"]
        self.p = 0.00005

        pd.DataFrame({
            "photo_id":["1","2","3","4"],
            "observation_id":["old-a","old-b","old-b","old-c"]
        }).to_parquet(self.main,index=False)
        self.main_receipt.write_text(json.dumps({
            "schema":"mois_sealed_actual_pixel_checkpoint_v1",
            "combined_unique_photo_attempts":4,
            "combined_four_hue_classified":2,
            "output_sha256":{self.main.name:ledger.sha(self.main)}}))
        pd.DataFrame({
            "photo_id":["5","6"],
            "observation_id":["prior-1","prior-2"]
        }).to_parquet(self.old,index=False)
        self.old_receipt.write_text(json.dumps({
            "schema":"mois_global93m_observation_baseline_real_RGB_photo_results_v1",
            "attempted":2,
            "classifiable":1,
            "measured_photo_parquet_sha256":ledger.sha(self.old)}))
        self.source.write_text(json.dumps({
            "schema":"mois_full_source_observation_sampling_cohort_v1",
            "full_frame_photo_ids":93209351,
            "source_tier":6,
            "source_group":"all16",
            "new_photo_attempts_planned":4,
            "compute_bucket_count":512,
            "admitted_bucket_count":4,
            "not_probability_complete_global_source":True,
            "colour_measured":False,
            "selection_probability":self.p,
            "source_arm":"synthetic-tier06-broad",
            "work_units":8,"group_ids":[0,1],
        }))
        for i,dir in enumerate(self.groups):
            dir.mkdir()
            rows=pd.DataFrame({
                "photo_id":[str(7+2*i),str(8+2*i)],
                "observation_id":[f"broad-{2*i}",f"broad-{2*i+1}"],
                "morph":["white","mixed_uncertain"],
                "measurement_status":["classified_four_state_morph",
                                      "not_evaluable_roi_or_flip_gate"],
                "source_arm":["synthetic-tier06-broad"]*2,
                "selection_probability":[self.p,self.p]
            })
            parquet=dir/"original_photo_ROIs.parquet"
            colors=dir/"group_original_photo_hue_counts.csv"
            receipt=dir/"original_photo_group_receipt.json"
            rows.to_parquet(parquet,index=False)
            colors.write_text("morph,original_photo_attempts\nwhite,1\nmixed_uncertain,1\n")
            receipt.write_text(json.dumps({
                "schema":"mois_tier06_all16_append_only_original_RGB_group_v1",
                "original_acquisition_workflow_run":ledger.BROAD_RUN,
                "source_original_photo_manifest_sha256":ledger.sha(self.source),
                "complete_expected_source_photo_attempts":4,
                "original_source_warm_group":i,
                "original_work_units":list(range(i*4,(i+1)*4)),
                "distinct_genuinely_completed_photo_ids":2,
                "distinct_original_observation_ids":2,
                "four_hue_model_classified":1,
                "acquisition_failed":0,
                "acquired_but_unclassifiable":1,
                "hue_counts":{"white":1,"mixed_uncertain":1},
                "all_source_photo_and_original_terminal_SHA_verified":True,
                "frozen_original_FCP_model_and_output_verified":True,
                "image_bytes_reacquired_in_sealing":False,
                "ROI_or_hue_reinference_in_sealing":False,
                "globally_probability_complete":False,
                "saved_file_sha256":{
                    parquet.name:ledger.sha(parquet),
                    colors.name:ledger.sha(colors)
                }
            }))

    def run_audit(self,groups=None):
        return ledger.register(
            self.main,self.main_receipt,self.old,self.old_receipt,
            self.source,self.groups if groups is None else groups,self.output,
            expected_main=(4,2),expected_old=(2,1),
            expected_source_sha=ledger.sha(self.source),expected_broad_total=4)

    def test_verified_partial_count_and_persist_sha_only_inventory(self):
        r=self.run_audit([self.groups[0]])
        self.assertEqual(r["broad_completed_photo_attempts"],2)
        self.assertEqual(r["broad_pending_photo_attempts"],2)
        self.assertEqual(r["original_photo_attempts_distinct_across_all_verified_arms"],8)
        self.assertEqual(r["four_hue_classifications_across_all_verified_arms"],4)
        self.assertTrue(r["all_original_photo_per_partition_SHA_receipts_verified"])
        self.assertFalse(r["old_pixel_Parquets_rewritten"])
        self.assertFalse(r["global_flora_white_flower_rate_estimable_by_naive_pooled_rows"])
        self.assertEqual(json.loads(self.output.read_text())["broad_sealed_groups"][0]["group"],0)

    def test_full_original_groups_distinct_photo_and_observation(self):
        r=self.run_audit()
        self.assertEqual(r["broad_completed_photo_attempts"],4)
        self.assertEqual(r["broad_pending_photo_attempts"],0)
        self.assertEqual(r["original_photo_attempts_distinct_across_all_verified_arms"],10)
        self.assertEqual(r["four_hue_classifications_across_all_verified_arms"],5)
        self.assertEqual(len(r["broad_sealed_groups"]),2)

    def test_group_parquet_sha_tamper_refused(self):
        parquet=self.groups[0]/"original_photo_ROIs.parquet"
        parquet.write_bytes(parquet.read_bytes()+b"tampered")
        with self.assertRaisesRegex(ValueError,"bytes changed"):
            self.run_audit()
        self.assertFalse(self.output.exists())

    def test_same_original_observation_across_two_compute_groups_is_legitimate(self):
        # Source grouping is a 4*256-PHOTO scheduling boundary, not a partition
        # of observation IDs. Same original observation with two DIFFERENT
        # photo IDs may straddle the boundary of two completed work groups.
        path=self.groups[1]/"original_photo_ROIs.parquet"
        data=pd.read_parquet(path)
        data.loc[0,"observation_id"]="broad-0"
        data.to_parquet(path,index=False)
        receipt=self.groups[1]/"original_photo_group_receipt.json"
        record=json.loads(receipt.read_text())
        record["saved_file_sha256"][path.name]=ledger.sha(path)
        receipt.write_text(json.dumps(record))
        out=self.run_audit()
        self.assertEqual(out["same_broad_arm_observations_spanning_original_compute_groups"],1)
        self.assertEqual(out["broad_arm_distinct_original_observations"],3)
        self.assertEqual(out["original_photo_attempts_distinct_across_all_verified_arms"],10)
        self.assertEqual(out["original_photo_or_observation_ID_cross_arm_overlap"],0)

    def test_same_observation_reused_from_historical_source_arm_is_rejected(self):
        path=self.groups[1]/"original_photo_ROIs.parquet"
        data=pd.read_parquet(path)
        data.loc[0,"observation_id"]="old-a"
        data.to_parquet(path,index=False)
        receipt=self.groups[1]/"original_photo_group_receipt.json"
        record=json.loads(receipt.read_text())
        record["saved_file_sha256"][path.name]=ledger.sha(path)
        receipt.write_text(json.dumps(record))
        with self.assertRaisesRegex(ValueError,"historical source-arm observation reused"):
            self.run_audit()
        self.assertFalse(self.output.exists())


    def test_group_reused_twice_is_not_double_counted(self):
        with self.assertRaisesRegex(ValueError,"work units missing or duplicated"):
            self.run_audit([self.groups[0],self.groups[0]])
        self.assertFalse(self.output.exists())

    def test_inconsistent_real_hue_and_label_rejected(self):
        path=self.groups[0]/"original_photo_ROIs.parquet"
        df=pd.read_parquet(path)
        df.loc[1,"morph"]="white"
        df.to_parquet(path,index=False)
        receipt=self.groups[0]/"original_photo_group_receipt.json"
        record=json.loads(receipt.read_text())
        record["saved_file_sha256"][path.name]=ledger.sha(path)
        receipt.write_text(json.dumps(record))
        with self.assertRaisesRegex(ValueError,"classification or source-arm"):
            self.run_audit()


if __name__ == "__main__":
    unittest.main()
