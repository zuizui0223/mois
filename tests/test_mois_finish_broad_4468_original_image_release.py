"""Tests for photo-by-photo immutable broad all16 original RGB release recovery."""
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import unittest

import pandas as pd

from scripts.acquisition import mois_finish_broad_4468_original_image_release as release
from scripts.acquisition.mois_completed_pixel_checkpoint import sha
from scripts.analysis.mois_reassemble_original_new4096_measurements import METHOD


class VerifiedFullOriginalColourRelease(unittest.TestCase):
    def setUp(self):
        self.tmp=TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.finished=self.root/"finished"
        self.finished.mkdir()
        self.manifest=self.root/"manifest.json"
        self.crossarm=self.root/"crossarm.json"
        self.groups=[self.root/f"group{i}" for i in range(5)]
        self.source_arm="global93m_tier06_all16_first4of512_broad_warm_R1"
        self.prob=0.00005
        self.manifest.write_text(json.dumps({
            "schema":"mois_full_source_observation_sampling_cohort_v1",
            "source_tier":6,
            "source_group":"all16",
            "admitted_bucket_count":4,
            "new_photo_attempts_planned":10,
            "source_arm":self.source_arm,
            "selection_probability":self.prob,
            "not_probability_complete_global_source":True,
            "colour_measured":False,
        }))
        self.rows=[]
        groups=[]
        for i,folder in enumerate(self.groups):
            folder.mkdir()
            rows=[]
            for j in range(2):
                p=str(80000+i*2+j)
                is_colour=(j==0)
                is_failed=(i==4 and j==1)
                rows.append({
                    "photo_id":p,
                    "measurement_id":"M"+p,
                    "observation_id":"obs"+p,
                    "inat_taxon_id":str(1+i),
                    "source_arm":self.source_arm,
                    "selection_probability":str(self.prob),
                    "morph":"white" if is_colour else "mixed_uncertain",
                    "measurement_status":("classified_four_state_morph" if is_colour
                                          else "image_acquisition_failed" if is_failed
                                          else "not_evaluable_roi_or_flip_gate"),
                    "image_sha256":"" if is_failed else "e"*64,
                    "latitude":"35" if j==0 else "-3",
                    "observed_on":"2025-06-14" if j==0 else "2025-02-28",
                })
            self.rows+=rows
            frame=pd.DataFrame(rows)
            parq=folder/"original_photo_ROIs.parquet"
            csv=folder/"group_original_photo_hue_counts.csv"
            frame.to_parquet(parq,index=False)
            csv.write_text("morph,original_photo_attempts\nwhite,1\nmixed_uncertain,1\n")
            group_receipt=folder/"original_photo_group_receipt.json"
            group_receipt.write_text(json.dumps({
                "schema":"mois_tier06_all16_append_only_original_RGB_group_v1",
                "original_acquisition_workflow_run":release.RUN,
                "source_original_photo_manifest_sha256":sha(self.manifest),
                "all_source_photo_and_original_terminal_SHA_verified":True,
                "frozen_original_FCP_model_and_output_verified":True,
                "image_bytes_reacquired_in_sealing":False,
                "original_source_warm_group":i,
                "saved_file_sha256":{parq.name:sha(parq),csv.name:sha(csv)},
            }))
            groups.append({
                "group":i,"photo_attempts":2,"four_hue_classified":1,
                "acquisition_failed":1 if i==4 else 0,
                "original_worker_receipt_sha256":sha(group_receipt),
                "original_photo_parquet_sha256":sha(parq),
                "original_hue_count_csv_sha256":sha(csv),
            })
        full=self.finished/"measured_photos.parquet"
        pd.DataFrame(self.rows).to_parquet(full,index=False)
        (self.finished/"source.json").write_bytes(self.manifest.read_bytes())
        (self.finished/"receipt.json").write_text(json.dumps({
            "schema":"mois_global93m_observation_baseline_real_RGB_photo_results_v1",
            "attempted":10,"classifiable":5,"acquisition_failed":1,
            "unclassifiable_acquired":4,"downloaded":9,
            "measured_photo_parquet_sha256":sha(full),
            "source_plan_sha256":sha(self.manifest),
            "frozen_original_ROI_method":METHOD,
            "raw_licensed_photo_pixels_redistributed":False,
        }))
        self.crossarm.write_text(json.dumps({
            "schema":"mois_original_pixel_immutable_cross_arm_group_registry_v1",
            "main_checkpoint_run":release.MAIN_CHECKPOINT,
            "main_checkpoint_photo_attempts":100,
            "main_checkpoint_four_hue_classified":70,
            "previous_tier06_source_photo_attempts":3,
            "previous_tier06_four_hue_classified":2,
            "broad_completed_photo_attempts":10,
            "broad_completed_four_hue_classified":5,
            "broad_pending_photo_attempts":0,
            "original_photo_attempts_distinct_across_all_verified_arms":113,
            "four_hue_classifications_across_all_verified_arms":77,
            "original_photo_or_observation_ID_cross_arm_overlap":0,
            "all_original_photo_per_partition_SHA_receipts_verified":True,
            "broad_sealed_groups":groups,
        }))

    def run_recovery(self, group_files=None):
        return release.validate_collector_against_append_groups(
            self.finished,self.manifest,self.groups if group_files is None else group_files,self.crossarm,
            source_sha=sha(self.manifest),
            expected_original=10,expected_classified=5,
            expected_acq_failed=1,expected_unknown=4,expected_groups=5,
            current_run=release.RUN,cross_arm_total=113,cross_arm_hues=77,
            main_run=release.MAIN_CHECKPOINT,main_photos=100,main_classified=70,
            older_tier06_photos=3,older_tier06_classified=2)

    def test_source_and_original_image_recovery_exact_same_photo_rows(self):
        proof,descriptive=self.run_recovery()
        self.assertEqual(proof["original_photo_attempts"],10)
        self.assertEqual(proof["genuine_four_hue_classifications"],5)
        self.assertEqual(proof["acquisition_failed"],1)
        self.assertTrue(proof["same_source_five_original_fcp_worker_groups_photo_by_photo_verified"])
        self.assertEqual(descriptive["source_taxon_ids_attempted"],5)
        self.assertEqual(descriptive["source_taxon_ids_with_automatic_four_hue"],5)
        self.assertEqual(descriptive["observed_month_by_automatic_ROI_morph"]["06"]["white"],5)
        self.assertEqual(descriptive["latitude_band_by_automatic_ROI_morph"]["tropical"]["mixed_uncertain"],5)
        self.assertFalse(descriptive["flowering_anthesis_confirmed"])
        self.assertEqual(descriptive["technical_acquisition_failed_photo_ids"],1)

    def test_wrong_original_worker_parquet_is_rejected(self):
        f=self.groups[0]/"original_photo_ROIs.parquet"
        f.write_bytes(f.read_bytes()+b"corrupted")
        with self.assertRaisesRegex(ValueError,"altered"):
            self.run_recovery()

    def test_same_source_sha_but_swapped_colour_fails_exact_row_parity(self):
        f=self.groups[0]/"original_photo_ROIs.parquet"
        x=pd.read_parquet(f)
        x.loc[0,"image_sha256"]="d"*64
        x.to_parquet(f,index=False)
        receipt=self.groups[0]/"original_photo_group_receipt.json"
        z=json.loads(receipt.read_text())
        z["saved_file_sha256"][f.name]=sha(f)
        receipt.write_text(json.dumps(z))
        cross=json.loads(self.crossarm.read_text())
        cross["broad_sealed_groups"][0]["original_photo_parquet_sha256"]=sha(f)
        cross["broad_sealed_groups"][0]["original_worker_receipt_sha256"]=sha(receipt)
        self.crossarm.write_text(json.dumps(cross))
        with self.assertRaises(AssertionError):
            self.run_recovery()

    def test_original_full_source_independent_record_missing_fails(self):
        (self.finished/"source.json").write_text("{}")
        with self.assertRaisesRegex(ValueError,"Collector source metadata"):
            self.run_recovery()

    def test_one_missing_worker_group_is_not_accepted_as_full_study(self):
        with self.assertRaisesRegex(ValueError,"Not every original model worker"):
            self.run_recovery(group_files=self.groups[:4])

    def test_original_acquisition_failed_not_petal_colour(self):
        f=self.finished/"measured_photos.parquet"
        x=pd.read_parquet(f)
        x.loc[9,"morph"]="white"
        x.to_parquet(f,index=False)
        rec=self.finished/"receipt.json"
        r=json.loads(rec.read_text())
        r["measured_photo_parquet_sha256"]=sha(f)
        rec.write_text(json.dumps(r))
        with self.assertRaisesRegex(ValueError,"unknown/white/colour status"):
            self.run_recovery()


if __name__=="__main__":
    unittest.main()
