"""Verify actual novel wave2 original photo+observation nonoverlap vs signed old sampling arms."""
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import pandas as pd

from scripts.acquisition import mois_crossarm_tier06_wave2_new_original_ROI_ledger as reg
from scripts.acquisition.mois_completed_pixel_checkpoint import sha, HUES


class Wave2TrueOriginalCrossArmTests(unittest.TestCase):
    def setUp(self):
        self.tmp=TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.prior=self.root/"prior.json"
        self.source=self.root/"source.json"
        self.early=self.root/"early.parquet"
        self.broad=self.root/"broad.parquet"
        self.delta=self.root/"delta.parquet"
        self.output=self.root/"audit.json"
        self.groups=[self.root/"group0",self.root/"group1"]
        self.sample_p=0.0001
        for filename,rows in [
            (self.early,[("1","early-a"),("2","early-b")]),
            (self.broad,[("3","broad-a"),("4","broad-b"),("5","broad-c")]),
            (self.delta,[("6","delta-a"),("11","delta-b")]),
        ]:
            pd.DataFrame(rows,columns=["photo_id","observation_id"]).to_parquet(filename,index=False)
        self.prior.write_text(json.dumps({
            "schema":reg.PRIOR_SCHEMA,
            "prior_historical_main_checkpoint_run":reg.PRIOR_RUN,
            "full_authenticated_distinct_original_photo_attempts":15,
            "full_authenticated_four_hue_original_image_classifications":8,
            "photo_ID_overlap_new_p0032_vs_other_original_sampling_arms":0,
            "observation_ID_overlap_new_p0032_vs_other_original_sampling_arms":0,
            "separately_source_arm_audited":True,
            "signed_old_tier06_first_group_parquet_SHA256":sha(self.early),
            "signed_old_tier06_all16_first4_parquet_SHA256":sha(self.broad),
            "new_p0032_original_photo_delta_parquet_SHA256":sha(self.delta),
        }))
        self.source.write_text(json.dumps({
            "schema":"mois_full_source_observation_sampling_cohort_v1",
            "full_frame_photo_ids":93_209_351,
            "complete_tier06_source_photo_IDs":597_070,
            "source_group":"all16","source_tier":6,"source_arm":reg.SOURCE_ARM,
            "new_photo_attempts_planned":4,
            "previous_broad_original_photo_attempts":3,
            "previous_broad_source_parquet_sha256":sha(self.broad),
            "metadata_original_photos_in_admitted_buckets":9306,
            "previously_measured_photo_ids_skipped":541,
            "compute_bucket_count":512,"admitted_buckets":list(range(4,12)),
            "prior_photo_and_observation_exclusion_verified":True,
            "not_probability_complete_global_source":True,
            "colour_measured":False,
            "selection_probability":self.sample_p,
            "group_ids":[0,1],"work_units":8,
        }))
        source_sha=sha(self.source)
        for i,folder in enumerate(self.groups):
            folder.mkdir()
            rows=pd.DataFrame({
                "photo_id":[str(7+i*2),str(8+i*2)],
                "observation_id":["shared-new" if i==0 else "shared-new",
                                  f"second-observation-{i}"],
                "morph":["white","mixed_uncertain"],
                "measurement_status":["classified_four_state_morph",
                                      "not_evaluable_roi_or_flip_gate"],
                "selection_probability":[self.sample_p]*2,
                "source_arm":[reg.SOURCE_ARM]*2,
                "inat_taxon_id":["100","101"],
            })
            parq=folder/"original_photo_ROIs.parquet"
            rows.to_parquet(parq,index=False)
            coverage=folder/"original_photo_source_coverage.json"
            coverage.write_text(json.dumps({
                "schema":"mois_wave2_real_photo_source_taxonomy_time_latitude_coverage_v1",
                "source_taxon_ids_are_not_species_validated":True,
                "individual_photo_open_flowers_independently_verified":False,
                "original_photo_attempts":2,
                "distinct_nominal_inat_taxon_ids_four_hue":1,
            }))
            hues=folder/"group_original_photo_hue_counts.csv"
            hues.write_text("morph,original_photo_attempts\nwhite,1\nmixed_uncertain,1\n")
            (folder/"original_photo_group_receipt.json").write_text(json.dumps({
                "schema":reg.R2_SCHEMA,
                "original_source_warm_group":i,
                "original_acquisition_workflow_run":reg.R2_RUN,
                "source_original_photo_manifest_sha256":source_sha,
                "complete_expected_source_photo_attempts":4,
                "all_source_photo_and_original_terminal_SHA_verified":True,
                "frozen_original_FCP_model_and_output_verified":True,
                "image_bytes_reacquired_in_sealing":False,
                "ROI_or_hue_reinference_in_sealing":False,
                "globally_probability_complete":False,
                "original_work_units":list(range(i*4,i*4+4)),
                "distinct_genuinely_completed_photo_ids":2,
                "four_hue_model_classified":1,
                "distinct_original_observation_ids":2,
                "acquisition_failed":0,
                "acquired_but_unclassifiable":1,
                "hue_counts":{**{h:1 if h=="white" else 0 for h in HUES},
                              "mixed_uncertain":1},
                "saved_file_sha256":{
                    parq.name:sha(parq),hues.name:sha(hues),coverage.name:sha(coverage)},
            }))

    def run_verify(self,folders=None):
        return reg.validate(
            self.prior,self.source,self.early,self.broad,self.delta,
            self.groups if folders is None else folders,self.output,
            prior_total=15,prior_hues=8,expected_source_sha=sha(self.source),
            expected_new=4,prior_early=(2,1),prior_broad=(3,2),
            prior_delta=(2,1))

    def test_two_complete_groups_and_same_arm_observation_crossing_group_boundary(self):
        d=self.run_verify()
        self.assertEqual(d["wave2_groups_completed"],[0,1])
        self.assertEqual(d["wave2_original_photo_attempts"],4)
        self.assertEqual(d["wave2_original_four_hue_classified"],2)
        self.assertEqual(d["wave2_original_distinct_observations"],3)
        self.assertEqual(d["wave2_same_original_observations_spanning_compute_groups"],1)
        self.assertEqual(d["full_distinct_original_photo_attempts"],19)
        self.assertEqual(d["full_automatically_classified_four_hue_original_photo_ROIs"],10)
        self.assertTrue(d["prior_separate_arms_and_new_wave2_no_duplicate_photo_or_observation"])
        self.assertFalse(d["photo_hues_not_botanically_verified_as_flowering"] is False)

    def test_one_sealed_group_can_be_published_without_counting_another_queued_group(self):
        d=self.run_verify([self.groups[0]])
        self.assertEqual(d["wave2_original_photo_attempts"],2)
        self.assertEqual(d["wave2_remaining_original_photo_attempts"],2)
        self.assertEqual(d["full_distinct_original_photo_attempts"],17)

    def test_same_photo_id_as_previous_tier06_arm_is_rejected(self):
        path=self.groups[0]/"original_photo_ROIs.parquet"
        tab=pd.read_parquet(path)
        tab.loc[0,"photo_id"]="3"
        tab.to_parquet(path,index=False)
        self.update_source_group_hash(0,path)
        with self.assertRaisesRegex(ValueError,"collision"):
            self.run_verify()

    def test_same_original_observation_as_p0032_delta_is_rejected(self):
        path=self.groups[0]/"original_photo_ROIs.parquet"
        tab=pd.read_parquet(path)
        tab.loc[0,"observation_id"]="delta-a"
        tab.to_parquet(path,index=False)
        self.update_source_group_hash(0,path)
        with self.assertRaisesRegex(ValueError,"collision"):
            self.run_verify()

    def test_technical_unknown_relabelled_as_white_is_rejected(self):
        path=self.groups[0]/"original_photo_ROIs.parquet"
        tab=pd.read_parquet(path)
        tab.loc[1,"morph"]="white"
        tab.to_parquet(path,index=False)
        self.update_source_group_hash(0,path)
        with self.assertRaisesRegex(ValueError,"changed"):
            self.run_verify()

    def test_forged_source_SHA_or_tampered_QA_refused(self):
        cov=self.groups[1]/"original_photo_source_coverage.json"
        cov.write_text(cov.read_text()+"\n")
        with self.assertRaisesRegex(ValueError,"byte SHA"):
            self.run_verify()
        self.assertFalse(self.output.exists())

    def test_second_execution_cannot_rewrite_signed_ledger(self):
        self.run_verify()
        with self.assertRaisesRegex(ValueError,"Never mutate"):
            self.run_verify()

    def update_source_group_hash(self,index,path):
        name=self.groups[index]/"original_photo_group_receipt.json"
        content=json.loads(name.read_text())
        content["saved_file_sha256"][path.name]=sha(path)
        name.write_text(json.dumps(content))


if __name__=="__main__":
    unittest.main()
