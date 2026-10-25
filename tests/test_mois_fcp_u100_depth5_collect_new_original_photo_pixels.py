"""No precolour leakage and original ROI-v4 terminal audit for next FCP 34898 wave."""
from __future__ import annotations

import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from scripts.acquisition import mois_fcp_u100_depth5_collect_new_original_photo_pixels as m
from scripts.acquisition.mois_completed_pixel_checkpoint import sha


class NewFCPNestedOriginalColourTerminalTests(unittest.TestCase):
    def setUp(self):
        self.temp=TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.source=self.root/"src"
        self.source.mkdir()
        self.output=self.root/"out"
        self.rows=[
            dict(photo_id="301",observation_id="O1",inat_taxon_id="1",
                 source_arm=m.ARM,morph="white",
                 measurement_status="classified_four_state_morph",
                 image_sha256="a"*64,selection_probability="0.05"),
            dict(photo_id="302",observation_id="O2",inat_taxon_id="2",
                 source_arm=m.ARM,morph="mixed_uncertain",
                 measurement_status="image_acquisition_failed",
                 image_sha256="",selection_probability="0.03"),
        ]
        self.audit={
            "wave_complete":True,"attempted":2,"planned_photo_attempts":2,
            "completed_units":[0],"classifiable":1,"acquisition_failed":1,
            "unclassifiable_acquired":0,
            "completed_unit_audit":[{"unit":0,"receipt_sha256":"a"*64}]
        }
        self.manifest=self.source/"manifest.json"
        self.write_manifest()

    def write_manifest(self):
        self.manifest.write_text(json.dumps(dict(
            schema=m.OUT_SCHEMA,original_sampling_schema=m.SOURCE_STAGE_SCHEMA,
            full_frame_photo_ids=93209351,source_arm=m.ARM,
            source_group="FCP_U100_nested5_rank2to5_only_NOMINAL_ELIGIBLE",
            colour_measured=False,new_photo_attempts_planned=2,
            previously_measurement_outcomes_reused_without_inference=3,
            original_first5_photo_source_IDs=5,
            original_incremental_rank2to5_candidate_photo_IDs=4,
            source_original_eligible_taxa=20242,
            original_photo_nested5_source_receipt_SHA256=
                "037e74dd236d04c44f0d84c7b9d0b3e436e72fe6ddd470cc315c1dc56b0e197f",
            previously_measured_historical_K1_FCP_1009032_photo_IDs_accounted=1009032,
            selection_probability=None,
            individual_original_photo_nested_prefix_inclusion_probabilities_preserved=True,
            pending_image_work_is_only_source_photo_IDs_not_previously_inferred=True,
            old_unknown_ROI_not_treated_as_flower_colour=True,
            new_image_downloads_or_ROI_inferences_by_this_source_step=0,
            source_metadata_wild_native_bloom_not_independently_verified=True,
            technical_unknowns_preserved=True,
            work_units=1,worker_ids=[0],group_ids=[0],
            units_per_group=4,original_model_batch_size=256,
            full_old_K1_FCP_reuse_photo_count_conservation=5,
            photo_dispositions={
                "new_original_photo_ROI_to_infer":2,
                "prior_original_photo_exact_ID_reuse":3,
                "original_photo_unsupported_extension_unknown":0},
            photo_shard_sha256={"0":"a"*64}
        )),encoding="utf-8")

    def invoke(self):
        with patch.object(m,"collect_wave",return_value=(self.rows,self.audit)):
            return m.run(self.source,self.root/"workers",self.output,
                source_sha=sha(self.manifest),expected_new=2,
                expected_reused=3,work_units=1,groups=1,
                total_first5=5,delta=4)

    def test_true_original_FCP5_new_images_only_with_unknown(self):
        d=self.invoke()
        self.assertEqual(d["genuinely_new_original_photo_pixels_attempted"],2)
        self.assertEqual(d["genuinely_new_automatic_four_hue_classified"],1)
        self.assertEqual(d["genuinely_new_image_access_failed"],1)
        self.assertEqual(d["previous_real_original_photo_ID_reuses_NOT_REINFERRED"],3)
        self.assertEqual(d["source_new_photo_Parquet_SHA256"],
                         sha(self.output/"new_FCP_U100_nested5_original_photo_ROIs.parquet"))
        self.assertTrue(d["original_photo_full_prefix5_5_over_N_over_M_selection_probability_retained"])
        self.assertEqual(d["original_image_downloads_or_ROI_reinferences_in_this_collector"],0)
        self.assertFalse(d["botanically_confirmed_wild_native_open_blooming_individuals"])

    def test_numerical_technical_unknown_not_mislabelled_as_white(self):
        self.rows[1]["morph"]="white"
        with self.assertRaisesRegex(ValueError,"terminal contradiction"):
            self.invoke()
        self.assertFalse(self.output.exists())

    def test_doubled_photo_ID_and_observation_are_rejected(self):
        self.rows[1]["photo_id"]=self.rows[0]["photo_id"]
        with self.assertRaisesRegex(ValueError,"duplicate"):
            self.invoke()
        self.assertFalse(self.output.exists())

    def test_source_probability_not_omitted(self):
        self.rows[0]["selection_probability"]="0"
        with self.assertRaisesRegex(ValueError,"inclusion probability absent"):
            self.invoke()
        self.assertFalse(self.output.exists())

    def test_missing_256_photo_unit_never_accepted_as_finished(self):
        self.audit["wave_complete"]=False
        with self.assertRaisesRegex(ValueError,"Missing finished"):
            self.invoke()

    def test_forged_original_source_SHA_fails_before_releasing_pixels(self):
        with self.assertRaisesRegex(ValueError,"SHA changed"):
            with patch.object(m,"collect_wave",return_value=(self.rows,self.audit)):
                m.run(self.source,self.root/"workers",self.output,
                    source_sha="0"*64,expected_new=2,
                    expected_reused=3,work_units=1,groups=1,total_first5=5,delta=4)

    def test_immutable_source_photo_ROIs_not_recomputed(self):
        self.invoke()
        with self.assertRaisesRegex(ValueError,"Never overwrite"):
            self.invoke()


if __name__=="__main__":
    unittest.main()
