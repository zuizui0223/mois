"""K1 original image full 75-group collector source contract regression."""
from __future__ import annotations
import copy
import unittest

from scripts.acquisition.mois_full_index_real_original_photo_measure import (
    SCHEMA, source_prior_dedup_receipt,
)
from scripts.acquisition.mois_fcp_k1_append_completed_original_photo_group import (
    PARENT_SOURCE_SHA,SOURCE_ARM,SOURCE_PLANNED,
    SELECTED_SOURCE_OBSERVATIONS,REUSED_OLD_ORIGINALS,
)

class FullK1CompletedSourceCollectorTests(unittest.TestCase):
    def data(self):
        return dict(
            schema=SCHEMA,
            source_tier=6,source_group="balanced_taxon_month_lat_K1",
            source_arm=SOURCE_ARM,full_frame_photo_ids=93209351,
            source_precolour_balanced_manifest_SHA256=PARENT_SOURCE_SHA,
            new_photo_attempts_planned=SOURCE_PLANNED,
            all_balanced_selected_original_observations=SELECTED_SOURCE_OBSERVATIONS,
            reused_original_photo_ids=REUSED_OLD_ORIGINALS,
            technical_holdout_images_not_redownloaded=2,
            unsupported_original_formats_technical_unknown=0,
            prior_original_crossarm_photo_overlap=dict(
                unique_old_photo_IDs=929138,signed_arm_photo_rows=929138),
            prior_original_pixel_arms=[{"source_arm":name}
                                     for name in ("main","prior_group00","prior_R1","prior_R2")],
            selection_probability=None,sampling_nested_joint_probability_per_photo=True,
            original_photo_source_probabilities_unchanged_by_prior_pixel_reuse=True,
            known_previous_full_photo_and_observation_IDs_source_SHA_verified=True,
            previously_counted_photo_ID_duplicates_are_NOT_new_inferences=True,
            all_original_photo_ids_in_master_probability_frame_preserved=True,
            technical_unknown_not_other_colour=True,colour_measured=False,
            images_downloaded_by_this_source_bridge=0,
            new_colour_ROIs_inferred_by_this_source_bridge=0,
            photo_specific_anthesis_native_wild_not_independently_reviewed=True,
        )

    def test_actual_FCP_warmK1_selected_observation_sample_conservation(self):
        s=source_prior_dedup_receipt(self.data())
        self.assertEqual(s["new_original_photo_IDs_measured"],76753)
        self.assertEqual(s["matching_prior_measured_photos_reused"],3280)
        self.assertEqual(s["technical_photo_ID_holdouts_not_redownloaded"],2)
        self.assertEqual(s["historical_unique_source_photo_IDs_unchanged"],929138)
        self.assertEqual(76753+3280+2,80035)
        self.assertTrue(s["K1_original_joint_photo_selection_probability_per_individual_photo"])
        self.assertTrue(s["botanical_flower_native_status_unreviewed"])
        self.assertTrue(s["source_precolour_photo_dispositions_reconciled"])

    def test_never_assume_global_scalar_probability_or_old_old_source(self):
        s=self.data()
        s["selection_probability"]=0.01
        with self.assertRaisesRegex(ValueError,"denominator changed"):
            source_prior_dedup_receipt(s)
        s=self.data()
        s["source_precolour_balanced_manifest_SHA256"]="0"*64
        with self.assertRaisesRegex(ValueError,"denominator changed"):
            source_prior_dedup_receipt(s)

    def test_missing_reuse_OR_technical_nonresponse_must_fail(self):
        for attr in ("reused_original_photo_ids",
                     "technical_holdout_images_not_redownloaded",
                     "new_photo_attempts_planned",
                     "all_balanced_selected_original_observations"):
            s=self.data()
            s[attr]+=1
            with self.subTest(attr=attr),self.assertRaisesRegex(ValueError,"denominator changed"):
                source_prior_dedup_receipt(s)

    def test_non_original_precolour_source_cannot_masquerade_as_K1(self):
        for attr,value in (
            ("source_arm","made-up"),
            ("source_tier",3),
            ("full_frame_photo_ids",80035),
            ("original_photo_source_probabilities_unchanged_by_prior_pixel_reuse",False),
            ("colour_measured",True),
            ("prior_original_crossarm_photo_overlap",{"unique_old_photo_IDs":99}),
        ):
            s=self.data()
            s[attr]=value
            with self.subTest(attr=attr),self.assertRaisesRegex(ValueError,"denominator changed"):
                source_prior_dedup_receipt(s)

    def test_legacy_source_previous_photo_reuse_path_untouched(self):
        r=source_prior_dedup_receipt({"previously_measured_photo_ids_skipped":3})
        self.assertEqual(r["matching_prior_measured_photos_reused"],3)
        self.assertEqual(r["source_deduplication_reporting_grain"],"prior_original_photo_id")


if __name__=="__main__":
    unittest.main()
