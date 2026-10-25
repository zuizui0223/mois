"""Real source-only FCP5 original ROI group sums: zero hues, missingness and SHAs."""
from __future__ import annotations

import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from scripts.acquisition import mois_fcp_u100_nested5_cumulative_signed_photos as m
from scripts.acquisition.mois_completed_pixel_checkpoint import sha


class TestFCP5OriginalPhotoPartialCumulative(unittest.TestCase):
    def setUp(self):
        self.temp=TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.previous=self.root/"old.json"
        self.previous.write_text(json.dumps({
            "schema":m.PRIOR_SCHEME,
            "exact_global_unique_original_photo_ID_attempts":1009032,
            "exact_global_automatic_four_hue_photo_ROIs":445805,
        }),encoding="utf-8")
        self.output=self.root/"new"
        self.docs=[]

    def group(self,g,*,classified=500,failed=4):
        n=min(1024,m.NEW-g*1024)
        extra=n-classified-failed
        self.assertGreaterEqual(extra,0)
        units=list(range(g*4,min((g+1)*4,m.WORK_UNITS)))
        row={
            "schema":m.GROUP_SCHEMA,
            "frozen_original_model_run":m.MODEL_RUN,
            "source_SHA256":m.SOURCE_MANIFEST_SHA,
            "source_original_FCP_U100_first5_nominal_taxa":20242,
            "source_101210_first5_original_photos":101210,
            "source_80968_new_nested_rank2to5_photo_candidates":80968,
            "prior_complete_original_photo_IDs_reused_without_inference":66312,
            "prior_global_actual_photo_IDs_SHA_audited_before_new_selection":1009032,
            "all_source_genuinely_new_original_photo_attempts":34898,
            "original_worker_group":g,
            "original_completed_source_work_units":units,
            "original_terminal_source_receipts_SHA256":[{
                "unit":u,"receipt_sha256":"a"*64,
                "photos_sha256":"b"*64,"terminal_sha256":"c"*64} for u in units],
            "new_actual_original_photo_ID_attempts":n,
            "new_actual_four_hue_classified":classified,
            "new_original_image_access_failed":failed,
            "new_original_acquired_hue_unclassifiable":extra,
            "new_distinct_original_observations":n,
            "new_original_hue_counts":{
                "white":min(200,classified),"yellow_orange":0,
                "red_pink":classified-min(200,classified),
                "blue_purple":0,"mixed_uncertain":n-classified},
            "distinct_original_per_photo_joint_probabilities":17,
            "frozen_original_FCP_yolo_efficientSAM_roi_v4_method_attested":True,
            "per_original_photo_source_inclusion_probability_preserved":True,
            "old_photo_ID_reuses_not_relabelled":True,
            "unknown_technical_ROIs_never_cast_to_white":True,
            "original_image_bytes_downloaded_OR_colour_ROI_inferred_by_sealer":0,
            "botanical_wild_native_open_flower_status_independently_verified":False,
            "global_photo_observation_probability_complete_cohort":False,
            "worldwide_unbiased_white_flora_frequency_estimated":False,
            "model_group_wall_seconds":3100.0,
            "original_photo_prefetch_stats":{
                "download_calls":n,"consumed":n,"permanent_failures":failed},
            "saved_file_SHA256":{
                "original_photo_ROIs.parquet":"d"*64,
                "original_photo_hue_counts.csv":"e"*64,
                "original_photo_source_coverage.json":"f"*64,
            },
        }
        file=self.root/f"group-{g}.json"
        file.write_text(json.dumps(row),encoding="utf-8")
        self.docs.append(file)
        return file

    def accumulate(self,docs=None):
        return m.aggregate(
            docs if docs is not None else self.docs,
            self.previous,self.output,previous_sha=sha(self.previous))

    def test_two_signed_source_group_hue_counts_and_global_total(self):
        a=self.group(0,classified=500)
        b=self.group(3,classified=550,failed=5)
        d=self.accumulate([b,a])
        self.assertEqual(d["new_finished_original_photo_group_IDs"],[0,3])
        self.assertEqual(d["new_genuine_original_photo_pixel_attempts"],2048)
        self.assertEqual(d["new_genuine_automatic_four_hue_ROIs"],1050)
        self.assertEqual(d["new_original_photo_acquisition_failed"],9)
        self.assertEqual(d["new_acquired_photo_hue_unclassifiable"],989)
        self.assertEqual(d["exact_precolour_disjoint_source_global_original_photo_attempts"],1011080)
        self.assertEqual(d["exact_precolour_disjoint_source_global_four_hue_ROIs"],446855)
        self.assertEqual(d["new_original_photo_morph_counts"]["yellow_orange"],0)
        self.assertEqual(d["new_pending_model_worker_groups"][0],1)
        self.assertTrue(d["counts_are_partial_model_completions_not_a_complete_probability_cohort"])
        self.assertTrue(d["new_source_34898_photo_IDs_disjoint_from_1009032_prior_by_frozen_manifest"])

    def test_last_group_has_only_82_real_original_photos(self):
        a=self.group(34,classified=40,failed=0)
        d=self.accumulate([a])
        self.assertEqual(d["new_genuine_original_photo_pixel_attempts"],82)
        self.assertEqual(d["new_genuine_automatic_four_hue_ROIs"],40)

    def test_missing_original_image_SHA_proof_cannot_publish(self):
        a=self.group(0)
        d=json.loads(a.read_text())
        d["original_terminal_source_receipts_SHA256"][0]["photos_sha256"]="a"
        a.write_text(json.dumps(d))
        with self.assertRaisesRegex(ValueError,"evidence absent"):
            self.accumulate()
        self.assertFalse(self.output.exists())

    def test_colourless_unknown_must_remain_technical_nonresponse(self):
        a=self.group(0)
        d=json.loads(a.read_text())
        d["new_original_hue_counts"]["white"]+=1
        a.write_text(json.dumps(d))
        with self.assertRaisesRegex(ValueError,"Zero-occurrence"):
            self.accumulate()

    def test_synthetic_old_global_cannot_be_forced_to_new_count(self):
        self.group(0)
        d=json.loads(self.previous.read_text())
        d["exact_global_unique_original_photo_ID_attempts"]=50000000
        self.previous.write_text(json.dumps(d))
        with self.assertRaisesRegex(ValueError,"not authenticated"):
            self.accumulate()

    def test_duplicate_model_worker_group_not_silently_added(self):
        a=self.group(0)
        with self.assertRaisesRegex(ValueError,"Duplicate signed"):
            self.accumulate([a,a])

    def test_cumulative_artifacts_are_content_deterministic_and_immutable(self):
        self.group(0)
        first=self.accumulate()
        self.assertTrue((self.output/"mois_FCP5_partial_genuine_original_pixel_cumulative.json").exists())
        with self.assertRaisesRegex(ValueError,"immutable"):
            self.accumulate()
        with TemporaryDirectory() as temp:
            second=m.aggregate(self.docs,self.previous,Path(temp)/"other",
                               previous_sha=sha(self.previous))
        self.assertEqual(first["original_group_receipt_set_SHA256"],
                         second["original_group_receipt_set_SHA256"])


if __name__=="__main__":
    unittest.main()
