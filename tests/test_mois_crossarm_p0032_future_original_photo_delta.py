"""Original p0032 signed worker 133+ cross-arm incremental photo ID/colour tests."""
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import pandas as pd

from scripts.acquisition import mois_crossarm_p0032_future_original_photo_delta as m
from scripts.acquisition.mois_completed_pixel_checkpoint import sha, HUES


class P0032FutureCrossarmTests(unittest.TestCase):
    def setUp(self):
        self.tmp=TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.prior=self.root/"prior.json"
        self.early=self.root/"early.parquet"
        self.broad=self.root/"broad.parquet"
        self.delta=self.root/"delta.parquet"
        self.wave2=self.root/"wave2"
        self.new=self.root/"group133"
        self.output=self.root/"registry.json"
        self.wave2.mkdir()
        self.new.mkdir()
        for path,rows in [
            (self.early,[("1","ear-one"),("2","ear-two")]),
            (self.broad,[("3","broad-a"),("4","broad-b"),("5","broad-c")]),
            (self.delta,[("6","delta-a"),("7","delta-b")]),
        ]:
            pd.DataFrame(rows,columns=["photo_id","observation_id"]).to_parquet(path,index=False)
        p=self.wave2/"original_photo_ROIs.parquet"
        pd.DataFrame([
            ("8","w2-a"),("9","w2-b")
        ],columns=["photo_id","observation_id"]).to_parquet(p,index=False)
        r=self.wave2/"original_photo_group_receipt.json"
        r.write_text(json.dumps({"original_source_warm_group":0}))
        self.prior.write_text(json.dumps({
            "schema":m.OLD_SCHEMA,
            "full_distinct_original_photo_attempts":15,
            "full_automatically_classified_four_hue_original_photo_ROIs":8,
            "prior_separate_arms_and_new_wave2_no_duplicate_photo_or_observation":True,
            "original_images_fetched_and_ROI_inferred_during_this_ledger":0,
            "all_wave2_group_photo_and_source_coverage_SHA_verified":True,
            "wave2_groups_completed":[0],
            "wave2_original_photo_attempts":2,
            "per_group_original_photo_source_audit":[{
                "group":0,
                "group_source_SHA256":sha(r),
                "group_photo_Parquet_SHA256":sha(p),
                "original_photo_attempts":2,
            }],
        }))
        self.create_new_group()

    def create_new_group(self):
        rows=pd.DataFrame([
            {"photo_id":"10","observation_id":"new-p32",
             "morph":"white","measurement_status":"classified_four_state_morph",
             "source_arm":"p0032-source","selection_probability":0.0032},
            {"photo_id":"11","observation_id":"delta-b",
             "morph":"mixed_uncertain",
             "measurement_status":"not_evaluable_roi_or_flip_gate",
             "source_arm":"p0032-source","selection_probability":0.0032},
        ])
        photo=self.new/"original_photo_ROIs.parquet"
        rows.to_parquet(photo,index=False)
        counts=self.new/"original_photo_hue_counts.csv"
        counts.write_text("morph,true_original_photo_attempts\nwhite,1\nmixed_uncertain,1\n")
        receipt=self.new/"p0032_original_group_receipt.json"
        receipt.write_text(json.dumps({
            "schema":m.P32_SCHEMA,
            "signed_original_source_run":m.P32_RUN,
            "source_manifest_sha256":m.P32_SOURCE_SHA,
            "original_source_photo_opportunities":148071,
            "original_frozen_source_worker_group":133,
            "original_completed_source_units":list(range(532,536)),
            "frozen_original_fcp_roi_v4_attested":True,
            "original_photo_and_terminal_SHA_certified":True,
            "new_S3_photo_downloads_or_model_reinferences":0,
            "other_prior_original_worker_ZIPs_downloaded":0,
            "original_true_photo_attempts":2,
            "original_FCP_roi_four_hue_classified":1,
            "original_photo_image_acquisition_failed":0,
            "original_acquired_unclassifiable":1,
            "original_distinct_observations":2,
            "hue_counts":{**{h:1 if h=="white" else 0 for h in HUES},
                          "mixed_uncertain":1},
            "actual_output_sha256":{
                photo.name:sha(photo),counts.name:sha(counts)},
        }))
        (self.new/"original_worker_artifact_zip_identity.json").write_text(json.dumps({
            "original_p0032_group":133,
            "source_run":m.P32_RUN,
            "original_source_SHA256":m.P32_SOURCE_SHA,
            "original_terminal_group_receipt_SHA256":sha(receipt),
            "artifact_name":"mois-93m-p0032-RGB-group-133",
            "artifact_zip_sha256":"sha256:"+"a"*64,
        }))

    def audit(self):
        return m.validate(
            self.prior,self.delta,self.early,self.broad,
            [self.wave2],[self.new],self.output,
            prior_total=15,prior_hues=8,
            p32_delta_n=2,expected_early=2,expected_broad=3,
        )

    def update_photo_sha(self):
        photo=self.new/"original_photo_ROIs.parquet"
        rp=self.new/"p0032_original_group_receipt.json"
        x=json.loads(rp.read_text())
        x["actual_output_sha256"][photo.name]=sha(photo)
        rp.write_text(json.dumps(x))
        a=self.new/"original_worker_artifact_zip_identity.json"
        y=json.loads(a.read_text())
        y["original_terminal_group_receipt_SHA256"]=sha(rp)
        a.write_text(json.dumps(y))

    def test_success_and_same_observation_within_same_p0032_source_arm_allowed(self):
        r=self.audit()
        self.assertEqual(r["new_p0032_postdelta_original_worker_groups"],[133])
        self.assertEqual(r["new_p0032_postdelta_authentic_original_photo_attempts"],2)
        self.assertEqual(r["new_p0032_postdelta_four_hue_classified"],1)
        self.assertEqual(r["historical_and_new_genuine_original_photo_attempts"],17)
        self.assertEqual(r["historical_and_new_four_hue_classifications"],9)
        self.assertEqual(r["new_original_observation_ID_overlap_with_independent_tier06_arms"],0)
        self.assertFalse(r["new_image_downloads_or_ROI_colour_inference_performed"])

    def test_same_photo_in_independent_wave2_source_is_rejected(self):
        p=self.new/"original_photo_ROIs.parquet"
        df=pd.read_parquet(p)
        df.loc[0,"photo_id"]="8"
        df.to_parquet(p,index=False)
        self.update_photo_sha()
        with self.assertRaisesRegex(ValueError,"clash"):
            self.audit()
        self.assertFalse(self.output.exists())

    def test_same_observation_in_independent_old_source_is_rejected(self):
        p=self.new/"original_photo_ROIs.parquet"
        df=pd.read_parquet(p)
        df.loc[0,"observation_id"]="broad-b"
        df.to_parquet(p,index=False)
        self.update_photo_sha()
        with self.assertRaisesRegex(ValueError,"clash"):
            self.audit()
        self.assertFalse(self.output.exists())

    def test_hue_relabelled_without_frozen_status_is_rejected(self):
        p=self.new/"original_photo_ROIs.parquet"
        df=pd.read_parquet(p)
        df.loc[1,"morph"]="white"
        df.to_parquet(p,index=False)
        self.update_photo_sha()
        with self.assertRaisesRegex(ValueError,"clash"):
            self.audit()

    def test_terminal_artifact_SHA_is_immutable(self):
        p=self.new/"original_worker_artifact_zip_identity.json"
        obj=json.loads(p.read_text())
        obj["artifact_zip_sha256"]="x"
        p.write_text(json.dumps(obj))
        with self.assertRaisesRegex(ValueError,"SHA"):
            self.audit()

    def test_duplicate_run_is_not_allowed_to_replace_existing_signed_union(self):
        self.audit()
        with self.assertRaisesRegex(ValueError,"cannot be overwritten"):
            self.audit()


if __name__=="__main__":
    unittest.main()
