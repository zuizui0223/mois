"""Photo-ID-immutable 93m balanced sample reuse and never-downloaded source bridge."""
from __future__ import annotations

import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import pandas as pd

from scripts.acquisition import mois_93m_balanced_prior_ROI_reuse_bridge as bridge
from scripts.acquisition.mois_completed_pixel_checkpoint import sha
from scripts.acquisition.mois_queue_from_full93m_target_tiers import WORKER_COLUMNS


class BalancingExistingOriginalPixelsTests(unittest.TestCase):
    def setUp(self):
        self.temp=TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.source=self.root/"fullsource.csv"
        self.manifest=self.root/"source.json"
        self.technical=self.root/"technical.csv"
        self.output=self.root/"new-photo-shards"
        self.priors=[]
        self.prior_counts={}
        specs=[
            ("main",[("1","obsA","white"),("2","obsM","mixed_uncertain")]),
            ("prior_group00",[("3","obsC","white")]),
            ("prior_R1",[("4","obsD","yellow_orange")]),
            ("prior_R2",[("5","obsE","mixed_uncertain")])
        ]
        for arm,items in specs:
            frame=pd.DataFrame([{
                "photo_id":pid, "observation_id":obs, "morph":hue,
                "measurement_status":"classified_four_state_morph" if hue!="mixed_uncertain"
                      else "not_evaluable_roi_or_flip_gate",
                "image_sha256":"d"*64
            } for pid,obs,hue in items])
            p=self.root/(arm+".parquet")
            rec=self.root/(arm+".json")
            frame.to_parquet(p,index=False)
            count=len(frame);hue=sum(z!="mixed_uncertain" for _,_,z in items)
            if arm=="main":
                rec.write_text(json.dumps({
                    "schema":"mois_sealed_actual_pixel_checkpoint_v1",
                    "combined_unique_photo_attempts":count,
                    "combined_four_hue_classified":hue,
                    "output_sha256":{p.name:sha(p)},
                }))
            else:
                rec.write_text(json.dumps({
                    "schema":"mois_global93m_observation_baseline_real_RGB_photo_results_v1",
                    "attempted":count,"classifiable":hue,
                    "measured_photo_parquet_sha256":sha(p),
                }))
            self.priors.append((arm,p,rec))
            self.prior_counts[arm]=(count,hue)
        self.technical.write_text(
            "photo_id,observation_id,acquired,status\n8,obsT,true,technical_cohort\n"
        )
        self.rows=[]
        fixtures=[
            ("1","obsA","supported",2,3),
            ("6","obsB","supported",4,2),
            ("7","obsM","supported",1,1),  # observation already measured, distinct photo -> STILL selected
            ("8","obsT","supported",1,1),  # technical benchmark -> no re-download
            ("9","obsI","unsupported_original_extension",1,1),
        ]
        for pid,obs,status,m,n in fixtures:
            self.rows.append({
                "photo_id":pid,"observation_id":obs,
                "inat_taxon_id":"1","species":"Campanula test",
                "observer_id":"photographer","observed_on":"2025-06-01",
                "latitude":"35","longitude":"140",
                "positional_accuracy":"5","photo_license":"cc-by",
                "photo_url_large":f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.jpg",
                "photo_source_disposition":status,
                "source_arm":bridge.SOURCE_ARM,
                "photos_in_original_observation":str(m),
                "observations_in_original_source_stratum":str(n),
                "joint_photo_selection_probability":format(.0064/m/n,".17g"),
            })
        with self.source.open("w",encoding="utf-8",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=self.rows[0].keys(),lineterminator="\n")
            writer.writeheader();writer.writerows(self.rows)
        self.manifest.write_text(json.dumps({
            "schema":"mois_full_source_observation_sampling_cohort_v1",
            "sampling_schema":bridge.SOURCE_SCHEMA,
            "source_arm":bridge.SOURCE_ARM,
            "all_upstream_original_photo_IDs":93_209_351,
            "original_tier06_source_photo_IDs_accounted":597070,
            "sampled_distinct_original_observations":5,
            "sampled_per_stratum_K":1,
            "source_sample_decided_without_photo_pixels":True,
            "selection_probability_is_per_photo_in_signed_shards":True,
            "original_source_sample_full_table_SHA256":sha(self.source),
        }))

    def bridge(self,**kwargs):
        return bridge.build(
            self.manifest,self.source,self.priors,self.technical,self.output,
            expected_sample=5,
            expected_sample_sha=sha(self.manifest),
            expected_selected_sha=sha(self.source),
            expected_prior_counts=self.prior_counts,
            expected_technical_sha=sha(self.technical),**kwargs)

    def test_every_selected_photo_disposition_conserved_and_old_pixels_not_refetched(self):
        r=self.bridge()
        self.assertEqual(r["all_balanced_selected_original_observations"],5)
        self.assertEqual(r["new_photo_attempts_planned"],2)
        self.assertEqual(r["reused_original_photo_ids"],1)
        self.assertEqual(r["technical_holdout_images_not_redownloaded"],1)
        self.assertEqual(r["unsupported_original_formats_technical_unknown"],1)
        self.assertEqual(r["previously_measured_original_observations_with_different_sample_photo"],1)
        self.assertEqual(r["work_units"],1)
        self.assertEqual(r["group_ids"],[0])
        self.assertFalse(r["colour_measured"])
        self.assertTrue(r["original_photo_source_probabilities_unchanged_by_prior_pixel_reuse"])
        self.assertEqual(len(r["prior_original_pixel_arms"]),4)
        with (self.output/"shards/shard-000.csv").open(newline="") as f:
            kept=list(csv.DictReader(f))
        self.assertEqual(tuple(kept[0].keys()),WORKER_COLUMNS)
        self.assertEqual({v["photo_id"] for v in kept},{"6","7"})
        self.assertEqual(set(v["observation_id"] for v in kept),{"obsB","obsM"})
        self.assertEqual(
            {v["photo_id"]:float(v["selection_probability"]) for v in kept},
            {"6":.0064/8,"7":.0064})
        self.assertEqual(sha(self.output/"shards/shard-000.csv"),
                         r["photo_shard_sha256"]["0"])
        with (self.output/"already_measured_original_photo_reuse.csv").open(newline="") as f:
            old=list(csv.DictReader(f))
        self.assertEqual(len(old),1)
        self.assertEqual(old[0]["photo_id"],"1")
        self.assertEqual(old[0]["prior_automatic_ROI_morph"],"white")

    def test_old_sealed_cross_arm_photo_overlap_is_refused(self):
        path=self.priors[1][1]
        d=pd.read_parquet(path)
        d.loc[0,"photo_id"]="1"
        d.to_parquet(path,index=False)
        receipt=self.priors[1][2]
        r=json.loads(receipt.read_text())
        r["measured_photo_parquet_sha256"]=sha(path)
        receipt.write_text(json.dumps(r))
        with self.assertRaisesRegex(ValueError,"Overlapping original photo"):
            self.bridge()
        self.assertFalse(self.output.exists())

    def test_unchanged_selected_photo_probability_is_mandatory(self):
        d=list(self.rows)
        d[1]=dict(d[1],joint_photo_selection_probability="0.5")
        with self.source.open("w",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=d[0],lineterminator="\n")
            writer.writeheader();writer.writerows(d)
        manifest=json.loads(self.manifest.read_text())
        manifest["original_source_sample_full_table_SHA256"]=sha(self.source)
        self.manifest.write_text(json.dumps(manifest))
        with self.assertRaisesRegex(ValueError,"probability changed"):
            self.bridge()
        self.assertFalse(self.output.exists())

    def test_byte_tampering_of_original_93m_sample_is_refused(self):
        with self.source.open("ab") as f:
            f.write(b"fake-extra-source")
        with self.assertRaisesRegex(ValueError,"manifest or selected roster SHA changed"):
            self.bridge()

    def test_sample_is_immutable_after_success(self):
        self.bridge()
        with self.assertRaisesRegex(ValueError,"immutable"):
            self.bridge()


if __name__=="__main__":
    unittest.main()
