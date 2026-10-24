"""True source-only original photo/observation-ID provenance diagnosis tests."""
from pathlib import Path
from tempfile import TemporaryDirectory
import csv
import json
import unittest

import pandas as pd

from scripts.acquisition import mois_diagnose_prior_photo_uuid_crosswalk as diag
from scripts.acquisition.mois_onepass_full93m_target_tiers import sha256


class OriginalPhotoIDCrosswalkDiagnosisTests(unittest.TestCase):
    def setUp(self):
        self.tmp=TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.base=Path(self.tmp.name)
        self.sample=self.base/"sample.csv"
        self.manifest=self.base/"manifest.json"
        self.output=self.base/"audit.json"
        self.uuid={
            i:f"12345678-1234-1234-1234-{i:012d}" for i in range(1,5)
        }
        with self.sample.open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=["photo_id","observation_id","inat_taxon_id"])
            w.writeheader()
            for i in range(1,5):
                w.writerow({"photo_id":str(i),"observation_id":self.uuid[i],
                            "inat_taxon_id":str(10+i)})
        self.manifest.write_text(json.dumps({
            "schema":"mois_full_source_observation_sampling_cohort_v1",
            "sampled_distinct_original_observations":4,
            "all_upstream_original_photo_IDs":93_209_351,
            "source_sample_decided_without_photo_pixels":True,
        }))
        arms={
            "main":[("1","998712"),("2",self.uuid[2])],
            "group00":[("3","nan")],
            "R1":[("4",self.uuid[4])],
            "R2":[("5","missing-prior-not-in-candidates")],
        }
        self.sources={}
        self.counts={}
        for arm,items in arms.items():
            p=self.base/f"{arm}.parquet"
            r=self.base/f"{arm}.json"
            pd.DataFrame(items,columns=["photo_id","observation_id"]).to_parquet(p,index=False)
            n=len(items)
            h=1 if arm in ("main","R1") else 0
            if arm=="main":
                r.write_text(json.dumps({
                    "schema":"mois_sealed_actual_pixel_checkpoint_v1",
                    "combined_unique_photo_attempts":n,
                    "combined_four_hue_classified":h,
                    "output_sha256":{p.name:sha256(p)},
                }))
            else:
                r.write_text(json.dumps({
                    "schema":"mois_global93m_observation_baseline_real_RGB_photo_results_v1",
                    "attempted":n,"classifiable":h,
                    "measured_photo_parquet_sha256":sha256(p),
                }))
            self.sources[arm]=(p,r)
            self.counts[arm]=(n,h)

    def execute(self):
        return diag.probe(self.sample,self.manifest,self.sources,self.output,
                          expected_source=4,expected_counts=self.counts,
                          expected_sample_sha=sha256(self.sample),
                          expected_manifest_sha=sha256(self.manifest))

    def test_distinguishes_numeric_vs_uuid_from_truly_missing_source_field(self):
        d=self.execute()
        self.assertEqual(d["schema"],diag.SCHEMA)
        self.assertEqual(d["total_per_arm_overlapping_photo_rows_including_any_cross_arm_duplicates"],4)
        main=d["per_arm_diagnosis"]["main"]
        self.assertEqual(main["photo_ID_overlap_with_new_probability_sample"],2)
        self.assertEqual(main["historical_actual_ID_disagreement"],1)
        self.assertEqual(main["exact_same_observation_ID"],1)
        self.assertEqual(main["identity_format_counts"]["different_numeric_versus_uuid"],1)
        self.assertEqual(d["per_arm_diagnosis"]["group00"]["historical_missing_observation_ID"],1)
        self.assertEqual(d["per_arm_diagnosis"]["R1"]["exact_same_observation_ID"],1)
        self.assertEqual(d["per_arm_diagnosis"]["R2"]["photo_ID_overlap_with_new_probability_sample"],0)
        self.assertEqual(d["hashed_format_examples"][0]["prior_observation_format"],"numeric")
        self.assertEqual(d["hashed_format_examples"][0]["source_observation_format"],"uuid")
        self.assertNotIn("original_photo_id", d["hashed_format_examples"][0])
        self.assertTrue(d["ambiguity_requires_canonical_observation_crosswalk_or_photo_only_technical_unknown"])
        self.assertFalse(d["no_original_photo_images_downloaded"] is False)

    def test_altered_old_pixel_parquet_bytes_refused(self):
        p=self.sources["main"][0]
        p.write_bytes(p.read_bytes()+b"tamper")
        with self.assertRaisesRegex(ValueError,"SHA changed"):
            self.execute()
        self.assertFalse(self.output.exists())

    def test_earlier_photo_frame_cannot_be_reissued_as_another_diagnosis(self):
        self.execute()
        with self.assertRaisesRegex(ValueError,"immutable"):
            self.execute()

    def test_original_selected_photo_ID_uniqueness_mandatory(self):
        df=pd.read_csv(self.sample,dtype=str)
        df.loc[1,"photo_id"]="1"
        df.to_csv(self.sample,index=False)
        with self.assertRaisesRegex(ValueError,"duplicated or empty"):
            self.execute()

    def test_observation_format_categories_do_not_imply_independence(self):
        self.assertEqual(diag.field_type("none"),"missing")
        self.assertEqual(diag.field_type("12345"),"numeric")
        self.assertEqual(diag.field_type(self.uuid[1]),"uuid")
        self.assertEqual(diag.field_type("abcdef"),"other")


if __name__=="__main__":
    unittest.main()
