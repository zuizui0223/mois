"""Data-only preflight of actual FCP 3141 photo pipeline before K1 completion."""
from __future__ import annotations

import csv
import gzip
import json
import math
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from scripts.acquisition import mois_fcp_u100_depth1_preflight as m
from scripts.acquisition.mois_completed_pixel_checkpoint import sha


class OriginalFCPPrecolourPreflightTests(unittest.TestCase):
    def setUp(self):
        self.tmp=TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.source=self.root/"source"
        (self.source/"shards").mkdir(parents=True)
        self.prefix=self.root/"prefix.csv"
        self.disp=self.root/"dispositions.csv"
        self.cert=self.root/"k1_cert.json"
        self.k1=self.root/"k1_roster.csv.gz"
        self.output=self.root/"out"
        self.elig=5
        self.reuse=2
        self.reserved=1
        self.new=2
        self.k1_photo_ids=["102","900","901"]
        self.rows=[]
        for index in range(self.elig):
            pid=str(100+index)
            n=100+index
            p=1/n/2
            self.rows.append(dict(
                photo_id=pid,observation_id=f"OBS-{pid}",inat_taxon_id=str(1+index),
                photo_url_large=f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.jpg",
                photo_license="cc-by",original_image_extension="jpg",
                photos_in_original_observation="2",
                full_93m_distinct_original_observations_of_taxon=str(n),
                depth1_joint_original_photo_inclusion_probability=format(p,".17g"),
                parent_U100_joint_original_photo_inclusion_probability=format(100*p,".17g"),
                progressive_source_arm=m.PREFIX_ARM
            ))
        self.prepare()

    @staticmethod
    def write_csv(path,fields,items):
        with path.open("w",newline="",encoding="utf-8") as stream:
            w=csv.DictWriter(stream,fieldnames=fields,lineterminator="\n")
            w.writeheader()
            w.writerows(items)

    def prepare(self):
        self.write_csv(self.prefix,list(self.rows[0]),self.rows)
        assignments=["prior_original_photo_ROI_reuse"]*2+[
            "reserved_by_inflight_K1_original_photo",
            "new_original_photo_ROI_queued","new_original_photo_ROI_queued"]
        self.dispositions=[]
        for index,row in enumerate(self.rows):
            action=assignments[index]
            self.dispositions.append(dict(
                photo_id=row["photo_id"],observation_id=row["observation_id"],
                inat_taxon_id=row["inat_taxon_id"],
                depth1_joint_original_photo_inclusion_probability=
                    row["depth1_joint_original_photo_inclusion_probability"],
                source_fcp_u100_parent_joint_photo_probability=
                    row["parent_U100_joint_original_photo_inclusion_probability"],
                source_arm=m.PREFIX_ARM,reuse_or_new_photo_outcome=action,
                prior_original_ROI_morph="mixed_uncertain" if index<2 else "",
                prior_photo_arm="old" if index<2 else "",
            ))
        self.write_csv(self.disp,list(self.dispositions[0]),self.dispositions)
        k1rows=[dict(photo_id=pid,observation_id=f"OBS-{pid}",
                     inat_taxon_id=str(3 if pid=="102" else 100),
                     selection_probability="0.1",
                     original_worker_group="0") for pid in self.k1_photo_ids]
        with gzip.open(self.k1,"wt",newline="",encoding="utf-8") as stream:
            w=csv.DictWriter(stream,fieldnames=list(k1rows[0]),lineterminator="\n")
            w.writeheader()
            w.writerows(k1rows)
        self.cert.write_text(json.dumps({
            "schema":m.K1_CERT_SCHEMA,
            "K1_all_new_precolour_source_original_photo_IDs":len(self.k1_photo_ids),
            "K1_exact_historical_photo_ID_collisions":0,
            "K1_new_source_photo_roster_sha256":sha(self.k1),
            "original_photo_joint_selection_probability_retained":True,
            "actual_photo_bytes_downloaded_or_model_inference_in_this_certificate":0,
        }),encoding="utf-8")
        expected_new=self.rows[3:]
        self.shard=self.source/"shards/shard-000.csv"
        self.workers=[]
        for row in expected_new:
            self.workers.append(dict(
                photo_id=row["photo_id"],observation_id=row["observation_id"],
                inat_taxon_id=row["inat_taxon_id"],
                source_arm=m.PREFIX_ARM,shard="0",
                selection_probability=row["depth1_joint_original_photo_inclusion_probability"],
                photo_url_large=row["photo_url_large"],
                measurement_id=f"MEAS-{row['photo_id']}",
            ))
        self.write_csv(self.shard,list(self.workers[0]),self.workers)
        self.manifest=self.source/"manifest.json"
        self.manifest.write_text(json.dumps({
            "schema":"mois_full_source_observation_sampling_cohort_v1",
            "full_frame_photo_ids":93209351,
            "source_arm":m.PREFIX_ARM,
            "selected_FCP_U100_eligible_source_taxa":self.elig,
            "selected_original_photo_observations":self.elig,
            "previously_measured_original_photo_ids_reused":self.reuse,
            "inflight_K1_photo_IDs_reserved_without_redownload":self.reserved,
            "new_photo_attempts_planned":self.new,
            "work_units":1,"worker_ids":[0],"group_ids":[0],
            "units_per_group":4,"original_model_batch_size":256,
            "technical_pilot_photo_IDs_not_refetched":0,
            "unsupported_original_photo_format_selected_as_unknown":0,
            "source_photo_inclusion_probability_is_per_original_photo":True,
            "photo_already_reserved_by_K1_never_fetched_twice":True,
            "FCP_top100_photo_prefix_inclusion_probability_corrected_to_1_over_N_M":True,
            "historic_pixels_reused_at_exact_original_photo_ID_only":True,
            "native_flowering_botanical_and_syndrome_status_verified":False,
            "colour_measured":False,
            "original_images_downloaded_by_source_bridge":0,
            "no_original_image_ROI_inference_by_source_bridge":True,
            "new_FCP_photo_ROI_jobs_not_launched_until_separate_resource_and_duplicate_gate":True,
            "frozen_all_selected_FCP_photo_dispositions_SHA256":sha(self.disp),
            "source_full93m_prefix1_selected_csv_SHA256":sha(self.prefix),
            "photo_shard_sha256":{"0":sha(self.shard)},
            "photo_ID_technical_unknown_and_reuse_dispositions":{
                "prior_original_photo_ROI_reuse":self.reuse,
                "reserved_by_inflight_K1_original_photo":self.reserved,
                "prior_independent_technical_photo_holdout":0,
                "unsupported_original_photo_format_technical_unknown":0,
                "new_original_photo_ROI_queued":self.new,
            }
        }),encoding="utf-8")

    def run_preflight(self):
        return m.preflight(
            self.source,self.prefix,self.disp,self.cert,self.k1,self.output,
            source_sha=sha(self.manifest),prefix_sha=sha(self.prefix),
            dispositions_sha=sha(self.disp),
            k1_certificate_sha=sha(self.cert),k1_roster_sha=sha(self.k1),
            eligible=self.elig,reuse=self.reuse,reserved=self.reserved,
            new=self.new,expected_work_units=1,expected_groups=1,
            expected_K1_source=len(self.k1_photo_ids))

    def test_eligible_taxa_and_exact_first_photo_ID_decisions_preflight(self):
        d=self.run_preflight()
        self.assertEqual(d["eligible_nominal_FCP_U100_source_taxa"],5)
        self.assertEqual(d["never_before_measured_FCP_U100_original_photo_ROI_target_IDs"],2)
        self.assertEqual(d["already_running_K1_new_photo_IDs_to_reuse_NOT_RECOMPUTED"],1)
        self.assertEqual(d["already_measured_old_photo_IDs_to_reuse_NOT_RECOMPUTED"],2)
        self.assertEqual(d["new_original_image_bytes_downloaded_here"],0)
        self.assertEqual(d["new_original_photo_colour_results_generated_here"],0)
        self.assertFalse((self.root/"source"/"original_photo_ROIs.parquet").exists())
        self.assertTrue((self.output/"fcp_u100_3141_photo_only_source_preflight_receipt.json").exists())

    def test_reserved_K1_photo_must_be_same_original_observation_and_taxon(self):
        rows=list(self.dispositions)
        rows[2]["observation_id"]="WRONG_UUID"
        self.write_csv(self.disp,list(rows[0]),rows)
        self.manifest.write_text(
            self.manifest.read_text().replace(
                self.dispositions[0]["photo_id"],self.dispositions[0]["photo_id"]))
        with self.assertRaisesRegex(ValueError,"photo/observation/taxon"):
            self.run_preflight()

    def test_new_photo_cannot_be_K1_already_reserved_and_redownloaded(self):
        rows=list(self.dispositions)
        rows[2]["reuse_or_new_photo_outcome"]="new_original_photo_ROI_queued"
        self.write_csv(self.disp,list(rows[0]),rows)
        with self.assertRaisesRegex(ValueError,"already active K1|photo"):
            self.run_preflight()

    def test_a_shard_byte_flip_cannot_be_published_as_new_source(self):
        self.shard.write_text(self.shard.read_text()+"\nforged")
        with self.assertRaisesRegex(ValueError,"SHA changed"):
            self.run_preflight()

    def test_wrong_1_over_N_over_M_probability_is_rejected(self):
        rows=list(self.rows)
        rows[0]=dict(rows[0])
        rows[0]["depth1_joint_original_photo_inclusion_probability"]="0.5"
        self.write_csv(self.prefix,list(rows[0]),rows)
        with self.assertRaisesRegex(ValueError,"1/N/M changed"):
            self.run_preflight()

    def test_immutable_original_target_preflight_output(self):
        self.run_preflight()
        with self.assertRaisesRegex(ValueError,"Refuse to overwrite"):
            self.run_preflight()


if __name__=="__main__":
    unittest.main()
