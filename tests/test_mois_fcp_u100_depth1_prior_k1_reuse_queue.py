"""No-new-image test: FCP U100 depth1 source vs old pixels and active K1.

A source photo belongs to exactly one of five dispositions: reuse old actual
ROI, await K1 in-flight photo, existing independent technical photo, unknown
unsupported loader format, or NEW original photo for later frozen ROI. All
original selection probabilities remain 1/N_taxon/M_original_photos.
"""
from __future__ import annotations
import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from scripts.acquisition import mois_fcp_u100_depth1_prior_k1_reuse_queue as bridge
from scripts.acquisition.mois_fcp_u100_progressive_depth1_photo_source import ARM as FCP_ARM
from scripts.acquisition.mois_completed_pixel_checkpoint import sha


class FCPFull93mFirstPhotoReuseTests(unittest.TestCase):
    def setUp(self):
        self.tmp=TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.prefix=self.root/"prefix.csv"
        self.receipt=self.root/"prefix.json"
        self.k1=self.root/"k1.csv"
        self.k1_receipt=self.root/"k1.json"
        self.technical=self.root/"technical.csv"
        self.out=self.root/"source"
        self.prefix_rows=[]
        specs=[
            ("1","obs-a","100","jpg"),
            ("2","obs-b","101","jpg"),
            ("3","obs-c","102","jpg"),
            ("4","obs-d","103","gif"),
            ("5","obs-e","104","jpg"),
            ("6","obs-b","105","jpg"), # different photo of K1-selected observation
        ]
        for pid,obs,taxon,ext in specs:
            p=1/100/2
            row={
                "photo_id":pid,
                "observation_id":obs,
                "inat_taxon_id":taxon,"species":"Flower synthetic",
                "observer_id":"observer","observed_on":"2025-06-10",
                "recorded_month":"06","latitude_band":"north_extratropical",
                "latitude":"35.5","longitude":"140",
                "positional_accuracy":"3","photo_license":"cc-by",
                "original_image_extension":ext,
                "photos_in_original_observation":"2",
                "full_93m_distinct_original_observations_of_taxon":"100",
                "U100_source_original_observation_rank":"1",
                "depth1_observation_inclusion_probability":".01",
                "depth1_one_photo_inclusion_probability":".5",
                "depth1_joint_original_photo_inclusion_probability":format(p,".17g"),
                "parent_U100_joint_original_photo_inclusion_probability":format(100*p,".17g"),
                "parent_original_photo_U100_source_arm":"mois_full93m_fcp_exact100_original_observation_per_taxon_v1",
                "progressive_source_arm":FCP_ARM,
                "photo_source_disposition":"supported" if ext=="jpg" else "unsupported_original_extension",
                "photo_url_large":f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.jpg" if ext=="jpg" else "",
                "independent_accepted_species_and_native_stage_verified":"false",
            }
            self.prefix_rows.append(row)
        with self.prefix.open("w",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=list(self.prefix_rows[0]),lineterminator="\n")
            writer.writeheader();writer.writerows(self.prefix_rows)
        self.receipt.write_text(json.dumps({
            "schema":"mois_full93m_fcp_U100_prefix1_complete_20242_nominal_taxa_v1",
            "first_original_photo_candidates":6,
            "full_93m_original_photo_id_opportunities":93_209_351,
            "selected_original_photo_candidates_SHA256":sha(self.prefix),
            "original_photo_inclusion_probability":"1/N_taxon/M_original_photos_per_observation, conditional on 100-obs source eligibility",
        }))
        self.k1.write_text("photo_id,observation_id,photo_source_disposition,joint_photo_selection_probability\n2,obs-b,supported,0.00001\n7,obs-z,supported,0.00001\n")
        self.k1_receipt.write_text(json.dumps({
            "schema":"mois_full_source_observation_sampling_cohort_v1",
            "sampling_schema":"mois_93m_tier06_taxon_month_lat_stratified_onephoto_v1",
            "sampled_distinct_original_observations":2,
            "source_sample_decided_without_photo_pixels":True,
        }))
        self.technical.write_text("photo_id,observation_id,acquired,status\n3,obs-c,true,technical\n")
        self.priors=[("main",self.root/"m.parquet",self.root/"m.json"),
                     ("prior_group00",self.root/"x.parquet",self.root/"x.json"),
                     ("prior_R1",self.root/"y.parquet",self.root/"y.json"),
                     ("prior_R2",self.root/"z.parquet",self.root/"z.json")]
        self.old={
            "1":("main","obs-a","white","classified_four_state_morph","a"*64)
        }

    def bridge(self):
        with patch.object(bridge,"prior_original_photo_index",return_value=(
              self.old,{"obs-a":("main","1")},[],{},[],[])):
            return bridge.make_work(
                self.prefix,self.receipt,self.priors,
                self.k1,self.k1_receipt,self.technical,self.out,
                expected_prefix=6,expected_prefix_receipt_sha=sha(self.receipt),
                expected_prefix_selected_sha=sha(self.prefix),
                expected_k1=2,expected_k1_selected_sha=sha(self.k1),
                expected_k1_manifest_sha=sha(self.k1_receipt),
                expected_old={},expected_technical_sha=sha(self.technical),
                max_new=6)

    def test_exact_source_dispositions_and_no_doubled_old_or_active_photo(self):
        d=self.bridge()
        self.assertEqual(d["selected_original_photo_observations"],6)
        self.assertEqual(d["previously_measured_original_photo_ids_reused"],1)
        self.assertEqual(d["inflight_K1_photo_IDs_reserved_without_redownload"],1)
        self.assertEqual(d["technical_pilot_photo_IDs_not_refetched"],1)
        self.assertEqual(d["unsupported_original_photo_format_selected_as_unknown"],1)
        self.assertEqual(d["new_photo_attempts_planned"],2)
        self.assertEqual(d["work_units"],1)
        self.assertEqual(d["other_photo_same_original_observation_kept_for_design"][
            "K1_observation_other_photo_preserved_as_separate_random_photo"],1)
        self.assertFalse(d["colour_measured"])
        self.assertEqual(d["original_images_downloaded_by_source_bridge"],0)
        with (self.out/"all_FCP_U100_prefix1_photo_source_dispositions.csv").open(newline="") as f:
            rows={r["photo_id"]:r for r in csv.DictReader(f)}
        self.assertEqual(rows["1"]["reuse_or_new_photo_outcome"],
                         "prior_original_photo_ROI_reuse")
        self.assertEqual(rows["1"]["prior_original_ROI_morph"],"white")
        self.assertEqual(rows["2"]["reuse_or_new_photo_outcome"],
                         "reserved_by_inflight_K1_original_photo")
        self.assertEqual(rows["3"]["reuse_or_new_photo_outcome"],
                         "prior_independent_technical_photo_holdout")
        self.assertEqual(rows["4"]["reuse_or_new_photo_outcome"],
                         "unsupported_original_photo_format_technical_unknown")
        self.assertEqual(rows["6"]["reuse_or_new_photo_outcome"],
                         "new_original_photo_ROI_queued")
        with (self.out/"shards/shard-000.csv").open(newline="") as f:
            work=list(csv.DictReader(f))
        self.assertEqual({r["photo_id"] for r in work},{"5","6"})
        self.assertTrue(all(float(r["selection_probability"])==.005 for r in work))
        self.assertTrue(all(r["source_arm"]==FCP_ARM for r in work))
        self.assertTrue(d["photo_already_reserved_by_K1_never_fetched_twice"])

    def test_source_photo_or_group_receipt_tamper_refused(self):
        self.prefix.write_text(self.prefix.read_text()+"\n")
        with self.assertRaisesRegex(ValueError,"SHA differs"):
            self.bridge()
        self.assertFalse(self.out.exists())

    def test_same_k1_photo_id_with_wrong_observation_UUID_fails(self):
        self.k1.write_text("photo_id,observation_id,photo_source_disposition,joint_photo_selection_probability\n2,obs-wrong,supported,0.00001\n7,obs-z,supported,0.00001\n")
        with self.assertRaisesRegex(ValueError,"different signed observation"):
            self.bridge()
        self.assertFalse(self.out.exists())

    def test_same_observation_other_original_photo_is_kept_not_swapped(self):
        d=self.bridge()
        self.assertEqual(d["new_photo_attempts_planned"],2)
        self.assertNotEqual(self.prefix_rows[1]["photo_id"],self.prefix_rows[5]["photo_id"])
        self.assertEqual(self.prefix_rows[1]["observation_id"],self.prefix_rows[5]["observation_id"])

    def test_past_source_output_is_immutable(self):
        self.bridge()
        with self.assertRaisesRegex(ValueError,"Never overwrite"):
            self.bridge()

    def test_photo_sampling_probabilities_do_not_silently_change(self):
        self.prefix_rows[-1]["depth1_joint_original_photo_inclusion_probability"]="0.8"
        with self.prefix.open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=list(self.prefix_rows[0]),lineterminator="\n")
            w.writeheader();w.writerows(self.prefix_rows)
        source=json.loads(self.receipt.read_text())
        source["selected_original_photo_candidates_SHA256"]=sha(self.prefix)
        self.receipt.write_text(json.dumps(source))
        with self.assertRaisesRegex(ValueError,"true FCP prefix|inclusion probability"):
            self.bridge()
        self.assertFalse(self.out.exists())


if __name__=="__main__":
    unittest.main()
