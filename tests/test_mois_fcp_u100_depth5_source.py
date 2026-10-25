"""Synthetic SOURCE-ONLY nested FCP rank1-to5/2-to5 exact photo ID sampling."""
from __future__ import annotations

import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from scripts.acquisition import mois_fcp_u100_depth5_source as m
from scripts.acquisition.mois_fcp_u100_progressive_depth1_photo_source import SOURCE_FIELDS
from scripts.acquisition.mois_fcp_u100_full93m_source_onepass_group import original_observation_priority
from scripts.acquisition.mois_onepass_full93m_target_tiers import sha256


class TestNestedSourceOriginalFCPDepth5(unittest.TestCase):
    def setUp(self):
        self.temp=TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.csv=self.root/"parent-top100.csv"
        self.group=self.root/"group.json"
        self.global_receipt=self.root/"global.json"
        self.output=self.root/"freeze"
        # Taxon 128 belongs to taxon modulo group 0; SOURCE rank is stable,
        # outcome-blind, and determined by the original 93m observation seed.
        ranked=sorted((f"orig-uuid-{i}" for i in range(100)),
            key=lambda x:(original_observation_priority(x),x))
        self.source_rows=[]
        for i,obs in enumerate(ranked):
            pid=str(100000+i)
            n=120
            photo_count=2
            self.source_rows.append({
                "photo_id":pid,"observation_id":obs,
                "inat_taxon_id":"128","species":"Campanula syntheticum",
                "observer_id":f"observer-{i}",
                "observed_on":"2024-07-16",
                "recorded_month":"7","latitude_band":"north",
                "latitude":"35.0","longitude":"135.0",
                "positional_accuracy":"100",
                "photo_license":"cc-by",
                "original_image_extension":"jpg",
                "photos_in_original_observation":str(photo_count),
                "full_93m_distinct_original_observations_of_taxon":str(n),
                "conditional_species_observation_probability":format(100/n,".17g"),
                "conditional_one_photo_probability":format(1/photo_count,".17g"),
                "conditional_full_93m_photo_selection_probability":
                    format(100/n/photo_count,".17g"),
                "photo_source_disposition":"supported",
                "photo_url_large":
                    f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.jpg",
                "source_arm":m.PARENT_ARM,
            })
        self.freeze_inputs()

    def freeze_inputs(self):
        with self.csv.open("w",encoding="utf-8",newline="") as file:
            writer=csv.DictWriter(file,fieldnames=SOURCE_FIELDS,lineterminator="\n")
            writer.writeheader()
            writer.writerows(self.source_rows)
        d={
            "schema":m.PARENT_SCHEMA,
            "original_taxon_modulo_source_group":0,
            "full_original_photo_source_IDs":100,
            "full_original_named_plant_taxa":1,
            "FCP_original_photo_observation_target_per_taxon":100,
            "no_original_image_downloads":True,
            "no_original_ROI_or_hue_inferences":True,
            "all_old_FCP_MOIS_original_photo_IDs_must_be_reused_before_any_new_ROI":True,
            "source_group_original_photo_ID_SHA":{
                "original_full_93m_source_registry_SHA256":"a"*64},
            "U100_complete_original_photo_candidate_SHA256":sha256(self.csv),
            "FCP_U100_source_eligible_nominal_taxa":1,
            "FCP_U100_source_original_photo_opportunities":100,
        }
        self.group.write_text(json.dumps(d),encoding="utf-8")
        g={
            "schema":m.GLOBAL_SCHEMA,
            "source_all_original_photo_metadata_IDs":100,
            "all_source_nominal_taxon_IDs":1,
            "source_FCP_style_exact100_photo_observation_taxa":1,
            "source_FCP_style_exact100_observation_photo_targets":100,
            "photo_ROI_pixels_opened_or_reinferred":0,
            "full93m_original_source_registry_sha256":"a"*64,
            "group_original_photo_SHA_proof":[
                {"group":j,"eligible_taxa":int(j==0),
                 "group_receipt_sha256":sha256(self.group) if j==0 else "f"*64}
                for j in range(16)],
        }
        self.global_receipt.write_text(json.dumps(g),encoding="utf-8")

    def freeze(self):
        return m.freeze_group(
            self.csv,self.group,self.global_receipt,self.output,0,
            expected_global_sha=sha256(self.global_receipt),
            expected_global_taxa=1,expected_all_taxa=1,expected_photos=100,
        )

    @staticmethod
    def read_csv(path):
        with Path(path).open(newline="",encoding="utf-8") as f:
            return list(csv.DictReader(f))

    def test_nested_rank5_and_first1_exact_original_observation_photo_source(self):
        d=self.freeze()
        self.assertEqual(d["source_group_eligible_nominal_taxa"],1)
        self.assertEqual(d["full_100_observation_photo_rank_candidates_audited"],100)
        self.assertEqual(d["first5_observation_independent_original_photo_candidates"],5)
        self.assertEqual(d["new_candidate_ranks2to5_NOT_yet_measured"],4)
        a=self.read_csv(self.output/"original_FCP_U100_prefix5_rank1to5_candidates.csv")
        b=self.read_csv(self.output/
            "original_FCP_U100_rank2to5_INCREMENTAL_CANDIDATES_not_new_ROIs.csv")
        self.assertEqual([x["photo_id"] for x in a],
                         [x["photo_id"] for x in self.source_rows[:5]])
        self.assertEqual([x["photo_id"] for x in b],
                         [x["photo_id"] for x in self.source_rows[1:5]])
        self.assertEqual([int(x["U100_source_original_observation_rank"]) for x in a],
                         [1,2,3,4,5])
        self.assertTrue(all(x["source_rank1_original_photo_id"]==a[0]["photo_id"]
                            for x in a))
        for i,x in enumerate(a):
            self.assertAlmostEqual(float(x["cumulative_prefix5_joint_original_photo_probability"]),
                                   5/120/2)
            self.assertAlmostEqual(float(x["rank2to5_incremental_joint_original_photo_probability"]),
                                   4/120/2)
            self.assertAlmostEqual(float(x["depth1_joint_original_photo_inclusion_probability"]),
                                   1/120/2)
            self.assertAlmostEqual(float(x["parent_U100_joint_original_photo_inclusion_probability"]),
                                   100/120/2)
            self.assertEqual(x["incremental_since_depth1"],str(i>0).lower())
            self.assertEqual(x["independent_accepted_species_and_native_stage_verified"],"false")
        self.assertFalse(d["original_photo_IDs_already_measured_by_old_FCP_MOIS_and_K1_excluded"])
        self.assertEqual(d["new_image_download_or_frozen_ROI_inference_this_stage"],0)

    def test_source_exact_byte_sha_mutation_fails_before_output(self):
        self.csv.write_bytes(self.csv.read_bytes()+b"\nforged")
        with self.assertRaisesRegex(ValueError,"ranked 100-photo taxon source or SHA"):
            self.freeze()
        self.assertFalse(self.output.exists())

    def test_rank_changed_but_resigned_sha_still_rejected(self):
        self.source_rows[0],self.source_rows[1]=self.source_rows[1],self.source_rows[0]
        self.freeze_inputs()
        with self.assertRaisesRegex(ValueError,"observation priority changed"):
            self.freeze()
        self.assertFalse(self.output.exists())

    def test_probability_changed_but_resigned_sha_still_rejected(self):
        self.source_rows[1]["conditional_species_observation_probability"]="0.1"
        self.freeze_inputs()
        with self.assertRaisesRegex(ValueError,"inclusion probability"):
            self.freeze()
        self.assertFalse(self.output.exists())

    def test_repeated_photo_id_even_after_resign_rejected(self):
        self.source_rows[4]["photo_id"]=self.source_rows[3]["photo_id"]
        self.freeze_inputs()
        with self.assertRaisesRegex(ValueError,"photo ID duplicated"):
            self.freeze()
        self.assertFalse(self.output.exists())

    def test_cannot_overwrite_real_stage5_original_source(self):
        self.freeze()
        with self.assertRaisesRegex(ValueError,"Never mutate"):
            self.freeze()


if __name__=="__main__":
    unittest.main()
