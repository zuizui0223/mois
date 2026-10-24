"""Independent synthetic audit for exact FCP U100 depth-1 photo sample.

Tests no network/ROI. Source is 2 eligible nominal taxa each with exactly 100
different observation IDs, in a 16-group globally source-signed 93m analogue.
"""
from __future__ import annotations
import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from scripts.acquisition import mois_fcp_u100_progressive_depth1_photo_source as p
from scripts.acquisition.mois_fcp_u100_full93m_source_onepass_group import (
    SOURCE_ARM as PARENT_ARM, original_observation_priority,
)


class FCPProgressiveDepthOneTests(unittest.TestCase):
    def setUp(self):
        self.temp=TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.infile=self.root/"FCP_U100_original.csv"
        self.group_receipt=self.root/"group.json"
        self.global_receipt=self.root/"global.json"
        self.output=self.root/"stage1"
        self.frozen_photo_rows=[]
        self.first={}
        seen=set()
        for taxon,n,m in (("128",150,3),("256",250,2)):
            group=[]
            for j in range(100):
                uuid=f"observation-{taxon}-{j:04d}"
                photo_id=str(10_000_000+int(taxon)*1000+j)
                extension="gif" if taxon=="256" and j==0 else "jpg"
                group.append((original_observation_priority(uuid),uuid,photo_id,extension))
            group.sort()
            for index,(_,uuid,photo_id,ext) in enumerate(group):
                rank="supported" if ext=="jpg" else "unsupported_original_extension"
                raw={
                    "photo_id":photo_id,"observation_id":uuid,
                    "inat_taxon_id":taxon,"species":f"Plantus example {taxon}",
                    "observer_id":"tester","observed_on":"2025-07-02",
                    "recorded_month":"07","latitude_band":"north_extratropical",
                    "latitude":"35","longitude":"140",
                    "positional_accuracy":"5","photo_license":"cc-by",
                    "original_image_extension":ext,
                    "photos_in_original_observation":str(m),
                    "full_93m_distinct_original_observations_of_taxon":str(n),
                    "conditional_species_observation_probability":format(100/n,".17g"),
                    "conditional_one_photo_probability":format(1/m,".17g"),
                    "conditional_full_93m_photo_selection_probability":format(100/n/m,".17g"),
                    "photo_source_disposition":rank,
                    "photo_url_large":(
                        f"https://inaturalist-open-data.s3.amazonaws.com/photos/{photo_id}/large.jpg"
                        if ext=="jpg" else ""),
                    "source_arm":PARENT_ARM,
                }
                self.frozen_photo_rows.append(raw)
                self.assertNotIn(photo_id,seen)
                seen.add(photo_id)
                if index==0:self.first[taxon]=raw
        with self.infile.open("w",newline="",encoding="utf-8") as f:
            w=csv.DictWriter(f,fieldnames=p.SOURCE_FIELDS,lineterminator="\n")
            w.writeheader()
            w.writerows(self.frozen_photo_rows)
        self.parent={
            "schema":p.PARENT_SCHEMA,
            "original_taxon_modulo_source_group":0,
            "full_original_photo_source_IDs":50000,
            "full_original_named_plant_taxa":3,
            "FCP_original_photo_observation_target_per_taxon":100,
            "no_original_image_downloads":True,
            "no_original_ROI_or_hue_inferences":True,
            "all_old_FCP_MOIS_original_photo_IDs_must_be_reused_before_any_new_ROI":True,
            "source_group_original_photo_ID_SHA":{
                "original_full_93m_source_registry_SHA256":"f"*64,
            },
            "U100_complete_original_photo_candidate_SHA256":p.sha256(self.infile),
            "FCP_U100_source_original_photo_opportunities":200,
            "FCP_U100_source_eligible_nominal_taxa":2,
            "source_nominal_taxa":3,
        }
        self.group_receipt.write_text(json.dumps(self.parent))
        self.source_global={
            "schema":p.GLOBAL_SCHEMA,
            "source_all_original_photo_metadata_IDs":50000,
            "all_source_nominal_taxon_IDs":3,
            "source_FCP_style_exact100_photo_observation_taxa":2,
            "source_FCP_style_exact100_observation_photo_targets":200,
            "photo_ROI_pixels_opened_or_reinferred":0,
            "full93m_original_source_registry_sha256":"f"*64,
            "group_original_photo_SHA_proof":[
                {"group":g,
                 "eligible_taxa":2 if g==0 else 0,
                 "group_receipt_sha256":p.sha256(self.group_receipt) if g==0 else "0"*64}
                for g in range(16)
            ],
        }
        self.global_receipt.write_text(json.dumps(self.source_global))

    def freeze(self):
        return p.freeze_group(self.infile,self.group_receipt,self.global_receipt,
                              self.output,0,expected_taxa=3,
                              expected_photos=50000,expected_global_eligible=2)

    def test_frozen_progressive_prefix_one_photo_per_eligible_taxon(self):
        out=self.freeze()
        self.assertEqual(out["depth1_original_photo_observations_selected"],2)
        self.assertEqual(out["U100_original_source_candidates_accounted"],200)
        self.assertEqual(out["depth1_supported_original_image_formats"]+
                         out["depth1_unsupported_original_image_formats_unknown"],2)
        self.assertEqual(out["source_taxa_below100_observations_retained_in_42111_census"],1)
        self.assertEqual(out["original_images_downloaded"],0)
        self.assertEqual(out["original_colour_inferences_executed"],0)
        with (self.output/"FCP_U100_prefix1_original_photo_per_taxon.csv").open(
              newline="",encoding="utf-8") as f:
            rows=list(csv.DictReader(f))
        self.assertEqual(len(rows),2)
        self.assertEqual({r["inat_taxon_id"] for r in rows},{"128","256"})
        for row in rows:
            parent=self.first[row["inat_taxon_id"]]
            self.assertEqual(row["photo_id"],parent["photo_id"])
            self.assertEqual(row["observation_id"],parent["observation_id"])
            n=int(parent["full_93m_distinct_original_observations_of_taxon"])
            m=int(parent["photos_in_original_observation"])
            self.assertAlmostEqual(float(row["depth1_observation_inclusion_probability"]),1/n)
            self.assertAlmostEqual(float(row["depth1_one_photo_inclusion_probability"]),1/m)
            self.assertAlmostEqual(float(row["depth1_joint_original_photo_inclusion_probability"]),1/n/m)
            self.assertAlmostEqual(float(row["parent_U100_joint_original_photo_inclusion_probability"]),100/n/m)
            self.assertEqual(row["independent_accepted_species_and_native_stage_verified"],"false")
            self.assertEqual(row["progressive_source_arm"],p.ARM)
        self.assertEqual(out["depth1_original_photo_candidate_SHA256"],
                         p.sha256(self.output/"FCP_U100_prefix1_original_photo_per_taxon.csv"))
        self.assertTrue(out["source_probability_true_for_eligible_taxa_only"])
        self.assertFalse(out["native_wild_or_individual_flowering_photo_verified"])

    def test_signed_original_sample_is_immutable(self):
        self.freeze()
        with self.assertRaisesRegex(ValueError,"immutable"):
            self.freeze()

    def test_original_group_photo_source_digest_corruption_fails_closed(self):
        self.infile.write_text(self.infile.read_text()+"\ntruncated")
        with self.assertRaisesRegex(ValueError,"pinned FCP U100"):
            self.freeze()
        self.assertFalse(self.output.exists())

    def test_original_U100_observation_order_tamper_not_hidden(self):
        records=list(self.frozen_photo_rows)
        records[0],records[1]=records[1],records[0]
        with self.infile.open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=p.SOURCE_FIELDS,lineterminator="\n")
            w.writeheader();w.writerows(records)
        self.parent["U100_complete_original_photo_candidate_SHA256"]=p.sha256(self.infile)
        self.group_receipt.write_text(json.dumps(self.parent))
        self.source_global["group_original_photo_SHA_proof"][0][
            "group_receipt_sha256"]=p.sha256(self.group_receipt)
        self.global_receipt.write_text(json.dumps(self.source_global))
        with self.assertRaisesRegex(ValueError,"priority order changed"):
            self.freeze()

    def test_photo_formats_unsupported_are_accounted_not_replaced(self):
        self.freeze()
        with (self.output/"FCP_U100_prefix1_original_photo_per_taxon.csv").open() as f:
            rows=list(csv.DictReader(f))
        unsupported=[r for r in rows if r["photo_source_disposition"]!="supported"]
        self.assertEqual(len(unsupported),self.first["256"]["original_image_extension"]=="gif")
        for row in unsupported:
            self.assertEqual(row["photo_url_large"],"")
            self.assertGreater(float(row["depth1_joint_original_photo_inclusion_probability"]),0)

    def test_same_taxon_extra_missing_or_reordered_original_source_photo_fails(self):
        rows=self.frozen_photo_rows[:-1]
        with self.infile.open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=p.SOURCE_FIELDS,lineterminator="\n")
            w.writeheader();w.writerows(rows)
        self.parent["U100_complete_original_photo_candidate_SHA256"]=p.sha256(self.infile)
        self.group_receipt.write_text(json.dumps(self.parent))
        self.source_global["group_original_photo_SHA_proof"][0][
            "group_receipt_sha256"]=p.sha256(self.group_receipt)
        self.global_receipt.write_text(json.dumps(self.source_global))
        with self.assertRaisesRegex(ValueError,"Final source taxon lost"):
            self.freeze()
        self.assertFalse(self.output.exists())


if __name__=="__main__":
    unittest.main()
