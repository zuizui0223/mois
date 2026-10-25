"""Source-only global U100 original observation photo ranks 1..5 SHA union tests."""
from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from scripts.acquisition import mois_fcp_u100_depth5_global_source as g
from scripts.acquisition.mois_fcp_u100_depth5_source import source_row
from scripts.acquisition.mois_onepass_full93m_target_tiers import sha256


class TestFullGlobalNestedU100Source(unittest.TestCase):
    def setUp(self):
        self.tmp=TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.first=self.root/"first.csv"
        self.folders=[self.root/"group0",self.root/"group1"]
        self.output=self.root/"global"
        self.source={}
        for group,taxon in ((0,"128"),(1,"136")):
            folder=self.folders[group];folder.mkdir()
            source=[]
            for rank in range(1,6):
                pid=str(12000+group*100+rank)
                uuid=f"uuid-{group}-{rank}"
                parent=dict(
                    photo_id=pid,observation_id=uuid,inat_taxon_id=taxon,
                    species="synthetic taxon",observer_id="person",
                    observed_on="2024-08-07",recorded_month="8",
                    latitude_band="north",latitude="35",longitude="139",
                    positional_accuracy="50",photo_license="cc-by",
                    original_image_extension="jpg",
                    photos_in_original_observation="1",
                    full_93m_distinct_original_observations_of_taxon="100",
                    conditional_species_observation_probability="1",
                    conditional_one_photo_probability="1",
                    conditional_full_93m_photo_selection_probability="1",
                    photo_source_disposition="supported",
                    photo_url_large=f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.jpg",
                    source_arm=g.__dict__.get("PARENT_ARM","mois_full93m_fcp_exact100_original_observation_per_taxon_v1"),
                )
                source.append(source_row(parent,rank,str(12000+group*100+1)))
            self.source[group]=source
        self.write_all()

    def make_csv(self,path,header,rows):
        with path.open("w",encoding="utf-8",newline="") as f:
            w=csv.DictWriter(f,fieldnames=header,lineterminator="\n")
            w.writeheader();w.writerows(rows)

    def write_all(self):
        self.make_csv(self.first,
            ["photo_id","observation_id","inat_taxon_id",
             "depth1_joint_original_photo_inclusion_probability"],
            [dict(photo_id=self.source[g][0]["photo_id"],
                  observation_id=self.source[g][0]["observation_id"],
                  inat_taxon_id=self.source[g][0]["inat_taxon_id"],
                  depth1_joint_original_photo_inclusion_probability=self.source[g][0][
                      "depth1_joint_original_photo_inclusion_probability"])
             for g in (0,1)])
        for group in (0,1):
            folder=self.folders[group]
            allfile=folder/g.PREFIX_FILENAME
            delta=folder/g.DELTA_FILENAME
            self.make_csv(allfile,g.FIELDS,self.source[group])
            self.make_csv(delta,g.FIELDS,self.source[group][1:])
            receipt=dict(
                schema=g.GROUP_SCHEMA,
                original_taxon_modulo_source_group=group,
                global_U100_complete_42111_source_receipt_SHA256=g.GLOBAL_RECEIPT_SHA,
                nested_original_photo_depth5_source_arm=g.ARM,
                global_source_original_photo_IDs_metadata_only=g.EXPECTED_PHOTOS,
                global_nominal_source_taxa=g.EXPECTED_TAXA,
                global_U100_source_eligible_taxa=2,
                source_group_original_photo_IDs_distinct=True,
                source_group_original_observations_distinct=True,
                source_group_candidate_photo_rank_unchanged=True,
                original_photo_IDs_already_measured_by_old_FCP_MOIS_and_K1_excluded=False,
                new_image_download_or_frozen_ROI_inference_this_stage=0,
                first5_csv_SHA256=sha256(allfile),
                ranks2to5_csv_SHA256=sha256(delta),
                first5_observation_independent_original_photo_candidates=5,
                source_group_eligible_nominal_taxa=1,
                new_candidate_ranks2to5_NOT_yet_measured=4,
            )
            (folder/"nested5_source_receipt.json").write_text(json.dumps(receipt))

    def assemble(self):
        return g.global_source(self.folders,self.first,self.output,
            expected_groups=2,expected_eligible=2,
            expected_prefix5=10,expected_incremental=8,
            expected_first1_sha=sha256(self.first))

    def test_all_groups_nested_photo_id_source_conserved_and_gz_sha(self):
        receipt=self.assemble()
        self.assertEqual(receipt["global_first5_nested_source_photo_IDs"],10)
        self.assertEqual(receipt["global_incremental_rank2to5_candidate_photo_IDs_NOT_NEW_INFERENCES"],8)
        self.assertEqual(receipt["global_first1_original_photo_IDs_unchanged"],2)
        self.assertTrue(receipt["ranks2to5_genuine_never_measured_photo_IDs_not_yet_deduplicated"])
        self.assertEqual(receipt["new_image_downloads_and_original_ROI_inferences_performed"],0)
        for filename,count in (
            ("global_FCP_U100_nested_first5_original_photo_candidates.csv.gz",10),
            ("global_FCP_U100_incremental_rank2to5_candidates_NEED_DEDUP.csv.gz",8)):
            file=self.output/filename
            self.assertEqual(sha256(file),receipt["csv_outputs_SHA256"][filename])
            with gzip.open(file,"rt",encoding="utf-8",newline="") as stream:
                rows=list(csv.DictReader(stream))
            self.assertEqual(len(rows),count)
            self.assertEqual(len({x["photo_id"] for x in rows}),count)

    def test_group_original_byte_tamper_fails_before_publication(self):
        file=self.folders[0]/g.DELTA_FILENAME
        file.write_bytes(file.read_bytes()+b"mutated")
        with self.assertRaisesRegex(ValueError,"CSV file bytes"):
            self.assemble()
        self.assertFalse(self.output.exists())

    def test_repeated_crossgroup_photo_ID_rejected_even_if_resigned(self):
        self.source[1][2]["photo_id"]=self.source[0][2]["photo_id"]
        self.write_all()
        with self.assertRaisesRegex(ValueError,"duplicated"):
            self.assemble()
        self.assertFalse(self.output.exists())

    def test_rank1_must_reproduce_original_staged_prefix(self):
        first=csv.DictReader(self.first.open(encoding="utf-8",newline=""))
        row=list(first)
        row[0]["photo_id"]="999999"
        self.make_csv(self.first,list(row[0]),row)
        with self.assertRaisesRegex(ValueError,"not the immutable original"):
            self.assemble()

    def test_immutable_previous_worldwide_stage5_cannot_be_rewritten(self):
        self.assemble()
        with self.assertRaisesRegex(ValueError,"immutable"):
            self.assemble()


if __name__=="__main__":
    unittest.main()
