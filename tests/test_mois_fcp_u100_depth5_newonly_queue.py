"""Photo-ID and 1/N/M safe original ROI work queue after nested U100 rank2to5."""
from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from scripts.acquisition import mois_fcp_u100_depth5_newonly_queue as m
from scripts.acquisition.mois_fcp_u100_depth5_source import source_row
from scripts.acquisition.mois_completed_pixel_checkpoint import sha


class FCPNested5OriginalSourceDedupTests(unittest.TestCase):
    def setUp(self):
        self.temp=TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.prefix=self.root/"first5.csv.gz"
        self.delta=self.root/"delta.csv.gz"
        self.receipt=self.root/"receipt.json"
        self.output=self.root/"to-download"
        self.rows=[]
        for j,taxon in enumerate(("128","136")):
            first_id=str(90000+j*100+1)
            for rank in range(1,6):
                pid=str(90000+j*100+rank)
                obs=f"FCP-UUID-{j}-{rank}"
                supported=not (j==0 and rank==4)
                source=dict(
                    photo_id=pid,observation_id=obs,inat_taxon_id=taxon,
                    species="Test specimen",observer_id="person",
                    observed_on="2024-06-01",recorded_month="6",latitude_band="north",
                    latitude="32",longitude="133",positional_accuracy="100",
                    photo_license="cc-by",
                    original_image_extension="jpg" if supported else "gif",
                    photos_in_original_observation="1",
                    full_93m_distinct_original_observations_of_taxon="100",
                    conditional_species_observation_probability="1",
                    conditional_one_photo_probability="1",
                    conditional_full_93m_photo_selection_probability="1",
                    photo_source_disposition=(
                        "supported" if supported else "unsupported_original_extension"),
                    photo_url_large=(
                        f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.jpg"
                        if supported else ""),
                    source_arm="mois_full93m_fcp_exact100_original_observation_per_taxon_v1",
                )
                self.rows.append(source_row(source,rank,first_id))
        prior_indices=[0,5,1,7]
        self.known={}
        for i in prior_indices:
            row=self.rows[i]
            h="mixed_uncertain" if i==5 else "white"
            status="unclassifiable_original_photo" if i==5 else "classified_four_state_morph"
            self.known[row["photo_id"]]=(
                "signed_historical" if i in (0,1) else "K1_finished",
                row["observation_id"],h,status,"",
            )
        self.write()

    @staticmethod
    def gzip_csv(path,fields,items):
        with gzip.open(path,"wt",encoding="utf-8",newline="") as file:
            w=csv.DictWriter(file,fieldnames=fields,lineterminator="\n")
            w.writeheader();w.writerows(items)

    def write(self):
        self.gzip_csv(self.prefix,m.STAGE5_FIELDS,self.rows)
        self.gzip_csv(self.delta,m.STAGE5_FIELDS,[
            x for x in self.rows if int(x["U100_source_original_observation_rank"])>1])
        receipt={
            "schema":m.GLOBAL_STAGE5_SCHEMA,
            "global_first5_nested_source_photo_IDs":10,
            "global_incremental_rank2to5_candidate_photo_IDs_NOT_NEW_INFERENCES":8,
            "eligible_taxa_with_at_least100_original_observations":2,
            "ranks2to5_genuine_never_measured_photo_IDs_not_yet_deduplicated":True,
            "new_image_downloads_and_original_ROI_inferences_performed":0,
            "csv_outputs_SHA256":{
                "global_FCP_U100_nested_first5_original_photo_candidates.csv.gz":sha(self.prefix),
                "global_FCP_U100_incremental_rank2to5_candidates_NEED_DEDUP.csv.gz":sha(self.delta),
            },
        }
        self.receipt.write_text(json.dumps(receipt),encoding="utf-8")

    def reconcile(self):
        return m.reconcile(self.prefix,self.delta,self.receipt,self.known,self.output,
            expected_source_sha=sha(self.prefix),
            expected_delta_sha=sha(self.delta),
            expected_receipt_sha=sha(self.receipt),
            expected_taxa=2,expected_first5=10,expected_delta=8,
            expected_prior=len(self.known))

    def test_distinct_original_photo_reuse_unknown_and_fresh_work_shards(self):
        doc=self.reconcile()
        self.assertEqual(doc["new_photo_attempts_planned"],5)
        self.assertEqual(doc["photo_dispositions"]["prior_original_photo_exact_ID_reuse"],4)
        self.assertEqual(doc["photo_dispositions"]["original_photo_unsupported_extension_unknown"],1)
        self.assertEqual(doc["photo_dispositions"]["new_original_photo_ROI_to_infer"],5)
        self.assertEqual(doc["work_units"],1)
        self.assertEqual(doc["group_ids"],[0])
        self.assertEqual(doc["new_image_downloads_or_ROI_inferences_by_this_source_step"],0)
        self.assertTrue(doc["old_unknown_ROI_not_treated_as_flower_colour"])
        with (self.output/"shards/shard-000.csv").open(newline="",encoding="utf-8") as f:
            rows=list(csv.DictReader(f))
        self.assertEqual(len(rows),5)
        self.assertTrue({x["photo_id"] for x in rows}.isdisjoint(self.known))
        self.assertTrue(all(abs(float(x["selection_probability"])-0.05)<1e-15
                            for x in rows))
        self.assertEqual(sha(self.output/"shards/shard-000.csv"),
                         doc["photo_shard_sha256"]["0"])
        with gzip.open(
            self.output/"original_FCP5_all_first5_photo_reuse_and_new_dispositions.csv.gz",
            "rt",encoding="utf-8",newline="") as file:
            results=list(csv.DictReader(file))
        self.assertEqual(len(results),10)
        old_unknown=[x for x in results if x["photo_id"]==self.rows[5]["photo_id"]]
        self.assertEqual(old_unknown[0]["already_measured_hue_or_unknown"],"mixed_uncertain")

    def test_later_source_duplicate_K1_original_id_skips_reinference(self):
        self.known[self.rows[4]["photo_id"]]=(
            "source_K1",self.rows[4]["observation_id"],"mixed_uncertain","not_evaluable","")
        doc=self.reconcile()
        self.assertEqual(doc["new_photo_attempts_planned"],4)
        self.assertEqual(doc["photo_dispositions"]["prior_original_photo_exact_ID_reuse"],5)

    def test_same_original_photo_disposition_source_mutation_rejected(self):
        changed=list(self.rows)
        changed[1]=dict(changed[1])
        changed[1]["cumulative_prefix5_joint_original_photo_probability"]="0.2"
        self.gzip_csv(self.prefix,m.STAGE5_FIELDS,changed)
        self.gzip_csv(self.delta,m.STAGE5_FIELDS,changed[1:5]+changed[6:])
        receipt=json.loads(self.receipt.read_text())
        receipt["csv_outputs_SHA256"]={
            "global_FCP_U100_nested_first5_original_photo_candidates.csv.gz":sha(self.prefix),
            "global_FCP_U100_incremental_rank2to5_candidates_NEED_DEDUP.csv.gz":sha(self.delta),
        }
        self.receipt.write_text(json.dumps(receipt))
        with self.assertRaisesRegex(ValueError,"probability"):
            self.reconcile()
        self.assertFalse(self.output.exists())

    def test_unmeasured_original_first1_cannot_be_misreported_as_done(self):
        self.known.pop(self.rows[0]["photo_id"])
        with self.assertRaisesRegex(ValueError,"first photo unexpectedly unavailable"):
            self.reconcile()
        self.assertFalse(self.output.exists())

    def test_nested5_photo_duplicated_in_rank2to5_rejected(self):
        changed=[dict(x) for x in self.rows if int(x["U100_source_original_observation_rank"])>1]
        changed[1]["photo_id"]=changed[0]["photo_id"]
        self.gzip_csv(self.delta,m.STAGE5_FIELDS,changed)
        receipt=json.loads(self.receipt.read_text())
        receipt["csv_outputs_SHA256"][
            "global_FCP_U100_incremental_rank2to5_candidates_NEED_DEDUP.csv.gz"]=sha(self.delta)
        self.receipt.write_text(json.dumps(receipt))
        with self.assertRaisesRegex(ValueError,"Wrong nested"):
            self.reconcile()

    def test_immutable_original_pixel_queue_is_not_rewritten(self):
        self.reconcile()
        with self.assertRaisesRegex(ValueError,"Cannot replace"):
            self.reconcile()


if __name__=="__main__":
    unittest.main()
