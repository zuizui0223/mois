"""Nested original photo ranks 2..5 exact prior-ID reuse and source-only SHA tests."""
from __future__ import annotations
import csv
import gzip
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from scripts.acquisition import mois_fcp_u100_depth5_new_original_photo_queue as q
from scripts.acquisition.mois_completed_pixel_checkpoint import sha


class NestedFiveNewOriginalPixelQueue(unittest.TestCase):
    def setUp(self):
        self.temp=TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.csv=self.root/"source.csv.gz"
        self.receipt=self.root/"source.json"
        self.output=self.root/"output"
        self.already={
            "1002":("old","OBS-2","mixed_uncertain","unclassified",""),
        }
        self.rows=[]
        for rank in range(2,6):
            pid=str(1000+rank)
            source_ext="heic" if rank==3 else "jpg"
            supported=rank!=3
            data={k:"" for k in q.SOURCE_FIELDS}
            data.update(
                photo_id=pid,observation_id=f"OBS-{rank}",
                inat_taxon_id="1",species="Synthetic plant",
                observer_id="photographer",observed_on="2024-08-07",
                recorded_month="8",latitude_band="n",latitude="35",
                longitude="139",positional_accuracy="40",photo_license="cc-by",
                original_image_extension=source_ext,
                photos_in_original_observation="1",
                full_93m_distinct_original_observations_of_taxon="100",
                U100_source_original_observation_rank=str(rank),
                source_rank1_original_photo_id="1001",
                cumulative_prefix5_joint_original_photo_probability=".05",
                rank2to5_incremental_joint_original_photo_probability=".04",
                depth1_joint_original_photo_inclusion_probability=".01",
                parent_U100_joint_original_photo_inclusion_probability="1",
                prefix_depth="5",incremental_since_depth1="true",
                source_arm=q.NESTED_ARM,
                parent_original_photo_U100_source_arm="mois_full93m_fcp_exact100_original_observation_per_taxon_v1",
                photo_source_disposition="supported" if supported else "unsupported_original_extension",
                photo_url_large=(
                    f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.jpg"
                    if supported else ""),
                independent_accepted_species_and_native_stage_verified="false",
                cumulative_prefix5_observation_probability=".05",
                rank2to5_incremental_observation_probability=".04",
            )
            self.rows.append(data)
        self.write()

    def write(self):
        with gzip.open(self.csv,"wt",encoding="utf-8",newline="") as stream:
            writer=csv.DictWriter(stream,fieldnames=q.SOURCE_FIELDS,lineterminator="\n")
            writer.writeheader();writer.writerows(self.rows)
        receipt={
            "schema":q.GLOBAL_SCHEMA,
            "eligible_taxa_with_at_least100_original_observations":1,
            "global_incremental_rank2to5_candidate_photo_IDs_NOT_NEW_INFERENCES":4,
            "global_first5_nested_source_photo_IDs":5,
            "global_first1_original_photo_IDs_unchanged":1,
            "ranks2to5_genuine_never_measured_photo_IDs_not_yet_deduplicated":True,
            "source_probabilities_original_5_over_N_over_M_and_delta_4_over_N_over_M":True,
            "new_image_downloads_and_original_ROI_inferences_performed":0,
            "csv_outputs_SHA256":{
                "global_FCP_U100_incremental_rank2to5_candidates_NEED_DEDUP.csv.gz":sha(self.csv)
            },
        }
        self.receipt.write_text(json.dumps(receipt),encoding="utf-8")

    def freeze(self,**kwargs):
        cfg=dict(expected_delta=4,expected_taxa=1,
                 expected_source_sha=sha(self.csv),
                 expected_receipt_sha=sha(self.receipt))
        cfg.update(kwargs)
        return q.freeze_new_rank2to5(self.csv,self.receipt,self.already,
                                     self.output,**cfg)

    def test_keep_only_two_genuinely_new_pixels_and_4_over_N_M(self):
        m=self.freeze()
        self.assertEqual(m["new_rank2to5_source_photo_candidates"],4)
        self.assertEqual(m["new_original_photo_IDs_selected_never_previously_measured"],2)
        self.assertEqual(m["original_photo_IDs_with_existing_ROI_reused"],1)
        self.assertEqual(m["original_photo_unsupported_formats_technical_unknown"],1)
        self.assertEqual(m["old_existing_global_real_original_photo_attempts"],1)
        self.assertTrue(m["source_U100_incremental_rank2to5_photo_probability_is_4_over_N_over_M"])
        self.assertFalse(m["true_wild_native_blooming_plant_individuals_not_verified"] is False)
        self.assertEqual(m["original_images_downloaded_by_this_source_process"],0)
        self.assertEqual(m["new_four_hue_ROIs_inferred_by_this_source_process"],0)
        part=self.output/"shards/shard-000.csv"
        with part.open(newline="",encoding="utf-8") as handle:
            rows=list(csv.DictReader(handle))
        self.assertEqual({r["photo_id"] for r in rows},{"1004","1005"})
        self.assertEqual({r["selection_probability"] for r in rows},{".04"})
        self.assertEqual(sha(part),m["worker_photo_shard_SHA256"]["0"])

    def test_old_technically_unclassifiable_pixel_does_not_get_redownloaded(self):
        self.freeze()
        file=self.output/"previous_original_image_photo_ID_reuse.csv"
        with file.open(newline="",encoding="utf-8") as f:
            reuse=list(csv.DictReader(f))
        self.assertEqual(len(reuse),1)
        self.assertEqual(reuse[0]["prior_automatic_four_hue_morph"],"mixed_uncertain")
        self.assertEqual(reuse[0]["prior_image_arm"],"old")

    def test_source_probability_mutation_rejected(self):
        self.rows[2]["rank2to5_incremental_joint_original_photo_probability"]=".2"
        self.write()
        with self.assertRaisesRegex(ValueError,"selection probability changed"):
            self.freeze()

    def test_duplicate_original_observation_UUID_rejected(self):
        self.rows[3]["observation_id"]=self.rows[2]["observation_id"]
        self.write()
        with self.assertRaisesRegex(ValueError,"Duplicate source photograph"):
            self.freeze()

    def test_unsupported_extension_becomes_unknown_not_a_new_photo_job(self):
        self.rows[1]["photo_source_disposition"]="supported"
        self.write()
        with self.assertRaisesRegex(ValueError,"technical extension state changed"):
            self.freeze()

    def test_mutated_signed_source_fails_before_new_work(self):
        self.write()
        with self.assertRaisesRegex(ValueError,"checksum changed"):
            self.freeze(expected_source_sha="0"*64)
        self.assertFalse(self.output.exists())

    def test_signed_genuine_ROI_worker_source_is_immutable(self):
        self.freeze()
        with self.assertRaisesRegex(ValueError,"immutable"):
            self.freeze()


if __name__=="__main__":
    unittest.main()
