"""Exact K1 global photo-ID/hue counts from compact old-vs-source SHA proof."""
from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from scripts.acquisition import mois_k1_compact_global_pixel_totals as m
from scripts.acquisition.mois_completed_pixel_checkpoint import sha


class TestCompactSourceAuthenticatedGlobalOriginalROI(unittest.TestCase):
    def setUp(self):
        self.temp=TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.cert=self.root/"sourcecert.json"
        self.gz=self.root/"source.csv.gz"
        self.cum=self.root/"cumulative.json"
        self.ref=self.root/"old16.json"
        self.out=self.root/"result"
        self.rows=[
            dict(photo_id="21",observation_id="new-A",inat_taxon_id="1",
                 measurement_id="A",selection_probability="0.02",original_worker_group="0"),
            dict(photo_id="22",observation_id="new-B",inat_taxon_id="2",
                 measurement_id="B",selection_probability="0.005",original_worker_group="0"),
        ]
        self.make_files()

    def make_files(self):
        with self.gz.open("wb") as binary:
            with gzip.GzipFile(fileobj=binary,filename="",mode="wb",mtime=0) as gz:
                with io.TextIOWrapper(gz,encoding="utf-8",newline="") as t:
                    w=csv.DictWriter(t,fieldnames=m.COLUMNS,lineterminator="\n")
                    w.writeheader();w.writerows(self.rows)
        cert={
            "schema":m.CERT_SCHEMA,
            "frozen_K1_original_source_manifest_SHA256":m.SOURCE_SHA,
            "historical_verified_unique_original_photo_IDs":1,
            "historical_verified_four_hue_original_photo_ROIs":1,
            "K1_all_new_precolour_source_original_photo_IDs":2,
            "K1_all_new_source_original_observations":2,
            "K1_exact_historical_photo_ID_collisions":0,
            "K1_source_group_count":1,
            "K1_new_source_photo_roster_sha256":sha(self.gz),
            "K1_new_source_photo_counts_by_group":{"0":2},
            "source_arm":m.SOURCE_ARM,
            "source_selected_before_any_K1_colour_inference":True,
            "original_photo_joint_selection_probability_retained":True,
            "prior_original_pixels_recomputed_or_refetched":False,
            "actual_photo_bytes_downloaded_or_model_inference_in_this_certificate":0,
            "independently_verified_open_flower_native_wild_status":False,
            "global_white_flora_fraction_unweighted_estimated":False,
        }
        self.cert.write_text(json.dumps(cert),encoding="utf-8")
        sha0="a"*64
        fingerprint=hashlib.sha256(f"0:{sha0}".encode()).hexdigest()
        cumul={
            "schema":m.CUM_SCHEMA,"source_run":m.SOURCE_RUN,
            "new_original_photo_target_total":2,
            "previously_measured_K1_photo_IDs_reused_elsewhere_NOT_NEW":m.REUSED_OLD_ORIGINALS,
            "total_scheduled_original_worker_groups":1,
            "source_precolour_probability_photo_sampling_preserved":True,
            "within_K1_original_source_group_units_disjoint":True,
            "new_original_images_downloaded_or_colour_inferred_in_THIS_ledger":0,
            "this_is_NOT_a_crossarm_photo_ID_deduplicated_global_total":True,
            "global_wild_flora_or_blooming_petal_status_not_independently_verified":True,
            "completed_original_worker_group_IDs":[0],
            "completed_and_SHA_sealed_original_worker_groups":1,
            "completed_real_new_original_photo_attempts":2,
            "completed_automatic_four_hue_classifications":1,
            "completed_original_acquisition_failures":1,
            "completed_acquired_but_automatic_four_hue_unclassifiable":0,
            "source_original_group_receipt_aggregate_SHA256":fingerprint,
            "per_group":[{
                "original_worker_group":0,"new_original_photo_attempts":2,
                "four_hue_classified":1,"image_acquisition_failed":1,
                "acquired_but_unclassifiable":0,
                "original_hue_counts":{
                    "white":1,"yellow_orange":0,"red_pink":0,
                    "blue_purple":0,"mixed_uncertain":1},
                "original_group_receipt_SHA256":sha0,
                "distinct_conditional_joint_photo_probabilities":2,
            }]
        }
        self.cum.write_text(json.dumps(cumul),encoding="utf-8")
        reference={
            "K1_source_finished_original_worker_groups":[0],
            "historical_actual_unique_original_photo_IDs":1,
            "historical_automatic_four_hue_photo_ROIs":1,
            "K1_vs_old_crossarm_exact_photo_ID_overlaps":0,
            "new_K1_sha_verified_original_photo_IDs_disjoint_from_history":2,
            "new_K1_frozen_four_hue_ROIs":1,
            "exact_global_historical_plus_K1_distinct_original_photo_IDs":3,
            "exact_global_historical_plus_K1_automatic_four_hue_photo_ROIs":2,
        }
        self.ref.write_text(json.dumps(reference),encoding="utf-8")

    def verify(self,**overrides):
        args=dict(
            certificate_sha=sha(self.cert),
            roster_sha=sha(self.gz),
            reference_sha=sha(self.ref),
            expected_historical=1,expected_hues=1,
            expected_new_source=2,expected_groups=1,reference_groups=1,
        )
        args.update(overrides)
        return m.aggregate(self.cert,self.gz,self.cum,self.ref,self.out,**args)

    def test_compact_authentic_original_source_reproduces_full_parquet_reference(self):
        d=self.verify()
        self.assertEqual(d["exact_worldwide_historical_plus_completed_K1_distinct_original_photo_IDs"],3)
        self.assertEqual(d["exact_worldwide_historical_plus_completed_K1_four_hue_ROIs"],2)
        self.assertEqual(d["K1_original_photo_acquisition_failed"],1)
        self.assertEqual(d["K1_original_acquired_unclassifiable"],0)
        self.assertTrue(d["zero_all_source_K1_vs_historical_photo_ID_overlap"])
        self.assertTrue(d["counts_reproduce_independent_full_929138_ID_union"])
        self.assertFalse(d["historical_929138_photo_Parquets_re_downloaded_for_this_update"])
        self.assertEqual(d["additional_raw_image_downloads_or_ROI_inference_during_this_update"],0)
        self.assertFalse(d["independent_botanical_open_blooming_wild_native_photo_status_verified"])

    def test_missing_or_forged_cumulative_classifiable_hue_rejected(self):
        d=json.loads(self.cum.read_text())
        d["per_group"][0]["four_hue_classified"]=2
        self.cum.write_text(json.dumps(d))
        with self.assertRaisesRegex(ValueError,"terminal/hue/unknown"):
            self.verify()
        self.assertFalse(self.out.exists())

    def test_new_photo_group_cannot_be_counted_twice(self):
        d=json.loads(self.cum.read_text())
        d["per_group"].append(d["per_group"][0])
        d["completed_and_SHA_sealed_original_worker_groups"]=2
        d["completed_original_worker_group_IDs"]=[0,0]
        self.cum.write_text(json.dumps(d))
        with self.assertRaisesRegex(ValueError,"repeated"):
            self.verify()

    def test_photo_roster_reorder_or_duplicate_is_refused_even_if_resigned(self):
        self.rows.reverse()
        self.make_files()
        with self.assertRaisesRegex(ValueError,"not original sorted unique"):
            self.verify()
        self.assertFalse(self.out.exists())

    def test_photo_gzip_bit_tampering_fails_before_any_publication(self):
        data=self.gz.read_bytes()
        self.gz.write_bytes(data[:-1]+bytes([data[-1]^1]))
        with self.assertRaisesRegex(ValueError,"Pinned independently published"):
            self.verify(roster_sha="0"*64)

    def test_reference_full_16group_independent_photocount_must_agree(self):
        d=json.loads(self.ref.read_text())
        d["exact_global_historical_plus_K1_distinct_original_photo_IDs"]=9
        self.ref.write_text(json.dumps(d))
        with self.assertRaisesRegex(ValueError,"independent full-ID audit"):
            self.verify()

    def test_immutable_publication_refuses_rewrite(self):
        self.verify()
        with self.assertRaisesRegex(ValueError,"Never overwrite"):
            self.verify()


if __name__=="__main__":
    unittest.main()
