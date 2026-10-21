"""Bounded-memory signed source group router: no repeated 93m scan, no hue."""
import csv
import gzip
import json
from pathlib import Path
import tarfile
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from scripts.acquisition import mois_singlepass_signed_tier_bucket_router as router
from scripts.acquisition.mois_onepass_full93m_target_tiers import COLUMNS


class SinglepassTierRoutingTests(unittest.TestCase):
    def setUp(self):
        self.temp=TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.root.mkdir(exist_ok=True)
        self.registry=self.root/"registry.json"
        self.receipt=self.root/"group-00-receipt.json"
        self.archive=self.root/"tier-06-group-00.tar"
        self.output=self.root/"bucketed"
        self.inputs=self.root/"input/tier-06"
        self.inputs.mkdir(parents=True)
        self.source_rows=[
            ("1","O1","jpg","2025-05-01"),
            ("2","O2","jpg","2025-05-01"),
            ("3","O1","gif",""),
            ("4","O3","jpg","2025-09-01"),
        ]
        part_receipts=[]
        for part,subset in enumerate((self.source_rows[:3],self.source_rows[3:])):
            path=self.inputs/f"part-{part:03d}.csv.gz"
            with gzip.open(path,"wt",encoding="utf-8",newline="") as f:
                w=csv.DictWriter(f,fieldnames=COLUMNS,lineterminator="\n")
                w.writeheader()
                for pid,obs,ext,date in subset:
                    row={col:"" for col in COLUMNS}
                    row.update(photo_id=pid,observation_uuid=obs,inat_taxon_id="1",
                               image_extension=ext,observed_on=date,
                               observation_inclusion_score="0.007",tier="6",
                               frozen_loader_supported="false" if ext=="gif" else "true")
                    w.writerow(row)
            part_receipts.append({
                "part":part,
                "counts_by_disjoint_tier":{"6":len(subset)},
                "ranked_target_csv_sha256":{"6":router.sha256(path)}
            })
        with tarfile.open(self.archive,"w") as t:
            for i in range(2):
                t.add(self.inputs/f"part-{i:03d}.csv.gz",
                      arcname=f"tier-06/part-{i:03d}.csv.gz")
        self.receipt.write_text(json.dumps({
            "schema":router.RESULT_SCHEMA,
            "source_snapshot_SHA256":"a"*64,
            "all_source_photo_IDs_retained_in_one_and_only_one_tier":True,
            "source_original_photo_IDs_uniqueness_proved_by_upstream_primary_key":True,
            "same_observation_uses_same_selection_score":True,
            "images_downloaded":0,
            "colour_pixels_measured":0,
            "counts_by_disjoint_tier":{"6":4},
            "parts":part_receipts
        }))
        self.registry.write_text(json.dumps({
            "schema":"mois_complete_full93m_original_photo_target_registry_v1",
            "exact_original_photo_IDs_in_all_disjoint_tiers":4,
            "original_named_taxa_from_full_source":2,
            "source_only_NO_photo_pixels_or_hue_measured":True,
            "native_vs_introduced_not_hard_filtered":True,
            "full_source_manifest_sha256":"a"*64,
            "published_tier_archive_SHA256":{self.archive.name:router.sha256(self.archive)}
        }))

    def route(self):
        with patch.object(router,"obs_draw",lambda _:0.007):
            return router.route_signed_group(
                self.registry,self.receipt,self.archive,self.output,
                bucket_count=16,max_open=1,expected_total=4,
                expected_taxa=2,parts_per_group=2)

    def test_singlepass_every_photo_id_retained_observations_share_bucket(self):
        report=self.route()
        self.assertEqual(report["original_source_photo_IDs_routed"],4)
        self.assertTrue(report["source_group_read_exactly_once"])
        self.assertTrue(report["no_repeated_source_tier_scan_per_bucket"])
        self.assertFalse(report["native_status_not_hard_filtered"] is False)
        ids=[]
        locations={}
        for info in report["routed_buckets"]:
            file=self.output/"bucket_csv"/info["source_file"]
            self.assertEqual(router.sha256(file),info["source_csv_gzip_sha256"])
            with gzip.open(file,"rt",encoding="utf-8",newline="") as f:
                rows=list(csv.DictReader(f))
            self.assertEqual(len(rows),info["original_photo_IDs"])
            for row in rows:
                ids.append(row["photo_id"])
                locations[row["photo_id"]]=info["bucket"]
        self.assertEqual(sorted(ids),["1","2","3","4"])
        self.assertEqual(locations["1"],locations["3"])
        # Ext / missing day are metadata uncertainty, never auto exclusion.
        self.assertEqual(report["images_downloaded"],0)
        self.assertEqual(report["flower_pixels_inferred"],0)

    def test_source_tar_tamper_rejected_without_partial_routing(self):
        with self.archive.open("ab") as f:f.write(b"damage")
        with self.assertRaisesRegex(ValueError,"archive SHA changed"):
            self.route()
        self.assertFalse(self.output.exists())

    def test_same_output_is_immutable(self):
        self.route()
        with self.assertRaisesRegex(ValueError,"not be overwritten"):
            self.route()

    def test_incomplete_source_part_signatures_fail(self):
        body=json.loads(self.receipt.read_text())
        body["parts"]=body["parts"][:1]
        self.receipt.write_text(json.dumps(body))
        with self.assertRaisesRegex(ValueError,"missing fixed original"):
            self.route()
        self.assertFalse(self.output.exists())


if __name__ == "__main__":
    unittest.main()
