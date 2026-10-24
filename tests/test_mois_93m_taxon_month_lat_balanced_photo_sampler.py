"""Full source-photo 93m tier observation-balanced selector, entirely metadata-only.

Real SHA source router + SQLite multiple photos/observation are tested as a
nested two-stage design. These are synthetic source observations, not ecological
training/evaluation photos. No provider URLs are ever fetched.
"""
from __future__ import annotations

import csv
import gzip
import io
import json
import math
from pathlib import Path
import tarfile
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from scripts.acquisition import mois_93m_taxon_month_lat_balanced_photo_sampler as s
from scripts.acquisition.mois_onepass_full93m_target_tiers import COLUMNS
from scripts.acquisition.mois_queue_from_full93m_target_tiers import WORKER_COLUMNS
from scripts.acquisition.mois_singlepass_signed_tier_bucket_router import ROUTE_SCHEMA


class BalancedSignedOriginalPhotoSampleTests(unittest.TestCase):
    def setUp(self):
        self.temp=TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.route=self.root/"route"
        self.registry=self.root/"93m.json"
        self.taxa=self.root/"taxa.csv"
        self.taxa.write_text("inat_taxon_id,species\n1,Flower A\n2,Flower B\n3,Flower C\n")
        self.output=self.root/"sample"
        self.taxa_map={"1":"Flower A","2":"Flower B","3":"Flower C"}
        self.rows_by_group={
            0:[
                ("1","obsA","1","2025-06-01","35","140","jpg"),
                ("2","obsA","1","2025-06-01","35","140","png"),
                ("3","obsA","1","2025-06-01","35","140","jpeg"),
            ],
            1:[("4","obsB","1","2025-06-05","40","141","jpg")],
            2:[
                ("5","obsC","1","2025-06-09","10","130","jpg"),
                ("6","obsC","1","2025-06-09","10","130","jpg"),
            ],
            3:[("7","obsD","2","2025-12-01","-35","150","jpeg")],
            4:[("8","obsE","3","","","","gif")],
            5:[("9","obsF","3","","","","jpg")],
        }
        self.pinned={}
        self.total=sum(len(x) for x in self.rows_by_group.values())
        for group in range(16):
            path=self.route/f"group{group:02d}"
            path.mkdir(parents=True)
            records=self.rows_by_group.get(group,[])
            buckets={}
            for pid,oid,taxon,observed,lat,lon,ext in records:
                row={key:"" for key in COLUMNS}
                row.update(photo_id=pid,observation_uuid=oid,inat_taxon_id=taxon,
                           observed_on=observed,latitude=lat,longitude=lon,
                           observer_id="test-obs-"+oid,
                           positional_accuracy="5",photo_license="cc-by",
                           image_extension=ext,observation_inclusion_score="0.01",
                           tier="6")
                buckets.setdefault(10,[]).append(row)
            files={}
            bucket_receipts=[]
            for bucket,items in buckets.items():
                mem=io.BytesIO()
                with gzip.GzipFile(fileobj=mem,mode="wb",mtime=0) as out:
                    with io.TextIOWrapper(out,encoding="utf-8",newline="") as f:
                        writer=csv.DictWriter(f,fieldnames=COLUMNS,lineterminator="\n")
                        writer.writeheader();writer.writerows(items)
                filename=f"bucket_csv/bucket-{bucket:04d}.csv.gz"
                payload=mem.getvalue()
                files[filename]=payload
                bucket_receipts.append({
                    "bucket":bucket,
                    "source_file":Path(filename).name,
                    "source_csv_gzip_sha256":s.hashlib.sha256(payload).hexdigest(),
                    "original_photo_IDs":len(items),
                })
            receipt={
                "schema":ROUTE_SCHEMA,
                "source_tier":6,"source_group":group,
                "compute_bucket_count":512,
                "full_original_source_photo_IDs":self.total,
                "native_status_not_hard_filtered":True,
                "images_downloaded":0,"flower_pixels_inferred":0,
                "source_group_read_exactly_once":True,
                "no_repeated_source_tier_scan_per_bucket":True,
                "original_source_photo_IDs_routed":len(records),
                "routed_buckets":bucket_receipts,
            }
            rp=path/"routing_receipt.json"
            rp.write_text(json.dumps(receipt,sort_keys=True))
            arc=path/f"tier06-group{group:02d}-512-source-buckets.tar"
            with tarfile.open(arc,"w") as tar:
                b=rp.read_bytes()
                header=tarfile.TarInfo("./routing_receipt.json")
                header.size=len(b)
                tar.addfile(header,io.BytesIO(b))
                for filename,payload in files.items():
                    entry=tarfile.TarInfo("./"+filename)
                    entry.size=len(payload)
                    tar.addfile(entry,io.BytesIO(payload))
            self.pinned[f"{group:02d}"]=(s.sha256(rp),s.sha256(arc))
        self.registry.write_text(json.dumps({
            "schema":s.TIER_REGISTRY_SCHEMA,
            "exact_original_photo_IDs_in_all_disjoint_tiers":self.total,
            "original_named_taxa_from_full_source":3,
            "new_target_counts_by_tier":{"6":self.total},
            "native_vs_introduced_not_hard_filtered":True,
            "source_only_NO_photo_pixels_or_hue_measured":True,
        }))

    def run_plan(self,k=1,pinned=None,max_candidates=100):
        with patch.object(s,"read_taxa",return_value=self.taxa_map),\
             patch.object(s,"obs_draw",return_value=0.01):
            return s.plan(
                self.registry,self.route,self.taxa,self.output,
                k=k,max_candidates=max_candidates,expected_photos=self.total,
                expected_taxa=3,expected_tier=self.total,
                pinned=self.pinned if pinned is None else pinned,
                check_taxa_sha=False,
            )

    def all_rows(self):
        with (self.output/"selected_observation_photo_probability_sample.csv").open(
             encoding="utf-8",newline="") as f:
            return list(csv.DictReader(f))

    def test_complete_signed_source_counted_once_and_scientific_probabilities(self):
        summary=self.run_plan(k=1)
        self.assertEqual(summary["original_tier06_source_photo_IDs_accounted"],9)
        self.assertEqual(summary["original_tier06_source_observations"],6)
        self.assertEqual(summary["original_tier06_taxon_month_latitude_strata"],4)
        self.assertEqual(summary["sampled_distinct_original_observations"],4)
        self.assertEqual(summary["unsupported_original_extension_selected_as_technical_unknown"]+
                         summary["supported_photo_attempt_targets"],4)
        rows=self.all_rows()
        self.assertEqual(len(rows),4)
        self.assertEqual(len({r["observation_id"] for r in rows}),4)
        self.assertTrue(all(float(r["joint_photo_selection_probability"])>0 for r in rows))
        for row in rows:
            n=int(row["observations_in_original_source_stratum"])
            m=int(row["photos_in_original_observation"])
            q=.0064*(1/n)*(1/m)
            self.assertAlmostEqual(float(row["joint_photo_selection_probability"]),q)
        from_earlier={"1","2","3","4","5","6","7","8","9"}
        self.assertTrue(set(r["photo_id"] for r in rows)<=from_earlier)
        self.assertTrue(summary["all_source_observations_have_strictly_positive_phase2_probability"])
        self.assertFalse(summary["no_botanical_anthesis_or_individual_wild_native_validation"])
        self.assertEqual(summary["source_sample_decided_without_photo_pixels"],True)
        self.assertEqual(summary["source_sampling_schema"] if "source_sampling_schema" in summary else
                         summary["sampling_schema"],s.SCHEMA)

        with (self.output/"manifest.json").open() as f:
            manifest=json.load(f)
        self.assertEqual(manifest["schema"],"mois_full_source_observation_sampling_cohort_v1")
        self.assertFalse(manifest["colour_measured"])
        self.assertTrue(manifest["selection_probability_is_per_photo_in_signed_shards"])
        self.assertEqual(manifest["new_photo_attempts_planned"],
                         summary["supported_photo_attempt_targets"])
        ids=[]
        for shard in sorted((self.output/"shards").glob("*.csv")):
            with shard.open(newline="") as f:
                indexed=list(csv.DictReader(f))
            self.assertEqual(tuple(indexed[0].keys()),WORKER_COLUMNS)
            for r in indexed:
                self.assertGreater(float(r["selection_probability"]),0)
                ids.append(r["photo_id"])
        self.assertEqual(len(ids),summary["supported_photo_attempt_targets"])
        self.assertEqual(len(set(ids)),len(ids))

    def test_k2_keeps_both_original_observations_in_each_populated_stratum(self):
        receipt=self.run_plan(k=2)
        self.assertEqual(receipt["sampled_distinct_original_observations"],6)
        rows=self.all_rows()
        self.assertEqual({r["observation_id"] for r in rows},
                         {"obsA","obsB","obsC","obsD","obsE","obsF"})
        self.assertTrue(all(r["source_arm"].endswith("_K2_v1") for r in rows))
        for r in rows:
            n=int(r["observations_in_original_source_stratum"])
            m=int(r["photos_in_original_observation"])
            self.assertAlmostEqual(
                float(r["joint_photo_selection_probability"]),.0064/m)

    def test_source_routing_byte_sha_tamper_fails_closed(self):
        x=self.route/"group08"/"tier06-group08-512-source-buckets.tar"
        with x.open("ab") as f:f.write(b"corruption")
        with self.assertRaisesRegex(ValueError,"routing group SHA changed"):
            self.run_plan()
        self.assertFalse(self.output.exists())

    def test_conflicting_source_observation_taxon_or_date_fails_before_publication(self):
        # Internal source-group SHA pins updated after source tampering to
        # demonstrate that repeated original observation metadata MUST agree.
        path=self.route/"group01"/"tier06-group01-512-source-buckets.tar"
        # Instead of fake SHA tampering we check exact obs contract directly.
        import sqlite3
        with sqlite3.connect(":memory:") as conn:
            s.create_database(conn)
            one={key:"" for key in COLUMNS}
            one.update(photo_id="10",observation_uuid="A",inat_taxon_id="1",
                       observed_on="2025-06-01",latitude="35",longitude="140",
                       observer_id="z",positional_accuracy="5",
                       photo_license="cc-by",image_extension="jpg",
                       observation_inclusion_score="0.01",tier="6")
            with patch.object(s,"obs_draw",return_value=0.01):
                s.import_photo(conn,one)
                changed=dict(one,photo_id="11",observed_on="2025-07-01")
                with self.assertRaisesRegex(ValueError,"Conflicting source metadata"):
                    s.import_photo(conn,changed)

    def test_stop_on_overcapacity_before_any_new_original_photo_download(self):
        with self.assertRaisesRegex(ValueError,"Too many stratified original image candidates"):
            self.run_plan(k=1,max_candidates=2)
        self.assertFalse(self.output.exists())

    def test_output_immutable_and_one_observation_has_only_one_selected_photo(self):
        self.run_plan()
        self.assertEqual(len({r["observation_id"] for r in self.all_rows()}),4)
        with self.assertRaisesRegex(ValueError,"may not be overwritten"):
            self.run_plan()


if __name__ == "__main__":
    unittest.main()
