"""All-16 signed original tier06 photo-source queue synthetic identity tests."""
from __future__ import annotations
import csv
import gzip
import io
import json
from pathlib import Path
import tarfile
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from scripts.acquisition import mois_tier06_all16_broad_warm_photo_queue as q
from scripts.acquisition.mois_onepass_full93m_target_tiers import COLUMNS, obs_draw
from scripts.acquisition.mois_queue_from_full93m_target_tiers import observation_compute_bucket


class BroadAll16WarmSourceTests(unittest.TestCase):
    def setUp(self):
        self.temp=TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.registry=self.root/"all.json"
        self.routed=self.root/"route"
        self.routed.mkdir()
        self.output=self.root/"queue"
        self.main=self.root/"main.parquet";self.main.write_bytes(b"synthetic-main")
        self.mainreceipt=self.root/"main.json";self.mainreceipt.write_text("{}")
        self.group00=self.root/"group00.parquet";self.group00.write_bytes(b"synthetic-group00")
        self.group00receipt=self.root/"group00.json";self.group00receipt.write_text("{}")
        self.tech=self.root/"technical.csv";self.tech.write_text("photo_id,observation_id,acquired,status\n")
        self.taxnames=self.root/"taxnames.csv";self.taxnames.write_text("inat_taxon_id,species\n1,Plant one\n2,Plant two\n")
        self.bucket_obs={}
        self.pins={}
        self.registry.write_text(json.dumps({
            "schema":"mois_complete_full93m_original_photo_target_registry_v1",
            "exact_original_photo_IDs_in_all_disjoint_tiers":48,
            "new_target_counts_by_tier":{"6":48},
            "native_vs_introduced_not_hard_filtered":True
        }))
        original_sha=q.sha256(self.registry)
        for n in range(16):
            folder=self.routed/f"group{n:02d}";folder.mkdir()
            buckets={}
            for desired in (0,1,6):
                i=0
                while True:
                    obs=f"O-{n}-{desired}-{i}"
                    if observation_compute_bucket(obs,512)==desired:
                        break
                    i+=1
                self.bucket_obs[(n,desired)]=obs
                pid=10_000+n*10+desired
                row={c:"" for c in COLUMNS}
                row.update(photo_id=str(pid),observation_uuid=obs,
                           inat_taxon_id=str(1+(n%2)),tier="6",
                           image_extension="gif" if n==4 and desired==0 else "jpg",
                           observation_inclusion_score=format(obs_draw(obs),".17g"),
                           observed_on="2025-07-01",photo_license="cc-by")
                stream=io.BytesIO()
                with gzip.GzipFile(fileobj=stream,mode="wb",mtime=0) as gz:
                    with io.TextIOWrapper(gz,encoding="utf-8",newline="") as out:
                        writer=csv.DictWriter(out,fieldnames=COLUMNS,lineterminator="\n")
                        writer.writeheader();writer.writerow(row)
                data=stream.getvalue()
                buckets[desired]=(f"bucket-{desired:04d}.csv.gz",data)
            entries=[{
                "bucket":b,"original_photo_IDs":1,
                "source_csv_gzip_sha256":q.hashlib.sha256(contents).hexdigest(),
                "source_file":file}
                for b,(file,contents) in sorted(buckets.items())]
            receipt={
                "schema":q.ROUTE_SCHEMA,
                "source_tier":6,"source_group":n,"compute_bucket_count":512,
                "full_original_source_photo_IDs":48,
                "full_source_registry_sha256":original_sha,
                "source_group_read_exactly_once":True,
                "images_downloaded":0,"flower_pixels_inferred":0,
                "native_status_not_hard_filtered":True,
                "original_source_photo_IDs_routed":3,
                "routed_buckets":entries
            }
            receiptfile=folder/"routing_receipt.json"
            receiptfile.write_text(json.dumps(receipt,sort_keys=True))
            archive=folder/f"tier06-group{n:02d}-512-source-buckets.tar"
            with tarfile.open(archive,"w") as tf:
                raw=receiptfile.read_bytes()
                h=tarfile.TarInfo("routing_receipt.json");h.size=len(raw)
                tf.addfile(h,io.BytesIO(raw))
                for index,(name,data) in sorted(buckets.items()):
                    h=tarfile.TarInfo("bucket_csv/"+name);h.size=len(data)
                    tf.addfile(h,io.BytesIO(data))
            self.pins[f"{n:02d}"]=(q.sha256(receiptfile),q.sha256(archive))

    def build(self,pins=None):
        main_photo={"10020"}
        main_obs={self.bucket_obs[(2,0)]}
        group00_photo={"10000","10001"}
        group00_obs={self.bucket_obs[(0,0)],self.bucket_obs[(0,1)]}
        with patch.object(q,"read_taxa",return_value={"1":"Plant one","2":"Plant two"}),\
             patch.object(q,"identity_from_parquet",side_effect=[
                 (main_photo,main_obs),(group00_photo,group00_obs)]),\
             patch.object(q,"technical_reservations",return_value=(
                 {"10031"},{self.bucket_obs[(3,1)]})):
            return q.build(self.registry,self.routed,self.main,self.mainreceipt,
                           self.group00,self.group00receipt,self.tech,self.taxnames,self.output,
                           pinned=pins or self.pins,expected_full=48,expected_taxa=2,
                           expected_tier_total=48,first_buckets=4,
                           check_taxa_sha=False,max_new=100)

    def test_every_signed_taxon_group_and_original_photo_unit_accounted(self):
        m=self.build()
        self.assertEqual(m["complete_tier06_source_photo_IDs"],48)
        self.assertEqual(len(m["all16_immutable_signed_source_group_counts"]),16)
        self.assertEqual(m["metadata_original_photos_in_admitted_buckets"],32)
        self.assertEqual(m["new_photo_attempts_planned"],27)
        self.assertEqual(m["work_units"],1)
        self.assertEqual(m["group_ids"],[0])
        self.assertTrue(m["not_probability_complete_global_source"])
        self.assertFalse(m["colour_measured"])
        path=self.output/"shards"/"shard-000.csv"
        self.assertEqual(q.sha256(path),m["photo_shard_sha256"]["0"])
        photos=list(csv.DictReader(path.open(newline="")))
        self.assertEqual(len(photos),27)
        self.assertEqual(len({p["photo_id"] for p in photos}),27)
        self.assertNotIn("10020",{p["photo_id"] for p in photos})
        self.assertNotIn("10031",{p["photo_id"] for p in photos})
        self.assertNotIn("10040",{p["photo_id"] for p in photos})
        self.assertNotIn("10000",{p["photo_id"] for p in photos})
        self.assertEqual({p["selection_probability"] for p in photos},
                         {format((.0128-.0064)*4/512,".17g")})
        self.assertEqual(m["selected_new_photos_unsupported_by_frozen_loader"],1)

    def test_original_signed_group_archive_tampering_fails_closed(self):
        damaged=self.pins.copy()
        a,b=damaged["07"];damaged["07"]=(a,"0"*64)
        with self.assertRaisesRegex(ValueError,"Release digest mismatch"):
            self.build(pins=damaged)
        self.assertFalse(self.output.exists())

    def test_immutable_repeat_queue_refuses_to_reselect_photos(self):
        self.build()
        with self.assertRaisesRegex(ValueError,"Cannot change original"):
            self.build()


if __name__=="__main__":
    unittest.main()
