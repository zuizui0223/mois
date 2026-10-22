"""Source-only real 93m routed-bucket -> fixed FCP warm queue contracts."""
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

from scripts.acquisition import mois_routed_tier06_warm_photo_queue as q
from scripts.acquisition.mois_onepass_full93m_target_tiers import COLUMNS, obs_draw
from scripts.acquisition.mois_queue_from_full93m_target_tiers import observation_compute_bucket


class Tier06WarmQueueTests(unittest.TestCase):
    def setUp(self):
        self.tmp=TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.receipt=self.root/"routing_receipt.json"
        self.archive=self.root/"bucketed.tar"
        self.old=self.root/"old.parquet"
        self.old.write_bytes(b"fake synthetic old test data")
        self.oldreceipt=self.root/"old.json"
        self.oldreceipt.write_text("{}")
        self.technical=self.root/"technical.csv"
        self.technical.write_text("photo_id,observation_id,acquired,status\n")
        self.names=self.root/"names.csv"
        self.names.write_text("inat_taxon_id,species\n1,Planta test\n2,Planta second\n")
        self.out=self.root/"current"

        obs={};found=set()
        for i in range(100000):
            o=f"OBS-{i}"
            b=observation_compute_bucket(o,512)
            if len(found)<4 and b<64 and b not in found:
                obs[len(found)]=(o,b)
                found.add(b)
            elif b>=64 and "outside" not in obs:
                obs["outside"]=(o,b)
            if len(found)==4 and "outside" in obs:
                break
        self.obs=obs
        specs=[
            (1,obs[0][0],"jpg"),(2,obs[0][0],"jpeg"),
            (3,obs[1][0],"jpg"),(4,obs[2][0],"png"),
            (5,obs[3][0],"gif"),(6,obs["outside"][0],"jpg"),
        ]
        by_bucket={}
        for pid,oid,ext in specs:
            bucket=observation_compute_bucket(oid,512)
            row={key:"" for key in COLUMNS}
            row.update(photo_id=str(pid),observation_uuid=oid,
                       inat_taxon_id="1",image_extension=ext,
                       tier="6",observed_on="2025-07-01",
                       observation_inclusion_score=format(obs_draw(oid),".17g"),
                       photo_license="cc-by",observer_id="observer1")
            by_bucket.setdefault(bucket,[]).append(row)
        zipped={}
        entries=[]
        for bucket,rows in sorted(by_bucket.items()):
            buf=io.BytesIO()
            with gzip.GzipFile(fileobj=buf,mode="wb",mtime=0) as g:
                with io.TextIOWrapper(g,encoding="utf-8",newline="") as f:
                    w=csv.DictWriter(f,fieldnames=COLUMNS,lineterminator="\n")
                    w.writeheader();w.writerows(rows)
            filename=f"bucket-{bucket:04d}.csv.gz"
            contents=buf.getvalue()
            zipped["bucket_csv/"+filename]=contents
            entries.append({"bucket":bucket,"original_photo_IDs":len(rows),
                            "source_csv_gzip_sha256":q.sha_bytes(contents),
                            "source_file":filename})
        report={
            "schema":q.ROUTE_SCHEMA,"source_tier":6,"source_group":0,
            "compute_bucket_count":512,"full_original_source_photo_IDs":6,
            "native_status_not_hard_filtered":True,"images_downloaded":0,
            "flower_pixels_inferred":0,"routed_buckets":entries,
            "original_source_photo_IDs_routed":len(specs),
            "full_source_registry_sha256":"a"*64
        }
        self.receipt.write_text(json.dumps(report,sort_keys=True))
        with tarfile.open(self.archive,"w") as t:
            contents=self.receipt.read_bytes()
            info=tarfile.TarInfo("routing_receipt.json")
            info.size=len(contents);t.addfile(info,io.BytesIO(contents))
            for filename,contents in sorted(zipped.items()):
                info=tarfile.TarInfo(filename)
                info.size=len(contents)
                t.addfile(info,io.BytesIO(contents))

    def go(self):
        with patch.object(q,"read_taxa",return_value={"1":"Planta test","2":"Planta second"}),\
             patch.object(q,"previous_ids",return_value=({"3"},{self.obs[1][0]})),\
             patch.object(q,"technical_reservations",return_value=({"4"},{self.obs[2][0]})):
            return q.build_queue(
                self.receipt,self.archive,self.names,
                self.old,self.oldreceipt,self.technical,self.out,
                expected_route_sha=None,expected_archive_sha=None,
                expected_prior_receipt_sha=None,expected_technical_sha=None,
                expected_total=6,expected_taxa=2,check_taxa_sha=False,
                first_buckets=64,max_new=10)

    def test_every_source_ID_accounted_and_only_new_supported_photos_queued(self):
        r=self.go()
        self.assertEqual(r["source_group_original_photo_IDs"],6)
        self.assertEqual(r["selected_bucket_photo_IDs_before_dedup"],5)
        self.assertEqual(r["new_photo_attempts_planned"],2)
        self.assertEqual(r["selected_new_photos_unsupported_by_frozen_loader"],1)
        self.assertEqual(r["prior_and_technical_excluded_source_photo_IDs"]["prior_actual_photo_or_observation"],1)
        self.assertEqual(r["prior_and_technical_excluded_source_photo_IDs"]["prior_technical_photo_or_observation"],1)
        self.assertEqual(r["group_ids"],[0])
        self.assertFalse(r["colour_measured"])
        self.assertTrue(r["not_probability_complete_global_source"])
        f=self.out/"shards/shard-000.csv"
        self.assertEqual(q.sha256(f),r["photo_shard_sha256"]["0"])
        rows=list(csv.DictReader(f.open(newline="")))
        self.assertEqual({x["photo_id"] for x in rows},{"1","2"})
        self.assertEqual(len({x["observation_id"] for x in rows}),1)
        self.assertEqual({x["selection_probability"] for x in rows},
                         {format(.0064*64/512,".17g")})
        self.assertEqual(json.loads((self.out/"unsupported_selected_original_photo_formats.jsonl").read_text().splitlines()[0])["photo_id"],"5")

    def test_repeated_same_source_cannot_rewrite_frozen_output(self):
        self.go()
        with self.assertRaisesRegex(ValueError,"immutable"):
            self.go()

    def test_bucket_sha_drift_fails_before_creating_queue(self):
        with tarfile.open(self.archive,"a") as t:
            # Another copy of a trusted file is an invalid duplicate tar member.
            path=self.root/"duplicate"
            path.write_bytes(b"replaced bucket data")
            t.add(path,arcname="bucket_csv/"+self.receipt.name)
        with self.assertRaisesRegex(ValueError,"Unknown/duplicate"):
            self.go()
        self.assertFalse(self.out.exists())

    def test_corrupt_source_manifest_counts_rejected(self):
        r=json.loads(self.receipt.read_text())
        r["original_source_photo_IDs_routed"]=99
        self.receipt.write_text(json.dumps(r))
        with self.assertRaisesRegex(ValueError,"manifest inconsistent"):
            self.go()
        self.assertFalse(self.out.exists())


if __name__=="__main__":
    unittest.main()
