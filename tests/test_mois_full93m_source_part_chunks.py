"""Prove source-complete, bounded 4096-photo chunks before large 93m runs."""
from pathlib import Path
from tempfile import TemporaryDirectory
import csv
import gzip
import json
import unittest

from scripts.acquisition import mois_full93m_source_part_chunks as q
from scripts.acquisition.mois_onepass_full93m_target_tiers import COLUMNS

class SplitSignedOriginalPhotoTests(unittest.TestCase):
    def setUp(self):
        self.tmp = TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.source = self.root / "signed_part.csv.gz"
        self.output = self.root / "split"
        self.records=[]
        for i in range(15):
            row={k:"" for k in COLUMNS}
            row.update(photo_id=str(500+i), observation_uuid="obs-"+str(i),
                       inat_taxon_id="123", observed_on="2025-07-01",
                       image_extension="jpg", tier="3")
            self.records.append(row)
        self.write_source()

    def write_source(self):
        with gzip.open(self.source, "wt", newline="", encoding="utf-8") as fh:
            w=csv.DictWriter(fh,fieldnames=COLUMNS,lineterminator="\n")
            w.writeheader()
            w.writerows(self.records)

    def execute(self, **options):
        return q.run(self.source, self.output,
                     expected_source_sha256=q.sha(self.source),
                     expected_photo_count=15, tier=3, part=6,
                     chunk_size=4, **options)

    def test_all_unique_original_ids_exact_once_and_sha_proof(self):
        d=self.execute()
        self.assertEqual(d["chunks"],4)
        self.assertEqual(d["actual_original_photo_IDs"],15)
        self.assertEqual(d["image_pixels_downloaded"],0)
        self.assertFalse(d["floral_hue_or_bloom_inferred"])
        self.assertTrue(d["source_fully_partitioned_once_into_nonoverlapping_photo_chains"])
        with (self.output/"chunk_manifest.csv").open() as fh:
            manifest=list(csv.DictReader(fh))
        ids=[]
        for c, expected in zip(manifest,[4,4,4,3]):
            self.assertEqual(int(c["expected_photo_ids"]),expected)
            p=self.output/c["chunk_relative_path"]
            self.assertEqual(q.sha(p),c["chunk_gzip_sha256"])
            with gzip.open(p,"rt") as handle:
                rows=list(csv.DictReader(handle))
            self.assertEqual(len(rows),expected)
            self.assertEqual(c["first_photo_id"],rows[0]["photo_id"])
            self.assertEqual(c["last_photo_id"],rows[-1]["photo_id"])
            ids.extend(r["photo_id"] for r in rows)
        self.assertEqual(ids,[str(500+i) for i in range(15)])
        self.assertEqual(len(set(ids)),15)

    def test_refuse_duplicate_photo_id(self):
        self.records[-1]["photo_id"]=self.records[0]["photo_id"]
        self.write_source()
        with self.assertRaisesRegex(ValueError,"Repeated original photo"):
            self.execute()
        self.assertFalse(self.output.exists())

    def test_refuse_unmatching_source_count(self):
        with self.assertRaisesRegex(ValueError,"does not reconcile"):
            q.run(self.source,self.output,expected_source_sha256=q.sha(self.source),
                  expected_photo_count=18,tier=3,part=6,chunk_size=4)
        self.assertFalse(self.output.exists())

    def test_refuse_tier_or_sha_drift(self):
        with self.assertRaisesRegex(ValueError,"SHA256 differs"):
            q.run(self.source,self.output,expected_source_sha256="0"*64,
                  expected_photo_count=15,tier=3,part=6,chunk_size=4)
        self.assertFalse(self.output.exists())
        with self.assertRaisesRegex(ValueError,"tier invalid"):
            q.run(self.source,self.output,expected_source_sha256=q.sha(self.source),
                  expected_photo_count=15,tier=4,part=6,chunk_size=4)
        self.assertFalse(self.output.exists())

if __name__=="__main__":
    unittest.main()
