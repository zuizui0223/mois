"""Exhaustive 93m workload estimate uses *actual signed photo counts*, no RNG."""
import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from scripts.acquisition import mois_full93m_4096_chunk_warm_job_plan as q

class Signed93mWarmJobs(unittest.TestCase):
    def setUp(self):
        self.temp=TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.csv=self.root/"original1792.csv"
        self.out=self.root/"out"
        self.create([7,0,9,5],parts=2,tiers=2)

    def create(self,counts,*,parts,tiers):
        with self.csv.open("w",newline="",encoding="utf-8") as f:
            w=csv.DictWriter(f,fieldnames=q.COLUMNS,lineterminator="\n")
            w.writeheader()
            for t in range(tiers):
                for p in range(parts):
                    w.writerow({"tier":t,"part":p,
                     "relative_source_gzip_path":f"tier-{t:02d}/part-{p:03d}.csv.gz",
                     "source_gzip_sha256":"d"*64,
                     "expected_photo_rows":counts[t*parts+p],
                     "photo_level_sampling_probability":"1"})

    def compute(self,**kwargs):
        return q.run(self.csv,self.out,expected_sha=q.sha(self.csv),
                     photo_count=21,parts=2,tiers=2,chunk_size=4,
                     chunks_per_warm=2,**kwargs)

    def test_complete_source_photo_count_exactly_once(self):
        r=self.compute()
        self.assertEqual(r["original_source_photo_IDs_entire_universe"],21)
        self.assertEqual(r["planned_child_photo_chunk_tasks"],7)
        self.assertEqual(r["planned_one_warm_model_worker_batches"],4)
        self.assertEqual(r["images_downloaded"],0)
        self.assertEqual(r["actual_photo_scans_completed_by_this_plan"],0)
        with (self.out/"all_photo_source_4096_chunk_tasks.csv").open() as f:
            child=list(csv.DictReader(f))
        with (self.out/"all_one_warm_model_8chunk_batches.csv").open() as f:
            warm=list(csv.DictReader(f))
        self.assertEqual(len(child),7)
        self.assertEqual(len(warm),4)
        self.assertEqual(sum(int(x["planned_original_photo_ids"]) for x in child),21)
        self.assertEqual(sum(int(x["expected_original_photo_ids"]) for x in warm),21)
        self.assertTrue(all(x["child_chunk_sha256_available_only_after_real_source_split"]=="true" for x in child))
        self.assertEqual([int(x["planned_original_photo_ids"]) for x in child],[4,3,4,4,1,4,1])
        self.assertTrue(all(int(x["expected_original_photo_ids"])<=8 for x in warm))
        receipt=json.loads((self.out/"all93m_chunk_workplan_receipt.json").read_text())
        self.assertEqual(receipt["chunk_task_csv_sha256"],q.sha(self.out/"all_photo_source_4096_chunk_tasks.csv"))

    def test_missing_photo_universe_refused(self):
        with self.assertRaisesRegex(ValueError,"Whole original"):
            q.run(self.csv,self.out,expected_sha=q.sha(self.csv),
                  photo_count=22,parts=2,tiers=2,chunk_size=4,
                  chunks_per_warm=2)
        self.assertFalse(self.out.exists())

    def test_missing_partition_or_zero_probability_refused(self):
        self.create([7,0,9],parts=1,tiers=3)
        with self.assertRaisesRegex(ValueError,"Source part"):
            q.run(self.csv,self.out,expected_sha=q.sha(self.csv),
                  photo_count=16,parts=2,tiers=2,chunk_size=4,chunks_per_warm=2)
        self.assertFalse(self.out.exists())
        self.create([7,0,9,5],parts=2,tiers=2)
        source=self.csv.read_text().replace(",1\n",",0.5\n",1)
        self.csv.write_text(source)
        with self.assertRaisesRegex(ValueError,"Do not silently drop"):
            self.compute()

if __name__=="__main__":
    unittest.main()
