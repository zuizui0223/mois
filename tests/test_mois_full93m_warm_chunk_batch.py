"""Multiple complete child-source pixel chunks share one warm detection object."""
from pathlib import Path
from tempfile import TemporaryDirectory
import csv
import gzip
import io
import json
import unittest

from PIL import Image
from scripts.acquisition import mois_full93m_source_part_chunks as chunk
from scripts.acquisition import mois_full93m_warm_chunk_batch as warm
from scripts.acquisition.mois_onepass_full93m_target_tiers import COLUMNS

class WarmSignedChunkBatchTests(unittest.TestCase):
    def setUp(self):
        self.temp=TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        source=self.root/"source.csv.gz"
        with gzip.open(source,"wt",newline="",encoding="utf-8") as f:
            writer=csv.DictWriter(f,fieldnames=COLUMNS,lineterminator="\n")
            writer.writeheader()
            for n in range(7):
                row={k:"" for k in COLUMNS}
                row.update(photo_id=str(900+n),observation_uuid="O"+str(n),
                           inat_taxon_id="73",tier="12",image_extension="jpg",
                           photo_license="cc0",observed_on="2022-07-13")
                writer.writerow(row)
        self.signed=self.root/"signed"
        chunk.run(source,self.signed,expected_source_sha256=chunk.sha(source),
                  expected_photo_count=7,tier=12,part=3,chunk_size=4)
        self.output=self.root/"stage1"
        x=Image.new("RGB",(40,40),(200,50,120));b=io.BytesIO();x.save(b,format="JPEG")
        self.photo=b.getvalue()
        self.downloads=[]
        class Detector:
            def __init__(self):
                self.calls=0
            def __call__(self, imgs):
                self.calls+=1
                return [None]*len(imgs)
        self.detector=Detector()

    def download(self,url):
        self.downloads.append(url)
        return self.photo

    def test_all_source_child_images_in_order_same_model_without_duplicate_fetch(self):
        x=warm.run(self.signed,self.output,first_chunk=0,chunk_count=2,
                   detector=self.detector,download=self.download,
                   max_real_photo_ids=7,workers=1,batch_size=2)
        self.assertEqual(x["photo_IDs_fully_pixel_attempted_in_this_batch"],7)
        self.assertEqual(x["source_chunks_in_batch"],2)
        self.assertEqual(self.detector.calls,4)
        self.assertEqual(len(self.downloads),7)
        self.assertFalse(x["original_photo_part_complete"])
        self.assertFalse(x["source_93m_all_pixels_complete"])
        self.assertEqual(len(x["completed_chunk_receipts"]),2)
        for child in x["completed_chunk_receipts"]:
            self.assertEqual(child["pixel_status_counts"],{"negative":child["original_photo_ids"]})
            self.assertTrue((self.output/f"chunk-{child['chunk_index']:05d}"/"complete_receipt.json").exists())
        # Complete rerun reads existing packet SHA and source identity; it NEVER
        # reaches either photo downloader or detector.
        old=self.detector.calls
        def forbid(_):
            self.fail("Previously completed source original photo redownloaded")
        again=warm.run(self.signed,self.output,first_chunk=0,chunk_count=2,
                       detector=self.detector,download=forbid,
                       max_real_photo_ids=7,workers=1,batch_size=2)
        self.assertEqual(x,again)
        self.assertEqual(self.detector.calls,old)

    def test_never_exceed_bounded_photo_cap(self):
        with self.assertRaisesRegex(ValueError,"explicitly bounded"):
            warm.run(self.signed,self.output,first_chunk=0,chunk_count=2,
                     detector=self.detector,download=self.download,
                     max_real_photo_ids=6)
        self.assertFalse(self.output.exists())
        self.assertFalse(self.downloads)

    def test_parent_never_complete_by_one_chunk(self):
        summary=warm.run(self.signed,self.output,first_chunk=1,chunk_count=1,
                         detector=self.detector,download=self.download,
                         max_real_photo_ids=4,workers=1,batch_size=2)
        self.assertEqual(summary["photo_IDs_fully_pixel_attempted_in_this_batch"],3)
        self.assertEqual(summary["source_chunks_in_batch"],1)
        self.assertFalse(summary["original_photo_part_complete"])

if __name__=="__main__":
    unittest.main()
