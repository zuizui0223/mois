"""Bounded full93m source part chunks are NOT false completed parent source parts."""
from pathlib import Path
from tempfile import TemporaryDirectory
import csv
import gzip
import io
import json
import unittest

from PIL import Image
from scripts.acquisition import mois_full93m_source_part_chunks as chunks
from scripts.acquisition import mois_full93m_signed_chunk_pixel_scan as pixel
from scripts.acquisition.mois_onepass_full93m_target_tiers import COLUMNS
from scripts.acquisition.mois_full93m_medium_pixel_scan import sha_file

class SignedOriginalPhotoChunkPixelTests(unittest.TestCase):
    def setUp(self):
        self.tmp=TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        source=self.root/"source.csv.gz"
        with gzip.open(source,"wt",newline="",encoding="utf-8") as fh:
            writer=csv.DictWriter(fh,fieldnames=COLUMNS,lineterminator="\n")
            writer.writeheader()
            for i in range(6):
                row={x:"" for x in COLUMNS}
                row.update(photo_id=str(105+i),observation_uuid="o-"+str(i),
                           inat_taxon_id="11",image_extension="jpg",
                           observed_on="2025-06-22",photo_license="cc-by",
                           latitude="45",tier="12")
                writer.writerow(row)
        self.chunks=self.root/"source_chunks"
        chunks.run(source,self.chunks,expected_source_sha256=sha_file(source),
                   expected_photo_count=6,tier=12,part=4,chunk_size=4)
        self.output=self.root/"original_medium"
        im=Image.new("RGB",(80,60),(80,100,30))
        f=io.BytesIO(); im.save(f,format="JPEG")
        self.bytes=f.getvalue()

    def test_resume_real_pixel_chunk_without_second_download_and_not_complete_parent(self):
        hits=[]
        def retrieve(url):
            hits.append(url)
            return self.bytes
        first=pixel.run(self.chunks,self.output,0,
                        detector=None,download=retrieve,
                        workers=1,batch_size=2,max_packets=1)
        self.assertEqual(first["actual_new_or_resumed_medium_photo_attempts"],2)
        self.assertFalse(first["source_chunk_all_photo_ids_scanned"])
        self.assertFalse((self.output/"complete_receipt.json").exists())
        self.assertEqual(len(hits),2)
        final=pixel.run(self.chunks,self.output,0,
                        detector=None,download=retrieve,
                        workers=1,batch_size=2)
        self.assertEqual(final["source_chunk_expected_photo_ids"],4)
        self.assertEqual(final["actual_new_or_resumed_medium_photo_attempts"],4)
        self.assertTrue(final["source_chunk_all_photo_ids_scanned"])
        self.assertFalse(final["entire_parent_original_photo_part_scanned"])
        self.assertFalse(final["source_93m_all_images_scanned"])
        self.assertEqual(len(hits),4)
        completed=json.loads((self.output/"complete_receipt.json").read_text())
        self.assertEqual(completed["source_role"],"full93m_chunk")
        self.assertEqual(completed["source_chunk_index"],0)
        self.assertFalse(completed["entire_original_source_part_completed"])
        self.assertEqual(completed["medium_image_pixels_decoded"],4)
        self.assertEqual(completed["screen_counts"],{"uncertain":4})

    def test_invalid_child_bytes_fail_without_image_download(self):
        path=self.chunks/"chunks/chunk-00000.csv.gz"
        path.write_bytes(b"tampered")
        with self.assertRaisesRegex(ValueError,"SHA verified"):
            pixel.run(self.chunks,self.output,0,detector=None,
                      download=lambda url: self.fail("no fetch"))
        self.assertFalse(self.output.exists())

    def test_parent_manifest_sha_drift_fails_before_image_download(self):
        p=self.chunks/"chunk_manifest.csv"
        p.write_text(p.read_text()+"\n")
        with self.assertRaisesRegex(ValueError,"manifest SHA"):
            pixel.run(self.chunks,self.output,0,detector=None,
                      download=lambda url:self.fail("no fetch"))
        self.assertFalse(self.output.exists())

    def test_second_child_can_complete_without_completing_parent(self):
        result=pixel.run(self.chunks,self.output,1,
                         detector=None,download=lambda _: self.bytes,
                         workers=1,batch_size=2)
        self.assertEqual(result["source_chunk_index"],1)
        self.assertEqual(result["actual_new_or_resumed_medium_photo_attempts"],2)
        self.assertFalse(result["entire_parent_original_photo_part_scanned"])

if __name__=="__main__":
    unittest.main()
