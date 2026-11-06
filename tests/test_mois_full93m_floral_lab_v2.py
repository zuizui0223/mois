"""CI for a genuinely new sRGB/Lab Stage1 mode; old v1 packets immutable."""
import csv
import gzip
import io
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from PIL import Image

from scripts.acquisition import mois_full93m_medium_pixel_scan as scanner
from scripts.acquisition.mois_onepass_full93m_target_tiers import COLUMNS
from scripts.acquisition.mois_flower_box_lab_colour import PROTOCOL

class SignedFloralLabV2Tests(unittest.TestCase):
    def setUp(self):
        self.temp=TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.source=self.root/"part.csv.gz"
        with gzip.open(self.source,"wt",newline="",encoding="utf-8") as f:
            w=csv.DictWriter(f,fieldnames=COLUMNS,lineterminator="\n")
            w.writeheader()
            for pid in range(4):
                row={c:"" for c in COLUMNS}
                row.update(photo_id=str(6000+pid),observation_uuid="obs"+str(pid),
                           inat_taxon_id="17",latitude="20",observed_on="2025-07-04",
                           image_extension="png",tier="0")
                w.writerow(row)
        self.bytes={}
        for i,rgb in enumerate(((20,160,35),(115,67,27),
                                (255,255,255),(255,255,255))):
            m=io.BytesIO()
            Image.new("RGB",(80,80),rgb).save(m,format="PNG")
            self.bytes[str(6000+i)]=m.getvalue()
        self.output=self.root/"v2"

    class Detector:
        def __call__(self,images):
            return ([{"score":.92,"xyxy":(0,0,80,80)} for _ in images[:-1]]
                    +[None])

    def download(self,url):
        photo_id=url.split("/")[-2]
        return self.bytes[photo_id]

    def scan(self,output,mode):
        return scanner.scan_source(
            self.source,output,
            expected_source_sha256=scanner.sha_file(self.source),
            expected_count=4,part=0,tier=0,
            detector=self.Detector(),download=self.download,
            batch_size=4,workers=2,feature_mode=mode)

    def test_floral_box_lab_preserves_green_brown_white_and_nonflower_unknown(self):
        proof=self.scan(self.output,"floral_box_lab_v2")
        self.assertEqual(proof["schema"],scanner.LAB_SCHEMA)
        self.assertEqual(proof["CIELAB_box_source_protocol"],PROTOCOL)
        self.assertTrue(proof["floral_region_is_not_petals_semantically_separated"])
        self.assertEqual(proof["screen_counts"],{"candidate":3,"negative":1})
        self.assertEqual(proof["actual_unique_source_photo_ids_attempted"],4)
        self.assertEqual(proof["medium_image_pixels_decoded"],4)
        with gzip.open(self.output/"packet-0000000.csv.gz","rt") as f:
            reader=csv.DictReader(f)
            self.assertEqual(tuple(reader.fieldnames),tuple(scanner.COLUMNS_LAB))
            r=list(reader)
        self.assertEqual([x["screen_status"] for x in r],
                         ["candidate","candidate","candidate","negative"])
        self.assertGreater(float(r[0]["floral_box_fraction_green_like"]),.95)
        self.assertLess(float(r[0]["floral_box_fraction_white_like"]),.01)
        self.assertGreater(float(r[1]["floral_box_fraction_brown_like"]),.95)
        self.assertLess(float(r[1]["floral_box_main_L_star"]),55)
        self.assertEqual(r[2]["floral_box_main_rgb_hex"],"#FFFFFF")
        self.assertEqual(r[2]["floral_box_main_hue_degrees"],"")
        self.assertEqual(r[3]["floral_box_main_L_star"],"")
        self.assertEqual(r[3]["technical_reason"],"model_screen_no_flower_candidate")
        with gzip.open(self.output/"packet-0000000.csv.gz","rb") as f:
            self.assertGreater(len(f.read()),0)
        # Exact replay does not download any original photo pixels.
        old=self.download
        self.download=lambda u: self.fail("Versioned source pixels were fetched again")
        self.assertEqual(self.scan(self.output,"floral_box_lab_v2"),proof)
        self.download=old

    def test_legacy_v1_never_mistaken_for_photo_lab_v2(self):
        out=self.root/"original_v1"
        legacy=self.scan(out,"legacy_hsv")
        self.assertEqual(legacy["schema"],scanner.SCHEMA)
        with gzip.open(out/"packet-0000000.csv.gz","rt") as f:
            self.assertEqual(tuple(csv.DictReader(f).fieldnames),tuple(scanner.COLUMNS_OUT))
        with self.assertRaisesRegex(ValueError,"Cannot resume"):
            self.scan(out,"floral_box_lab_v2")
        self.assertTrue((out/"complete_receipt.json").exists())

    def test_bad_feature_mode_rejected_before_fetch(self):
        with self.assertRaisesRegex(ValueError,"supported immutable"):
            self.scan(self.output,"invent_organ_lab")
        self.assertFalse(self.output.exists())

if __name__=="__main__":
    unittest.main()
