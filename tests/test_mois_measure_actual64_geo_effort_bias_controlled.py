"""Source ID, inference and missing source inclusion probability QA of actual-pilot code."""
import csv,gzip,io,json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
import numpy as np
import pandas as pd
from PIL import Image

from scripts.acquisition.mois_full93m_medium_pixel_scan import (
    scan_source,sha_file,LAB_SCHEMA
)
from scripts.acquisition.mois_onepass_full93m_target_tiers import COLUMNS
from scripts.acquisition.mois_measure_actual64_geo_effort_bias_controlled import (
    describe_scan
)

class BiasAwarePhotoPixelPilotTests(unittest.TestCase):
    class Detector:
        def __call__(self,images):
            return [{"score":.87,"xyxy":[0,0,80,80]} for x in images]
    def setUp(self):
        t=TemporaryDirectory()
        self.addCleanup(t.cleanup)
        self.root=Path(t.name)
        self.gzip=self.root/"true_source.csv.gz"
        with gzip.open(self.gzip,"wt",encoding="utf-8",newline="") as f:
            w=csv.DictWriter(f,fieldnames=COLUMNS,lineterminator="\n")
            w.writeheader()
            for i in range(4):
                v={key:"" for key in COLUMNS}
                v.update({
                    "photo_id":str(7000+i),
                    "observation_uuid":"independent"+str(i),
                    "inat_taxon_id":str(8000+i),
                    "tier":"7",
                    "observed_on":"2026-05-06",
                    "latitude": ["70","-37","17","-19"][i],
                    "longitude":"120",
                    "image_extension":"png","photo_license":"cc-by",
                    "frozen_loader_supported":"true",
                })
                w.writerow(v)
        self.payload={}
        for i,rgb in enumerate([(255,255,255),(20,160,35),
                                (115,67,27),(215,120,190)]):
            f=io.BytesIO()
            Image.new("RGB",(80,80),rgb).save(f,format="PNG")
            self.payload[str(7000+i)]=f.getvalue()
    def test_all_pixels_including_green_and_brown_one_original_obs_each(self):
        out=self.root/"measure"
        proof=scan_source(self.gzip,out,expected_source_sha256=sha_file(self.gzip),
                 expected_count=4,part=0,tier=7,
                 detector=self.Detector(),
                 download=lambda url:self.payload[url.split("/")[-2]],
                 workers=2,batch_size=4,max_packets=1,
                 source_role="bias_controlled_pilot",feature_mode="floral_box_lab_v2")
        self.assertEqual(proof["schema"],LAB_SCHEMA)
        self.assertEqual(proof["source_role"],"bias_controlled_pilot")
        self.assertTrue(proof["outcome_conditioned_supplement_not_population_unbiased"])
        self.assertFalse(proof["all_original_input_photos_attempted"])
        self.assertEqual(proof["medium_image_pixels_decoded"],4)
        with gzip.open(out/"packet-0000000.csv.gz","rt",encoding="utf-8") as f:
            d=list(csv.DictReader(f))
        self.assertEqual(len(d),4)
        self.assertTrue(all(not r["source_photo_inclusion_probability"] for r in d))
        self.assertTrue(all(r["source_arm"].startswith("bias_controlled_") for r in d))
        self.assertEqual([r["observation_id"] for r in d],
             ["independent"+str(i) for i in range(4)])
        self.assertGreater(float(d[1]["floral_box_fraction_green_like"]),.9)
        self.assertGreater(float(d[2]["floral_box_fraction_brown_like"]),.9)
        d=pd.DataFrame(d)
        d["source_hemisphere_band"]=["north_extratropical","south_extratropical",
                                     "north_tropical","south_tropical"]
        d["sampling_prior_attempt_class"]=["one_attempt"]*4
        summary=describe_scan(d)
        self.assertEqual(summary[0]["original_source_photo_IDs_selected"],4)
        self.assertEqual(summary[0]["photo_RGB_successfully_downloaded"],4)
        self.assertEqual(summary[0]["flower_boxes_with_actual_Lab_proxy"],4)
        self.assertEqual(len(summary),1+4+4)
        self.assertEqual(sum(x["original_source_photo_IDs_selected"] for x in summary
                        if x["design_stratum"].startswith("latitude:")),4)

    def test_resume_refuses_to_change_first_frozen_pilot_photo_source(self):
        out=self.root/"once"
        method=dict(expected_source_sha256=sha_file(self.gzip),
             expected_count=4,part=0,tier=7,detector=self.Detector(),
             download=lambda url:self.payload[url.split("/")[-2]],
             workers=2,batch_size=4,max_packets=1,
             source_role="bias_controlled_pilot",feature_mode="floral_box_lab_v2")
        scan_source(self.gzip,out,**method)
        # Immutable packet SHA and model identity remain stable; no network
        # requests are made when the exact same source is replayed.
        method["download"]=lambda url:self.fail("Already signed source was re-fetched")
        scan_source(self.gzip,out,**method)

if __name__=="__main__":
    unittest.main()
