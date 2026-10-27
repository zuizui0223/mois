"""Source-signed old ROI unknown pixels get box diagnostic, never pigment rescue."""
import csv,hashlib,io,json,tempfile,unittest
from pathlib import Path

import numpy as np
from PIL import Image

from scripts.acquisition.mois_existing_roi_unknown_342_medium_flower_diagnostic import (
    run,sha,photo_medium,classify_box_predictions,RequestPacer,
    MODEL_SHA,DIAGNOSTIC_CLASSES
)

class Boxes:
    def __init__(self,xyxy,conf):
        self._xy=np.asarray(xyxy,dtype=float)
        self._conf=np.asarray(conf,dtype=float)
    @property
    def xyxy(self):
        class Tensor:
            def __init__(self,values):self.value=values
            def cpu(self):return self
            def numpy(self):return self.value
        return Tensor(self._xy)
    @property
    def conf(self):
        class Tensor:
            def __init__(self,values):self.value=values
            def cpu(self):return self
            def numpy(self):return self.value
        return Tensor(self._conf)
    def __len__(self):return len(self._conf)

class Detector:
    def predict(self,source,**opts):
        assert opts["conf"]==0.02
        return [type("Prediction",(),{
            "boxes":Boxes([[5,5,40,40]],[.9])})()
            for _ in source]

class RealOldPhotoBoxDiagnosisTests(unittest.TestCase):
    def prepared(self,p):
        root=Path(p)
        (root/"workers").mkdir(parents=True)
        rows=[]
        for i in range(43):
            rows.append({
                "photo_id":str(50000000+i),
                "observation_id":f"already-observed-{i}",
                "inat_taxon_id":str(1500+i),
                "source_design":"historical_K1_FCP_source_heterogeneous",
                "latitude_zone":"N_temperate",
                "local_season":"summer",
                "postoutcome_source_gap":"P2_source_has_photos_no_colour_10plus",
                "previous_original_photo_attempts":"5",
                "source_photo_opportunities":"250",
                "reason":"FOUR_HUE_UNKNOWN_AND_NO_STRICT9",
            })
        dest=root/"workers"/"shard_00.csv"
        with dest.open("w",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]))
            writer.writeheader();writer.writerows(rows)
        m={"schema":"mois_512_prior_strictROI9_unknown_existing_medium500_reuse_and_bounded_screen_v1",
           "source_review_photo_rows":512,
           "existing_500px_photo_screen_results_reused":170,
           "existing_unknown_photos_needing_small_bounded_500px_recheck":342,
           "new_photo_IDs_selected_for_world_census":0,
           "eight_worker_photo_counts":[43,43,43,43,43,43,42,42],
           "worker_source_SHA256":{"0":sha(dest)}}
        (root/"manifest.json").write_text(json.dumps(m))
        image=Image.new("RGB",(100,100),(245,120,40))
        buf=io.BytesIO();image.save(buf,format="JPEG")
        return buf.getvalue()

    def test_43_previously_measured_images_still_9palette_unknown(self):
        with tempfile.TemporaryDirectory() as t:
            source=Path(t)/"source";raw=self.prepared(source)
            out=Path(t)/"screens"
            record=run(source,Path(t)/"placeholder-model.pt",0,out,
                detector=Detector(),verify_model=False,requests_per_second=1,
                fetch=lambda url,pace:raw)
            self.assertEqual(record["original_old_photo_ids_attempted_for_first_medium_diagnostic"],43)
            self.assertEqual(record["decoded_medium500"],43)
            self.assertEqual(record["statuses"],{"candidate":43})
            self.assertEqual(record["technical_box_classes"],
                             {"detectable_large_flower_box_ROI9_failed":43})
            self.assertEqual(record["old_source_ROI9_unknown_statuses_preserved"],43)
            with (out/"shard_00_original_medium_diagnostic.csv").open(newline="") as f:
                rows=list(csv.DictReader(f))
            self.assertTrue(all(x["strict_original_roi9_remains_unknown"]=="True" for x in rows))
            self.assertTrue(all(float(x["largest_box_area_fraction"])>.06 for x in rows))
            self.assertTrue(all(x["screen_status"]=="candidate" for x in rows))
            self.assertFalse(any("flower_colour" in str(x) for x in rows))

    def test_detector_negative_does_not_become_no_flower_truth(self):
        class Empty:
            boxes=Boxes([],[])
        d=classify_box_predictions(Image.new("RGB",(400,400)),Empty())
        self.assertEqual(d["screen_status"],"negative")
        self.assertIn("requires_botanical_review",d["technical_diagnostic_class"])
        self.assertNotIn("white",d["technical_diagnostic_class"])
        self.assertIn(d["technical_diagnostic_class"],DIAGNOSTIC_CLASSES)

if __name__=="__main__":
    unittest.main()
