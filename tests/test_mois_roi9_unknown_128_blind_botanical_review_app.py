"""Botanical review must be independent image-case selection, never automatic label."""
import csv,json,tempfile,unittest
from pathlib import Path
from collections import Counter

from scripts.analysis.mois_roi9_unknown_128_blind_botanical_review_app import (
    choose,create
)

class Human128PhotoReviewTests(unittest.TestCase):
    def data(self):
        result=[]
        for i in range(512):
            status="negative" if i<14 else ("error" if i<17 else "candidate")
            diag=("detector_negative_requires_botanical_review" if status=="negative" else
                 "fetch_or_decode_error_without_botanical_inference" if status=="error" else
                 "detectable_large_flower_box_ROI9_failed" if i<400 else
                 "existing_cached_flower_box_ROI9_failed")
            result.append(dict(
                photo_id=str(9000000+i),observation_id="obs-"+str(i),
                inat_taxon_id=str(20000+i),old_ROI9_missing="True",
                source_design=["historical_K1_FCP_source_heterogeneous",
                               "p2_first_4096","bias_refill_2048",
                               "bulk128_131072"][i%4],
                latitude_zone=["N_temperate","S_temperate","N_tropical",
                               "S_tropical"][i%4],
                local_season=["winter","spring","summer","autumn"][i%4],
                postoutcome_source_gap=(
                    "P2_source_has_photos_no_colour_10plus" if i%2 else
                    "P1_source_rare_no_colour_1to9"),
                screen_status=status,technical_diagnostic_class=diag
            ))
        return result

    def test_detector_negatives_kept_but_not_claimed_biological_absence(self):
        a=choose(self.data())
        self.assertEqual(len(a),128)
        self.assertEqual(len({r["photo_id"] for r in a}),128)
        self.assertEqual(len({r["inat_taxon_id"] for r in a}),128)
        self.assertEqual(sum(r["screen_status"]!="candidate" for r in a),17)
        self.assertGreaterEqual(sum(r["technical_diagnostic_class"]==
                             "detectable_large_flower_box_ROI9_failed" for r in a),64)

    def test_real_annotation_UI_labels_begin_empty(self):
        with tempfile.TemporaryDirectory() as tmp:
            source=Path(tmp)/"all512.csv"
            data=self.data()
            with source.open("w",newline="") as f:
                writer=csv.DictWriter(f,fieldnames=list(data[0]))
                writer.writeheader();writer.writerows(data)
            output=Path(tmp)/"review"
            result=create(source,output,validate_sha=False)
            self.assertEqual(result["botanical_review_original_photo_IDs_selected"],128)
            self.assertTrue(result["all_botanical_open_flower_or_petal_labels_still_blank"])
            html=(output/"mois_128_botanical_flower_ROI9_failure_annotation.html").read_text()
            self.assertIn("localStorage",html)
            self.assertIn("human_flower_visible",html)
            self.assertNotIn("__SOURCE_ROWS__",html)
            self.assertNotIn("__SOURCE_SHA__",html)
            with (output/"manual128_botanical_annotation_BLANK_not_yet_reviewed.csv").open(newline="") as f:
                rows=list(csv.DictReader(f))
            self.assertEqual(len(rows),128)
            self.assertTrue(all(not r["human_flower_visible"] and
                                not r["human_failure_reason"] and
                                r["human_review_completed"]=="False" for r in rows))
            self.assertEqual(len(set(r["inat_taxon_id"] for r in rows)),128)

    def test_duplicate_photo_or_source_taxon_rejected(self):
        cases=self.data()
        cases[1]["photo_id"]=cases[0]["photo_id"]
        with self.assertRaisesRegex(ValueError,"unique original source identities"):
            choose(cases)
        cases=self.data()
        cases[2]["inat_taxon_id"]=cases[0]["inat_taxon_id"]
        with self.assertRaisesRegex(ValueError,"unique original source identities"):
            choose(cases)

if __name__=="__main__":
    unittest.main()
