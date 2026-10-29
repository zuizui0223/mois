"""Independent blind botanical human review files contain only blank labels."""
import csv
import json
import tempfile
import unittest
from pathlib import Path

from scripts.analysis.mois_observation_probability_audit_20261029 import (
    SCHEMA as SOURCE_SCHEMA, sha_file
)
from scripts.analysis.mois_1200_blind_human_review_ui_20261029 import (
    create_review_bundle, REVIEW_SCHEMA
)

class BlindReviewUIContractTests(unittest.TestCase):
    def test_1200_primary_and_independent_240_overlap(self):
        with tempfile.TemporaryDirectory() as tmp:
            base=Path(tmp)
            roster=base/"audit1200_source_observation_roster.csv"
            fields=("photo_id","observation_id","inat_taxon_id","source_design",
                    "latitude_zone","coarse_zone","local_season","roi_status",
                    "stratum","stratum_N","stratum_n","audit_inclusion_q",
                    "review_inverse_q")
            with roster.open("w",newline="") as f:
                writer=csv.DictWriter(f,fieldnames=fields)
                writer.writeheader()
                for i in range(1200):
                    status="success" if i<400 else "unknown"
                    stratum="arm|zone|season|"+status+"|"+str(i//20)
                    writer.writerow({
                        "photo_id":str(10000000+i),"observation_id":"obs-"+str(i),
                        "inat_taxon_id":str(i),"source_design":"arm",
                        "latitude_zone":"N_temperate","coarse_zone":"N_extratropical",
                        "local_season":"summer","roi_status":status,
                        "stratum":stratum,"stratum_N":"200",
                        "stratum_n":"20","audit_inclusion_q":"0.1",
                        "review_inverse_q":"10"
                    })
            source={"schema":SOURCE_SCHEMA,"sample_total":1200,
                    "sample_success_count":400,"sample_unknown_count":800,
                    "output_sha256":{"audit1200_source_observation_roster.csv":sha_file(roster)}}
            (base/"audit1200_frozen_design_receipt.json").write_text(json.dumps(source))
            result=create_review_bundle(base)
            self.assertEqual(result["schema"],REVIEW_SCHEMA)
            self.assertEqual(result["n_second_reviewer_photo_IDs"],240)
            self.assertEqual(result["double_review_success"],80)
            self.assertEqual(result["double_review_unknown"],160)
            with (base/"audit1200_secondary240_review_BLANK.csv").open() as f:
                second=list(csv.DictReader(f))
            self.assertEqual(len({r["photo_id"] for r in second}),240)
            self.assertTrue(all(r["human_open_corolla_visible"]=="" for r in second))
            for suffix in ("audit1200_primary_blind_review.html",
                           "audit1200_secondary240_blind_review.html"):
                html=(base/suffix).read_text()
                self.assertNotIn("__ROWS__",html)
                self.assertNotIn("__ROLE__",html)
                self.assertNotIn('"roi_status":',html)
                self.assertNotIn('"latitude_zone":',html)
                self.assertNotIn('"source_design":',html)
                self.assertIn("inaturalist-open-data.s3.amazonaws.com/photos/",html)
            # Any untracked roster edit invalidates the precommitted audit.
            with roster.open("a") as f:f.write("tamper")
            with self.assertRaisesRegex(ValueError,"SHA"):
                create_review_bundle(base)

if __name__=="__main__":
    unittest.main()
