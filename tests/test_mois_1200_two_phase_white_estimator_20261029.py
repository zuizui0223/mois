"""Human review blanks never become white, uncertainty yields partial-ID bounds."""
from __future__ import annotations
import csv
import json
import tempfile
import unittest
from pathlib import Path

from scripts.analysis.mois_observation_probability_audit_20261029 import (
    SCHEMA as AUDIT_SCHEMA, sha_file
)
from scripts.analysis.mois_1200_two_phase_white_estimator_20261029 import (
    compute
)

class TwoPhaseHumanEstimatorTests(unittest.TestCase):
    def prepare(self, folder):
        root=Path(folder)
        roster=root/"audit1200_source_observation_roster.csv"
        annotations=root/"human.csv"
        rs=["photo_id","observation_id","stratum","stratum_N","stratum_n",
            "audit_inclusion_q","review_inverse_q","roi_status"]
        anns=["photo_id","observation_id","review_role","reviewer_id",
              "human_review_completed","human_open_corolla_visible",
              "human_corolla_colour"]
        all_r=[];all_a=[]
        for i in range(1200):
            success=i<400
            group=(i//25) if success else (i-400)//50
            stratum=("S|" if success else "U|")+str(group)
            n=25 if success else 50
            N=100
            p=str(123000000+i)
            o="obs-"+str(i)
            all_r.append({
                "photo_id":p,"observation_id":o,"stratum":stratum,
                "stratum_N":N,"stratum_n":n,"audit_inclusion_q":n/N,
                "review_inverse_q":N/n,"roi_status":"success" if success else "unknown"
            })
            all_a.append({
                "photo_id":p,"observation_id":o,"review_role":"R1",
                "reviewer_id":"observer-X","human_review_completed":"False",
                "human_open_corolla_visible":"","human_corolla_colour":""
            })
        with roster.open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=rs);w.writeheader();w.writerows(all_r)
        with annotations.open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=anns);w.writeheader();w.writerows(all_a)
        source=root/"receipt.json"
        source.write_text(json.dumps({
            "schema":AUDIT_SCHEMA,"sample_total":1200,
            "output_sha256":{"audit1200_source_observation_roster.csv":sha_file(roster)}
        }))
        return root,roster,source,annotations,all_a,anns

    def rewrite(self,path,rows,fields):
        with path.open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)

    def test_blanks_block_prevalence_until_all_reviewed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root,r,s,a,rows,fields=self.prepare(tmp)
            out=root/"blocked.json"
            result=compute(r,s,a,out)
            self.assertEqual(result["status"],"BLOCKED_PENDING_HUMAN_REVIEW")
            self.assertEqual(result["n_review_pending"],1200)
            self.assertIsNone(result["design_corrected_white_fraction"])

    def test_full_known_labels_produce_design_based_ratio_and_se(self):
        with tempfile.TemporaryDirectory() as tmp:
            root,r,s,a,rows,fields=self.prepare(tmp)
            for i,row in enumerate(rows):
                row["human_review_completed"]="True"
                row["human_open_corolla_visible"]="yes"
                row["human_corolla_colour"]="white" if i%2==0 else "red_pink"
            self.rewrite(a,rows,fields)
            result=compute(r,s,a,root/"complete.json")
            self.assertEqual(result["status"],"DESIGN_BASED_FULLY_OBSERVED_HUMAN_E1")
            self.assertAlmostEqual(result["design_corrected_white_fraction"],0.5,places=8)
            self.assertGreater(result["taylor_linearization_standard_error"],0)
            self.assertEqual(len(result["approximate_normal_95_interval"]),2)
            self.assertFalse(result["natural_flora_white_frequency_identified"])

    def test_manual_uncertainty_yields_bounds_and_never_a_fabricated_point(self):
        with tempfile.TemporaryDirectory() as tmp:
            root,r,s,a,rows,fields=self.prepare(tmp)
            for i,row in enumerate(rows):
                row["human_review_completed"]="True"
                row["human_open_corolla_visible"]="yes"
                row["human_corolla_colour"]="white" if i%2==0 else "blue_purple"
            rows[5]["human_open_corolla_visible"]="inaccessible"
            rows[5]["human_corolla_colour"]=""
            rows[7]["human_corolla_colour"]="mixed"
            self.rewrite(a,rows,fields)
            result=compute(r,s,a,root/"bounds.json")
            self.assertEqual(result["status"],"PARTIALLY_IDENTIFIED_AFTER_COMPLETE_REVIEW")
            self.assertIsNone(result["design_corrected_white_fraction"])
            bounds=result["white_frequency_partial_bounds"]
            self.assertLess(bounds["lower"],bounds["upper"])
            self.assertGreater(bounds["lower"],0)
            self.assertLess(bounds["upper"],1)

    def test_mismatched_identity_and_incoherent_colour_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root,r,s,a,rows,fields=self.prepare(tmp)
            rows[0]["observation_id"]="not-the-frozen-observation"
            self.rewrite(a,rows,fields)
            with self.assertRaisesRegex(ValueError,"exactly 1200|mismatched|join|Annotation"):
                compute(r,s,a,root/"bad.json")
            self.assertFalse((root/"bad.json").exists())

if __name__=="__main__":
    unittest.main()
