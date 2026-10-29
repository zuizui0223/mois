"""Automatic-only ROI9 bias bounds: no human truth, no unknowable colour imputation."""
from __future__ import annotations
import json
import math
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from scripts.analysis.mois_observation_roi9_auto_bias_20261029 import (
    SCHEMA, logit_tilt, group_summary, source_groups, run,
)
from scripts.analysis.mois_roi9_observation_frame_20261029 import (
    sha_file,
)

class AutomaticBiasTests(unittest.TestCase):
    def fake_frame(self):
        arm="historical_K1_FCP_source_heterogeneous"
        rec=[]
        # Two distinct source strata: 10 each, 2 vs 8 valid colours.
        for zone,k,white_n in (("N_tropical",2,2),("S_tropical",8,0)):
            for i in range(10):
                valid = i<k
                rec.append({
                    "photo_id":str(89000+len(rec)),
                    "observation_id":"obs-"+str(len(rec)),
                    "inat_taxon_id":str(1000+i),
                    "source_design":arm, "coarse_zone":zone,
                    "local_season":"summer",
                    "strict9":valid,
                    "hue":("white" if i<white_n else "red_pink") if valid else "mixed_uncertain",
                })
        return pd.DataFrame(rec)

    def test_MAR_reweighting_has_distinct_target_and_unrestricted_bounds(self):
        t=source_groups(self.fake_frame())
        result=group_summary(t,scope="synthetic")
        self.assertEqual(result["observations"],20)
        self.assertEqual(result["roi9_success"],10)
        self.assertAlmostEqual(result["complete_case_white_rate"],0.2)
        self.assertAlmostEqual(result["MAR_stratum_poststandardized_photo_white_rate"],0.5)
        self.assertEqual(result["nonparametric_white_fraction_identification_bounds"],[0.1,0.6])
        self.assertEqual(result["strata_without_any_valid_roi9"],0)
        self.assertTrue(result["half_white_possible_with_arbitrary_missing_colours"])
        minus=result["MNAR_delta_log_odds_scenarios"]["-2.0"]["scenario_white_fraction_if_all_strata_supported"]
        plus=result["MNAR_delta_log_odds_scenarios"]["+2.0"]["scenario_white_fraction_if_all_strata_supported"]
        self.assertLess(minus,plus)

    def test_missing_labels_not_smuggled_in_from_four_hue(self):
        d=self.fake_frame()
        groups=source_groups(d)
        original=group_summary(groups,scope="toy")
        d.loc[~d.strict9,"hue"]="white"  # technical failure remains absent
        revised=group_summary(source_groups(d),scope="toy")
        self.assertEqual(original["complete_case_white_rate"],
                         revised["complete_case_white_rate"])
        self.assertEqual(original["MAR_stratum_poststandardized_photo_white_rate"],
                         revised["MAR_stratum_poststandardized_photo_white_rate"])
        self.assertEqual(original["nonparametric_white_fraction_identification_bounds"],
                         revised["nonparametric_white_fraction_identification_bounds"])
        d.loc[0,"hue"]="not_a_colour"
        with self.assertRaisesRegex(ValueError,"invalid recorded"):
            source_groups(d)

    def test_no_support_rejects_false_MAR_point(self):
        d=self.fake_frame()
        d.loc[d.coarse_zone.eq("N_tropical"),"strict9"]=False
        result=group_summary(source_groups(d),scope="unsupported")
        self.assertIsNone(result["MAR_stratum_poststandardized_photo_white_rate"])
        self.assertEqual(result["strata_without_any_valid_roi9"],1)
        self.assertEqual(result["observations_without_stratum_colour_support"],10)
        zero=result["MNAR_delta_log_odds_scenarios"]["+0.0"]
        self.assertIsNone(zero["scenario_white_fraction_if_all_strata_supported"])
        self.assertGreater(zero["interval_if_no_colour_support"][1],
                           zero["interval_if_no_colour_support"][0])

    def test_odds_scenarios_do_not_replace_empirical_rates(self):
        self.assertAlmostEqual(logit_tilt(0.5,0),0.5)
        self.assertAlmostEqual(logit_tilt(0.5,math.log(3)),0.75)
        self.assertLess(logit_tilt(0.3,-2),0.3)
        self.assertGreater(logit_tilt(0.3,2),0.3)
        with self.assertRaises(ValueError):
            logit_tilt(-0.1,1)

    def test_actual_pipeline_contract_and_no_human_data(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            p=root/"toy.parquet"
            table=[]
            for i in range(96):
                good=i%3==0
                table.append({
                    "photo_id":str(800000+i),
                    "observation_id":"obs-id-"+str(i),
                    "inat_taxon_id":str(i//3+1),
                    "source_design":("historical_K1_FCP_source_heterogeneous" if i%2==0 else
                                     "bulk128_131072"),
                    "observed_on":"2025-07-15" if i%4 else "2025-01-12",
                    "latitude":-15.0 if i%3 else 44.0,
                    "longitude":127.2,
                    "strict9":good,
                    "original_roi_four_hue":("white" if i%6==0 else "red_pink") if good
                        else "mixed_uncertain",
                })
            pd.DataFrame(table).to_parquet(p,index=False)
            receipt=root/"receipt.json"
            receipt.write_text(json.dumps({
                "schema":"mois_signed1181146_trueROI9_bulk128_original_metadata_georestore_v2",
                "same_original_image_ID_and_original_9palette_unmodified":True,
                "source_georestored_original_parquet_SHA256":sha_file(p),
            }))
            result=run(p,receipt,root/"automatic",verify_counts=False)
            self.assertEqual(result["schema"],SCHEMA)
            self.assertEqual(result["source_distinct_observations"],96)
            self.assertEqual(result["frozen_roi9_unknown_observations"],64)
            self.assertEqual(result["strict9_unknowns_upgraded_to_known"],0)
            self.assertEqual(result["actual_human_labels"],0)
            self.assertTrue(result["MNAR_true_missing_colour_not_identified"])
            self.assertIsNotNone(result["overall"]["complete_case_white_rate"])
            for name,digest in result["output_SHA256"].items():
                self.assertEqual(sha_file(root/"automatic"/name),digest)
            with self.assertRaisesRegex(ValueError,"replace"):
                run(p,receipt,root/"automatic",verify_counts=False)
            receipt.write_text(json.dumps({
                "schema":"wrong",
                "source_georestored_original_parquet_SHA256":sha_file(p)
            }))
            with self.assertRaisesRegex(ValueError,"corrected v2"):
                run(p,receipt,root/"wrong",verify_counts=False)

if __name__=="__main__":
    unittest.main()
