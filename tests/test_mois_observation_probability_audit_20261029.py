"""No outcome leakage, no invented source probabilities, reproducible 1200-photo design."""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from scripts.analysis.mois_observation_probability_audit_20261029 import (
    allocate, draw_audit, freeze, observation_frame, sha_file
)

class ObservationProbabilityAuditTests(unittest.TestCase):
    def toy(self):
        records = []
        arms = ("historical_K1_FCP_source_heterogeneous", "bulk128_131072")
        zones = (42.0, -14.0, -33.0, 13.0)
        for status in (False, True):
            for arm in arms:
                for zi, lat in enumerate(zones):
                    for season, month in enumerate((1, 4, 7, 10)):
                        for i in range(18):
                            oid = "test-obs-%s-%s-%d-%d-%d" % (
                                int(status), arm, zi, season, i
                            )
                            pid = str(10000000 + len(records))
                            records.append({
                                "photo_id": pid, "observation_id": oid,
                                "inat_taxon_id": "taxon-" + str(zi) + str(i),
                                "source_design": arm,
                                "observed_on": "2025-%02d-10" % month,
                                "latitude": lat, "longitude": 123.0,
                                "strict9": status,
                            })
        # Duplicate photo of first observation, with the opposite status.
        extra = records[0].copy()
        extra["photo_id"] = "999999999999"
        extra["strict9"] = True
        records.append(extra)
        return pd.DataFrame(records)

    def test_positive_probability_allocation_and_exact_budget(self):
        out = allocate({"a": 3, "b": 100, "c": 1000}, 50)
        self.assertEqual(sum(out.values()), 50)
        self.assertTrue(all(0 < out[k] <= n for k, n in {
            "a": 3, "b": 100, "c": 1000}.items()))
        self.assertEqual(allocate({"a": 1, "b": 1, "c": 9}, 3),
                         {"a": 1, "b": 1, "c": 1})
        with self.assertRaisesRegex(ValueError, "budget"):
            allocate({"a": 100, "b": 100, "c": 100}, 2)

    def test_canonical_observation_choice_is_outcome_and_order_blind(self):
        d = self.toy()
        f1, report = observation_frame(d)
        self.assertEqual(report["source_original_photo_attempts"], len(d))
        self.assertEqual(report["additional_photos_of_already_represented_observations"], 1)
        self.assertTrue(f1.observation_id.is_unique)
        f2, _ = observation_frame(d.sample(frac=1, random_state=42).reset_index(drop=True))
        a = dict(zip(f1.observation_id, f1.photo_id))
        self.assertEqual(a, dict(zip(f2.observation_id, f2.photo_id)))
        d["strict9"] = ~d["strict9"]
        f3, _ = observation_frame(d)
        self.assertEqual(a, dict(zip(f3.observation_id, f3.photo_id)))
        d.loc[0, "photo_id"] = d.loc[1, "photo_id"]
        with self.assertRaisesRegex(ValueError, "unique"):
            observation_frame(d)

    def test_review_sample_inverse_q_recovers_frame_stratum_sizes(self):
        frame, _ = observation_frame(self.toy())
        s, strata = draw_audit(frame, n_success=80, n_unknown=120)
        self.assertEqual(len(s), 200)
        self.assertEqual(s.roi_status.value_counts().to_dict(),
                         {"success": 80, "unknown": 120})
        self.assertTrue(s.photo_id.is_unique)
        self.assertTrue(s.observation_id.is_unique)
        self.assertTrue((s.audit_inclusion_q > 0).all())
        self.assertTrue((s.audit_inclusion_q <= 1).all())
        for key, chunk in s.groupby("_stratum"):
            actual = next(row for row in strata.itertuples(index=False) if row.stratum == key)
            self.assertEqual(len(chunk), actual.sampled_observations)
            self.assertAlmostEqual(float(chunk.review_inverse_q.sum()),
                                   actual.frame_observations, places=8)
        self.assertEqual(int(strata.frame_observations.sum()), len(frame))

    def test_real_metadata_contract_and_blank_annotation_csv(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / "toy.parquet"
            receipt = root / "geo.json"
            self.toy().to_parquet(path, index=False)
            receipt.write_text(json.dumps({
                "schema": "mois_signed1181146_trueROI9_bulk128_original_metadata_georestore_v2",
                "same_original_image_ID_and_original_9palette_unmodified": True,
                "source_georestored_original_parquet_SHA256": sha_file(path)
            }))
            out = root / "frozen"
            done = freeze(path, receipt, out, verify_frozen=False,
                          n_success=80, n_unknown=120)
            self.assertEqual(done["sample_total"], 200)
            self.assertTrue(done["audit_status_all_human_labels_blank"])
            self.assertFalse(done["native_wild_flora_population_probability_identified"])
            roster = pd.read_csv(out / "audit1200_source_observation_roster.csv")
            self.assertEqual(len(roster), 200)
            self.assertNotIn("strict9_white", roster.columns)
            self.assertNotIn("original_roi_four_hue", roster.columns)
            self.assertNotIn("latitude", roster.columns)
            self.assertNotIn("longitude", roster.columns)
            blank = pd.read_csv(out / "audit1200_human_annotation_BLANK.csv")
            self.assertTrue(blank.human_corolla_colour.isna().all())
            self.assertTrue(blank.human_open_corolla_visible.isna().all())
            self.assertEqual(len(blank), 200)
            with self.assertRaisesRegex(ValueError, "overwrite"):
                freeze(path, receipt, out, verify_frozen=False,
                       n_success=80, n_unknown=120)
            receipt.write_text(json.dumps({
                "source_georestored_original_parquet_SHA256": "badSHA"
            }))
            with self.assertRaisesRegex(ValueError, "signed georestoration"):
                freeze(path, receipt, root / "bad", verify_frozen=False,
                       n_success=80, n_unknown=120)

if __name__ == "__main__":
    unittest.main()
