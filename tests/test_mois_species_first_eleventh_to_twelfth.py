"""Eleventh authentic SHA provenance and twelfth image selection regressions."""
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis import mois_species_first_eleventh_to_twelfth as w
from scripts.analysis.mois_species_first_palette9 import PALETTE, sha256


class TestEleventhToTwelfth(unittest.TestCase):
    def make_workers(self, root):
        source = root / "source"
        (source / "shards").mkdir(parents=True)
        workers = root / "workers"
        receipts = {}
        sourcehashes = {}
        actual = []
        for i in range(8):
            n = np.arange(i * 256, (i + 1) * 256)
            photo_id = [str(40000000 + int(j)) for j in n]
            obs = [f"obs-independent-{j}" for j in n]
            taxa = [str(50000 + int(j)) for j in n]
            good = n < 1000
            raw = pd.DataFrame({
                "photo_id": photo_id,
                "observation_id": obs,
                "inat_taxon_id": taxa,
                "morph": np.where(good, "white", "unknown"),
                "measurement_status": np.where(good, "measured", "not_detected"),
                "image_sha256": ["original-pixel"]*256,
            })
            for hue in PALETTE:
                raw["flower_fraction_" + hue] = np.where(
                    good, "1.0" if hue == "white" else "0.0", "")
            selected = raw[["photo_id", "observation_id", "inat_taxon_id"]].copy()
            selected["shard"] = str(i)
            p = source / "shards" / f"shard-{i:03d}.csv"
            selected.to_csv(p, index=False)
            sourcehashes[str(i)] = sha256(p)
            folder = workers / f"mois-eleventh4096-real-original-worker-{i//4}" / f"shard_{i:03d}"
            folder.mkdir(parents=True)
            terminal = folder / "photos.csv"
            raw.to_csv(terminal, index=False)
            receipt = folder / "receipt.json"
            receipt.write_text(json.dumps({
                "schema": "mois_full_original93m_real_RGB_photo_unit_v1",
                "shard": i, "attempted": 256, "photo_colours_measured": True,
                "source_sha256": sha256(p), "photos_sha256": sha256(terminal),
            }))
            receipts[str(i)] = sha256(receipt)
            actual.append(raw)
        (source / "manifest.json").write_text(json.dumps({
            "schema": "mois_full_source_observation_sampling_cohort_v1",
            "source_arm": w.SOURCE_ARM,
            "colour_measured": False,
            "source_previous_attempts": w.PRIOR_PHOTOS,
            "source_prior_no_reused_original_photo_or_observation": True,
            "work_units": 8, "worker_ids": list(range(8)),
            "new_photo_attempts_planned": 2048,
            "photo_shard_sha256": sourcehashes,
        }))
        raw = pd.concat(actual, ignore_index=True)
        raw["strict9_original_photo_numeric"] = raw.morph.eq("white").astype(str)
        raw["dominant_hue9_label"] = np.where(raw.morph.eq("white"), "white", "")
        raw.loc[~raw.morph.eq("white"),
                ["flower_fraction_" + p for p in PALETTE]] = ""
        audited = root / "audit"
        audited.mkdir()
        p = audited / "mois_eleventh_actual_original_nine_roi.csv.gz"
        raw.to_csv(p, index=False, compression="gzip")
        (audited / "receipt.json").write_text(json.dumps({
            "schema": w.AUDIT_SCHEMA,
            "actual_independent_original_photo_attempts": 2048,
            "actual_genuine_strict9_numeric_success": 1000,
            "actual_clear_nine_palette_hues": 1000,
            "actual_technical_unknowns": 1048,
            "photo_observation_source_taxon_ID_collisions": 0,
            "frozen_original_manifest_sha256": sha256(source / "manifest.json"),
            "audited_derivative_sha256": sha256(p),
            "original_photo_worker_receipt_sha256": receipts,
        }))
        return source, workers, audited

    def test_signed_eleventh_original_numeric_and_unknown_integrity(self):
        with tempfile.TemporaryDirectory() as tmp:
            source, workers, audited = self.make_workers(Path(tmp))
            frame, proof = w.read_eleventh(source, workers, audited)
            self.assertEqual(len(frame), 2048)
            self.assertEqual(int(frame.strict9.sum()), 1000)
            self.assertEqual(int(frame.dominant_label9.ne("").sum()), 1000)
            self.assertTrue(frame.loc[~frame.strict9,
                         ["strict9_" + p for p in PALETTE]].isna().all().all())
            self.assertEqual(len(proof["eight_original_worker_receipt_sha256"]), 8)
            receipt = workers / "mois-eleventh4096-real-original-worker-1" / "shard_005" / "receipt.json"
            receipt.write_text(receipt.read_text() + " ")
            with self.assertRaises(ValueError):
                w.read_eleventh(source, workers, audited)

    def test_mixed_retry_preserves_independent_photos_and_first_coverage(self):
        taxa = [str(k) for k in range(1, 42112)]
        supported, labelled = w.PRIOR_SUPPORTED_TAXA, w.PRIOR_LABELLED_TAXA
        atlas = pd.DataFrame({
            "inat_taxon_id": taxa,
            "after_strict9": [1] * supported + [0] * (42111-supported),
            "after_photo_attempts": [1] * supported + [0] * (42111-supported),
            "independent_roi9_observations": [1] * supported + [0] * (42111-supported),
            "dominant_hue9": ["white"] * labelled + [""] * (42111-labelled),
            "observed_hue9_set": ["white"] * labelled + [""] * (42111-labelled),
            "hue_status": ["provisional_single"] * labelled +
                          ["mixed_image_unresolved"] * (supported-labelled) +
                          ["no_strict9"] * (42111-supported),
            "distinct_palette_votes": [1] * labelled + [0] * (42111-labelled),
            "ambiguous_observation_rois": [0] * labelled +
                                         [1] * (supported-labelled) + [0] * (42111-supported),
            "next_measure_priority": [2] * labelled +
                                      [1] * (supported-labelled) + [0] * (42111-supported),
        })
        for hue in PALETTE:
            atlas["observed_mean_" + hue] = (
                [1.0 if hue == "white" else 0.0] * supported
                + [np.nan] * (42111-supported))
        first = taxa[supported:supported + 1800]
        retry = taxa[labelled:labelled + 248]
        selected = first + retry
        assert len(selected) == 2048
        good = np.arange(2048) < 700
        good[1800:1900] = True  # 100 successful formerly numeric-mixed retries
        labelled_now = np.arange(2048) < 600
        labelled_now[1800:1880] = True  # 80 become clear after retry
        eleventh = pd.DataFrame({
            "photo_id": [str(7000000 + i) for i in range(2048)],
            "observation_id": [f"eleventh-unique-observation-{i}" for i in range(2048)],
            "inat_taxon_id": selected,
            "strict9": good,
            "dominant_label9": np.where(labelled_now, "yellow", ""),
        })
        for hue in PALETTE:
            eleventh["strict9_" + hue] = np.where(
                good, (1.0 if hue == "yellow" else 0.0), np.nan)
        # For real numerical strict9 photos with no dominant hue, assign
        # a valid ambiguous 50/50 mixture.
        inconclusive = good & ~labelled_now
        for hue in PALETTE:
            eleventh.loc[inconclusive, "strict9_" + hue] = (
                0.5 if hue in ("white", "yellow") else 0.0)
        old_ids = pd.DataFrame({"photo_id": ["5","6"],
                                "observation_id": ["past-1","past-2"],
                                "inat_taxon_id": ["1","2"]})
        possible = taxa[supported + 1800:supported + 1812]
        options = pd.DataFrame({
            "photo_id": [str(8000000 + j) for j in range(12)],
            "observation_id": [f"twelfth-obs-{j}" for j in range(12)],
            "inat_taxon_id": possible,
            "source_observation_rank": [1] * 12,
        })
        updated, queue, new_covered = w.append_with_retries(
            atlas, eleventh, old_ids, options, quota=8, expected_previous=2)
        self.assertEqual(new_covered, 700)
        self.assertEqual(int(updated.after_strict9.gt(0).sum()), supported + 700)
        self.assertEqual(int(updated.dominant_hue9.ne("").sum()), labelled + 680)
        self.assertTrue(queue.photo_id.is_unique)
        self.assertTrue(queue.acquisition_priority.eq("first_unmeasured").all())
        old_row = updated.set_index("inat_taxon_id").loc[retry[0]]
        self.assertEqual(int(old_row.independent_roi9_observations), 2)
        self.assertEqual(int(old_row.after_photo_attempts), 2)
        self.assertAlmostEqual(float(old_row.observed_mean_yellow), 0.5)
        self.assertAlmostEqual(float(old_row.observed_mean_white), 0.5)
        self.assertEqual(old_row.dominant_hue9, "yellow")


if __name__ == "__main__":
    unittest.main()
