"""Original tier06 genuine image worker groups become append-only photo Parquets."""
import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import pandas as pd

from scripts.acquisition import mois_append_only_second8_original_worker_group as m
from scripts.analysis.mois_reassemble_original_new4096_measurements import METHOD


class AppendOnlySecondEightOriginalTests(unittest.TestCase):
    def setUp(self):
        self.tmp = TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.source = self.root / "source"
        self.worker = self.root / "workers"
        (self.source / "shards").mkdir(parents=True)
        target = self.worker / "shard_000"
        (target / "sealed").mkdir(parents=True)
        self.target = target
        self.output = self.root / "sealed"
        self.rows = []
        for i in range(2):
            photo = str(1_234_500 + i)
            self.rows.append({
                "measurement_id": "MOIS-TEST" + str(i),
                "photo_id": photo,
                "inat_taxon_id": "1",
                "species": "Asteraceae test",
                "observation_id": "obs" + str(i),
                "observer_id": "observer",
                "source_arm": "global93m_tier06_all16_buckets004to011_broad_warm_R2",
                "selection_probability": format(0.0064 * 4 / 512, ".17g"),
                "shard": "0",
            })
        def csv_write(file, rows):
            with file.open("w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
        source_file = self.source / "shards" / "shard-000.csv"
        csv_write(source_file, self.rows)
        manifest = {
            "schema": "mois_full_source_observation_sampling_cohort_v1",
            "full_frame_photo_ids": 93_209_351,
            "source_tier": 6,
            "source_group": "all16",
            "source_arm": self.rows[0]["source_arm"],
            "complete_tier06_source_photo_IDs": 597070,
            "compute_bucket_count": 512,
            "admitted_bucket_count": 8,
            "admitted_bucket_start": 4,
            "admitted_buckets": list(range(4,12)),
            "previous_broad_original_photo_attempts": 4468,
            "not_probability_complete_global_source": True,
            "native_photo_status_not_inferred_or_hard_filtered": True,
            "colour_measured": False,
            "new_photo_attempts_planned": 2,
            "work_units": 1,
            "worker_ids": [0],
            "group_ids": [0],
            "photo_shard_sha256": {"0": m.sha(source_file)},
            "selection_probability": 0.0064 * 4 / 512,
        }
        (self.source / "manifest.json").write_text(json.dumps(manifest))
        terminal = []
        photos = []
        for i, row in enumerate(self.rows):
            known = i == 0
            extra = {
                "measurement_id": row["measurement_id"],
                "morph": "white" if known else "mixed_uncertain",
                "measurement_status": "classified_four_state_morph" if known else "not_evaluable_roi_or_flip_gate",
                "image_sha256": "f" * 64,
                "flower_fraction_white": "1" if known else "0",
            }
            terminal.append(extra)
            photos.append({
                **row, **extra,
                "classified_four_hue": str(known).lower(),
            })
        term_file = target / "sealed/terminal_p000_b0000000.csv"
        csv_write(term_file, terminal)
        photos_file = target / "photos.csv"
        csv_write(photos_file, photos)
        original_receipt = {
            "schema": "mois_full_original93m_real_RGB_photo_unit_v1",
            "shard": 0,
            "attempted": 2,
            "downloaded": 2,
            "classifiable": 1,
            "acquisition_failed": 0,
            "unclassifiable_acquired": 1,
            "source_sha256": m.sha(source_file),
            "method": METHOD,
            "terminal_sha256": m.sha(term_file),
            "photos_sha256": m.sha(photos_file),
            "photo_colours_measured": True,
            "raw_licensed_image_bytes_persisted": False,
        }
        (target / "receipt.json").write_text(json.dumps(original_receipt))

    def seal(self, **kwargs):
        return m.seal(self.source, self.worker, 0, self.output,
                      expected_source_sha=m.sha(self.source / "manifest.json"),
                      expected_attempts=2, **kwargs)

    def test_terminal_source_colour_and_unknowns_reconciled(self):
        receipt = self.seal()
        self.assertEqual(receipt["distinct_genuinely_completed_photo_ids"], 2)
        self.assertEqual(receipt["four_hue_model_classified"], 1)
        self.assertEqual(receipt["acquired_but_unclassifiable"], 1)
        self.assertFalse(receipt["image_bytes_reacquired_in_sealing"])
        self.assertFalse(receipt["globally_probability_complete"])
        file = self.output / "original_photo_ROIs.parquet"
        self.assertEqual(receipt["saved_file_sha256"][file.name], m.sha(file))
        photos = pd.read_parquet(file)
        self.assertEqual(set(photos["morph"]), {"white", "mixed_uncertain"})
        self.assertEqual(set(photos["photo_id"]), {"1234500", "1234501"})
        self.assertEqual(len(receipt["original_work_units"]), 1)

    def test_previously_sealed_partition_cannot_be_overwritten(self):
        self.seal()
        with self.assertRaisesRegex(ValueError, "overwrite"):
            self.seal()

    def test_unknown_group_source_or_sha_cannot_be_sealed(self):
        with self.assertRaisesRegex(ValueError, "source SHA"):
            m.seal(self.source, self.worker, 0, self.output,
                   expected_source_sha="0" * 64, expected_attempts=2)
        with self.assertRaisesRegex(ValueError, "not scheduled"):
            m.seal(self.source, self.worker, 1, self.output,
                   expected_source_sha=m.sha(self.source / "manifest.json"),
                   expected_attempts=2)
        self.assertFalse(self.output.exists())

    def test_tampered_terminal_photo_file_is_rejected(self):
        path = self.target / "photos.csv"
        path.write_text(path.read_text() + "\ncorrupted")
        with self.assertRaisesRegex(ValueError, "SHA mismatch"):
            self.seal()
        self.assertFalse(self.output.exists())


if __name__ == "__main__":
    unittest.main()
