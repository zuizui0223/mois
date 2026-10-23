"""Sealed original p0032 workers133..144: source hash, history ZIP identity, ROI photo terminal parity."""
import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import pandas as pd

from scripts.acquisition import mois_append_p0032_postcheckpoint_source_worker as m
from scripts.analysis.mois_reassemble_original_new4096_measurements import METHOD


class OriginalP0032WorkerTests(unittest.TestCase):
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
                "source_arm": "original_p0032_real_photo_source",
                "selection_probability": format(0.0032, ".17g"),
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
                "new_photo_attempts_planned": 2,
            "work_units": 1,
            "worker_ids": [0],
            "group_ids": list(range(145)),
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
        with patch.object(m,"BEGIN_GROUP",0):
            return m.seal(self.source, self.worker, 0, self.output,
                          expected_sha=m.sha(self.source / "manifest.json"),
                          expected_attempts=2, **kwargs)

    def test_terminal_source_colour_and_unknowns_reconciled(self):
        receipt = self.seal()
        self.assertEqual(receipt["original_true_photo_attempts"], 2)
        self.assertEqual(receipt["original_FCP_roi_four_hue_classified"], 1)
        self.assertEqual(receipt["original_acquired_unclassifiable"], 1)
        self.assertFalse(receipt["new_S3_photo_downloads_or_model_reinferences"])
        self.assertTrue(receipt["globally_biological_flower_hue_rate_unestimated"])
        file = self.output / "original_photo_ROIs.parquet"
        self.assertEqual(receipt["actual_output_sha256"][file.name], m.sha(file))
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
            with patch.object(m,"BEGIN_GROUP",0):
                m.seal(self.source, self.worker, 0, self.output,
                       expected_sha="0" * 64, expected_attempts=2)
        with self.assertRaisesRegex(ValueError, "Missing source photo units"):
            with patch.object(m,"BEGIN_GROUP",0):
                m.seal(self.source, self.worker, 1, self.output,
                       expected_sha=m.sha(self.source / "manifest.json"),
                       expected_attempts=2)
        self.assertFalse(self.output.exists())

    def test_tampered_terminal_photo_file_is_rejected(self):
        path = self.target / "photos.csv"
        path.write_text(path.read_text() + "\ncorrupted")
        with self.assertRaisesRegex(ValueError, "SHA mismatch"):
            self.seal()
        self.assertFalse(self.output.exists())


class FrozenFutureP0032ArtifactTests(unittest.TestCase):
    def setUp(self):
        self.tmp=TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.p=self.root/"prior.json"
        self.now=self.root/"live.json"
        self.previous=[]
        for i in range(m.BEGIN_GROUP):
            self.previous.append({
                "id":i+1000,
                "name":f"mois-93m-p0032-RGB-group-{i}",
                "digest":"sha256:"+f"{i+12345:064x}",
                "size_in_bytes":250000,
            })
        self.p.write_text(json.dumps({
            "source_run_id":m.RUN,
            "artifacts":self.previous,
        }))
        self.now.write_text(json.dumps({
            "source_run_id":m.RUN,
            "artifacts":self.previous+[{
                "id":3001,
                "name":f"mois-93m-p0032-RGB-group-{m.BEGIN_GROUP}",
                "digest":"sha256:"+"a"*64,
                "size_in_bytes":270000,
            }],
        }))

    def test_only_next_unseen_group_is_selected_and_old133_group_hashes_preserved(self):
        record=m.verify_frozen_original_artifacts(self.p,self.now)
        self.assertEqual(record["new_completed_original_photo_work_groups"],[133])
        self.assertEqual(len(record["prior_completed_original_photo_work_groups"]),133)
        self.assertTrue(record["prior_source_worker_SHA_records_unchanged"])
        self.assertEqual(record["new_original_photo_downloads_or_roi_reinferences"],0)

    def test_earlier_completed_original_group_zip_tamper_is_rejected(self):
        live=json.loads(self.now.read_text())
        live["artifacts"][5]["digest"]="sha256:"+"b"*64
        self.now.write_text(json.dumps(live))
        with self.assertRaisesRegex(ValueError,"changed or disappeared"):
            m.verify_frozen_original_artifacts(self.p,self.now)

    def test_no_fresh_worker_artifacts_is_valid_empty_gate(self):
        self.now.write_text(self.p.read_text())
        record=m.verify_frozen_original_artifacts(self.p,self.now)
        self.assertEqual(record["new_completed_original_photo_work_groups"],[])


if __name__ == "__main__":
    unittest.main()
