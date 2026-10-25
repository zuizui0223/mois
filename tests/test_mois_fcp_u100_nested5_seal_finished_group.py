"""Synthetic original photo IDs and authentic signed K1-style terminal outcomes."""
from __future__ import annotations

import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import pandas as pd
from scripts.acquisition import mois_fcp_u100_nested5_seal_finished_group as m
from scripts.acquisition.mois_completed_pixel_checkpoint import sha
from scripts.analysis.mois_reassemble_original_new4096_measurements import METHOD


class FCP5FinishedGroupOriginalPhotoUnitTests(unittest.TestCase):
    def setUp(self):
        self.temp=TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.source=self.root/"source"
        self.workers=self.root/"worker"
        (self.source/"shards").mkdir(parents=True)
        (self.workers/"shard_000"/"sealed").mkdir(parents=True)
        self.output=self.root/"new-original-group"
        self.rows=[]
        for n in range(2):
            self.rows.append(dict(
                photo_id=str(700000+n),measurement_id=f"MOIS-FCP5-TEST-{n}",
                inat_taxon_id=str(n+100),observation_id=f"ORIGINAL-UUID-{n}",
                source_arm=m.ARM,shard="0",
                selection_probability="0.024" if n==0 else "0.0098",
                observed_on="2026-05-04" if n==0 else "",
                latitude="15.1" if n==0 else "-45.3",
                image_filename=f"MOIS-FCP5-TEST-{n}.jpg",
                photo_license="cc-by",
            ))
        self.source_shard=self.source/"shards/shard-000.csv"
        self.csv_write(self.source_shard,self.rows)
        self.source_manifest=self.source/"manifest.json"
        self.source_manifest.write_text(json.dumps({
            "schema":"mois_full_source_observation_sampling_cohort_v1",
            "full_frame_photo_ids":93209351,"colour_measured":False,
            "source_arm":m.ARM,
            "new_photo_attempts_planned":2,
            "work_units":1,"worker_ids":[0],"group_ids":[0],
            "photo_shard_sha256":{"0":sha(self.source_shard)}
        }))
        terminal=[]
        for n,row in enumerate(self.rows):
            terminal.append({
                "measurement_id":row["measurement_id"],
                "morph":"white" if n==0 else "mixed_uncertain",
                "measurement_status":("classified_four_state_morph"
                    if n==0 else "not_evaluable_roi_or_flip_gate"),
                "image_sha256":("a" if n==0 else "b")*64,
            })
        self.terminal=self.workers/"shard_000/sealed/terminal_p000_b0000000.csv"
        self.csv_write(self.terminal,terminal)
        self.photo=self.workers/"shard_000/photos.csv"
        self.csv_write(self.photo,[
            {**row,**end,
             "classified_four_hue":str(end["morph"]=="white").lower()}
            for row,end in zip(self.rows,terminal)
        ])
        self.unit=self.workers/"shard_000/receipt.json"
        self.unit.write_text(json.dumps({
            "schema":"mois_full_original93m_real_RGB_photo_unit_v1",
            "shard":0,"photo_colours_measured":True,
            "raw_licensed_image_bytes_persisted":False,"method":METHOD,
            "source_sha256":sha(self.source_shard),
            "photos_sha256":sha(self.photo),"terminal_sha256":sha(self.terminal),
            "attempted":2,"downloaded":2,"classifiable":1,
            "acquisition_failed":0,"unclassifiable_acquired":1,
        }))
        self.diagnostic=self.workers/"warm_roi_group_000_diagnostic.json"
        self.diagnostic.write_text(json.dumps({
            "group":0,"frozen_source_blob":m.FROZEN_SOURCE_BLOB,
            "frozen_yolo_and_onnx_initializations":1,
            "separate_original_main_calls":1,"new_photos_downloaded":2,
            "method_or_thresholds_modified":False,
            "old_photo_RGB_workers_modified":False,
            "artificial_flower_labels_created":False,
            "group_wall_seconds":12.0,"total_frozen_evaluator_seconds":8.0,
            "photo_prefetch":{"download_calls":2,"consumed":2,
                              "permanent_failures":0},
        }))

    @staticmethod
    def csv_write(path, rows):
        with path.open("w",newline="",encoding="utf-8") as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator="\n")
            writer.writeheader();writer.writerows(rows)

    def run_seal(self,**kw):
        # Source-wide 35 group SHA/selection validation is covered by
        # the real stage5 preflight. Synthetic fixture has one work unit only.
        with patch.object(m,"check_source",return_value={"schema":"tested"}):
            return m.seal(self.source,self.workers,0,self.output,
                expected_sha=sha(self.source_manifest),
                planned=2,units=1,groups=1,
                old_count=1009032,reused=3,candidates=4,all_first5=5,
                eligible=20242,**kw)

    def test_signed_group_preserves_real_flower_hue_and_technical_unknown(self):
        d=self.run_seal()
        self.assertEqual(d["new_actual_original_photo_ID_attempts"],2)
        self.assertEqual(d["new_actual_four_hue_classified"],1)
        self.assertEqual(d["new_original_image_access_failed"],0)
        self.assertEqual(d["new_original_acquired_hue_unclassifiable"],1)
        self.assertEqual(d["new_original_hue_counts"]["yellow_orange"],0)
        self.assertEqual(d["new_original_hue_counts"]["blue_purple"],0)
        self.assertEqual(d["distinct_original_per_photo_joint_probabilities"],2)
        self.assertEqual(d["original_image_bytes_downloaded_OR_colour_ROI_inferred_by_sealer"],0)
        self.assertFalse(d["botanical_wild_native_open_flower_status_independently_verified"])
        dest=self.output/"original_photo_ROIs.parquet"
        self.assertEqual(sha(dest),d["saved_file_SHA256"][dest.name])
        self.assertEqual(set(pd.read_parquet(dest).morph),{"white","mixed_uncertain"})
        cov=json.loads((self.output/"original_photo_source_coverage.json").read_text())
        self.assertFalse(cov["individual_photo_open_flowers_independently_verified"])

    def test_reject_actual_terminal_byte_change(self):
        self.terminal.write_text(self.terminal.read_text()+"\nforged")
        with self.assertRaisesRegex(ValueError,"Original terminal SHA mismatch"):
            self.run_seal()
        self.assertFalse(self.output.exists())

    def test_original_photo_probabilities_are_not_recomputed_or_scalar(self):
        data=list(csv.DictReader(self.photo.open()))
        data[0]["selection_probability"]="0"
        self.csv_write(self.photo,data)
        source_receipt=json.loads(self.unit.read_text())
        source_receipt["photos_sha256"]=sha(self.photo)
        self.unit.write_text(json.dumps(source_receipt))
        with self.assertRaisesRegex(ValueError,"source identity changed"):
            self.run_seal()

    def test_reject_swapped_undocumented_real_original_photo(self):
        data=list(csv.DictReader(self.photo.open()))
        data[0]["photo_id"]="88888"
        self.csv_write(self.photo,data)
        unit=json.loads(self.unit.read_text())
        unit["photos_sha256"]=sha(self.photo)
        self.unit.write_text(json.dumps(unit))
        with self.assertRaisesRegex(ValueError,"source identity changed"):
            self.run_seal()

    def test_reject_forged_model_source_weight_or_extra_session(self):
        diag=json.loads(self.diagnostic.read_text())
        diag["frozen_yolo_and_onnx_initializations"]=2
        self.diagnostic.write_text(json.dumps(diag))
        with self.assertRaisesRegex(ValueError,"worker/provenance timing changed"):
            self.run_seal()

    def test_permanent_image_acquisition_error_stays_unknown(self):
        terminal=list(csv.DictReader(self.terminal.open()))
        terminal[1]["measurement_status"]="image_acquisition_failed"
        terminal[1]["image_sha256"]=""
        self.csv_write(self.terminal,terminal)
        photos=list(csv.DictReader(self.photo.open()))
        photos[1]["measurement_status"]="image_acquisition_failed"
        photos[1]["image_sha256"]=""
        self.csv_write(self.photo,photos)
        unit=json.loads(self.unit.read_text())
        unit.update(photos_sha256=sha(self.photo),
                    terminal_sha256=sha(self.terminal),
                    downloaded=1,acquisition_failed=1,
                    unclassifiable_acquired=0)
        self.unit.write_text(json.dumps(unit))
        diag=json.loads(self.diagnostic.read_text())
        diag["photo_prefetch"]["permanent_failures"]=1
        self.diagnostic.write_text(json.dumps(diag))
        got=self.run_seal()
        self.assertEqual(got["new_original_image_access_failed"],1)
        self.assertEqual(got["new_original_hue_counts"]["mixed_uncertain"],1)

    def test_missing_original_work_unit_must_not_be_called_measured(self):
        self.unit.unlink()
        with self.assertRaisesRegex(ValueError,"Incomplete or duplicate"):
            self.run_seal()

    def test_released_finished_group_retries_transient_403_then_verifies_all_asset_SHAs(self):
        # The signed source pixels/terminal are produced once. If GitHub
        # release publishing briefly returns HTTP 403, retry metadata only,
        # and refuse partially created or mismatched immutable assets.
        import yaml
        workflow=(Path(__file__).resolve().parents[1] / ".github/workflows" /
                  "mois-fcp-u100-nested5-34898-append-finished-original-photo-groups.yml")
        flow=yaml.safe_load(workflow.read_text(encoding="utf-8"))
        step=next(x for x in flow["jobs"]["seal"]["steps"]
                  if x.get("name","").startswith("Publish immutable signed original"))
        run=step["run"]
        self.assertIn("for attempt in 1 2 3; do",run)
        self.assertIn("gh release create \"$TAG\"",run)
        self.assertIn("gh api \"repos/$REPO/releases/tags/$TAG\"",run)
        self.assertIn("'original_group_receipt.json'",run)
        self.assertIn("'original_photo_ROIs.parquet'",run)
        self.assertIn("'original_photo_hue_counts.csv'",run)
        self.assertIn("'original_photo_source_coverage.json'",run)
        self.assertIn("hashlib.sha256(source.read_bytes()).hexdigest()",run)
        self.assertIn("asset['digest']==digest",run)
        self.assertIn("asset['size']==source.stat().st_size",run)
        self.assertIn("Missing, duplicated or extra FCP5 original ROI release asset",run)
        self.assertIn("source_pixels_re_downloaded_or_re_inferred",run)

    def test_immutable_original_photo_group_no_overwrite(self):
        self.run_seal()
        with self.assertRaisesRegex(ValueError,"Never replace immutable"):
            self.run_seal()


if __name__=="__main__":
    unittest.main()
