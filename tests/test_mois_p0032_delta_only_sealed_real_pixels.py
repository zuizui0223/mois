"""Regression tests: only new p0032 frozen original-photo workers, never replay old."""
from __future__ import annotations

import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import pandas as pd

from scripts.acquisition import mois_p0032_delta_only_sealed_real_pixels as delta
from scripts.acquisition.mois_completed_pixel_checkpoint import sha
from scripts.analysis.mois_reassemble_original_new4096_measurements import METHOD


def asset(group, id):
    return {"id":id, "name":f"mois-93m-p0032-RGB-group-{group}",
            "digest":"sha256:"+f"{id:064x}", "size_in_bytes":200}


class OriginalIncrementalDeltaTests(unittest.TestCase):
    def setUp(self):
        self.tmp=TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.source=self.root/"source"
        self.worker=self.root/"workers"
        self.prev=self.root/"prior.parquet"
        self.receipt=self.root/"checkpoint_receipt.json"
        self.snapshot=self.root/"frozen_artifacts.json"
        self.output=self.root/"sealed_delta"
        (self.source/"shards").mkdir(parents=True)
        (self.worker/"RGB-group-0/shard_000/sealed").mkdir(parents=True)
        worker=self.worker/"RGB-group-0/shard_000"
        self.worker_dir=worker

        # One original photo was signed and classified in the existing main
        # checkpoint. Its actual original image bytes never enter this test.
        pd.DataFrame([{
            "photo_id":"100","morph":"white","checkpoint_wave":"p0032",
            "observation_id":"previous-observation",
        }]).to_parquet(self.prev,index=False)
        self.baseline_artifacts=[asset(1,1235)]
        source_rows=[
            {"measurement_id":"MOIS-001","photo_id":"101",
             "observation_id":"another-observation","source_arm":"p0032-synth",
             "selection_probability":"0.002","shard":"0"},
            {"measurement_id":"MOIS-002","photo_id":"102",
             "observation_id":"another-observation","source_arm":"p0032-synth",
             "selection_probability":"0.002","shard":"0"},
        ]

        def write_csv(file,rows):
            with file.open("w",newline="") as f:
                w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n")
                w.writeheader();w.writerows(rows)
        self.input_file=self.source/"shards/shard-000.csv"
        write_csv(self.input_file,source_rows)
        self.source_manifest={
            "schema":"mois_full_source_observation_sampling_cohort_v1",
            "new_photo_attempts_planned":2,
            "full_frame_photo_ids":93209351,
            "work_units":1,"worker_ids":[0],"group_ids":list(range(145)),
            "photo_shard_sha256":{"0":sha(self.input_file)},
            "colour_measured":False,
        }
        (self.source/"manifest.json").write_text(json.dumps(self.source_manifest))
        terminal=[
            {"measurement_id":"MOIS-001","morph":"white",
             "measurement_status":"classified_four_state_morph",
             "image_sha256":"f"*64},
            {"measurement_id":"MOIS-002","morph":"mixed_uncertain",
             "measurement_status":"not_evaluable_roi_or_flip_gate",
             "image_sha256":"e"*64},
        ]
        self.terminal=worker/"sealed/terminal_p000_b0000000.csv"
        write_csv(self.terminal,terminal)
        photos=[{**r,**out,"classified_four_hue":str(out["morph"]=="white").lower()}
                for r,out in zip(source_rows,terminal)]
        self.photos=worker/"photos.csv"
        write_csv(self.photos,photos)
        self.photo_receipt=worker/"receipt.json"
        self.photo_receipt.write_text(json.dumps({
            "schema":"mois_full_original93m_real_RGB_photo_unit_v1",
            "shard":0,"attempted":2,"downloaded":2,"classifiable":1,
            "acquisition_failed":0,"unclassifiable_acquired":1,
            "source_sha256":sha(self.input_file),
            "method":METHOD,"terminal_sha256":sha(self.terminal),
            "photos_sha256":sha(self.photos),
            "photo_colours_measured":True,"raw_licensed_image_bytes_persisted":False,
        }))
        self.old_receipt={
            "schema":"mois_sealed_actual_pixel_checkpoint_v1",
            "image_downloads_or_colour_reinference_performed":False,
            "combined_unique_photo_attempts":1,"combined_four_hue_classified":1,
            "output_sha256":{self.prev.name:sha(self.prev)},
            "source_artifact_snapshots":{
                "p0032":{"source_run_id":delta.SOURCE_RUN,"artifacts":self.baseline_artifacts}
            }
        }
        self.receipt.write_text(json.dumps(self.old_receipt))
        self.live_artifacts=[*self.baseline_artifacts,asset(0,3456)]
        self.snapshot.write_text(json.dumps({
            "source_run_id":delta.SOURCE_RUN,
            "artifacts":self.live_artifacts,
        }))

    def run_delta(self):
        return delta.verify_new_pixel_delta(
            self.prev,self.receipt,self.source,self.worker,self.snapshot,self.output,
            prior_attempts=1,prior_hues=1,
            source_hash=sha(self.source/"manifest.json"),
            expected_photo_population=2)

    def test_old_p0032_worker_is_not_downloaded_and_two_new_original_pixels_seal(self):
        r=self.run_delta()
        self.assertEqual(r["new_sealed_p0032_groups"],[0])
        self.assertEqual(r["completed_true_original_photo_attempts"],2)
        self.assertEqual(r["new_true_original_four_hue_classified"],1)
        self.assertEqual(r["combined_disjoint_main_plus_delta_photo_attempts"],3)
        self.assertTrue(r["previous_photo_ids_and_new_pixel_ids_exactly_disjoint"])
        self.assertEqual(r["other_original_wave_worker_archives_redownloaded"],0)
        self.assertFalse(r["previous_photo_dataset_rewritten"])
        self.assertEqual(r["new_unique_original_observations"],1)
        measured=pd.read_parquet(self.output/"new_p0032_original_photo_ROIs.parquet")
        self.assertEqual(set(measured.photo_id),{"101","102"})
        self.assertEqual(set(measured.morph),{"white","mixed_uncertain"})

    def test_same_snapshot_only_is_not_new_pixels(self):
        self.snapshot.write_text(json.dumps({
            "source_run_id":delta.SOURCE_RUN,
            "artifacts":self.baseline_artifacts,
        }))
        with self.assertRaisesRegex(ValueError,"No genuinely new"):
            self.run_delta()
        self.assertFalse(self.output.exists())

    def test_existing_historical_photo_artifact_digest_cannot_drift(self):
        changed=dict(self.baseline_artifacts[0])
        changed["digest"]="sha256:"+"a"*64
        self.snapshot.write_text(json.dumps({
            "source_run_id":delta.SOURCE_RUN,
            "artifacts":[changed,asset(0,3456)],
        }))
        with self.assertRaisesRegex(ValueError,"Previously signed"):
            self.run_delta()
        self.assertFalse(self.output.exists())

    def test_duplicate_new_photo_id_against_baseline_is_rejected(self):
        p=pd.read_parquet(self.prev)
        p.loc[0,"photo_id"]="101"
        p.to_parquet(self.prev,index=False)
        r=json.loads(self.receipt.read_text())
        r["output_sha256"][self.prev.name]=sha(self.prev)
        self.receipt.write_text(json.dumps(r))
        with self.assertRaisesRegex(ValueError,"processed twice"):
            self.run_delta()
        self.assertFalse(self.output.exists())

    def test_altered_original_photo_terminal_or_existing_output_fails(self):
        self.terminal.write_text(self.terminal.read_text()+"\nchanged")
        with self.assertRaisesRegex(ValueError,"terminal SHA mismatch"):
            self.run_delta()
        self.assertFalse(self.output.exists())

    def test_invalid_original_observation_probability_is_not_accepted(self):
        data=list(csv.DictReader(self.photos.open(newline="")))
        data[1]["selection_probability"]="0.009"
        with self.photos.open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=list(data[0]),lineterminator="\n")
            w.writeheader();w.writerows(data)
        p=json.loads(self.photo_receipt.read_text())
        p["photos_sha256"]=sha(self.photos)
        self.photo_receipt.write_text(json.dumps(p))
        with self.assertRaisesRegex(ValueError,"source identity changed|sampling probability differs"):
            self.run_delta()


if __name__ == "__main__":
    unittest.main()
