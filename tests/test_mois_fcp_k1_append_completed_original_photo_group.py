"""FCP K1 completed original photo SHA/method/propensity append-only regression."""
from __future__ import annotations
import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import pandas as pd
from scripts.acquisition import mois_fcp_k1_append_completed_original_photo_group as m
from scripts.acquisition.mois_completed_pixel_checkpoint import sha
from scripts.analysis.mois_reassemble_original_new4096_measurements import METHOD


class TrueFrozenK1GroupTests(unittest.TestCase):
    def setUp(self):
        self.tmp=TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.source=self.root/"source"
        self.worker=self.root/"original-workers"
        self.photo=self.worker/"shard_000"
        (self.photo/"sealed").mkdir(parents=True)
        (self.source/"shards").mkdir(parents=True)
        self.output=self.root/"append"
        self.rows=[]
        for i in range(2):
            self.rows.append(dict(
                photo_id=str(4001+i),measurement_id=f"FCP_K1_TEST_{i}",
                observation_id=f"OBS_{i}",inat_taxon_id="1" if i==0 else "2",
                observed_on="2025-06-11" if i==0 else "",
                latitude="42" if i==0 else "-11.2",
                image_filename=f"FCP_K1_TEST_{i}.jpg",
                photo_license="cc-by",source_arm=m.SOURCE_ARM,
                selection_probability="0.0008" if i==0 else "0.0002",
                shard="0"))
        def csv_write(path,rows):
            with path.open("w",newline="",encoding="utf-8") as f:
                w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n")
                w.writeheader();w.writerows(rows)
        self.csv_write=csv_write
        shard=self.source/"shards/shard-000.csv"
        csv_write(shard,self.rows)
        self.plan=dict(
            schema="mois_full_source_observation_sampling_cohort_v1",
            source_arm=m.SOURCE_ARM,full_frame_photo_ids=93209351,
            source_precolour_balanced_manifest_SHA256=m.PARENT_SOURCE_SHA,
            new_photo_attempts_planned=2,
            all_balanced_selected_original_observations=3,
            reused_original_photo_ids=1,
            original_photo_source_probabilities_unchanged_by_prior_pixel_reuse=True,
            images_downloaded_by_this_source_bridge=0,
            colour_measured=False,
            work_units=1,worker_ids=[0],group_ids=[0],
            photo_shard_sha256={"0":sha(shard)},
        )
        (self.source/"manifest.json").write_text(json.dumps(self.plan))
        terminal=[]
        for i in range(2):
            terminal.append(dict(
                measurement_id=self.rows[i]["measurement_id"],
                morph="white" if i==0 else "mixed_uncertain",
                measurement_status="classified_four_state_morph" if i==0 else "not_evaluable_roi_or_flip_gate",
                image_sha256="a"*64 if i==0 else "b"*64,
            ))
        self.terminal=self.photo/"sealed/terminal_p000_b0000000.csv"
        csv_write(self.terminal,terminal)
        measured=[]
        for src,end in zip(self.rows,terminal):
            measured.append({**src,**end,
                             "classified_four_hue":str(end["morph"]=="white").lower()})
        self.photos=self.photo/"photos.csv"
        csv_write(self.photos,measured)
        self.receipt=self.photo/"receipt.json"
        self.receipt.write_text(json.dumps(dict(
            schema="mois_full_original93m_real_RGB_photo_unit_v1",
            shard=0,photo_colours_measured=True,
            raw_licensed_image_bytes_persisted=False,
            method=METHOD,source_sha256=sha(shard),
            photos_sha256=sha(self.photos),
            terminal_sha256=sha(self.terminal),
            attempted=2,downloaded=2,classifiable=1,acquisition_failed=0,
            unclassifiable_acquired=1,
        )))
        self.diagnostic=self.worker/"warm_roi_group_000_diagnostic.json"
        self.diagnostic.write_text(json.dumps(dict(
            group=0,frozen_source_blob=m.FROZEN_SOURCE_BLOB,
            separate_original_main_calls=1,
            frozen_yolo_and_onnx_initializations=1,
            new_photos_downloaded=2,
            group_wall_seconds=8.0,total_frozen_evaluator_seconds=5.0,
            method_or_thresholds_modified=False,
            old_photo_RGB_workers_modified=False,
            photo_prefetch={"download_calls":2,"loader_seconds":2.7,
                            "consumer_wait_seconds":0.5},
        )))
    def seal(self,**kw):
        return m.seal(self.source,self.worker,0,self.output,
                      source_sha=sha(self.source/"manifest.json"),
                      expected_new=2,expected_selected=3,expected_reused=1,**kw)

    def test_fcp_new_photo_weighted_original_ROIs_attested_per_group(self):
        receipt=self.seal()
        self.assertEqual(receipt["actual_new_original_photo_attempts"],2)
        self.assertEqual(receipt["actual_automatic_four_hue_ROI_classifications"],1)
        self.assertEqual(receipt["original_warm_FCP_model_inits"],1)
        self.assertEqual(receipt["distinct_conditional_joint_photo_probabilities"],2)
        self.assertEqual(receipt["previously_measured_old_original_photo_IDs_reused_outside_this_group"],1)
        self.assertEqual(receipt["image_downloads_OR_ROI_reinferences_in_this_sealer"],0)
        f=self.output/"original_photo_ROIs.parquet"
        self.assertEqual(sha(f),receipt["artifact_SHA256"][f.name])
        data=pd.read_parquet(f)
        self.assertEqual(set(data["selection_probability"]),{"0.0008","0.0002"})
        self.assertEqual(set(data["morph"]),{"white","mixed_uncertain"})
        coverage=json.loads((self.output/"original_photo_source_coverage.json").read_text())
        self.assertEqual(coverage["distinct_nominal_inat_taxon_ids_four_hue"],1)
        self.assertEqual(coverage["latitude_band_by_ROI_state"]["tropical|mixed_uncertain"],1)
        self.assertFalse(coverage["individual_photo_open_flowers_independently_verified"])

    def test_bad_original_photo_terminal_SHA_fails_without_publication(self):
        self.terminal.write_text(self.terminal.read_text()+"\nmodified")
        with self.assertRaisesRegex(ValueError,"Original terminal SHA mismatch"):
            self.seal()
        self.assertFalse(self.output.exists())

    def test_no_scalar_photo_sampling_probability_allowed(self):
        with self.photos.open(newline="") as f:
            data=list(csv.DictReader(f))
        data[1]["selection_probability"]="0"
        self.csv_write(self.photos,data)
        r=json.loads(self.receipt.read_text())
        r["photos_sha256"]=sha(self.photos)
        self.receipt.write_text(json.dumps(r))
        with self.assertRaisesRegex(ValueError,"source identity changed"):
            self.seal()

    def test_frozen_FCP_model_initializes_once_and_cannot_drift(self):
        diag=json.loads(self.diagnostic.read_text())
        diag["frozen_yolo_and_onnx_initializations"]=4
        self.diagnostic.write_text(json.dumps(diag))
        with self.assertRaisesRegex(ValueError,"altered its frozen FCP model"):
            self.seal()

    def test_source_manifest_SHA_change_fails_closed(self):
        with self.assertRaisesRegex(ValueError,"manifest SHA changed"):
            m.seal(self.source,self.worker,0,self.output,source_sha="0"*64,
                   expected_new=2,expected_selected=3,expected_reused=1)

    def test_k1_idle_poll_skips_heavy_ci_but_live_originals_require_it(self):
        import yaml
        wf=(Path(__file__).resolve().parents[1]/".github/workflows"/
            "mois-fcp-K1-new-pixel-append-only.yml")
        jobs=yaml.safe_load(wf.read_text(encoding="utf-8"))["jobs"]
        self.assertNotIn("needs",jobs["inspect"])
        self.assertEqual(jobs["test"]["needs"],"inspect")
        self.assertIn("always()",jobs["test"]["if"])
        self.assertIn("github.event_name == 'pull_request'",jobs["test"]["if"])
        self.assertIn("github.event_name == 'push'",jobs["test"]["if"])
        self.assertIn("needs.inspect.outputs.ready == 'true'",jobs["test"]["if"])
        self.assertEqual(jobs["seal"]["needs"],["inspect","test"])
        self.assertIn("needs.test.result == 'success'",jobs["seal"]["if"])
        self.assertIn("needs.inspect.outputs.ready == 'true'",jobs["seal"]["if"])
        self.assertEqual(jobs["seal"]["strategy"]["max-parallel"],6)
        self.assertFalse(jobs["seal"]["strategy"]["fail-fast"])
        self.assertIn("'gh','api','--paginate'", wf.read_text(encoding="utf-8"))
        self.assertIn("releases?per_page=100", wf.read_text(encoding="utf-8"))
        gate=(Path(__file__).resolve().parents[1]/"scripts/acquisition"/
              "mois_fcp_k1_fast_completion_gate.py").read_text(encoding="utf-8")
        self.assertIn("Malformed or duplicate K1 original work group",gate)
        self.assertIn("Original completed photo ZIP expired or lacked authenticated SHA",gate)
        self.assertIn('artifact["size_in_bytes"]>0',gate)
        self.assertIn("signed_cumulative_archive_groups",wf.read_text(encoding="utf-8"))

    def test_sealing_is_immutable(self):
        self.seal()
        with self.assertRaisesRegex(ValueError,"overwrite"):
            self.seal()

if __name__=="__main__":
    unittest.main()
