"""Synthetic full-FCP source/photo terminal tests for first 3141 actual original ROIs."""
from __future__ import annotations

import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import pandas as pd
from scripts.acquisition import mois_fcp_u100_collect_new_first_photo_ROIs as m
from scripts.acquisition.mois_completed_pixel_checkpoint import sha
from scripts.analysis.mois_reassemble_original_new4096_measurements import METHOD


class OriginalFCPPrefixOnePhotoTests(unittest.TestCase):
    def setUp(self):
        self.tmp=TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.source=self.root/"source"
        self.worker=self.root/"workers"/"shard_000"
        self.output=self.root/"finished"
        (self.source/"shards").mkdir(parents=True)
        (self.worker/"sealed").mkdir(parents=True)
        self.prefix=self.root/"prefix.csv"
        self.prefixes=[]
        for index,pid in enumerate(("100","101","102")):
            n=100+index
            mphotos=2+index
            self.prefixes.append(dict(
                photo_id=pid,observation_id=f"UUID-{index}",
                inat_taxon_id=str(index+1),
                full_93m_distinct_original_observations_of_taxon=str(n),
                photos_in_original_observation=str(mphotos),
                depth1_joint_original_photo_inclusion_probability=format(1/n/mphotos,".17g")
            ))
        self.csv_write(self.prefix,self.prefixes)
        self.rows=[]
        for index,src in enumerate(self.prefixes[:2]):
            self.rows.append(dict(
                measurement_id=f"PREFIX_{index}",photo_id=src["photo_id"],
                observation_id=src["observation_id"],
                inat_taxon_id=src["inat_taxon_id"],
                source_arm=m.PREFIX_ARM,
                selection_probability=src["depth1_joint_original_photo_inclusion_probability"],
                observed_on="2025-06-01",latitude="40.1",photo_license="cc-by",
                image_filename=f"PREFIX_{index}.jpg",shard="0",
            ))
        self.work_csv=self.source/"shards/shard-000.csv"
        self.csv_write(self.work_csv,self.rows)
        self.manifest=self.source/"manifest.json"
        self.manifest.write_text(json.dumps(dict(
            schema="mois_full_source_observation_sampling_cohort_v1",
            full_frame_photo_ids=93_209_351,
            source_arm=m.PREFIX_ARM,
            selected_FCP_U100_eligible_source_taxa=3,
            selected_original_photo_observations=3,
            previously_measured_original_photo_ids_reused=1,
            inflight_K1_photo_IDs_reserved_without_redownload=0,
            new_photo_attempts_planned=2,
            technical_pilot_photo_IDs_not_refetched=0,
            unsupported_original_photo_format_selected_as_unknown=0,
            colour_measured=False,
            photo_already_reserved_by_K1_never_fetched_twice=True,
            FCP_top100_photo_prefix_inclusion_probability_corrected_to_1_over_N_M=True,
            work_units=1,worker_ids=[0],group_ids=[0],
            photo_shard_sha256={"0":sha(self.work_csv)}
        )))
        self.terminal=self.worker/"sealed/terminal_p000_b0000000.csv"
        ends=[
            dict(measurement_id=self.rows[0]["measurement_id"],
                 morph="white",measurement_status="classified_four_state_morph",
                 image_sha256="a"*64),
            dict(measurement_id=self.rows[1]["measurement_id"],
                 morph="mixed_uncertain",
                 measurement_status="not_evaluable_roi_or_flip_gate",
                 image_sha256="b"*64)
        ]
        self.csv_write(self.terminal,ends)
        self.photos=self.worker/"photos.csv"
        measured=[{**row,**end,
                   "classified_four_hue":str(end["morph"]=="white").lower()}
                  for row,end in zip(self.rows,ends)]
        self.csv_write(self.photos,measured)
        self.receipt=self.worker/"receipt.json"
        self.receipt.write_text(json.dumps(dict(
            schema="mois_full_original93m_real_RGB_photo_unit_v1",
            shard=0,attempted=2,downloaded=2,classifiable=1,
            acquisition_failed=0,unclassifiable_acquired=1,
            photo_colours_measured=True,raw_licensed_image_bytes_persisted=False,
            source_sha256=sha(self.work_csv),
            method=METHOD,photos_sha256=sha(self.photos),
            terminal_sha256=sha(self.terminal)
        )))

    def csv_write(self,p,rows):
        with p.open("w",newline="",encoding="utf-8") as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n")
            w.writeheader();w.writerows(rows)

    def run_collect(self):
        return m.collect(self.source,self.root/"workers",self.prefix,self.output,
                         expected_manifest_sha=sha(self.manifest),
                         prefix_sha=sha(self.prefix),
                         expected_new=2,expected_reuse=1,
                         expected_inflight=0,expected_eligible=3)

    def test_exact_one_per_species_original_photo_probability_and_unknowns(self):
        result=self.run_collect()
        self.assertEqual(result["newly_true_image_measured_original_photo_attempts"],2)
        self.assertEqual(result["newly_true_image_FCP_four_hue_classified"],1)
        self.assertEqual(result["newly_true_image_unclassifiable_acquired"],1)
        self.assertTrue(result["individual_photo_design_probabilities_1_over_N_over_M_preserved"])
        self.assertFalse(result["globally_unbiased_white_flora_fraction_estimated"])
        photo=self.output/"new_FCP_depth1_original_photo_ROIs.parquet"
        self.assertEqual(sha(photo),result["actual_new_image_parquet_SHA256"])
        df=pd.read_parquet(photo)
        self.assertEqual(set(df.morph),{"white","mixed_uncertain"})
        self.assertEqual(len(df),2)

    def test_wrong_original_FCP_photo_sampling_probability_is_rejected(self):
        with self.photos.open(newline="") as stream:
            d=list(csv.DictReader(stream))
        d[0]["selection_probability"]="0.1"
        self.csv_write(self.photos,d)
        rec=json.loads(self.receipt.read_text())
        rec["photos_sha256"]=sha(self.photos)
        self.receipt.write_text(json.dumps(rec))
        with self.assertRaisesRegex(ValueError,"source identity changed"):
            self.run_collect()
        self.assertFalse(self.output.exists())

    def test_two_photo_IDs_for_same_fcp_taxon_fail(self):
        with self.photos.open(newline="") as stream:
            d=list(csv.DictReader(stream))
        d[1]["inat_taxon_id"]="1"
        self.csv_write(self.photos,d)
        rec=json.loads(self.receipt.read_text())
        rec["photos_sha256"]=sha(self.photos)
        self.receipt.write_text(json.dumps(rec))
        with self.assertRaisesRegex(ValueError,"source identity changed"):
            self.run_collect()

    def test_original_image_terminal_modified_refused(self):
        self.terminal.write_text(self.terminal.read_text()+"\nforged")
        with self.assertRaisesRegex(ValueError,"Original terminal SHA mismatch"):
            self.run_collect()

    def test_missing_real_original_photo_work_unit_refused(self):
        self.receipt.unlink()
        with self.assertRaisesRegex(ValueError,"have not all been sealed"):
            self.run_collect()

    def test_fcp_source_manifest_tampering_fails(self):
        self.manifest.write_text("{}")
        with self.assertRaisesRegex(ValueError,"FCP first-photo source not the pinned"):
            self.run_collect()

    def test_frozen_photo_prefix_probability_1_over_N_M_gate(self):
        with self.prefix.open(newline="") as stream:
            d=list(csv.DictReader(stream))
        d[0]["depth1_joint_original_photo_inclusion_probability"]="0.1"
        self.csv_write(self.prefix,d)
        with self.assertRaisesRegex(ValueError,"1/N/M changed"):
            m.check_prefix_probabilities(self.prefix,prefix_sha=sha(self.prefix),
                                         expected_eligible=3)

    def test_four_parallel_new_photo_groups_only_after_all_original_K1_complete(self):
        # The FCP depth-1 3141-photo source has already excluded all K1 IDs.
        # Four groups may run concurrently only once the parent 8-worker wave
        # finished, preventing double S3 acquisition and accidental CPU overlap.
        workflow = (Path(__file__).resolve().parents[1] / ".github/workflows" /
                    "mois-fcp-u100-depth1-new3141-original-ROI.yml")
        text = workflow.read_text(encoding="utf-8")
        self.assertIn("K1_RUN: '40066261878'", text)
        self.assertIn('if [ "$completed" -lt 75 ]; then', text)
        self.assertIn('test "$completed" -eq 75', text)
        self.assertIn("max-parallel: 4", text)
        self.assertIn("groups=[0,1,2,3]", text)
        self.assertIn("cancel-in-progress: false", text)
        self.assertIn("new_photo_attempts_planned", text)
        self.assertIn("inflight_K1_photo_IDs_reserved_without_redownload", text)

    def test_signed_K1_group_release_gate_requires_canonical_complete_75(self):
        # Execute the actual embedded source gate, rather than only checking
        # text. The previous GitHub releases page contains >100 other tags.
        import re
        import subprocess
        import sys
        import yaml
        workflow=(Path(__file__).resolve().parents[1] / ".github/workflows" /
                  "mois-fcp-u100-depth1-new3141-original-ROI.yml")
        jobs=yaml.safe_load(workflow.read_text(encoding="utf-8"))["jobs"]
        gate=next(step["run"] for step in jobs["original-source-only"]["steps"]
                  if step.get("id")=="gate")
        self.assertIn("gh api --paginate",gate)
        self.assertIn("select(.draft == false)",gate)
        self.assertIn('if [ "$completed" -lt 75 ]; then',gate)
        self.assertIn('if [ "$sealed" -lt 75 ]; then',gate)
        self.assertIn('test "$sealed" -eq 75',gate)
        self.assertIn("original photo terminal groups have been SHA-published",gate)
        snippet=re.search(r"sealed=\$\(python - <<'PY'\n(.*?)\nPY\n\)",gate,re.S)
        self.assertIsNotNone(snippet)
        def run_group_gate(groups):
            with TemporaryDirectory() as directory:
                names=[f"mois-full93m-balanced-onephoto-K1-original-RGB-append-group-{g}-v1"
                       for g in groups]
                (Path(directory)/"k1-immutable-published-tags.txt").write_text(
                    "\n".join(names)+"\n",encoding="utf-8")
                return subprocess.run([sys.executable,"-c",snippet.group(1)],
                    cwd=directory,text=True,capture_output=True)
        for expected in (0,1,74,75):
            got=run_group_gate(range(expected))
            self.assertEqual(got.returncode,0,got.stderr)
            self.assertEqual(got.stdout.strip(),str(expected))
        wrong=run_group_gate(list(range(74))+[75])
        self.assertNotEqual(wrong.returncode,0,"Invalid K1 group 75 must fail")
        repeated=run_group_gate([0,0])
        self.assertNotEqual(repeated.returncode,0,"Duplicate signed K1 group must fail")

    def test_post_K1_sealed_release_triggers_FCP_without_idle_heavy_CI(self):
        # Final K1 model completion is not final publication; a source-linked
        # workflow_run from the last signed group must also wake FCP U100.
        # Neither unfinished nor failed K1 runs may start original photo work.
        import yaml
        wf=Path(__file__).resolve().parents[1] / ".github/workflows" / (
            "mois-fcp-u100-depth1-new3141-original-ROI.yml")
        d=yaml.safe_load(wf.read_text(encoding="utf-8"))
        on=d.get("on",d.get(True))
        triggers=on["workflow_run"]["workflows"]
        self.assertIn("MOIS new 93m balanced onephoto ROI-v4 warm 8-runner acquisition",triggers)
        self.assertIn(
            "MOIS FCP K1 immediately append completed source-verified original ROI groups",
            triggers)
        jobs=d["jobs"]
        self.assertNotIn("needs",jobs["original-source-only"])
        gate_condition=jobs["original-source-only"]["if"]
        self.assertIn("github.event.workflow_run.conclusion == 'success'",gate_condition)
        self.assertIn("github.event.workflow_run.head_branch == 'main'",gate_condition)
        self.assertIn("head_repository.full_name == 'zuizui0223/mois'",gate_condition)
        tests=jobs["frozen-source-probability-tests"]
        self.assertEqual(tests["needs"],"original-source-only")
        self.assertIn("always()",tests["if"])
        self.assertIn("github.event_name == 'pull_request'",tests["if"])
        self.assertIn("github.event_name == 'push'",tests["if"])
        self.assertIn("needs.original-source-only.outputs.ready == 'true'",tests["if"])
        infer=jobs["true-original-image-ROI"]
        self.assertEqual(infer["needs"],["original-source-only","frozen-source-probability-tests"])
        self.assertIn("needs.frozen-source-probability-tests.result == 'success'",infer["if"])
        seal=jobs["seal-new-photo-FCP-taxon-coverage"]
        self.assertEqual(seal["needs"],[
            "original-source-only","frozen-source-probability-tests",
            "true-original-image-ROI"])
        self.assertIn("needs.true-original-image-ROI.result == 'success'",seal["if"])
        self.assertIn("needs.frozen-source-probability-tests.result == 'success'",seal["if"])

    def test_expensive_FCP_image_workers_require_independent_3141_source_preflight(self):
        import yaml
        workflow=Path(__file__).resolve().parents[1]/".github/workflows"/(
            "mois-fcp-u100-depth1-new3141-original-ROI.yml")
        jobs=yaml.safe_load(workflow.read_text(encoding="utf-8"))["jobs"]
        source_gate=next(step["run"] for step in jobs["original-source-only"]["steps"]
                         if step.get("id")=="gate")
        self.assertIn("PREFLIGHT_TAG: mois-FCP-U100-20242-taxa-only3141-photo-K1-disjoint-preflight-v1",
                      workflow.read_text(encoding="utf-8"))
        self.assertIn("7920cdefb4dd1f5a35def3cdf3fb3c1d19e86e87d4f8db8dab08645e37542eef",
                      source_gate)
        self.assertIn("MOIS_FCP_U100_3141_K1_DISJOINT_SOURCE_SAFETY_PREFLIGHT_SHA_REUSED",
                      source_gate)
        self.assertIn("never_before_measured_FCP_U100_original_photo_ROI_target_IDs",source_gate)
        self.assertIn("new_original_photos_overlap_with_K1",source_gate)
        self.assertIn("all_photo_inclusion_probabilities_individually_source_checked",source_gate)
        self.assertTrue(source_gate.index('test "$sealed" -eq 75') <
                        source_gate.index('gh release download "$PREFLIGHT_TAG"') <
                        source_gate.index('echo ready=true >> "$GITHUB_OUTPUT"'))
        self.assertIn("needs.frozen-source-probability-tests.result == 'success'",
                      jobs["true-original-image-ROI"]["if"])

    def test_already_published_outcome_is_immutable(self):
        self.run_collect()
        with self.assertRaisesRegex(ValueError,"must not be overwritten"):
            self.run_collect()


if __name__=="__main__":
    unittest.main()
