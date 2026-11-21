#!/usr/bin/env python3
"""Prevent mismatched GitHub Actions worker artifact names across wave13 append and wave14 audit."""
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]


class TestWave14ArtifactPaths(unittest.TestCase):
    def test_wave13_append_uses_actual_wave13_worker_artifact(self):
        producer = (ROOT / ".github/workflows/mois-species-first-twelfth-append-thirteenth.yml").read_text()
        consumer = (ROOT / "scripts/analysis/mois_species_first_thirteenth_to_fourteenth.py").read_text()
        launcher = (ROOT / ".github/workflows/mois-species-first-thirteenth-append-fourteenth.yml").read_text()
        self.assertIn("name: mois-thirteenth13-real-original-worker-${{ matrix.group }}", producer)
        self.assertIn("mois-thirteenth13-real-original-worker-$g", launcher)
        self.assertIn('f"mois-thirteenth13-real-original-worker-{i//4}"', consumer)
        self.assertNotIn("mois-thirteenth11-real-original-worker", consumer)

    def test_wave14_audit_uses_actual_wave14_worker_artifact(self):
        producer = (ROOT / ".github/workflows/mois-species-first-thirteenth-append-fourteenth.yml").read_text()
        consumer = (ROOT / "scripts/analysis/mois_species_first_fourteenth_audit.py").read_text()
        self.assertIn("name: mois-fourteenth14-real-original-worker-${{ matrix.group }}", producer)
        self.assertIn("pattern: mois-fourteenth14-real-original-worker-*", producer)
        self.assertIn('f"mois-fourteenth14-real-original-worker-{i//4}"', consumer)
        self.assertNotIn("mois-fourteenth13-real-original-worker", consumer)


    def test_wave14_never_treats_skipped_success_as_pixel_source_completion(self):
        workflow = (ROOT / ".github/workflows/mois-species-first-thirteenth-append-fourteenth.yml").read_text()
        # A GHA job can be SUCCESS while all image work is skipped by a gate.
        # Only a signed first-stage source artifact is a real completed append.
        self.assertIn("real_artifact_count", workflow)
        self.assertIn("mois-thirteenth2793-verified-append-fourteenth-prioritized-original-source", workflow)
        self.assertIn("MOIS_PRIOR_WAVE14_SUCCESS_WITHOUT_SOURCE_ARTIFACT_NOT_A_COMPLETED_APPEND", workflow)
        self.assertIn("MOIS_KNOWN_GITHUB_QUEUE_ORPHAN_43095236635_ZERO_JOB_ZERO_ARTIFACT_IGNORED", workflow)
        self.assertIn('stale_jobs" = 0 ] && [ "$stale_artifacts" = 0', workflow)


if __name__ == "__main__":
    unittest.main()
