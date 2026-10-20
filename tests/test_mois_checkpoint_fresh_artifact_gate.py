"""Source-authenticated schedule gate: synthetic release artifacts, not ecology."""
import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from scripts.acquisition.mois_checkpoint_fresh_artifact_gate import (
    decide_checkpoint, latest_checkpoint_tag, validate_snapshots,
)

WAVES = {
    "p0004": 39353343403,
    "p0008": 39371199127,
    "p0016": 39374811084,
    "p0032": 39401484641,
    "shadow": 39369517433,
}
def artifact(name, id, digest="sha256:abc"):
    return {"name":name, "id":id, "digest":digest, "size_in_bytes":1452}
def snapshot(name, items):
    return {"source_run_id": WAVES[name], "artifacts": items}
def prior():
    return {w: snapshot(w, [artifact("signed-"+w+"-group-0", i+10)])
            for i,w in enumerate(WAVES)}

class TestScheduledCheckpointGate(unittest.TestCase):
    def test_no_new_complete_artifact_no_release_or_dispatch(self):
        old = prior()
        self.assertEqual(decide_checkpoint(old, old, WAVES)["new_worker_artifacts"], 0)
        self.assertFalse(decide_checkpoint(old, old, WAVES)["ready"])

    def test_real_new_artifact_triggers_once_with_provenance(self):
        old = prior()
        current = json.loads(json.dumps(old))
        current["p0016"]["artifacts"].append(artifact("signed-p0016-group-1",99))
        r = decide_checkpoint(old, current, WAVES)
        self.assertTrue(r["ready"])
        self.assertEqual(r["new_worker_artifacts"],1)
        self.assertEqual(r["new_artifacts_by_wave"]["p0016"],1)
        self.assertEqual(r["new_artifacts_by_wave"]["p0008"],0)

    def test_altered_old_worker_checksum_rejected_not_relabelled(self):
        old=prior();current=json.loads(json.dumps(old))
        current["p0008"]["artifacts"][0]["digest"]="sha256:tampered"
        with self.assertRaisesRegex(ValueError,"changed|missing"):
            decide_checkpoint(old,current,WAVES)

    def test_removed_worker_receipt_rejected(self):
        old=prior();current=json.loads(json.dumps(old))
        current["p0004"]["artifacts"]=[]
        with self.assertRaisesRegex(ValueError,"changed|missing"):
            decide_checkpoint(old,current,WAVES)

    def test_invalid_source_run_and_duplicate_artifact_rejected(self):
        old=prior();current=json.loads(json.dumps(old))
        current["p0032"]["source_run_id"]=18
        with self.assertRaisesRegex(ValueError,"source run"):
            decide_checkpoint(old,current,WAVES)
        current=prior()
        current["p0032"]["artifacts"].append(current["p0032"]["artifacts"][0])
        with self.assertRaisesRegex(ValueError,"duplicate"):
            validate_snapshots(current,WAVES)

    def test_unrecognized_wave_never_silently_skipped(self):
        old=prior();current=prior()
        current["rogue"]={"source_run_id":99,"artifacts":[]}
        with self.assertRaisesRegex(ValueError,"waves"):
            decide_checkpoint(old,current,WAVES)

    def test_latest_checkpoint_is_published_and_exactly_tagged(self):
        releases=[
            {"tag_name":"mois-real-pixel-checkpoint-88-v1","draft":False,"published_at":"2026-10-20T09:00:00Z"},
            {"tag_name":"mois-real-pixel-checkpoint-99-v1","draft":True,"published_at":"2026-10-20T12:00:00Z"},
            {"tag_name":"mois-complete-93209351-original-photo-target-tiers-v1","draft":False,"published_at":"2026-10-20T20:00:00Z"},
            {"tag_name":"mois-real-pixel-checkpoint-95-v1","draft":False,"published_at":"2026-10-20T10:00:00Z"},
        ]
        self.assertEqual(latest_checkpoint_tag(releases),"mois-real-pixel-checkpoint-95-v1")
        self.assertIsNone(latest_checkpoint_tag([]))

    def test_schedule_workflow_never_reacquires_or_overwrites_original_data(self):
        wf = (Path(__file__).resolve().parents[1]/".github/workflows"/
              "mois-auto-checkpoint-when-new-pixels.yml")
        s = wf.read_text()
        self.assertIn("  schedule:", s)
        self.assertIn("gh workflow run mois-completed-pixel-checkpoint.yml", s)
        self.assertIn("ready",s)
        self.assertNotIn("mois_full_index_real_original_photo_measure work",s)
        self.assertNotIn("bootstrap_fcp_original_roi_v4",s)
        self.assertIn("cancel-in-progress: false",s)

if __name__=="__main__":
    unittest.main()
