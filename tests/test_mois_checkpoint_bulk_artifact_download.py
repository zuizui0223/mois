"""The checkpoint must download exact frozen worker ZIP names in a single CLI call per wave.

The tests inspect the actual inline GitHub Actions downloader. This never invokes
photo image acquisition, YOLO/EfficientSAM inference, or real cloud downloads.
"""
import ast
from pathlib import Path
import re
import textwrap
import unittest

WORKFLOW=(Path(__file__).resolve().parents[1]/
          ".github/workflows/mois-completed-pixel-checkpoint.yml")


class BulkDownloadContract(unittest.TestCase):
    def setUp(self):
        raw=WORKFLOW.read_text(encoding="utf-8")
        blocks=re.findall(r"python - <<'PY'\n(.*?)\n\s+PY(?:\n|$)",raw,re.S)
        self.source=next(textwrap.dedent(x) for x in blocks
                         if "MOIS_STORED_PIXEL_ARTIFACTS_RECOVERED" in x)
        self.code=ast.parse(self.source)

    def test_all_sealed_names_passed_in_a_single_exact_name_download(self):
        self.assertIn("for a in selected:",self.source)
        self.assertIn("download_args.extend(['--name',a['name']])",self.source)
        self.assertIn("subprocess.run(download_args,check=True)",self.source)
        self.assertNotIn("subprocess.run(['gh','run','download',str(run),'--repo','zuizui0223/mois',\\n                      '--name',a['name']",self.source)
        self.assertNotIn("'--pattern',prefix+'*'",self.source)
        self.assertNotIn("'--pattern',prefix",self.source)

    def test_multiple_exact_archives_extract_to_their_individual_names(self):
        # CLI documents separate subdirectories for multi-artifact downloads,
        # and direct unpacking for a single artifact.
        for phrase in (
            "expected_names={a['name'] for a in selected}",
            "actual_names={p.name for p in (root/'workers').iterdir() if p.is_dir()}",
            "actual_names==expected_names",
            "if len(selected)>1:",
            "if selected:",
        ):
            self.assertIn(phrase,self.source)

    def test_snapshot_freeze_precedes_download_and_requires_sha_and_no_duplicates(self):
        self.assertLess(self.source.index("selected=[a for a in artifacts"),
                        self.source.index("download_args="))
        self.assertLess(self.source.index("'artifact_snapshot.json'"),
                        self.source.index("download_args="))
        self.assertIn("len({a['name'] for a in selected})==len(selected)",self.source)
        self.assertIn("all(not a['expired'] for a in selected)",self.source)
        self.assertIn("'digest':a.get('digest')",self.source)
        self.assertIn("repos/zuizui0223/mois/actions/runs/{run}/artifacts",self.source)

    def test_data_only_pipeline_remains_source_authenticated(self):
        s=WORKFLOW.read_text(encoding="utf-8")
        self.assertIn("mois_completed_pixel_checkpoint",s)
        self.assertNotIn("mois_full_index_real_original_photo_measure work",s)
        self.assertNotIn("bootstrap_fcp_original_roi_v4",s)
        self.assertIn("cancel-in-progress: false",s)


if __name__=="__main__":
    unittest.main()
