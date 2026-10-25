"""Smoke tests: real image bytes are decoded but NO network/API/model job is launched."""
import csv
import gzip
import io
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from PIL import Image
from scripts.acquisition import mois_full93m_medium_pixel_scan as scan
from scripts.acquisition import mois_full93m_every_photo_scan_manifest as plan
from scripts.acquisition.mois_onepass_full93m_target_tiers import COLUMNS


class OriginalMediumPixelScanTests(unittest.TestCase):
    def setUp(self):
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root / "tier-00-part-000.csv.gz"
        self.rows = []
        for pid, ext in [("11", "png"), ("12", "png"), ("13", "gif"), ("14", "jpg")]:
            row = {k: "" for k in COLUMNS}
            row.update(photo_id=pid, observation_uuid="O" + pid,
                       inat_taxon_id="1234", observer_id="test",
                       observed_on="2025-07-20", latitude="36.3",
                       longitude="140.1", photo_license="cc-by",
                       image_extension=ext, tier="0",
                       frozen_loader_supported="true")
            self.rows.append(row)
        with gzip.open(self.source, "wt", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=COLUMNS, lineterminator="\n")
            writer.writeheader()
            writer.writerows(self.rows)
        img = Image.new("RGB", (64, 64), (251, 251, 251))
        mem = io.BytesIO()
        img.save(mem, format="PNG")
        self.photo = mem.getvalue()
        self.output = self.root / "stage1"

    def test_source_image_pixels_and_screen_errors_conserved(self):
        seen = []
        def download(url):
            seen.append(url)
            if "/14/" in url:
                raise RuntimeError("synthetic offline HTTP error")
            return self.photo
        class Detector:
            def __call__(self, images):
                self.last_size = len(images)
                return [{"score": .92, "xyxy": [0, 0, 64, 64]}, None]
        det = Detector()
        receipt = scan.scan_source(
            self.source, self.output,
            expected_source_sha256=scan.sha_file(self.source),
            expected_count=4, part=0, tier=0, detector=det,
            download=download, batch_size=4, workers=2)
        self.assertTrue(receipt["all_original_input_photos_attempted"])
        self.assertEqual(receipt["actual_unique_source_photo_ids_attempted"], 4)
        self.assertEqual(det.last_size, 2)
        self.assertEqual(len(seen), 3)
        with gzip.open(self.output / "packet-0000000.csv.gz", "rt") as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual([r["screen_status"] for r in rows],
                         ["candidate", "negative", "error", "error"])
        self.assertEqual(rows[0]["image_pixel_sha256"], scan.sha_bytes(self.photo))
        self.assertGreater(float(rows[0]["box_proxy_white"]), 0.95)
        self.assertEqual(rows[2]["technical_reason"], "unsupported_image_extension")
        self.assertEqual(rows[1]["technical_reason"],
                         "model_screen_no_flower_candidate")
        self.assertEqual(rows[0]["source_photo_inclusion_probability"], "1")
        self.assertEqual(len({r["photo_id"] for r in rows}), 4)
        self.assertEqual(json.loads((self.output / "complete_receipt.json").read_text())[
            "screen_counts"]["error"], 2)
        # Genuine restart: no redownload, no repeated detection.
        def forbidden(url):
            self.fail("Same source image was fetched on restart")
        again = scan.scan_source(
            self.source, self.output,
            expected_source_sha256=scan.sha_file(self.source),
            expected_count=4, part=0, tier=0, detector=det,
            download=forbidden, batch_size=4, workers=2)
        self.assertEqual(again["screen_counts"], receipt["screen_counts"])

    def test_limited_pilot_never_claims_full_photo_census(self):
        doc = scan.scan_source(
            self.source, self.output,
            expected_source_sha256=scan.sha_file(self.source),
            expected_count=4, part=0, tier=0, detector=None,
            download=lambda url: self.photo, batch_size=2,
            max_packets=1, workers=1)
        self.assertFalse(doc["all_original_input_photos_attempted"])
        self.assertEqual(doc["actual_unique_source_photo_ids_attempted"], 2)
        self.assertFalse((self.output / "complete_receipt.json").exists())

    def test_provider_backpressure_aborts_without_sealing_packet(self):
        def throttled(_url):
            raise scan.ProviderBackpressure("simulated provider HTTP 429")
        with self.assertRaisesRegex(scan.ProviderBackpressure, "429"):
            scan.scan_source(
                self.source, self.output,
                expected_source_sha256=scan.sha_file(self.source),
                expected_count=4, part=0, tier=0,
                detector=None, download=throttled, batch_size=4, workers=2)
        self.assertFalse((self.output / "packet-0000000.receipt.json").exists())
        self.assertFalse((self.output / "complete_receipt.json").exists())

    def test_bounded_original_chunk_has_distinct_scope_and_strict_restart_identity(self):
        doc = scan.scan_source(
            self.source, self.output,
            expected_source_sha256=scan.sha_file(self.source),
            expected_count=4, tier=0, part=0, source_role="full93m_chunk",
            parent_source_sha256="a"*64, chunk_index=17,
            detector=None, download=lambda url: self.photo,
            batch_size=4, workers=1)
        self.assertTrue(doc["all_original_input_photos_attempted"])
        self.assertEqual(doc["source_role"], "full93m_chunk")
        self.assertEqual(doc["source_chunk_index"], 17)
        self.assertEqual(doc["source_parent_sha256"], "a"*64)
        self.assertFalse(doc["entire_original_source_part_completed"])
        with gzip.open(self.output/"packet-0000000.csv.gz", "rt") as fh:
            rows=list(csv.DictReader(fh))
        self.assertTrue(all(row["source_arm"]==scan.SOURCE_ARM for row in rows))
        self.assertTrue(all(row["source_photo_inclusion_probability"]=="1" for row in rows))
        with self.assertRaisesRegex(ValueError, "Cannot resume"):
            scan.scan_source(
                self.source, self.output,
                expected_source_sha256=scan.sha_file(self.source),
                expected_count=4, tier=0, part=0, source_role="full93m_chunk",
                parent_source_sha256="a"*64, chunk_index=18,
                detector=None, download=lambda url: self.photo, batch_size=4,
                workers=1)

    def test_chunk_cannot_masquerade_as_full_source_without_parent_sha(self):
        with self.assertRaisesRegex(ValueError, "signed original parent"):
            scan.scan_source(
                self.source, self.output,
                expected_source_sha256=scan.sha_file(self.source),
                expected_count=4, tier=0, part=0, source_role="full93m_chunk",
                detector=None, download=lambda url: self.photo,
                batch_size=4, workers=1)

    def test_bad_source_sha_and_malformed_photo_fail_closed(self):
        with self.assertRaisesRegex(ValueError, "SHA changed"):
            scan.scan_source(self.source, self.output,
                expected_source_sha256="0"*64, expected_count=4, part=0,
                tier=0, detector=None, download=lambda url: self.photo)
        self.rows[0]["photo_id"] = "not-an-id"
        with gzip.open(self.source, "wt", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=COLUMNS, lineterminator="\n")
            writer.writeheader();writer.writerows(self.rows)
        with self.assertRaisesRegex(ValueError, "Invalid original photo-ID"):
            scan.scan_source(self.source, self.output,
                expected_source_sha256=scan.sha_file(self.source),
                expected_count=4, part=0, tier=0, detector=None,
                download=lambda url: self.photo)

class CompleteUniverseManifestTests(unittest.TestCase):
    def setUp(self):
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.registry = self.root / "registry.json"
        self.group = self.root / "group.json"
        self.output = self.root / "plan"
        full = "f" * 64
        self.registry.write_text(json.dumps({
            "schema": plan.REGISTRY_SCHEMA,
            "exact_original_photo_IDs_in_all_disjoint_tiers": 4,
            "source_only_NO_photo_pixels_or_hue_measured": True,
            "native_vs_introduced_not_hard_filtered": True,
            "new_target_counts_by_tier": {"0": 2, "1": 2},
            "full_source_manifest_sha256": full,
        }))
        parts = []
        for p in range(2):
            parts.append({
                "part": p, "original_photo_IDs_in_part": 2,
                "counts_by_disjoint_tier": {"0": 1, "1": 1},
                "ranked_target_csv_sha256": {"0": "a"*64, "1": "b"*64},
            })
        self.group.write_text(json.dumps({
            "schema": plan.GROUP_SCHEMA,
            "source_snapshot_SHA256": full, "parts": parts,
        }))

    def test_exhaustive_original_photo_manifest_without_images(self):
        out = plan.build_jobs(self.registry, [self.group], self.output,
                              expected_photos=4, expected_parts=2,
                              expected_tiers=2)
        self.assertEqual(out["all_disjoint_scan_jobs_including_empty"], 4)
        self.assertEqual(out["source_original_photo_ID_count"], 4)
        self.assertEqual(out["actual_images_downloaded_or_scanned"], 0)
        self.assertEqual(len(out["signed_group_receipts_sha256"]), 1)
        self.assertFalse(out["scan_complete"])
        with (self.output / "all_original_photos_stage1_jobs.csv").open() as f:
            rows = list(csv.DictReader(f))
        self.assertEqual({(r["part"], r["tier"]) for r in rows},
                         {("0","0"), ("0","1"), ("1","0"), ("1","1")})
        self.assertEqual(sum(int(x["expected_photo_rows"]) for x in rows), 4)
        self.assertEqual({r["photo_level_sampling_probability"] for r in rows}, {"1"})

    def test_missing_partition_fails_full_universe_claim(self):
        with self.assertRaisesRegex(ValueError, "Overlapping"):
            plan.build_jobs(self.registry, [self.group, self.group],
                            self.output, expected_photos=4,
                            expected_parts=2, expected_tiers=2)
        self.assertFalse(self.output.exists())

if __name__ == "__main__":
    unittest.main()
