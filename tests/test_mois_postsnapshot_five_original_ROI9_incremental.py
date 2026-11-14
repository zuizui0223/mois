"""Synthetic source-SHA and exact-photo tests for five signed post-snapshot ROI9 deltas."""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis import mois_postsnapshot_five_original_ROI9_incremental as m
from scripts.analysis.mois_species_first_palette9 import PALETTE,sha256


class FiveTruePostSnapshotOriginalTests(unittest.TestCase):
    def fixture(self,root):
        roots={}
        count=0
        for offset,(n,positive,clear) in m.WAVES.items():
            d=Path(root)/str(offset)
            d.mkdir()
            entries=[]
            for i in range(n):
                tid=str(2_000_000+count)
                pid=str(750_000_000+count)
                row={
                    "inat_taxon_id":tid,"photo_id":pid,
                    "observation_id":f"12345678-1234-1234-1234-{count:012d}",
                    "photo_license":"cc-by",
                    "photo_url_large":f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.jpg",
                    "source_observation_created_at":"2026-11-01T02:03:04+00:00",
                    "morph":"white" if i<positive else "mixed_uncertain",
                    "strict9_original_photo_numeric":"True" if i<positive else "False",
                    "dominant_hue9_label":"white" if i<clear else "",
                    "image_sha256":f"{count+1:064x}",
                }
                for name in PALETTE:
                    if i>=positive:row["flower_fraction_"+name]=""
                    elif i<clear:row["flower_fraction_"+name]="1" if name=="white" else "0"
                    else:row["flower_fraction_"+name]="0.50" if name in ("white","yellow") else "0"
                entries.append(row)
                count+=1
            csv=d/"mois_new_inat_post20261028_actual_original_nine_roi.csv.gz"
            pd.DataFrame(entries).to_csv(csv,index=False,compression="gzip")
            report={
                "schema":m.AUDIT_SCHEMA,
                "original_licensed_photo_attempts":n,
                "actual_genuine_strict9_numeric_success":positive,
                "actual_clear_photo_nine_palette_labels":clear,
                "actual_technical_unknowns":n-positive,
                "photo_observation_taxon_collision_count":0,
                "new_first_strict9_nominal_taxa":positive,
                "full_historical_raw_image_SHA_collision_check_possible":False,
                "official_iNat_42111_original_frame_auto_appended":False,
                "audited_derivative_sha256":sha256(csv),
                "source_photo_manifest_sha256":f"{offset+3:064x}",
                "new_photo_image_SHA256":{r["photo_id"]:r["image_sha256"] for r in entries},
            }
            (d/"receipt.json").write_text(json.dumps(report)+"\n")
            roots[offset]=d
        return roots

    def test_true_three_measured_and_two_hues_across_eleven_new_photos(self):
        with tempfile.TemporaryDirectory() as td:
            data,cert=m.load_five_audits(self.fixture(td))
            self.assertEqual(len(data),11)
            self.assertEqual(int(data.strict9.sum()),3)
            self.assertEqual(int(data.dominant_label9.eq("white").sum()),2)
            self.assertEqual(data.photo_id.nunique(),11)
            self.assertEqual(data.observation_id.nunique(),11)
            self.assertEqual(data.inat_taxon_id.nunique(),11)
            self.assertEqual(int(data["strict9_white"].notna().sum()),3)
            self.assertEqual(len(cert),5)

    def test_photo_sha_or_new_taxon_reuse_cannot_create_false_species(self):
        with tempfile.TemporaryDirectory() as td:
            roots=self.fixture(td)
            d=roots[128]
            csv=d/"mois_new_inat_post20261028_actual_original_nine_roi.csv.gz"
            with csv.open("ab") as handle:handle.write(b"tamper")
            with self.assertRaisesRegex(ValueError,"SHA or numeric"):
                m.load_five_audits(roots)

    def test_source_requires_all_five_signed_original_photo_receipts(self):
        with tempfile.TemporaryDirectory() as td:
            roots=self.fixture(td)
            roots.pop(192)
            with self.assertRaisesRegex(ValueError,"five independent"):
                m.load_five_audits(roots)


if __name__=="__main__":
    unittest.main()
