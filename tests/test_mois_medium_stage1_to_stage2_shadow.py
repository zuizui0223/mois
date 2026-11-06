"""No-network source-exact packet integrity and Stage1-to-Stage2 SHADOW test."""
import csv
import gzip
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from scripts.acquisition.mois_medium_stage1_to_stage2_shadow import bridge
from scripts.acquisition.mois_full93m_medium_pixel_scan import SCHEMA, COLUMNS_OUT
from scripts.acquisition.mois_93m_stage2_probability_queue import sha256

class Stage1Stage2OriginalPhotoTests(unittest.TestCase):
    def setUp(self):
        self.tmp=TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.source=self.root/"original_pixels"
        self.source.mkdir()
        self.out=self.root/"bridge"
        self.rows=[
            ("14","O14","candidate","N1|JJA"),
            ("11","O11","negative","N1|JJA"),
            ("13","O13","error","N1|DJF"),
            ("12","O12","uncertain","S0|DJF"),
        ]
        self.save()

    def save(self):
        for path in self.source.iterdir():path.unlink()
        states={}
        for n, sub in enumerate((self.rows[:2],self.rows[2:])):
            file=self.source/f"packet-{n:07d}.csv.gz"
            with gzip.open(file,"wt",encoding="utf-8",newline="") as stream:
                writer=csv.DictWriter(stream,fieldnames=COLUMNS_OUT,lineterminator="\n")
                writer.writeheader()
                for pid,obs,status,stratum in sub:
                    row={k:"" for k in COLUMNS_OUT}
                    row.update(photo_id=pid,observation_id=obs,screen_status=status,
                        audit_stratum=stratum,
                        source_arm="full93m_observation_and_photo_complete_source",
                        source_photo_inclusion_probability="1")
                    writer.writerow(row)
                    states[status]=states.get(status,0)+1
            proof={
                "schema":SCHEMA,"source_role":"full93m",
                "source_gzip_sha256":"f"*64,
                "detector_identity":"frozen-yolo11-item",
                "output_sha256":sha256(file),
                "photo_ids":len(sub),
            }
            (self.source/f"packet-{n:07d}.receipt.json").write_text(json.dumps(proof))
        (self.source/"complete_receipt.json").write_text(json.dumps({
            "schema":SCHEMA,"source_role":"full93m",
            "all_original_input_photos_attempted":True,
            "calibration_subset_not_global_census":False,
            "expected_photo_ids":len(self.rows),
            "actual_unique_source_photo_ids_attempted":len(self.rows),
            "packets_sealed":2,"screen_counts":states,
            "source_gzip_sha256":"f"*64,"detector_identity":"frozen-yolo11-item",
        }))

    def test_signed_full_source_4_image_shard_stays_all_four_in_stage2(self):
        r=bridge(self.source,self.out,expected_count=4)
        self.assertEqual(r["source_original_photo_count"],4)
        self.assertEqual(r["stage2_full_ROI_candidates_including_screen_negatives"],4)
        self.assertEqual(r["screen_negatives_excluded"],0)
        decisions=self.out/"stage2_queue"/"all_photo_screen_decisions.csv"
        with decisions.open(newline="") as f:
            x=list(csv.DictReader(f))
        self.assertEqual([v["photo_id"] for v in x],["11","12","13","14"])
        self.assertEqual([v["full_roi_selected"] for v in x],["true"]*4)
        self.assertEqual([v["stage2_selection_probability"] for v in x],["1"]*4)

    def test_duplicate_original_photo_id_is_rejected(self):
        self.rows[-1]=("14","O12","uncertain","S0|DJF")
        self.save()
        with self.assertRaisesRegex(ValueError,"Duplicate photo ID"):
            bridge(self.source,self.out,expected_count=4)

    def test_fake_completed_global_source_rejected(self):
        proof=json.loads((self.source/"complete_receipt.json").read_text())
        proof["all_original_input_photos_attempted"]=False
        (self.source/"complete_receipt.json").write_text(json.dumps(proof))
        with self.assertRaisesRegex(ValueError,"Only completed real full93m"):
            bridge(self.source,self.out,expected_count=4)

if __name__=="__main__":
    unittest.main()
