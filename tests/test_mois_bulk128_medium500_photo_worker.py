"""Original-photo bulk128 medium photo acquisition without iNaturalist network."""
from __future__ import annotations

import csv,gzip,hashlib,io,json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np
from PIL import Image

from scripts.acquisition.mois_bulk128_medium500_photo_worker import (
    run,build_signed_input,reconcile_pixels,sha_bytes,sha_file,SCHEMA,
    SELECTED_ROLE,WEIGHT_SHA
)
from scripts.acquisition.mois_full93m_medium_pixel_scan import LAB_SCHEMA
from scripts.acquisition.mois_queue_from_full93m_target_tiers import WORKER_COLUMNS


class FastBulk128OriginalSourceTests(unittest.TestCase):
    def setUp(self):
        self.tmp=TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.src=self.root/"source"
        (self.src/"shards").mkdir(parents=True)
        (self.src/"provenance").mkdir()
        self.output=self.root/"result"
        self.ids=[]
        per_source=[]
        shas={}
        for unit in range(4):
            photo=str(800_000_000+unit)
            obs=f"test-uuid-of-source-{unit}"
            tax=str(128*(unit+1))
            role=("reference" if unit==0 else "targeted")
            self.ids.append((photo,obs,tax,role))
            row={k:"" for k in WORKER_COLUMNS}
            row.update(
                measurement_id="MOIS-"+str(unit),photo_id=photo,inat_taxon_id=tax,
                species="Taxon "+str(unit),observation_id=obs,observer_id="operator",
                observed_on="2026-04-15",latitude="-15",longitude="132",
                positional_accuracy="20",photo_license="cc-by",
                photo_url_large=f"https://inaturalist-open-data.s3.amazonaws.com/photos/{photo}/large.jpg",
                source_arm=("mois_bulk128_tier10_photo_reference_v1"
                            if unit==0 else "mois_bulk128_tier10_taxon_region_season_gap_v1"),
                selection_probability=".0002",source_snapshot_date="2026-09-27",
                shard=str(unit)
            )
            path=self.src/"shards"/f"shard-{unit:03d}.csv"
            with path.open("w",newline="") as fh:
                writer=csv.DictWriter(fh,fieldnames=WORKER_COLUMNS)
                writer.writeheader();writer.writerow(row)
            shas[str(unit)]=sha_file(path)
            per_source.append({
                "photo_id":photo,"observation_id":obs,"inat_taxon_id":tax,
                "acquisition_design_cohort":role,
                "original_source_priority":("P5_existing_10plus_strict9"
                    if role=="reference" else "P2_source_has_photos_no_colour_10plus"),
                "latitude_region":"S_tropical","recorded_local_season":"autumn",
                "previous_photo_ID_and_observation_overlap":False,
            })
        file=self.src/"provenance"/"part-000.csv"
        with file.open("w",newline="") as fh:
            w=csv.DictWriter(fh,fieldnames=list(per_source[0]))
            w.writeheader();w.writerows(per_source)
        proof={
            "selected_new_original_photo_IDs":4,
            "source_photo_provenance_SHA256":sha_file(file),
            "source_photo_worker_unit_SHA256":shas,
            "selected_photo_attempts_by_design":{"reference":1,"targeted":3},
        }
        proof_path=self.src/"part-000.receipt.json"
        proof_path.write_text(json.dumps(proof))
        manifest={
            "schema":"mois_full_source_observation_sampling_cohort_v1",
            "sampling_schema":"mois_original93m_tier10_bulk128_probability_reference_plus_gap_v1",
            "source_128_compute_shards":128,"source_tier":10,
            "new_photo_attempts_planned":131072,
            "selected_reference_photo_attempts":32768,
            "selected_targeted_photo_attempts":98304,
            "P5_targeted_photo_attempts":0,
            "all_upstream_128_source_taxon_shards_verified":True,
            "per_128_shard_receipts_sha256":{"0":sha_file(proof_path)},
            "photo_shard_sha256":shas,
        }
        (self.src/"manifest.json").write_text(json.dumps(manifest))
        image=Image.new("RGB",(96,96),(245,135,55))
        buffer=io.BytesIO();image.save(buffer,format="JPEG")
        self.photo_bytes=buffer.getvalue()

    def test_reference_and_targeted_original_source_ID_mapping(self):
        output=self.root/"source_out"
        proof=build_signed_input(self.src,0,output)
        self.assertEqual(proof["original_photo_ID_count"],4)
        self.assertEqual(set(proof["source_roles"].values()),{"reference","targeted"})
        with gzip.open(proof["source_file"],"rt",newline="") as f:
            rows=list(csv.DictReader(f))
        self.assertEqual(len(rows),4)
        self.assertEqual(rows[0]["tier"],"10")
        self.assertTrue(all(x["frozen_loader_supported"]=="true" for x in rows))
        self.assertEqual(len({x["observation_uuid"] for x in rows}),4)

    def test_genuine_medpixels_and_lab_proxy_without_strict9_fake_photos(self):
        class DummyFlowerBox:
            model_identity="item"
            model_sha=WEIGHT_SHA
            conf=.02
            def __call__(self, images):
                return [{"score":.89,"xyxy":[4.,3.,90.,95.]} for _ in images]
        report=run(self.src,self.output,0,self.root/"notloaded.pt",
            detector_override=DummyFlowerBox(),
            download_override=lambda url:self.photo_bytes,
            rps=.5,batch=256)
        self.assertEqual(report["schema"],SCHEMA)
        self.assertEqual(report["expected_original_photos"],4)
        self.assertEqual(report["actually_medium_photo_IDs_attempted"],4)
        self.assertEqual(report["real_medium_photo_pixels_decoded"],4)
        self.assertEqual(report["by_design"]["reference"]["photo_attempts"],1)
        self.assertEqual(report["by_design"]["targeted"]["photo_attempts"],3)
        self.assertEqual(report["by_design"]["targeted"]["Lab_modal_descriptor"],3)
        self.assertEqual(report["new_targeted_original_priority"],
                         {"P2_source_has_photos_no_colour_10plus":3})
        self.assertFalse(report["original_FCP_four_hue_full_image_ROI_inferred"])
        self.assertEqual(report["real_strict9_flower_pixel_vectors_generated"],0)
        packets=list((self.output/"medium500").glob("packet-*.csv.gz"))
        self.assertEqual(len(packets),1)
        with gzip.open(packets[0],"rt",newline="") as f:
            row=list(csv.DictReader(f))[0]
        self.assertEqual(row["source_photo_inclusion_probability"],"")
        self.assertEqual(row["source_arm"],
                         "bulk128_tier10_source_selected_reference_and_gap_SEPARATE")
        self.assertEqual(row["screen_status"],"candidate")
        self.assertTrue(bool(row["floral_box_main_L_star"]))

if __name__=="__main__":
    unittest.main()
