"""Strict original ROI-v4 source/colour reconciliation and no repeat inference."""
from pathlib import Path
from tempfile import TemporaryDirectory
import csv,json,unittest
from scripts.acquisition import mois_full_index_real_original_photo_measure as m
import pandas as pd

class TestFullSourceRealRGB(unittest.TestCase):
    def setUp(self):
        self.ctx=TemporaryDirectory();self.root=Path(self.ctx.name)
        self.plan=self.root/"source";(self.plan/"shards").mkdir(parents=True)
        self.data=[]
        for i in range(5):
            n=i//3;mid="MOIS-"+str(i).zfill(24)
            self.data.append(dict(measurement_id=mid,photo_id=str(10000+i),
                inat_taxon_id=str(i+10),species="Sample species"+str(i),
                observer_id=str(i+1),observation_id="obs"+str(i),
                observed_on="2025-06-21",latitude="37.1",longitude="136.8",
                positional_accuracy="100",photo_license="cc-by",
                photo_url_large=f"https://inaturalist-open-data.s3.amazonaws.com/photos/{10000+i}/large.jpg",
                image_filename=mid+".jpg",
                source_arm="probability_baseline",selection_probability="0.0002",
                source_snapshot_date="2026-09-27",shard=str(n)))
        files={}
        for n in range(2):
            filename=self.plan/"shards"/f"shard-{n:03d}.csv"
            with filename.open("w",newline="") as out:
                w=csv.DictWriter(out,fieldnames=list(self.data[0]))
                w.writeheader();w.writerows(x for x in self.data if int(x["shard"])==n)
            files[str(n)]=m.sha(filename)
        self.manifest={"schema":m.SCHEMA,"work_units":2,"worker_ids":[0,1],
            "photo_shard_sha256":files,"new_photo_attempts_planned":5,
            "selection_probability":0.0002,"full_frame_photo_ids":93209351,
            "previously_measured_photo_ids_skipped":0,"colour_measured":False}
        (self.plan/"manifest.json").write_text(json.dumps(self.manifest))
        self.calls=[]
    def tearDown(self):self.ctx.cleanup()
    def fake_runner(self,path,part,batch,size,dest,*,indexed_source,**kwargs):
        photos,urls,n=indexed_source.get_batch(part,batch,size)
        self.assertEqual((len(photos),len(urls),n),(size,size,1))
        self.calls.append(part)
        dest.mkdir(parents=True,exist_ok=True)
        file=dest/f"terminal_p{part:03d}_b0000000.csv"
        pd.DataFrame({"measurement_id":photos.measurement_id,
            "morph":["white"]*size,
            "measurement_status":["classified_four_state_morph"]*size,
            "image_sha256":["f"*64]*size}).to_csv(file,index=False)
        return {"status":"COMPLETE_WORK_UNIT",
            "n_source_photo_ids":size,"n_terminal_photo_ids":size,
            "number_batches_for_this_partition":1,
            "image_bytes_persisted":False,
            "reused_historical_measurement_csv":False,
            "terminal_csv_sha256":m.sha(file),
            "n_classifiable":size,"n_acquisition_failed":0,
            "method":{"model":"frozen-test"}}
    def test_all_original_pixels_once_and_signed_completion(self):
        folder=self.root/"measured"
        out=m.measure_group(self.plan,0,self.root/"runtime",self.root/"sam",
                             folder,runner=self.fake_runner)
        self.assertEqual(self.calls,[0,1])
        self.assertEqual(sum(r["attempted"] for r in out),5)
        def forbidden(*args,**kw):raise AssertionError("Old measured pixels were downloaded twice")
        again=m.measure_group(self.plan,0,self.root/"runtime",self.root/"sam",
                              folder,runner=forbidden)
        self.assertEqual(sum(r["attempted"] for r in again),5)
        receipt=m.collect(self.plan,folder,self.root/"final")
        self.assertEqual((receipt["attempted"],receipt["classifiable"]),(5,5))
        self.assertTrue(pd.read_parquet(self.root/"final"/"measured_photos.parquet").photo_id.is_unique)
    def test_tampered_worker_cannot_be_resumed(self):
        folder=self.root/"measured"
        m.measure_group(self.plan,0,self.root/"runtime",self.root/"sam",
                        folder,runner=self.fake_runner)
        (folder/"shard_000"/"photos.csv").write_text("tampered")
        with self.assertRaisesRegex(ValueError,"corrupted"):
            m.measure_group(self.plan,0,self.root/"runtime",self.root/"sam",
                            folder,runner=self.fake_runner)
    def test_non_frozen_image_extensions_become_permanent_technical_unknown(self):
        called=[]
        class FakePermanentOriginalPhotoError(ValueError):
            pass
        def original_loader(url):
            called.append(url)
            return b"authentic-photo-bytes"
        for ext in ("gif","webp","unsupported","avif","heic",""):
            url=f"https://inaturalist-open-data.s3.amazonaws.com/photos/1234/large.{ext}"
            self.assertFalse(m.original_source_loader_supported(url))
            with self.assertRaisesRegex(FakePermanentOriginalPhotoError,
                                        "technical unknown"):
                m.original_source_image_loader(url, original_loader,
                                               FakePermanentOriginalPhotoError)
        self.assertEqual(called,[])
        for ext in ("jpg","jpeg","png","JPG"):
            url=f"https://inaturalist-open-data.s3.amazonaws.com/photos/1234/large.{ext}"
            self.assertTrue(m.original_source_loader_supported(url))
            self.assertEqual(m.original_source_image_loader(
                url,original_loader,FakePermanentOriginalPhotoError),
                b"authentic-photo-bytes")
        self.assertEqual(len(called),4)

    def test_source_original_URL_cannot_escape_pinned_provider(self):
        for url in (
            "https://evil.example/photos/1234/large.jpg",
            "http://inaturalist-open-data.s3.amazonaws.com/photos/1234/large.jpg",
            "https://inaturalist-open-data.s3.amazonaws.com/photos/1234/large.jpg?redirect=x",
            "https://inaturalist-open-data.s3.amazonaws.com/photos/1234/large.jpg#fragment",
            "https://inaturalist-open-data.s3.amazonaws.com/photos/1234/../../../large.jpg"
        ):
            self.assertFalse(m.original_source_loader_supported(url))

    def test_precolour_source_id_drift_rejected(self):
        filename=self.plan/"shards"/"shard-000.csv"
        filename.write_text(filename.read_text()+"unlicensed,photo\n")
        with self.assertRaisesRegex(ValueError,"SHA"):
            m.input_rows(self.plan,0)

if __name__=="__main__":unittest.main()
