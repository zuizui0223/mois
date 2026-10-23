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
    def test_precolour_source_id_drift_rejected(self):
        filename=self.plan/"shards"/"shard-000.csv"
        filename.write_text(filename.read_text()+"unlicensed,photo\n")
        with self.assertRaisesRegex(ValueError,"SHA"):
            m.input_rows(self.plan,0)


    def test_signed_routed_tier06_old_source_recovery_without_new_image_inference(self):
        self.manifest.pop("previously_measured_photo_ids_skipped")
        self.manifest.update({
            "source_tier":6,"source_group":0,"compute_bucket_count":512,
            "admitted_bucket_count":64,"not_probability_complete_global_source":True,
            "previous_source_photo_IDs_not_redownloaded":True,
            "selected_bucket_photo_IDs_before_dedup":10,
            "selected_new_photos_unsupported_by_frozen_loader":0,
            "prior_and_technical_excluded_source_photo_IDs":{
                "prior_actual_photo_or_observation":3,
                "prior_technical_photo_or_observation":2},
        })
        (self.plan/"manifest.json").write_text(json.dumps(self.manifest))
        folder=self.root/"measured"
        m.measure_group(self.plan,0,self.root/"runtime",self.root/"sam",
                        folder,runner=self.fake_runner)
        def forbidden(*args,**kwargs):
            raise AssertionError("Collector unexpectedly ran original pixel inference")
        # Production recovery only consumes already signed unit files.
        data=m.collect(self.plan,folder,self.root/"final")
        self.assertEqual((data["attempted"],data["classifiable"]),(5,5))
        self.assertEqual(data["matching_prior_measured_photos_reused"],3)
        self.assertEqual(data["prior_photo_or_observation_source_rows_excluded"],3)
        self.assertEqual(data["independent_technical_photo_or_observation_rows_excluded"],2)
        self.assertEqual(data["source_deduplication_reporting_grain"],
                         "prior_original_photo_or_observation")
        self.assertTrue(data["source_precolour_photo_dispositions_reconciled"])
        self.assertTrue(pd.read_parquet(self.root/"final"/"measured_photos.parquet").photo_id.is_unique)

    def test_unknown_or_inconsistent_routed_prior_source_fails_closed(self):
        self.manifest.pop("previously_measured_photo_ids_skipped")
        self.manifest.update({
            "source_tier":6,"source_group":0,"compute_bucket_count":512,
            "admitted_bucket_count":64,"not_probability_complete_global_source":True,
            "previous_source_photo_IDs_not_redownloaded":True,
            "selected_bucket_photo_IDs_before_dedup":100,
            "selected_new_photos_unsupported_by_frozen_loader":0,
            "prior_and_technical_excluded_source_photo_IDs":{
                "prior_actual_photo_or_observation":3,
                "prior_technical_photo_or_observation":2},
        })
        with self.assertRaisesRegex(ValueError,"precolour prior/technical"):
            m.source_prior_dedup_receipt(self.manifest)
        self.manifest["selected_bucket_photo_IDs_before_dedup"]=10
        self.manifest["source_group"]=5
        with self.assertRaisesRegex(ValueError,"Unknown source metadata layout"):
            m.source_prior_dedup_receipt(self.manifest)

    def test_legacy_prior_photo_skip_contract_is_unchanged(self):
        s=m.source_prior_dedup_receipt(self.manifest)
        self.assertEqual(s["matching_prior_measured_photos_reused"],0)
        self.assertEqual(s["source_deduplication_reporting_grain"],"prior_original_photo_id")

if __name__=="__main__":unittest.main()
