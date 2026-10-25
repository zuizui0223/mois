"""Offline tests of long-lived S3 photo workers (zero public media requests)."""
from collections import Counter
import csv
import gzip
import io
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from PIL import Image
from scripts.acquisition import mois_93m_persistent_worker as p
from scripts.acquisition import mois_full93m_source_part_chunks as split
from scripts.acquisition.mois_onepass_full93m_target_tiers import COLUMNS
from scripts.acquisition.mois_full93m_medium_pixel_scan import ProviderBackpressure

class Response:
    def __init__(self, status, payload=b"image"):
        self.status=status
        self.payload=payload
        self.headers={"Content-Length":str(len(payload))}
        self.released=False
    def read(self,count):
        return self.payload[:count]
    def release_conn(self):
        self.released=True
    def close(self):
        self.released=True

class FakePool:
    def __init__(self,status=200,payload=b"hello"):
        self.status=status
        self.payload=payload
        self.calls=0
        self.responses=[]
    def request(self,method,url,**kwargs):
        assert method=="GET"
        assert kwargs["preload_content"] is False
        self.calls+=1
        response=Response(self.status,self.payload)
        self.responses.append(response)
        return response

class PooledOpenDataUnitTests(unittest.TestCase):
    def test_connection_pool_permits_multiple_url_bytes_with_one_pacer(self):
        pool=FakePool()
        getter=p.PooledOpenDataDownloader(100000,pool=pool,retries=0)
        urls=[p.HOST+"%s/medium.jpg"%i for i in [10,20,30]]
        self.assertEqual([getter(u) for u in urls],[b"hello"]*3)
        self.assertEqual(pool.calls,3)
        self.assertTrue(all(response.released for response in pool.responses))
        self.assertEqual(getter.received["downloaded_photo_count"],3)
        self.assertEqual(getter.received["downloaded_medium_bytes"],15)

    def test_403_429_503_stop_immediately_but_404_is_technical_unknown(self):
        for status in [403,429,503]:
            with self.subTest(status=status):
                pool=FakePool(status=status)
                getter=p.PooledOpenDataDownloader(100000,pool=pool,retries=3)
                with self.assertRaisesRegex(ProviderBackpressure,str(status)):
                    getter(p.HOST+"77/medium.jpg")
                self.assertEqual(pool.calls,1)
                self.assertTrue(pool.responses[0].released)
        pool=FakePool(status=404)
        getter=p.PooledOpenDataDownloader(100000,pool=pool,retries=3)
        with self.assertRaisesRegex(RuntimeError,"no_longer_available"):
            getter(p.HOST+"88/medium.jpg")
        self.assertEqual(pool.calls,1)

    def test_reject_any_media_source_outside_signed_open_data_bucket(self):
        getter=p.PooledOpenDataDownloader(100000,pool=FakePool(),retries=0)
        with self.assertRaisesRegex(ValueError,"Only signed"):
            getter("https://example.com/33/medium.jpg")

class LongLivedExactSourcePhotoTests(unittest.TestCase):
    def setUp(self):
        self.temp=TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.source=self.root/"source.csv.gz"
        with gzip.open(self.source,"wt",newline="",encoding="utf-8") as f:
            w=csv.DictWriter(f,fieldnames=COLUMNS,lineterminator="\n")
            w.writeheader()
            for i in range(7):
                item={k:"" for k in COLUMNS}
                item.update(photo_id=str(5000+i),
                            observation_uuid="obs-"+str(i),inat_taxon_id="13",
                            observed_on="2024-07-20",latitude="40",
                            longitude="140",photo_license="cc-by",
                            image_extension="jpg",tier="12")
                w.writerow(item)
        self.chunks=self.root/"signed_chunks"
        split.run(self.source,self.chunks,expected_source_sha256=p.sha_file(self.source),
                  expected_photo_count=7,tier=12,part=4,chunk_size=4)
        self.output=self.root/"scanned"
        im=Image.new("RGB",(30,30),(220,220,220))
        mem=io.BytesIO()
        im.save(mem,format="JPEG")
        self.picture=mem.getvalue()
        self.requested=[]
        self.detector_calls=0

    class Detector:
        def __call__(self,images):
            return [{"score":.99,"xyxy":[0,0,30,30]} for _ in images]

    def fetch(self,url):
        self.requested.append(url)
        return self.picture

    def test_source_slots_disjoint_and_sum_to_all_original_photos(self):
        plans=[p.plan([self.chunks],slot=i,slots=4) for i in range(4)]
        self.assertEqual(sum(x["assigned_original_photo_ids_to_this_slot"] for x in plans),7)
        self.assertEqual(sum(x["assigned_chunks_to_this_slot"] for x in plans),2)
        self.assertTrue(all(x["metadata_covers_all_93209351_photo_ids"] is False
                            for x in plans))
        assigned=[(str(root),idx) for x in plans for root,receipt,idx,n in x["assigned"]]
        self.assertEqual(len(set(assigned)),2)

    def test_one_warm_detector_and_executor_reused_across_chunks_then_resumed(self):
        a=p.execute([self.chunks],self.output,weights=self.root/"weight.pt",
                    weight_sha="a"*64,slot=0,slots=1,approved_aggregate_rps=50,
                    detector=self.Detector(),downloader=self.fetch,
                    download_workers=2,batch_size=16)
        self.assertEqual(a["actual_source_photo_attempts_in_verified_completed_chunks"],7)
        self.assertEqual(a["medium_photo_pixels_decoded"],7)
        self.assertEqual(a["chunks_source_sha_verified_complete"],2)
        self.assertEqual(a["technical_statuses"],{"candidate":7})
        self.assertEqual(len(self.requested),7)
        self.assertEqual(len(set(self.requested)),7)
        self.assertFalse(a["full_93m_pixel_census_complete"])
        b=p.execute([self.chunks],self.output,weights=self.root/"weight.pt",
                    weight_sha="a"*64,slot=0,slots=1,approved_aggregate_rps=50,
                    detector=self.Detector(),downloader=self.fetch,
                    download_workers=2,batch_size=32)
        self.assertEqual(b["actual_source_photo_attempts_in_verified_completed_chunks"],7)
        self.assertEqual(b["previously_finished_original_source_photo_pixels_reused"],7)
        self.assertEqual(b["new_source_photo_attempts_this_process"],0)
        self.assertEqual(len(self.requested),7,"Existing signed pixels were fetched again")

    def test_chunk_tampering_fails_before_any_photo_request(self):
        target=self.chunks/"chunks/chunk-00000.csv.gz"
        target.write_bytes(b"tampered")
        with self.assertRaisesRegex(ValueError,"SHA verified"):
            p.execute([self.chunks],self.output,weights=self.root/"weight.pt",
                      weight_sha="a"*64,approved_aggregate_rps=3,
                      detector=self.Detector(),downloader=self.fetch,
                      download_workers=1,batch_size=16)
        self.assertEqual(self.requested,[])

    def test_execute_requires_positive_explicit_budget(self):
        with self.assertRaisesRegex(ValueError,"explicit aggregate"):
            p.execute([self.chunks],self.output,weights=self.root/"weight.pt",
                      weight_sha="a"*64,approved_aggregate_rps=0,
                      detector=self.Detector(),downloader=self.fetch,
                      download_workers=1,batch_size=16)
        self.assertEqual(self.requested,[])

if __name__=="__main__":
    unittest.main()
