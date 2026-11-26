"""No-network tests for fresh post-93m-snapshot iNaturalist zero-ROI photo metadata."""
import unittest
import pandas as pd

from scripts.acquisition import mois_zero9_recent_inat_new_created_photo_source as m


class PostSnapshotZeroNinePhotoTests(unittest.TestCase):
    def test_select_only_zero_numeric_taxa_with_reproducible_disjoint_hash_windows(self):
        a=pd.DataFrame({
            "inat_taxon_id":[str(i+1) for i in range(42111)],
            "after_strict9":[0]*2229+[1]*(42111-2229),
        })
        q=m.choose_zero(a,sample=64,offset=0)
        p=m.choose_zero(a,sample=64,offset=64)
        self.assertEqual(len(q),64)
        self.assertEqual(len(p),64)
        self.assertFalse(set(q.inat_taxon_id)&set(p.inat_taxon_id))
        self.assertTrue(q.after_strict9.eq(0).all())
        self.assertTrue(p.after_strict9.eq(0).all())

    def test_new_creation_date_exact_taxon_and_licensed_original_photo_URL(self):
        good={"id":50,"uuid":"00000000-0000-4000-8000-000000000001",
              "created_at":"2026-11-01T10:00:00+00:00",
              "observed_on":"2026-07-01",
              "taxon":{"id":42,"name":"Plantus alba","rank":"species"},
              "captive":False,
              "photos":[{"id":555,
                  "url":"https://inaturalist-open-data.s3.amazonaws.com/photos/555/square.jpeg",
                  "license_code":"cc-by"}]}
        opts,st=m.source_rows([good],tid="42")
        self.assertEqual(len(opts),1)
        self.assertEqual(opts[0]["photo_id"],"555")
        self.assertEqual(opts[0]["photo_url_large"],
            "https://inaturalist-open-data.s3.amazonaws.com/photos/555/large.jpeg")
        self.assertEqual(opts[0]["observed_on"],"2026-07-01")
        self.assertFalse(opts[0]["previous_1_231_244_photo_IDs_verified_disjoint"])
        self.assertFalse(opts[0]["strict9_original_RGB9_measured"])
        self.assertEqual(st["license_eligible_new_photo"],1)

    def test_older_unknown_wrong_taxa_or_license_not_imputed(self):
        source={"id":1,"uuid":"00000000-0000-4000-8000-000000000001",
                "created_at":"2026-11-01T00:00:00Z","captive":False,
                "taxon":{"id":20},"photos":[{"id":300,
                "url":"https://inaturalist-open-data.s3.amazonaws.com/photos/300/square.jpg",
                "license_code":"cc-by-nc"}]}
        old=dict(source,created_at="2026-10-12T01:23:00Z")
        wrong=dict(source,taxon={"id":21})
        captive=dict(source,captive=True)
        unlicensed=dict(source,photos=[dict(source["photos"][0],license_code=None)])
        url_wrong=dict(source,photos=[dict(source["photos"][0],url="https://static.inaturalist.org/photos/300/square.jpg")])
        for x in (old,wrong,captive,unlicensed,url_wrong):
            p,_=m.source_rows([x],tid="20")
            self.assertEqual(len(p),0)
        self.assertFalse(m.created_after_freeze("2026-10-27T23:59:59Z"))
        self.assertTrue(m.created_after_freeze("2026-10-28T00:00:00Z"))

    def test_different_photo_id_and_host_must_not_be_manufactured(self):
        self.assertIsNone(m.strict_new_photo({"id":2,"license_code":"cc-by",
          "url":"https://inaturalist-open-data.s3.amazonaws.com/photos/3/square.jpg"},"1"))
        self.assertIsNone(m.strict_new_photo({"id":2,"license_code":"cc-by",
          "url":"https://evil.example.org/photos/2/square.jpg"},"1"))
        self.assertIsNone(m.strict_new_photo({"id":2,"license_code":"cc-by-nc-nd",
          "url":"https://inaturalist-open-data.s3.amazonaws.com/photos/2/square.jpg"},"1"))


if __name__=="__main__":
    unittest.main()
