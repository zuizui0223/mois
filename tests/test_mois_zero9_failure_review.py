import unittest
import pandas as pd
from scripts.analysis.mois_zero9_failure_review import (
    select_zero_taxa, shortlist_original_images, html_review
)

class TestMoisZero9FailureReview(unittest.TestCase):
    def test_42111_and_2200_source_taxa_preserved(self):
        rows=[]
        for n in range(42111):
            if n<90:
                attempts=13
            elif n<96:
                attempts=18
            else:
                attempts=3 if n<2200 else 1
            rows.append({"inat_taxon_id":str(500000+n),
                         "after_strict9":0 if n<2200 else 1,
                         "after_photo_attempts":attempts})
        src=pd.DataFrame(rows)
        zero,high=select_zero_taxa(src)
        self.assertEqual(len(zero),2200)
        self.assertEqual(len(high),96)
        self.assertEqual(int(high.after_photo_attempts.sum()),1278)
        self.assertEqual((zero.review_tier=="P0_10plus").sum(),96)

    def test_independent_photo_id_first_per_taxon(self):
        rows=[]
        first=[]
        for k in range(96):
            tax=str(10000+k)
            first.append({"inat_taxon_id":tax,"species":"Plant "+tax,
                          "after_photo_attempts":10,
                          "full_93m_original_photo_count":20,
                          "full_93m_original_observation_count":10})
            for j in range(2):
                rows.append({"inat_taxon_id":tax,"photo_id":str(100000+k*2+j),
                             "observation_id":"uuid"+str(k),
                             "strict9":False,"source_subarm":"FCP",
                             "source_design":"historical","observed_on":"2025-01-01",
                             "original_roi_four_hue":"unknown"})
        frame=pd.DataFrame(rows);first=pd.DataFrame(first)
        base=first[["inat_taxon_id"]].copy()
        base["after_photo_attempts"]=2
        queue,n=shortlist_original_images(frame,first,base,per_taxon=3)
        self.assertEqual(n,192)
        self.assertEqual(queue.inat_taxon_id.nunique(),96)
        self.assertEqual(len(queue),96)
        self.assertTrue(queue.review_order.eq(1).all())
        self.assertTrue(queue.manual_stage.eq("").all())
        self.assertTrue(queue.photo_page_url.str.startswith(
            "https://www.inaturalist.org/photos/").all())
        page=html_review(queue)
        self.assertIn("CSVを書き出す",page)
        self.assertIn("photo_id",page)

    def test_reject_positive_roi_in_zero_taxon(self):
        frame=pd.DataFrame([{"inat_taxon_id":str(i),"photo_id":str(1000+i),
                             "observation_id":"obs"+str(i),
                             "strict9":(i==7),"source_subarm":"FCP",
                             "source_design":"old","observed_on":"2025",
                             "original_roi_four_hue":"unknown"}
                            for i in range(96)])
        high=pd.DataFrame([{"inat_taxon_id":str(i),"species":"Plant",
                            "after_photo_attempts":10,
                            "full_93m_original_photo_count":20,
                            "full_93m_original_observation_count":10}
                           for i in range(96)])
        base=high[["inat_taxon_id"]].copy()
        base["after_photo_attempts"]=1
        with self.assertRaisesRegex(ValueError,"unexpected original true ROI9"):
            shortlist_original_images(frame,high,base)

if __name__=="__main__":
    unittest.main()
