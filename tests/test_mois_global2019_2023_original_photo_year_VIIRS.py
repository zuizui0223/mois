"""Source-year and cloudfree safeguards for multi-year true-year VIIRS."""
import unittest
from scripts.analysis import mois_global2019_2023_original_photo_year_VIIRS as p
from scripts.analysis.mois_global_multipoint_numeric_20261107 import extract_samples

class TrueYearPhotoAlansourceTests(unittest.TestCase):
    def rows(self,n,year):
        return [{"target_id":f"{year}|{i}_1|77","land_class":"island" if i%2 else "mainland",
           "query_longitude":str(i*.001-100),"query_latitude":str(20+i*.0001)} for i in range(n)]
    def test_2019_and_all_year_batch_counts(self):
        for y,n,plans in [(2019,27020,212),(2021,45966,360),
                          (2022,52075,408),(2023,59862,468)]:
            self.assertEqual(p.SITES[y],n)
            q=p.plans(y,self.rows(n,y))
            self.assertEqual(len(q),plans)
            self.assertEqual(sum(k[2] for k in q),n*4)
            self.assertEqual(max(k[2] for k in q),512)
    def test_no_2020_substitution(self):
        with self.assertRaisesRegex(ValueError,"Invalid"):
            p.params(2020,1,self.rows(3,2020))
    def test_year_month_params_vary(self):
        a=p.params(2019,1,self.rows(2,2019))
        b=p.params(2023,7,self.rows(2,2023))
        self.assertNotEqual(a["time"],b["time"])
    def test_crossyear_photo_ids_fail(self):
        x=self.rows(27020,2019)
        x[0]["target_id"]="2020|a"
        with self.assertRaisesRegex(ValueError,"different observation year"):
            p.verify_source(x,2019)
    def test_duplicate_API_raster_source_rejected(self):
        x=self.rows(27020,2019)
        x[0]["query_longitude"]=x[1]["query_longitude"]
        x[0]["query_latitude"]=x[1]["query_latitude"]
        with self.assertRaisesRegex(ValueError,"Duplicate original raster"):
            p.verify_source(x,2019)
    def test_real_2022_october_provider_error_is_not_light_zero(self):
        error={"error":{"code":400,"message":"Invalid or missing input parameters.",
                        "details":["Transformation is unavailable for the current image."]}}
        self.assertTrue(p.original_2022_october_provider_error(2022,10,error))
        self.assertFalse(p.original_2022_october_provider_error(2023,10,error))
        self.assertFalse(p.original_2022_october_provider_error(2022,7,error))
        self.assertFalse(p.original_2022_october_provider_error(2022,10,{"samples":[]}))
    def test_original_month_and_cloud_count(self):
        r=self.rows(2,2021)
        obj={"samples":[{"locationId":0,
          "location":{"x":float(r[0]["query_longitude"]),"y":float(r[0]["query_latitude"])},
          "attributes":{"year":2021,"month":4},"value":"0.0 5.0"}]}
        v=extract_samples("viirs",4,r,obj,year=2021)
        self.assertEqual(v[0]["status"],"OK")
        self.assertEqual(v[0]["value"],0.0)
        self.assertEqual(v[1]["status"],"NO_DATA")
        with self.assertRaisesRegex(ValueError,"wrong source year/month"):
            extract_samples("viirs",4,r,obj,year=2023)

if __name__=="__main__":
    unittest.main()
