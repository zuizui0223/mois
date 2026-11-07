"""Fail-closed checks for MOIS 2020 multipoint raw VIIRS / cloudfree / population."""
import unittest, tempfile, csv
from pathlib import Path
from scripts.analysis import mois_global_multipoint_numeric_20261107 as m

class MultipointSourceTests(unittest.TestCase):
    def rows(self):
        return [{"target_id":"2020|a","land_class":"island","query_longitude":"139.7","query_latitude":"35.7"},
                {"target_id":"2020|b","land_class":"mainland","query_longitude":"-78.3","query_latitude":"40.1"}]
    def sample(self,year=2020,month=1,cloud=5):
        return {"samples":[{"locationId":0,"location":{"x":139.7,"y":35.7},
                 "value":f"3.92 {cloud}","attributes":{"year":year,"month":month}}]}
    def test_monthly_source_has_radiance_and_cloudfree_bands(self):
        x=m.extract_samples("viirs",1,self.rows(),self.sample())
        self.assertEqual(x[0],{"status":"OK","value":3.92,"cloud":5})
        self.assertEqual(x[1]["status"],"NO_DATA")
    def test_cloudfree_one_is_NOT_usable_radiance(self):
        x=m.extract_samples("viirs",1,self.rows(),self.sample(cloud=1))
        self.assertEqual(x[0]["status"],"LOW_CLOUD_COUNT")
        self.assertIsNone(x[0]["value"])
        self.assertEqual(x[0]["cloud"],1)
    def test_wrong_month_rejected(self):
        with self.assertRaisesRegex(ValueError,"wrong source year/month"):
            m.extract_samples("viirs",1,self.rows(),self.sample(month=2))
    def test_nonfinite_radiance_rejected(self):
        x=self.sample();x["samples"][0]["value"]="NaN 5"
        with self.assertRaisesRegex(ValueError,"nonfinite"):
            m.extract_samples("viirs",1,self.rows(),x)
    def test_duplicate_location_id_rejected(self):
        s=self.sample()["samples"][0];x={"samples":[s,s]}
        with self.assertRaisesRegex(ValueError,"duplicate"):
            m.extract_samples("viirs",1,self.rows(),x)
    def test_moved_coordinate_rejected(self):
        x=self.sample();x["samples"][0]["location"]["x"]+=.1
        with self.assertRaisesRegex(ValueError,"location"):
            m.extract_samples("viirs",1,self.rows(),x)
    def test_population_zero_is_valid_but_NODATA_not_zero(self):
        p={"samples":[{"locationId":0,"location":{"x":139.7,"y":35.7},"value":"0.0"}]}
        v=m.extract_samples("population",0,self.rows(),p)
        self.assertEqual(v[0]["status"],"OK")
        self.assertEqual(v[0]["value"],0.)
        self.assertEqual(v[1]["status"],"NO_DATA")
        self.assertIsNone(v[1]["value"])
    def test_unsupport_year_and_area_latitude_not_selective(self):
        with self.assertRaisesRegex(ValueError,"year"):
            m.api_params("viirs",1,self.rows(),year=2024)
        params=m.api_params("viirs",1,self.rows())
        self.assertIn('"points"',params["geometry"])
        self.assertNotIn("area_filter",params)
    def test_target_identity_and_year_gate(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"manifest.csv"
            with p.open("w",newline="") as f:
                w=csv.DictWriter(f,fieldnames=list(self.rows()[0]));w.writeheader();w.writerows(self.rows())
            self.assertEqual(len(m.validate_targets(p)),2)
            d=self.rows();d[1]["target_id"]="2024|b"
            with p.open("w",newline="") as f:
                w=csv.DictWriter(f,fieldnames=list(d[0]));w.writeheader();w.writerows(d)
            with self.assertRaisesRegex(ValueError,"year"):
                m.validate_targets(p)
    def test_unexpected_source_dimension_rejected(self):
        with self.assertRaisesRegex(ValueError,"dimension"):
            m.extract_samples("viirs",1,self.rows(),
                {"samples":[{"locationId":0,"location":{"x":139.7,"y":35.7},
                "value":"3.92","attributes":{"year":2020,"month":1}}]})
if __name__=="__main__":unittest.main()
