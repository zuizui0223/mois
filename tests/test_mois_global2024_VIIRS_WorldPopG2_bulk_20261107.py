"""2024 source-type and no-2020-population mixing contracts."""
import unittest
from scripts.analysis import mois_global2024_VIIRS_WorldPopG2_bulk_20261107 as a

class SourceTests(unittest.TestCase):
    def source(self):
        return [
         dict(target_id="2024|1_2|1",land_class="mainland",query_longitude="139.7",
              query_latitude="35.7"),
         dict(target_id="2024|3_4|2",land_class="island",query_longitude="-123.4",
              query_latitude="48.6"),
        ]
    def test_four_2024_months_unique(self):
        for m in (1,4,7,10):
            p=a.get_samples_params(self.source(),m)
            self.assertIn('"points"',p["geometry"])
            self.assertIn("time",p)
        with self.assertRaisesRegex(ValueError,"Invalid"):
            a.get_samples_params(self.source(),2)
    def test_source_2020_cannot_be_mislabeled_2024(self):
        d=self.source();d[0]['target_id']="2020|1_2|1"
        with self.assertRaisesRegex(ValueError,"Wrong photo year"):
            a.validate_targets(d,2)
    def test_true_2024_raster_file_hash_pinned(self):
        self.assertEqual(len(a.RASTER_SHA),64)
        self.assertEqual(a.RASTER_BYTES,289380599)
        self.assertEqual(a.YEAR,2024)
        self.assertIn("R2025A",a.RASTER_NAME)
    def test_island_and_mainland_are_included_without_area_threshold(self):
        self.assertEqual(len(a.validate_targets(self.source(),2)),2)
    def test_duplicate_GIS_coordinates_are_not_false_independence(self):
        d=self.source();d[1]["query_longitude"]=d[0]["query_longitude"]
        d[1]["query_latitude"]=d[0]["query_latitude"]
        with self.assertRaisesRegex(ValueError,"Duplicate source raster"):
            a.validate_targets(d,2)
    def test_numeric_source_photos_never_reclassified(self):
        self.assertEqual(a.SELECT_SHA,"541b6134090b61a285ffa9c45ee3df8a382e66ae5628c394ca0c4d47af9207a0")
        self.assertEqual(a.REST_SHA,"a792602a05735c5049c960a80fc6a75f4a1241c0a3abaa9769f452330559bf2d")
if __name__=="__main__":unittest.main()
