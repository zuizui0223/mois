"""2025 source-type and no-2020-population mixing contracts."""
import unittest
from scripts.analysis import mois_global2025_VIIRS_WorldPopR2024B_UC_bulk_20261107 as a

class SourceTests(unittest.TestCase):
    def source(self):
        return [
         dict(target_id="2025|1_2|1",land_class="mainland",query_longitude="139.7",
              query_latitude="35.7"),
         dict(target_id="2025|3_4|2",land_class="island",query_longitude="-123.4",
              query_latitude="48.6"),
        ]
    def test_four_2025_months_unique(self):
        for m in (1,4,7,10):
            p=a.get_samples_params(self.source(),m)
            self.assertIn('"points"',p["geometry"])
            self.assertIn("time",p)
        with self.assertRaisesRegex(ValueError,"Invalid"):
            a.get_samples_params(self.source(),2)
    def test_source_2020_cannot_be_mislabeled_2025(self):
        d=self.source();d[0]['target_id']="2020|1_2|1"
        with self.assertRaisesRegex(ValueError,"Wrong photo year"):
            a.validate_targets(d,2)
    def test_true_2025_raster_file_hash_pinned(self):
        self.assertEqual(len(a.RASTER_SHA),64)
        self.assertEqual(a.RASTER_BYTES,865211479)
        self.assertEqual(a.YEAR,2025)
        self.assertIn("R2024B",a.RASTER_NAME)
    def test_island_and_mainland_are_included_without_area_threshold(self):
        self.assertEqual(len(a.validate_targets(self.source(),2)),2)
    def test_duplicate_GIS_coordinates_are_not_false_independence(self):
        d=self.source();d[1]["query_longitude"]=d[0]["query_longitude"]
        d[1]["query_latitude"]=d[0]["query_latitude"]
        with self.assertRaisesRegex(ValueError,"Duplicate source raster"):
            a.validate_targets(d,2)
    def test_numeric_source_photos_never_reclassified(self):
        self.assertEqual(a.SELECT_SHA,"45c5152e1edad997ff8749c645818996fb2b0bc03a2ab7c0753af446afaab143")
        self.assertEqual(a.REST_SHA,"68990f90ecc5157ffd21f0fb2269802dbd79be8d86de2e56ab80cfe8b7c72a65")
if __name__=="__main__":unittest.main()
