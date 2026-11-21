import unittest
from scripts.analysis.mois_existing38_remeasurement import large_url, classify

class TestExistingPhotoReplay(unittest.TestCase):
    def test_large_version_is_explicit_not_original(self):
        url=large_url('https://inaturalist-open-data.s3.amazonaws.com/photos/123/square.jpeg','123')
        self.assertEqual(url,'https://inaturalist-open-data.s3.amazonaws.com/photos/123/large.jpeg')
        self.assertNotIn('/original.',url)
    def test_wrong_photo_id_rejected(self):
        with self.assertRaises(ValueError):large_url('https://inaturalist-open-data.s3.amazonaws.com/photos/9/medium.jpg','123')
    def test_wrong_host_rejected(self):
        with self.assertRaises(ValueError):large_url('https://example.com/photos/123/medium.jpg','123')
    def test_complete_finite_composition(self):
        self.assertEqual(classify('white',[.8,.2,0,0,0,0,0,0,0]),(True,'white'))
    def test_numeric_mixed_is_not_forced_into_hue(self):
        self.assertEqual(classify('red_pink',[.3,0,0,.3,.4,0,0,0,0]),(True,''))
    def test_unknown_does_not_become_measured_colour(self):
        self.assertEqual(classify('mixed_uncertain',[1,0,0,0,0,0,0,0,0]),(False,''))
    def test_classified_but_invalid_numeric_stops(self):
        with self.assertRaises(ValueError):classify('white',[1,1,0,0,0,0,0,0,0])

if __name__=='__main__':unittest.main()
