import unittest
from scripts.analysis.mois_bounded_visual_review_sources import image_url
class TestSafeReviewSource(unittest.TestCase):
 def row(self,**kwargs):
  d=dict(photo_id='123',photo_url_from_api='https://inaturalist-open-data.s3.amazonaws.com/photos/123/square.jpeg',record_status='resolved',photo_membership='matched',taxon_agreement='exact',photo_license_code='cc-by-nc');d.update(kwargs);return d
 def test_only_size_changed(self):
  self.assertEqual(image_url(self.row()),'https://inaturalist-open-data.s3.amazonaws.com/photos/123/medium.jpeg')
 def test_wrong_photo_rejected(self):
  with self.assertRaises(ValueError):image_url(self.row(photo_id='456'))
 def test_wrong_host_rejected(self):
  with self.assertRaises(ValueError):image_url(self.row(photo_url_from_api='https://example.com/photos/123/square.jpeg'))
 def test_missing_licence_rejected(self):
  with self.assertRaises(ValueError):image_url(self.row(photo_license_code=''))
 def test_unverified_record_rejected(self):
  with self.assertRaises(ValueError):image_url(self.row(photo_membership='unverified'))
if __name__=='__main__':unittest.main()
