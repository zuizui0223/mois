import unittest
from scripts.analysis import mois_visual38_existing_v4_20261115 as m

class TestMoisRepeat38(unittest.TestCase):
    def test_signed_photo_url_uses_large_derivative_not_original(self):
        x=m.source_large_url("https://inaturalist-open-data.s3.amazonaws.com/photos/123/medium.jpeg","123")
        self.assertEqual(x,"https://inaturalist-open-data.s3.amazonaws.com/photos/123/large.jpeg")

    def test_refuses_wrong_photo_and_wrong_host(self):
        for url in ("https://evil.test/photos/123/medium.jpg",
                    "https://inaturalist-open-data.s3.amazonaws.com/photos/456/medium.jpg",
                    "https://inaturalist-open-data.s3.amazonaws.com/photos/123/original.jpg"):
            with self.assertRaises(ValueError):
                m.source_large_url(url,"123")

    def test_success_requires_finite_nine_fractions(self):
        r={"measurement_status":"classified_four_state_morph","morph":"white",
           **{"flower_fraction_"+k:("1" if k=="white" else "0") for k in m.PALETTE}}
        self.assertEqual(m.numeric_vector(r),[1.,0.,0.,0.,0.,0.,0.,0.,0.])
        r["flower_fraction_white"]="nan"
        with self.assertRaisesRegex(ValueError,"invalid real numerical"):
            m.numeric_vector(r)

    def test_unmeasured_images_remain_unknown_not_imputed(self):
        r={"measurement_status":"not_evaluable_roi_or_flip_gate","morph":"mixed_uncertain",
           **{"flower_fraction_"+k:("1" if k=="white" else "0") for k in m.PALETTE}}
        self.assertIsNone(m.numeric_vector(r))

    def test_disagreement_between_status_and_morph_fails_closed(self):
        r={"measurement_status":"classified_four_state_morph","morph":"mixed_uncertain"}
        with self.assertRaisesRegex(ValueError,"morph measurement status changed"):
            m.numeric_vector(r)

    def test_zero_nuisance_mass_is_not_replaced_by_guess(self):
        r={"measurement_status":"not_evaluable_no_biological_palette_mass","morph":"mixed_uncertain"}
        self.assertIsNone(m.numeric_vector(r))

if __name__=="__main__":
    unittest.main()
