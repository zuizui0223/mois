"""Offline device-sRGB -> Lab dominant flower BOX colour, no botany re-labelling."""
import math
import unittest
import numpy as np
from PIL import Image

from scripts.acquisition.mois_flower_box_lab_colour import (
    PROTOCOL, DESCRIPTOR_COLUMNS, flower_box_lab_descriptor, rgb_to_lab
)

class LabFloralBoxPixelTests(unittest.TestCase):
    def test_genuine_D65_reference_RGB_primaries_and_achromatic(self):
        rgb=np.array([[255,255,255],[0,0,0],[255,0,0],
                      [0,255,0],[0,0,255]],dtype=np.uint8)
        l=rgb_to_lab(rgb)
        self.assertAlmostEqual(float(l[0,0]),100,places=2)
        self.assertLess(np.hypot(l[0,1],l[0,2]),0.005)
        self.assertAlmostEqual(float(l[1,0]),0,places=3)
        self.assertAlmostEqual(float(l[2,0]),53.24,delta=.10)
        self.assertAlmostEqual(float(l[3,0]),87.73,delta=.10)
        self.assertAlmostEqual(float(l[4,0]),32.30,delta=.10)
        with self.assertRaisesRegex(ValueError,"sRGB values"):
            rgb_to_lab(np.array([[2.0,0,0]]))

    def test_white_roi_is_bright_neutral_without_fake_hue(self):
        im=Image.new("RGB",(80,80),(255,255,255))
        r=flower_box_lab_descriptor(im,(0,0,80,80))
        self.assertEqual(set(r),set(DESCRIPTOR_COLUMNS))
        self.assertEqual(r["floral_box_colour_protocol"],PROTOCOL)
        self.assertEqual(r["floral_box_main_hue_degrees"],"")
        self.assertEqual(r["floral_box_hue_circular_mean_degrees"],"")
        self.assertEqual(r["floral_box_main_rgb_hex"],"#FFFFFF")
        self.assertEqual(float(r["floral_box_main_mode_fraction"]),1)
        self.assertEqual(float(r["floral_box_fraction_white_like"]),1)
        self.assertEqual(float(r["floral_box_fraction_green_like"]),0)
        self.assertEqual(float(r["floral_box_fraction_brown_like"]),0)
        self.assertEqual(r["floral_box_may_include_nonfloral_pixels"],"true")

    def test_green_flower_box_is_never_forced_into_white(self):
        # Deliberately allow genuine green flower OR nuisance foliage signal:
        # no source image pixel may be dropped merely for being green.
        im=Image.new("RGB",(80,80),(20,160,35))
        r=flower_box_lab_descriptor(im,(0,0,80,80))
        self.assertGreater(float(r["floral_box_main_C_star"]),20)
        self.assertGreater(float(r["floral_box_fraction_green_like"]),.95)
        self.assertLess(float(r["floral_box_fraction_white_like"]),.01)
        self.assertGreater(float(r["floral_box_main_hue_degrees"]),100)
        self.assertLess(float(r["floral_box_main_hue_degrees"]),195)

    def test_brown_is_distinct_from_white_despite_lower_brightness(self):
        im=Image.new("RGB",(80,80),(115,67,27))
        r=flower_box_lab_descriptor(im,(0,0,80,80))
        self.assertGreater(float(r["floral_box_fraction_brown_like"]),.95)
        self.assertLess(float(r["floral_box_fraction_white_like"]),.01)
        self.assertLess(float(r["floral_box_main_L_star"]),55)
        self.assertGreater(float(r["floral_box_main_C_star"]),12)

    def test_two_colour_flower_region_stores_two_modes_not_forced_single_hue(self):
        pixel=np.zeros((80,80,3),dtype=np.uint8)
        pixel[:,:40,:]=[255,255,255]
        pixel[:,40:,:]=[42,95,195]
        im=Image.fromarray(pixel)
        r=flower_box_lab_descriptor(im,(0,0,80,80))
        # BILINEAR resizing is same size, so each colour occupies half the box.
        self.assertAlmostEqual(float(r["floral_box_main_mode_fraction"]),.5,delta=.01)
        self.assertAlmostEqual(float(r["floral_box_second_mode_fraction"]),.5,delta=.01)
        self.assertGreater(float(r["floral_box_fraction_white_like"]),.45)
        self.assertNotEqual(r["floral_box_may_include_nonfloral_pixels"],"false")

    def test_box_does_not_require_semantically_classified_petals(self):
        im=Image.new("RGB",(80,80),(123,125,129))
        r=flower_box_lab_descriptor(im,(15,20,70,66))
        self.assertEqual(r["floral_box_main_hue_degrees"],"")
        self.assertEqual(float(r["floral_box_fraction_low_chroma"]),1)
        with self.assertRaisesRegex(ValueError,"geometry"):
            flower_box_lab_descriptor(im,(1,float("nan"),70,80))

if __name__=="__main__":
    unittest.main()
