"""Guard mask QA diagnostic-only capture; never alter canonical FCP-v4."""
import unittest
import numpy as np
from scripts.analysis.mois_visual5_frozen_mask_qa_20261115 import (
    mask_geometry, require_exact_source_url, verify_replayed_counts, TARGETS
)

class TestFrozenMaskQA(unittest.TestCase):
    def test_five_post_hoc_frozen_targets(self):
        self.assertEqual(set(TARGETS),{"164665576","111641760","255188296",
                                       "446129531","487988985"})
    def test_mask_geometry_pixel_bbox_and_centroid(self):
        mask=np.zeros((8,9),dtype=bool)
        mask[2:5,3:7]=True
        g=mask_geometry(mask)
        self.assertEqual(g["mask_pixels"],12)
        self.assertEqual(g["bbox_xyxy"],[3,2,7,5])
        self.assertEqual(g["center_xy"],[4.5,3.0])
        self.assertAlmostEqual(g["fraction_of_frame"],12/72)
    def test_empty_mask_cannot_be_claimed_valid(self):
        with self.assertRaisesRegex(ValueError,"empty mask"):
            mask_geometry(np.zeros((6,8),dtype=bool))
    def test_mismatched_photo_or_host_rejected(self):
        require_exact_source_url("https://inaturalist-open-data.s3.amazonaws.com/photos/255188296/large.jpeg","255188296")
        for url in ("https://evil.tld/photos/255188296/large.jpeg",
                    "https://inaturalist-open-data.s3.amazonaws.com/photos/5/large.jpeg",
                    "http://inaturalist-open-data.s3.amazonaws.com/photos/255188296/large.jpeg",
                    "https://inaturalist-open-data.s3.amazonaws.com/photos/255188296/medium.jpeg"):
            with self.assertRaisesRegex(ValueError,"photo URL"):
                require_exact_source_url(url,"255188296")
    def test_count_comparison_is_exact(self):
        names=("white","yellow","orange","red","pink","magenta","purple","blue",
               "bronze","green","brown","black")
        expected={"palette_count_"+name:(2 if name=="white" else 0) for name in names}
        verify_replayed_counts({"white":2,**{n:0 for n in names if n!="white"}},expected)
        expected["palette_count_white"]=1
        with self.assertRaisesRegex(ValueError,"palette counts differ"):
            verify_replayed_counts({"white":2,**{n:0 for n in names if n!="white"}},expected)

if __name__=="__main__":
    unittest.main()
