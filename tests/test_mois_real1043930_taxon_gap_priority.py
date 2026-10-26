"""Decision support tests for the original 42,111-taxon photo-gap register."""
import unittest
import pandas as pd
import numpy as np

from scripts.analysis.mois_real1043930_taxon_gap_priority import (
    priority,source_bin,conservative_photo_slots,profile
)

class SourcePhotoPriorityTests(unittest.TestCase):
    def test_all_photo_opportunity_tiers_and_no_photos(self):
        cases=[
            (0,0,"P0_no_original_source_photo"),
            (3,0,"P1_source_rare_no_colour_1to9"),
            (9,0,"P1_source_rare_no_colour_1to9"),
            (10,0,"P2_source_has_photos_no_colour_10plus"),
            (10000,0,"P2_source_has_photos_no_colour_10plus"),
            (1,1,"P3_thin_colour_1to4"),
            (500,4,"P3_thin_colour_1to4"),
            (500,5,"P4_thin_colour_5to9"),
            (1500,9,"P4_thin_colour_5to9"),
            (50000,10,"P5_existing_10plus_strict9"),
        ]
        for n,v,want in cases:
            with self.subTest(n=n,v=v):
                self.assertEqual(priority(n,v),want)

    def test_frame_available_photo_slots_not_claimed_new_photo_download(self):
        self.assertEqual(source_bin(0),"00_no_source")
        self.assertEqual(source_bin(9),"01_1to9")
        self.assertEqual(source_bin(99),"02_10to99")
        self.assertEqual(source_bin(100),"03_100to999")
        self.assertEqual(source_bin(1000),"04_1000to9999")
        self.assertEqual(source_bin(10000),"05_10000plus")
        self.assertEqual(conservative_photo_slots(2,1),2)
        self.assertEqual(conservative_photo_slots(100,5),3)
        self.assertEqual(conservative_photo_slots(0,0),0)

    def test_priorities_sum_actual_taxa_and_photo_opportunities(self):
        d=pd.DataFrame({
            "priority":["P1_source_rare_no_colour_1to9",
                        "P2_source_has_photos_no_colour_10plus",
                        "P3_thin_colour_1to4",
                        "P5_existing_10plus_strict9"],
            "n_unique_source_photo_ids":[2,20,100,1000],
            "n_image_photo_attempts":[1,2,3,100],
            "n_strict_roi9":[0,0,1,25]
        })
        p=profile(d)
        self.assertEqual(sum(x["source_frame_taxa"] for x in p),4)
        self.assertEqual(sum(x["source_original_photo_opportunities"] for x in p),1122)
        self.assertEqual(sum(x["source_candidate_slots_first_three_per_taxon"] for x in p),11)
        self.assertTrue(all(x["not_actual_new_photo_IDs_unless_source_exactly_deduplicated"] for x in p))

if __name__=="__main__":
    unittest.main()
