"""True original-photo ID technical ROI review must reuse already scored images."""
import unittest
import pandas as pd
from scripts.analysis.mois_512_prior_unknown_technical_screen_plan import make_plan

class RoiUnknownReviewReuseTests(unittest.TestCase):
    def frame(self):
        rows=[]
        for i in range(12):
            rows.append(dict(
                photo_id=str(100000+i),
                observation_id=f"old-observation-{i}",
                inat_taxon_id=str(2000+i),
                source_design=("bulk128_131072" if i<4 else
                    ("historical_K1_FCP_source_heterogeneous" if i<8 else
                     "p2_first_4096")),
                latitude_zone=["N_temperate","S_tropical","S_temperate","N_tropical"][i%4],
                local_season=["spring","winter","autumn"][i%3],
                postoutcome_source_gap=(
                    "P2_source_has_photos_no_colour_10plus" if i%2 else
                    "P1_source_rare_no_colour_1to9"),
                previous_original_photo_attempts=2,
                source_photo_opportunities=100,
                reason="FOUR_HUE_UNKNOWN_AND_NO_STRICT9",
                botanical_organ_or_flower_stage_verified=False,
            ))
        return pd.DataFrame(rows)
    def bulk(self):
        rows=[]
        for i in range(4):
            rows.append(dict(
                photo_id=str(100000+i),
                observation_id=f"old-observation-{i}",
                inat_taxon_id=str(2000+i),
                strict_original_roi9=False,
                screen_status="candidate" if i%2 else "negative",
                medium_actual_RGB_decoded=True,
                floral_box_main_L_star="32.4" if i%2 else "",
            ))
        return pd.DataFrame(rows)
    def test_original_unknown_500px_reuse_excludes_screened_photoIDs(self):
        cached,queue,counts=make_plan(self.frame(),self.bulk(),
                                      workers=2,strict=False)
        self.assertEqual(counts["source_review_photo_rows"],12)
        self.assertEqual(len(cached),4)
        self.assertEqual(len(queue),8)
        self.assertEqual(set(cached.photo_id),set(str(100000+i) for i in range(4)))
        self.assertFalse(set(cached.photo_id)&set(queue.photo_id))
        self.assertEqual(counts["eight_worker_photo_counts"],[4,4])
        self.assertTrue(counts["prior_9palette_colour_outcomes_all_remain_UNKNOWN"])
    def test_mismatched_source_taxon_and_fake_existing_strict9_are_rejected(self):
        wrong=self.bulk()
        wrong.loc[0,"inat_taxon_id"]="999"
        with self.assertRaisesRegex(ValueError,"mismatched"):
            make_plan(self.frame(),wrong,workers=2,strict=False)
        wrong=self.bulk()
        wrong.loc[0,"strict_original_roi9"]=True
        with self.assertRaisesRegex(ValueError,"actually coloured"):
            make_plan(self.frame(),wrong,workers=2,strict=False)
    def test_duplicate_review_photo_and_already_botanical_truth_fail_closed(self):
        sample=self.frame()
        sample.loc[1,"photo_id"]=sample.loc[0,"photo_id"]
        with self.assertRaisesRegex(ValueError,"independent"):
            make_plan(sample,self.bulk(),workers=2,strict=False)
        sample=self.frame()
        sample.loc[4,"botanical_organ_or_flower_stage_verified"]=True
        with self.assertRaisesRegex(ValueError,"already botanically verified"):
            make_plan(sample,self.bulk(),workers=2,strict=False)

if __name__=="__main__":
    unittest.main()
