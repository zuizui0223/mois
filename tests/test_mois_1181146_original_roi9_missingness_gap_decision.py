"""ROI9 missingness decisions never infer white from failed image pixels."""
import tempfile,unittest
from pathlib import Path
import numpy as np
import pandas as pd

from scripts.analysis.mois_1181146_original_roi9_missingness_gap_decision import (
    latitude_band,recorded_local_season,source_action,audit,review_queue
)

class RealPhotoBiasTests(unittest.TestCase):
    def test_local_season_and_geography_do_not_impute_missing_values(self):
        lat=pd.Series([38.,-24.,-15.,70.,np.nan])
        dates=pd.Series(["2025-01-10","2025-01-10","2025-08-20",
                         "2025-07-01","unknown"])
        self.assertEqual(latitude_band(lat).tolist(),[
            "N_temperate","S_temperate","S_tropical","N_arctic","unknown"
        ])
        self.assertEqual(recorded_local_season(dates,lat).tolist(),[
            "winter","summer","winter","summer","unknown"
        ])

    def test_taxon_actions_use_number_of_old_failed_photo_attempts(self):
        f=lambda n,a,v,p:pd.Series({
            "n_unique_source_photo_ids":n,"after_photo_attempts":a,
            "after_strict9":v,"priority_after":p
        })
        p2="P2_source_has_photos_no_colour_10plus"
        self.assertEqual(source_action(f(800,13,0,p2)),
             "REVIEW_EXISTING_REPEATED_ROI9_FAILURES_BEFORE_FETCH")
        self.assertEqual(source_action(f(800,3,0,p2)),
             "REVIEW_FEW_OLD_UNKNOWN_PHOTOS_AND_SOURCE_IF_NEEDED")
        self.assertEqual(source_action(f(800,0,0,p2)),
             "PREFER_FIRST_SOURCE_PHOTO_IF_AVAILABLE")
        self.assertEqual(source_action(f(0,0,0,"P0_no_original_source_photo")),
             "P0_NO_SOURCE_OPPORTUNITY")

    def fake(self,root):
        original=[];taxa=[]
        for i in range(25):
            good=i>=21
            tax=str(100+i)
            original.append({
                "photo_id":str(90000+i),
                "observation_id":"uu-"+str(i),
                "inat_taxon_id":tax,
                "source_design":("original_1043930","bulk128_131072",
                  "bias_refill_2048","p2_first_4096")[i%4],
                "source_subarm":"targeted" if i%4 else "historical",
                "observed_on":"2025-07-12",
                "latitude":-11. if i%2 else 40.,
                "longitude":120.,
                "original_roi_four_hue":"white" if good else "mixed_uncertain",
                "strict9":good,
                "strict9_white":1.0 if good else np.nan
            })
            taxa.append({
                "inat_taxon_id":tax,"n_unique_source_photo_ids":100,
                "after_photo_attempts":1,"after_strict9":int(good),
                "priority_after":"P3_thin_colour_1to4" if good else
                                "P2_source_has_photos_no_colour_10plus"
            })
        d=Path(root)/"fake.parquet"
        q=Path(root)/"fake-taxa.csv.gz"
        pd.DataFrame(original).to_parquet(d,index=False)
        pd.DataFrame(taxa).to_csv(q,index=False,compression="gzip")
        return d,q

    def test_independent_technical_review_samples_do_not_promote_unknown(self):
        with tempfile.TemporaryDirectory() as folder:
            p,q=self.fake(folder)
            res=audit(p,q,Path(folder)/"results",
                      verify_sha=False,verify_census=False,review_size=16)
            self.assertEqual(res["exact_source_original_photo_IDs"],25)
            self.assertEqual(res["real_strict9_original_photos"],4)
            self.assertEqual(res["source_nominal_taxa_zero_strict9"],21)
            self.assertEqual(res["manual_unknown_original_photo_review_targets"],16)
            review=pd.read_csv(Path(folder)/"results"/
                        "manual_roi_failure_review_photoIDs.csv")
            self.assertEqual(review.photo_id.nunique(),16)
            self.assertEqual(review.inat_taxon_id.nunique(),16)
            self.assertTrue(review.reason.eq(
                        "FOUR_HUE_UNKNOWN_AND_NO_STRICT9").all())
            self.assertFalse(review.botanical_organ_or_flower_stage_verified.any())
            original=pd.read_parquet(p)
            original.loc[0,"photo_id"]=original.loc[1,"photo_id"]
            original.to_parquet(p,index=False)
            with self.assertRaisesRegex(ValueError,"Duplicate photo"):
                audit(p,q,Path(folder)/"invalid",
                      verify_sha=False,verify_census=False,review_size=16)

if __name__=="__main__":
    unittest.main()
