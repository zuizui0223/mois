import unittest
import pandas as pd
from scripts.analysis.mois_reconcile_strict9_counts import reconcile

class TestMoisSignedPhotoCountRepair(unittest.TestCase):
    def sample(self):
        b=pd.DataFrame({
          "inat_taxon_id":["a","b","c"],
          "after_strict9":[0,4,1],
          "independent_roi9_observations":[0,3,1],
          "after_photo_attempts":[2,5,1],
          "hue_status":["no_strict9","mixed_image_unresolved","provisional_single"]
        })
        c=b.copy()
        c["after_strict9"]=[1,1,1]
        c["independent_roi9_observations"]=[2,5,1]
        c["after_photo_attempts"]=[4,7,2]
        c["hue_status"]=["provisional_single","mixed_image_unresolved","provisional_single"]
        for key in ("dominant_hue9","observed_hue9_set",
                    "observed_mean_white","observed_mean_yellow","observed_mean_orange",
                    "observed_mean_red","observed_mean_pink","observed_mean_magenta",
                    "observed_mean_purple","observed_mean_blue","observed_mean_bronze"):
            b[key]=c[key]=("white" if key in ("dominant_hue9","observed_hue9_set") else 0.1)
        return b,c

    def test_realistic_counter_reset_restored_without_hue_change(self):
        b,c=self.sample()
        result=reconcile(b,c,new_attempts=5,new_numeric=4)
        self.assertEqual(result.after_strict9.tolist(),[2,6,1])
        self.assertEqual(result.after_strict9_original_signed.tolist(),[1,1,1])
        self.assertEqual(result.after_strict9_repaired.tolist(),[True,True,False])
        self.assertTrue(result.hue_status.equals(c.hue_status))
        self.assertTrue(result.dominant_hue9.equals(c.dominant_hue9))

    def test_refuses_wrong_signed_photo_delta(self):
        b,c=self.sample()
        with self.assertRaisesRegex(ValueError,"signed new original ROI success"):
            reconcile(b,c,new_attempts=5,new_numeric=5)

    def test_refuses_nonmonotone_independent_history(self):
        b,c=self.sample()
        c.loc[1,"independent_roi9_observations"]=0
        with self.assertRaisesRegex(ValueError,"nonmonotone"):
            reconcile(b,c,new_attempts=5,new_numeric=4)

if __name__=="__main__":
    unittest.main()
