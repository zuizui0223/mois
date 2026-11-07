"""No-network source identity tests for 2020 MOIS all-land island × nightlight."""
import unittest
import pandas as pd
from scripts.analysis import mois_2020_sameyear_island_alanspatial_refit_20261107 as m

class SameYearSourceTests(unittest.TestCase):
    def test_original_photo_target_landmass_is_preserved(self):
        d=pd.DataFrame([
            {"target_id":"2020|6011_2755|44","land_class":"island"},
            {"target_id":"2020|6011_2755|20","land_class":"mainland"},
            {"target_id":"2020|10_300|12","land_class":"island"}])
        v=m.parse_source_ids(d)
        self.assertEqual(v.cell_id.nunique(),2)
        self.assertEqual(v.land_id.nunique(),3)
        self.assertEqual(v.land_class.tolist(),["island","mainland","island"])
        self.assertEqual(len(v),3)
    def test_refuse_falsely_merged_island_and_mainland_target(self):
        d=pd.DataFrame([
            {"target_id":"2020|1_2|12","land_class":"island"},
            {"target_id":"2020|1_2|12","land_class":"mainland"}])
        with self.assertRaisesRegex(ValueError,"duplicate"):
            m.parse_source_ids(d)
    def test_unresolved_landmass_never_fabricated(self):
        d=pd.DataFrame([{"target_id":"2020|1_2|UNRESOLVED:shore",
                         "land_class":"island"}])
        with self.assertRaisesRegex(ValueError,"not a GSHHG"):
            m.parse_source_ids(d)
    def test_all_acquired_new2020_increment_counts_locked(self):
        self.assertEqual(sum(m.COUNTS.values()),36070)
        self.assertEqual(m.COUNTS["last30630"],30630)
        self.assertEqual(len(m.COAST_SHA),64)
        self.assertEqual(len(m.PINS["canonical"]),64)
    def test_island_light_interaction_refuses_eight_island_cells(self):
        rows=[]
        for k in range(8):
            rows.append({"region5":str(k),"land_class":"island",
                         "land_id":str(k),"white_pp":60+k})
        for k in range(39):
            rows.append({"region5":str(k),"land_class":"mainland",
                         "land_id":"100","white_pp":40+k/10})
        outcome=m.available_fit(pd.DataFrame(rows),"white_pp ~ 1",
                     "region5",require_island_support=True)
        self.assertFalse(outcome["estimable"])
        self.assertEqual(outcome["reason"],
                         "island_interaction_below_prespecified_support")
        self.assertEqual(outcome["island_cells"],8)

    def test_insufficient_clusters_not_inferred(self):
        d=pd.DataFrame({"region5":["a","b"],"white_pp":[20,40]})
        result=m.available_fit(d,"white_pp ~ 1","region5")
        self.assertFalse(result["estimable"])
        self.assertEqual(result["clusters"],2)
if __name__=="__main__":unittest.main()
