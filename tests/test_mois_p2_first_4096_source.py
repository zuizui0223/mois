"""P2-only taxon-once source quotas and signed true-photo identity sanity tests.

No provider photos are downloaded by any regression test.
"""
import sqlite3
import unittest

from scripts.acquisition.mois_balance1043930_tier07_2048_source import (
    create_db,UPSERT
)
from scripts.acquisition.mois_p2_first_4096_source import (
    ARM,SCHEMA,effort_bin,gap_quota,select_one_per_taxon,
    stable_rank,
)

class OriginalP2SelectionTests(unittest.TestCase):
    def test_rare_regions_are_supported_not_fabricated(self):
        counts={
            ("N_temperate","summer"):1500,
            ("S_tropical","winter"):200,
            ("N_tropical","autumn"):260,
            ("S_subpolar","spring"):2,
        }
        quota=gap_quota(counts,512)
        self.assertEqual(sum(quota.values()),512)
        self.assertGreater(quota[("S_tropical","winter")],0)
        self.assertEqual(quota[("S_subpolar","spring")],2)
        self.assertLessEqual(quota[("N_temperate","summer")],512)
        self.assertFalse(any(x not in counts for x in quota))

    def test_source_one_observation_one_photo_per_p2_taxon(self):
        db=sqlite3.connect(":memory:")
        try:
            create_db(db)
            allowed={}
            for i in range(640):
                taxon=str(i+10000)
                uuid=f"test-P2-{taxon}"
                selected_cell=i%3
                z,t=(("N_temperate","summer"),
                     ("N_tropical","spring"),
                     ("S_tropical","winter"))[selected_cell]
                lat=(38.0,10.0,-11.0)[selected_cell]
                date=("2024-07-10","2023-04-16","2024-07-22")[selected_cell]
                photo_id=str(800000000+i)
                db.execute(UPSERT,(uuid,taxon,date,str(lat),"100",
                    "15","observer"+str(i),"cc-by",photo_id,"jpg",
                    stable_rank("photo",uuid,photo_id),z,t,
                    stable_rank("reference",uuid),
                    stable_rank("anchor_observation",uuid),i%128))
                allowed[taxon]={
                    "species":"Taxon "+taxon,
                    "prior_image_attempts":0 if i%5<3 else (3 if i%5==3 else 6),
                    "frame_photo_opportunities":45,
                }
            selected,cells,effort,eligible,obs=select_one_per_taxon(
                db,allowed,target=512)
            self.assertEqual(len(selected),512)
            self.assertEqual(eligible,640)
            self.assertEqual(obs,640)
            self.assertEqual(len({x["taxon"] for x in selected}),512)
            self.assertEqual(len({x["uuid"] for x in selected}),512)
            self.assertEqual(len({x["photo_id"] for x in selected}),512)
            self.assertTrue(all(
                 x["conditional_design_probability"]>0 for x in selected))
            self.assertEqual(sum(x["selected_photo_quota"] for x in cells),512)
            self.assertEqual(sum(effort.values()),512)
            self.assertGreater(effort.get("E0_1",0),0)
            self.assertIn("P2",ARM)
            self.assertIn("taxon_once",SCHEMA)
        finally:
            db.close()

    def test_selection_effort_bucket_contract(self):
        self.assertEqual(effort_bin(0),"E0_1")
        self.assertEqual(effort_bin(1),"E0_1")
        self.assertEqual(effort_bin(2),"E2_4")
        self.assertEqual(effort_bin(4),"E2_4")
        self.assertEqual(effort_bin(5),"E5plus")

if __name__=="__main__":
    unittest.main()
