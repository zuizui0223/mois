"""Global-photo bias repair sampling must reduce coverage imbalance without hue selection."""
import sqlite3
import unittest
from tempfile import TemporaryDirectory

from scripts.acquisition.mois_balance1043930_tier07_2048_source import (
    zone,season,quotas,rank,create_db,UPSERT,sample,
    REFERENCE,AUGMENTATION,SOURCE_P,SCHEMA
)

class GenuineCurrentDataGapSourceTests(unittest.TestCase):
    def test_hemispheres_tropics_arctic_and_local_photographic_seasons(self):
        cases={
            85:"N_arctic",60:"N_subpolar",35:"N_temperate",
            4:"N_tropical",-5:"S_tropical",-37:"S_temperate",
            -58:"S_subpolar",-72:"S_antarctic",
        }
        for x,name in cases.items():
            self.assertEqual(zone(x),name)
        self.assertEqual(zone("wrong"),"unknown")
        self.assertEqual(season(35,"2025-07-15"),"summer")
        self.assertEqual(season(-35,"2025-07-15"),"winter")
        self.assertEqual(season(-35,"2025-01-15"),"summer")
        self.assertEqual(season(35,"2025-01-15"),"winter")
        self.assertEqual(season(35,"2025-02-30"),"unknown")
        self.assertEqual(season(35,""),"unknown")

    def test_source_support_nonresponse_and_previous_color_deficits(self):
        src={("N_temperate","summer"):500,
             ("S_tropical","summer"):120,
             ("N_tropical","winter"):170,
             ("S_temperate","spring"):300}
        old={("N_temperate","summer"):{"attempts":10000,"strict":4300},
             ("S_tropical","summer"):{"attempts":100,"strict":12},
             ("N_tropical","winter"):{"attempts":100,"strict":20},
             ("S_temperate","spring"):{"attempts":500,"strict":150}}
        quotas_out,weights=quotas(src,old,total=256,cap=120,floor=8)
        self.assertEqual(sum(quotas_out.values()),256)
        self.assertGreater(weights[("S_tropical","summer")],
                           weights[("N_temperate","summer")])
        self.assertGreaterEqual(quotas_out[("S_tropical","summer")],8)
        self.assertLessEqual(quotas_out[("S_tropical","summer")],120)
        self.assertGreaterEqual(quotas_out[("S_tropical","summer")],
                                quotas_out[("N_temperate","summer")])
        self.assertEqual(len(rank("sample","A")),64)
        self.assertEqual(len(SCHEMA)>0,True)
        self.assertAlmostEqual(SOURCE_P,.0128)

    def test_disjoint_reference_and_region_season_taxon_balanced_sampling(self):
        with TemporaryDirectory() as root:
            db=sqlite3.connect(":memory:")
            create_db(db)
            # Independent signed-source-style photos and observation UUIDs;
            # every taxon has one photo/observation. Across the two strata
            # species with broad distribution may be chosen again intentionally.
            for j in range(320):
                z="N_temperate" if j<160 else "S_tropical"
                y=40. if j<160 else -10.
                ss="summer" if j<160 else "winter"
                uuid=f"observe-{j}"
                vals=(uuid,str(j), "2025-07-01",str(y),"140.0","12",
                      f"observer-{j%7}","cc-by",
                      str(800000+j),"jpg",rank("photo",uuid,j),
                      z,ss,rank("reference",uuid),rank("independent_observation",uuid),0)
                db.execute(UPSERT,vals)
            old={
                ("N_temperate","summer"):{"attempts":500,"strict":300},
                ("S_tropical","winter"):{"attempts":500,"strict":40},
            }
            chosen,doc=sample(db,old,reference=32,targeted=96)
            self.assertEqual(len(chosen),128)
            self.assertEqual(len({v["uuid"] for v in chosen}),128)
            self.assertEqual(len({v["photo_id"] for v in chosen}),128)
            self.assertEqual(sum(v["cohort"]=="reference" for v in chosen),32)
            self.assertEqual(sum(v["cohort"]=="targeted" for v in chosen),96)
            self.assertGreater(
                sum(v["cohort"]=="targeted" and v["zone"]=="S_tropical" for v in chosen),
                sum(v["cohort"]=="targeted" and v["zone"]=="N_temperate" for v in chosen))
            self.assertEqual(doc["eligible_disjoint_original_photo_observations"],320)
            self.assertAlmostEqual(doc["reference_conditional_probability_per_eligible_observation"],.1)
            db.close()

if __name__=="__main__":
    unittest.main()
