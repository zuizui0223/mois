"""Source-independent ecological zone categorization and missingness contracts."""
import unittest
import numpy as np
import pandas as pd

from scripts.analysis.mois_real1043930_global_zonal_colour import (
    latitude_zone,season_label,elev_band,q5,
    describe_stratum,PALETTE,OUTCOMES
)

class GlobalZonalColourTests(unittest.TestCase):
    def test_arctic_and_antarctic_distinct_with_geographic_poles(self):
        self.assertEqual(latitude_zone(66.5622),"N_polar_66.56_90")
        self.assertEqual(latitude_zone(-66.5622),"S_polar_66.56_90")
        self.assertEqual(latitude_zone(-64),"S_subpolar_55_66.56")
        self.assertEqual(latitude_zone(67),"N_polar_66.56_90")
        self.assertEqual(latitude_zone(15),"N_tropical_0_23.5")
        self.assertEqual(latitude_zone(-39),"S_temperate_23.5_55")
        self.assertEqual(latitude_zone(96),"unknown")
        self.assertEqual(latitude_zone(np.nan),"unknown")

    def test_opposite_hemisphere_meteorological_season(self):
        for m,n,s in ((1,"winter","summer"),(2,"winter","summer"),
                      (4,"spring","autumn"),(7,"summer","winter"),
                      (10,"autumn","spring")):
            self.assertEqual(season_label(42,m),n)
            self.assertEqual(season_label(-42,m),s)
        self.assertEqual(season_label(75,7),"summer")
        self.assertEqual(season_label(-70,7),"winter")
        self.assertEqual(season_label(np.nan,7),"unknown")
        self.assertEqual(season_label(30,np.nan),"unknown")

    def test_real_elevation_bands_and_invalid(self):
        self.assertEqual(elev_band(0),"low_0_500m")
        self.assertEqual(elev_band(500),"mid_500_1500m")
        self.assertEqual(elev_band(1500),"high_1500_2500m")
        self.assertEqual(elev_band(2500),"alpine_2500m_plus")
        self.assertEqual(elev_band(-7),"below_sealevel")
        self.assertEqual(elev_band(np.nan),"unknown")

    def test_quintiles_source_complete_and_missing(self):
        out,cut=q5(np.array([1,2,3,4,5,np.nan],dtype=float))
        self.assertEqual(out.tolist(),["Q1","Q2","Q3","Q4","Q5","unknown"])
        self.assertEqual(len(cut),4)

    def test_roi_denominator_and_no_false_white_for_unknowns(self):
        rows=[]
        for pid,known,white,ent,dom in [
            ("1",True,1.,0.,"white"),
            ("2",True,.1,.9,"bronze"),
            ("3",False,np.nan,np.nan,"unknown")]:
            row={"photo_id":pid,"taxon":"13" if pid=="1" else "14",
                 "genus":"Campanula",
                 "strict":known,"ambiguous":False,"classified4":known,
                 "white_fraction":white,"entropy9":ent,"dominant_colour":dom}
            for k in OUTCOMES:
                row.setdefault(k,np.nan if not known else .2)
            rows.append(row)
        d=pd.DataFrame(rows)
        r=describe_stratum(d,"latitude_zone","N_polar_66.56_90")
        self.assertEqual(r["source_photo_attempts"],3)
        self.assertEqual(r["strict_original9_ROIs"],2)
        self.assertAlmostEqual(r["strict9_success_pct_of_all_attempts"],100*2/3,places=8)
        self.assertAlmostEqual(r["white_fraction_photo_mean_pct"],55)
        self.assertEqual(r["share_dominant_white_pct"],50)
        self.assertEqual(r["share_dominant_bronze_pct"],50)
        self.assertNotIn("share_dominant_green_pct",r)
        self.assertIn("n_dominant_bronze",r)

if __name__=="__main__":
    unittest.main()
