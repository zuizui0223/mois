import unittest
import numpy as np
import pandas as pd
from scripts.analysis.mois_gbif_continent_two_cohort_qa_20261115 import bh,chi_stat,analysis,HUES,REGIONS

class TestGBIFContinentReplicate(unittest.TestCase):
    def test_BH_monotonic(self):
        q=bh([0.01,0.10,0.02,0.80])
        self.assertTrue(np.allclose(q,[.04,.13333333333333333,.04,.80]))
    def test_mixture_has_positive_statistic_and_equal_has_zero(self):
        g=np.repeat(np.arange(9),20)
        y=np.tile(np.array([0,1]),90)
        self.assertAlmostEqual(chi_stat(g,y),0.0)
        y2=np.repeat([1]*3+[0]*6,20)
        self.assertGreater(chi_stat(g,y2),100.)
    def test_two_cohorts_retained_and_all_regions_tested(self):
        rows=[]
        for cohort in ("discovery","holdout"):
            for h in HUES:
                for i in range(20):
                    row={'cohort':cohort,'gbif_colxr_species_key':cohort+h+str(i),
                         'photo_hue9':h}
                    for c in REGIONS:row[c]=int(i%2==0)
                    rows.append(row)
        d=pd.DataFrame(rows)
        result=analysis(d,nperm=1000,seed=42)
        self.assertEqual(len(result),12)
        self.assertTrue(result.p_perm_10000.eq(1.0).all())
        self.assertTrue(result.q_BH_in_6_continents.eq(1.0).all())
    def test_duplicate_species_rejected(self):
        d=pd.DataFrame([{'gbif_colxr_species_key':'one','cohort':'discovery'}]*360)
        with self.assertRaisesRegex(ValueError,"two disjoint"):
            analysis(d,nperm=1000)

if __name__=="__main__":unittest.main()
