import unittest
import pandas as pd
from scripts.analysis.mois_experimental_stage2 import controlled_ratios

class ControlContrasts(unittest.TestCase):
    def test_negative_contrast_retained_and_missing_control_not_zero(self):
        raw,net=controlled_ratios(pd.Series([1.,1.,0.]),pd.Series([3.,3.,0.]),pd.Series([2.,float('nan'),0.]))
        self.assertEqual(raw.iloc[0],.25);self.assertEqual(net.iloc[0],-.25)
        self.assertTrue(pd.isna(raw.iloc[1]));self.assertTrue(pd.isna(net.iloc[2]))

if __name__=='__main__':unittest.main()
