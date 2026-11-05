import unittest
import numpy as np
import pandas as pd
from scripts.analysis.mois_macro_fit import fit,absorb

class RegressionContracts(unittest.TestCase):
    def test_absorbed_weighted_coefficients_match_explicit_dummies(self):
        rng=np.random.default_rng(123);n=600
        a=rng.integers(0,12,n);b=rng.integers(0,5,n);x=rng.normal(size=(n,2));w=rng.uniform(.2,2,n)
        y=x@np.array([2.,-3.])+a/2+b+rng.normal(size=n)
        dummy=pd.get_dummies(pd.DataFrame({'a':a.astype(str),'b':b.astype(str)}),drop_first=True,dtype=float).to_numpy()
        design=np.column_stack([x,np.ones(n),dummy])
        expected=np.linalg.lstsq(design*np.sqrt(w[:,None]),y*np.sqrt(w),rcond=None)[0][:2]
        actual=fit(y,x,[a,b],[a],w)[0]
        np.testing.assert_allclose(actual,expected,atol=1e-8)
    def test_constant_within_fixed_effect_is_unidentified(self):
        groups=np.repeat(np.arange(10),4)
        with self.assertRaises(ValueError):fit(np.arange(40),groups[:,None].astype(float),[groups],[groups])

if __name__=='__main__':unittest.main()
