import unittest
import numpy as np
import pandas as pd
from scripts.analysis.mois_moisture_bias_refit_20261030 import strata_weights,within_block_sums,solve_bootstrap,validate_frame

class MoistureRefitTests(unittest.TestCase):
    def frame(self):
        rows=[]
        for h,k,white in [('a',2,2),('b',8,0),('c',0,0)]:
            for i in range(10):
                rows.append({'photo_id':str(len(rows)+1),'observation_id':'obs'+str(len(rows)+1),'stratum':h,'block':h,'q':0,'strict9':i<k,'original_roi_four_hue':('white' if i<white else 'red_pink'),'strict9_white':float(i<white) if i<k else np.nan})
        return pd.DataFrame(rows)
    def test_weights_use_all_attempts_and_leave_unsupported_missing(self):
        d=self.frame();table,w=strata_weights(d)
        self.assertEqual(int(table.N.sum()),30)
        self.assertEqual(int(table.K.sum()),10)
        self.assertEqual(int(table.loc[table.K.eq(0),'N'].sum()),10)
        self.assertAlmostEqual(float(w[d.strict9].sum()),20)
        y=d.original_roi_four_hue.eq('white').to_numpy()
        self.assertAlmostEqual(np.average(y[d.strict9],weights=w[d.strict9]),.5)
        self.assertTrue(np.isnan(w[~d.strict9]).all())
    def test_unknown_colour_text_cannot_change_weights_or_observed_totals(self):
        d=self.frame();a,w=strata_weights(d)
        d.loc[~d.strict9,'original_roi_four_hue']='white'
        b,v=strata_weights(d)
        pd.testing.assert_frame_equal(a,b)
        np.testing.assert_allclose(w,v,equal_nan=True)
    def test_duplicate_source_ids_fail_closed(self):
        d=self.frame();validate_frame(d)
        d.loc[0,'photo_id']=d.loc[1,'photo_id']
        with self.assertRaises(ValueError):validate_frame(d)
    def test_nonboolean_roi_status_fails_closed(self):
        d=self.frame();d.strict9=d.strict9.astype(str)
        with self.assertRaises(ValueError):validate_frame(d)
    def test_known_within_slope_survives_group_intercepts_and_unequal_weights(self):
        x=np.tile(np.array([-1.,0.,1.,2.]),6)[:,None]
        groups=np.repeat(np.arange(6),4);blocks=groups//2
        y=(2*x[:,0]+groups*10)[:,None]
        A,b=within_block_sums(x,y,np.arange(len(x))+1,groups,blocks)
        slope=np.linalg.solve(A.sum(axis=0),b.sum(axis=0))
        np.testing.assert_allclose(slope,[[2.]],atol=1e-10)
        draws=np.array([[3,0,0],[0,3,0],[0,0,3],[1,1,1]])
        boot=solve_bootstrap(A,b,draws)
        np.testing.assert_allclose(boot,2.,atol=1e-10)
    def test_weighted_within_matches_explicit_dummy_regression(self):
        rng=np.random.default_rng(77);groups=np.repeat(np.arange(8),5)
        x=rng.normal(size=(40,3));y=(x@np.array([.2,-.3,.5])+groups+rng.normal(size=40)*.1)[:,None]
        w=rng.uniform(.5,8,size=40);blocks=groups//2
        A,b=within_block_sums(x,y,w,groups,blocks)
        actual=np.linalg.solve(A.sum(0),b.sum(0))[:,0]
        design=np.c_[x,np.eye(8)[groups]]
        ref=np.linalg.lstsq(design*np.sqrt(w)[:,None],y*np.sqrt(w)[:,None],rcond=None)[0][:3,0]
        np.testing.assert_allclose(actual,ref,atol=1e-10)

if __name__=='__main__':unittest.main()
