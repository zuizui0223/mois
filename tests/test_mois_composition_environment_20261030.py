import unittest
import numpy as np
import pandas as pd
from scripts.analysis import mois_composition_environment_20261030 as m

class CompositionTests(unittest.TestCase):
    def test_palette_unknown_is_not_a_dull_or_white_flower(self):
        d=pd.DataFrame({'strict9':[True,False]})
        for c in m.PALETTE:d['strict9_'+c]=[float(c=='white'),np.nan]
        out=m.palette_traits(d)
        self.assertEqual(out.loc[0,'palette_entropy'],0)
        self.assertEqual(out.loc[0,'palette_dominance'],1)
        self.assertTrue(out.loc[1,list(m.METRICS)].isna().all())
        self.assertNotIn('inconspicuousness',out.columns)
    def test_nonmissing_colour_in_unknown_fails(self):
        d=pd.DataFrame({'strict9':[False]})
        for c in m.PALETTE:d['strict9_'+c]=0.
        with self.assertRaises(ValueError):m.palette_traits(d)
    def test_composition_must_be_closed(self):
        d=pd.DataFrame({'strict9':[True]})
        for c in m.PALETTE:d['strict9_'+c]=.2
        with self.assertRaises(ValueError):m.palette_traits(d)
    def test_exact_fourterm_decomposition(self):
        a=pd.DataFrame({'taxon':['x','y'],'value':[.8,.4],'n':[1,3]}).set_index('taxon')
        b=pd.DataFrame({'taxon':['x','z'],'value':[.6,.2],'n':[3,1]}).set_index('taxon')
        r=m.decompose(a,b,'value')
        self.assertAlmostEqual(r['taxon_equal_difference_pp'],20)
        self.assertAlmostEqual(r['shared_taxon_outcome_pp'],10)
        self.assertAlmostEqual(r['shared_taxon_weight_pp'],0)
        self.assertAlmostEqual(r['focal_only_taxa_pp'],20)
        self.assertAlmostEqual(r['reference_only_taxa_pp'],-10)
        self.assertAlmostEqual(r['common_taxa_equal_difference_pp'],20)
        self.assertEqual(r['shared_taxa'],1)
    def test_no_overlap_does_not_invent_shared_traits(self):
        a=pd.DataFrame({'value':[.8],'n':[1]},index=['x'])
        b=pd.DataFrame({'value':[.2],'n':[1]},index=['z'])
        r=m.decompose(a,b,'value')
        self.assertIsNone(r['common_taxa_equal_difference_pp'])
        self.assertEqual(r['shared_taxa'],0)
        self.assertAlmostEqual(r['taxon_equal_difference_pp'],60)
    def test_join_rejects_taxon_drift(self):
        a=pd.DataFrame({'photo_id':['1'],'observation_id':['o'],'inat_taxon_id':['t'],'strict9':[True],'strict9_white':[1.]})
        b=a.copy();b.loc[0,'inat_taxon_id']='wrong'
        for c in m.PALETTE:b['strict9_'+c]=float(c=='white')
        with self.assertRaises(ValueError):m.join_palette(a,b)
    def test_join_rejects_duplicate_photo(self):
        a=pd.DataFrame({'photo_id':['1'],'observation_id':['o'],'inat_taxon_id':['t'],'strict9':[True],'strict9_white':[1.]})
        b=pd.concat([a,a],ignore_index=True)
        for c in m.PALETTE:b['strict9_'+c]=float(c=='white')
        with self.assertRaises(ValueError):m.join_palette(a,b)
    def test_join_preserves_photo_order_and_unknown(self):
        a=pd.DataFrame({'photo_id':['2','1'],'observation_id':['p','o'],'inat_taxon_id':['t','s'],'strict9':[False,True],'strict9_white':[np.nan,1.]})
        b=a.iloc[::-1].copy()
        for c in m.PALETTE:b['strict9_'+c]=[float(c=='white'),np.nan]
        out=m.join_palette(a,b)
        self.assertEqual(out.photo_id.tolist(),['2','1'])
        self.assertTrue(np.isnan(out.loc[0,'strict9_bronze']))
    def test_within_regression_matches_explicit_dummy_fit(self):
        rng=np.random.default_rng(20261030);n=200
        g=np.repeat(np.arange(10),20);blocks=np.repeat(np.arange(5),40)
        x=rng.normal(size=(n,2));y=np.c_[x@np.array([2.,-.3])+g, x@np.array([-2.,.3])-g]
        w=rng.uniform(.2,2.,n)
        r=m.within_fit(x,y,w,g,blocks)
        dummy=np.eye(10)[g];X=np.c_[x,dummy]
        expected=np.linalg.lstsq(X*np.sqrt(w[:,None]),y*np.sqrt(w[:,None]),rcond=None)[0][:2]
        np.testing.assert_allclose(r['beta'],expected,atol=1e-10)
        np.testing.assert_allclose(r['beta'].sum(1),0,atol=1e-10)
    def test_geographic_definition_not_south_polar_pooling(self):
        lat=pd.Series([66.5622,66.5,55.,54.9,-55.,-40.,-10.,np.nan])
        masks=m.region_masks(pd.DataFrame({'latitude':lat,'land_type':['unassigned']*8}))
        self.assertEqual(masks['arctic'].tolist(),[True,False,False,False,False,False,False,False])
        self.assertTrue(masks['south_temperate'].iloc[5])
        self.assertFalse(masks['south_temperate'].iloc[6])

if __name__=='__main__':unittest.main()
