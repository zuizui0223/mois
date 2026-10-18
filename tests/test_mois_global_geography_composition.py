import unittest
import pandas as pd
import geopandas as gpd
from shapely.geometry import box
from scripts.analysis.mois_global_geography_composition import classify_land, community_table, hues

class GeoTests(unittest.TestCase):
    def test_dateline_siblings_not_mislabelled_islands(self):
        land=gpd.GeoDataFrame({'sibling_id':[1,1,2], 'area':[100,5,8]},geometry=[box(-10,-10,10,10),box(170,-10,180,10),box(20,0,22,2)],crs=4326)
        points=pd.DataFrame({'latitude':[0,0,1,50], 'longitude':[0,175,21,50]})
        got=classify_land(points,land,{1})
        self.assertEqual(got.land_type.tolist(),['mainland','mainland','island','unassigned'])
    def test_unknown_never_nonwhite_and_photos_not_species(self):
        d=pd.DataFrame({'photo_id':list('abcdef'),'inat_taxon_id':['1']*4+['2','3'],'species':['A a']*4+['B b','C c'],'morph':['white']*4+['red_pink','mixed_uncertain'],'latitude':[5.1]*6,'longitude':[10.1]*6,'land_type':['mainland']*6,'landmass':['1']*6,'area_km2':[10.]*6})
        for c in ('heat','rain','solar','elevation'):d[c]=20.
        q=community_table(d)
        self.assertEqual(len(q),1)
        self.assertEqual(q.n_photos.iloc[0],5)
        self.assertEqual(q.n_taxa.iloc[0],2)
        self.assertEqual(q.p_white.iloc[0],.5)
        self.assertAlmostEqual(q[['p_'+h for h in hues]].sum(axis=1).iloc[0],1)
if __name__=='__main__':unittest.main()
