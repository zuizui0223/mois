import unittest
import pandas as pd
from scripts.acquisition import mois_global_acquisition_20261106 as m

class ScopeTests(unittest.TestCase):
 def test_no_lower_island_area(self):
  self.assertEqual(m.land_scope(.0001),'island')
  self.assertEqual(m.land_scope(4.99),'island')
  self.assertEqual(m.land_scope(5),'island')
  self.assertEqual(m.land_scope(8e6),'mainland')
  self.assertEqual(m.land_scope(float('nan')),'unresolved')
 def test_acquisition_retains_northern_small_island_unknown_colours(self):
  f=pd.DataFrame({'photo_id':['1','2','3'],'observation_id':['a','b','c'],
   'observed_on':['2024-01-01']*3,'latitude':[78,10,None],'longitude':[11,20,None],
   'strict9':[False,True,False]})
  t,r=m.global_cells(f,2024)
  self.assertEqual(len(t),2)
  self.assertEqual(r['invalid_coordinates_retained'],1)
  self.assertEqual(t.n_attempts.sum(),2)
 def test_duplicate_canonical_photo_refused(self):
  f=pd.DataFrame({'photo_id':['1','1'],'observation_id':['a','b'],
   'observed_on':['2024-01-01']*2,'latitude':[10,10],'longitude':[10,10]})
  with self.assertRaises(ValueError):m.global_cells(f,2024)
 def test_geography_not_an_eligibility_gate(self):
  self.assertEqual(m.CONFIG['minimum_island_area_km2'],None)
  self.assertTrue(m.CONFIG['include_mainland'])
  self.assertTrue(m.CONFIG['retain_unresolved_geography'])
 def test_zero_radiance_not_missing(self):
  j={'value':0,'properties':{'Values':['0 4']},
   'catalogItems':{'features':[{'attributes':{'year':2024,'month':1}}]}}
  self.assertEqual(m.parse_viirs(j,2024,1),(0.,4.,'ok'))
 def test_wrong_month_refused(self):
  j={'value':2,'properties':{'Values':['2 4']},
   'catalogItems':{'features':[{'attributes':{'year':2020,'month':1}}]}}
  self.assertEqual(m.parse_viirs(j,2024,1)[2],'wrong_scene_date')
 def test_photo_dedup_and_license_gate(self):
  rows=[{'photo_id':'1','observation_uuid':'a','license_code':'cc-by','original_url':'https://inaturalist-open-data.s3.amazonaws.com/photos/1/original.jpg'},
   {'photo_id':'2','observation_uuid':'b','license_code':'cc-by','original_url':'https://inaturalist-open-data.s3.amazonaws.com/photos/2/original.jpg'},
   {'photo_id':'3','observation_uuid':'c','license_code':None,'original_url':'https://static.inaturalist.org/photos/3/original.jpg'}]
  good=m.unseen_photos(rows,{'1'},{'a'},set())
  self.assertEqual([r['photo_id'] for r in good],['2'])
 def test_no_extra_photo_retries(self):
  self.assertEqual(m.CONFIG['photo_get_attempts'],1)
  self.assertEqual(m.CONFIG['new_strict9_automatically_claimed'],0)

if __name__=='__main__':unittest.main()
