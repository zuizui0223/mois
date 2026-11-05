import unittest
from unittest.mock import patch
from scripts.acquisition import mois_diel_environment as m

class EnvironmentalContracts(unittest.TestCase):
    def test_population_duplicate_cannot_replace_missing_year(self):
        response={'features':[{'id':'713_568','area_km2':100,'years':[{'data_year':2020,'population_density':4},{'data_year':2020,'population_density':5}]}]}
        with patch.object(m,'cached_request',return_value=response):
            with self.assertRaises(ValueError):m.population_batch(['713_568'],[2020,2024],'unused')
    def test_population_exact_requested_keys(self):
        response={'features':[{'id':'713_568','area_km2':100,'years':[{'data_year':2024,'population_density':0},{'data_year':2020,'population_density':5}]}]}
        with patch.object(m,'cached_request',return_value=response):
            rows=m.population_batch(['713_568'],[2020,2024],'unused')
        self.assertEqual({(r['cell025'],r['year']) for r in rows},{('713_568',2020),('713_568',2024)})
        self.assertEqual(rows[0]['population_density'],0)
    def sample(self,year=2019,month=7,cloud=5):
        return {'locationId':0,'location':{'x':10.,'y':20.},'value':f'2 {cloud}','attributes':{'year':year,'month':month}}
    def test_wrong_viirs_year_fails(self):
        with self.assertRaises(ValueError):m.parse_sample(self.sample(2018),2019,7,'viirs',10.,20.)
    def test_low_cloud_count_stays_missing(self):
        d=m.parse_sample(self.sample(cloud=1),2019,7,'viirs',10.,20.)
        self.assertEqual(d['status'],'LOW_CLOUD_FREE');self.assertIsNone(d['value'])
    def test_wrong_point_fails(self):
        with self.assertRaises(ValueError):m.parse_sample(self.sample(),2019,7,'viirs',11.,20.)
    def test_zero_radiance_with_valid_cloud_is_real_zero(self):
        d=self.sample();d['value']='0 5'
        self.assertEqual(m.parse_sample(d,2019,7,'viirs',10.,20.)['value'],0.)
    def test_population_returned_slice_must_match(self):
        d=self.sample();d['value']='20';d['attributes']={'StdTime':m.year_ms(2020)}
        with self.assertRaises(ValueError):m.parse_sample(d,2019,None,'population',10.,20.)

if __name__=='__main__':unittest.main()
