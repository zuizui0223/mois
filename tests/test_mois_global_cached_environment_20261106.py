import json,tempfile,unittest
from pathlib import Path
from scripts.analysis import mois_global_cached_environment_20261106 as c

def viirs(value=0,year=2020,month=1,cloud=4):
    return {'value':value,'catalogItems':{'features':[{'attributes':{'year':year,'month':month}}]},'properties':{'Values':[f'{value} {cloud}']}}

class CollectorTests(unittest.TestCase):
    def test_zero_is_a_valid_radiance(self):
        self.assertEqual(c.parse(viirs(0),'viirs',2020,1),(0.0,'OK'))
    def test_negative_raw_radiance_is_preserved(self):
        self.assertEqual(c.parse(viirs(-.1),'viirs',2020,1),(-.1,'OK'))
    def test_wrong_source_month_is_not_measurement(self):
        self.assertEqual(c.parse(viirs(month=4),'viirs',2020,1)[1],'WRONG_YEAR_MONTH')
    def test_too_few_cloud_free_obs_remains_missing(self):
        self.assertEqual(c.parse(viirs(cloud=1),'viirs',2020,1),(None,'LOW_CLOUD_COUNT'))
    def test_pop_nodata_is_not_zero(self):
        self.assertEqual(c.parse({'value':'NoData'},'population_global1',2020,0),(None,'NO_DATA'))
    def test_population_zero_is_valid(self):
        self.assertEqual(c.parse({'value':0},'population_global1',2020,0),(0.0,'OK'))
    def test_wrong_population_year_is_not_substituted(self):
        with self.assertRaises(ValueError):c.spec('population_global1',2024,0,0,0)
    def test_request_identity_includes_source_time_and_position(self):
        a=c.spec('viirs',2020,1,140,65)
        self.assertNotEqual(c.key(a),c.key(c.spec('viirs',2020,4,140,65)))
        self.assertNotEqual(c.key(a),c.key(c.spec('viirs',2020,1,141,65)))
        self.assertNotEqual(c.key(a),c.key(c.spec('population_global1',2020,0,140,65)))
    def test_polar_coordinates_are_retained(self):
        self.assertEqual(c.spec('viirs',2020,1,0,82)['latitude'],82.0)
    def test_cache_reuses_missing_and_valid_without_http(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d);s=c.spec('viirs',2020,1,0,0)
            c.save_response(folder,s,json.dumps({'value':'NoData'}).encode(),200)
            value=c.read_cache(folder,s)
            self.assertEqual(value['status'],'NO_DATA');self.assertIsNone(value['value'])
    def test_cache_checks_exact_original_bytes(self):
        with tempfile.TemporaryDirectory() as d:
            folder=Path(d);s=c.spec('viirs',2020,1,0,0)
            c.save_response(folder,s,json.dumps(viirs()).encode(),200)
            (folder/(c.key(s)+'.response')).write_text('tampered')
            with self.assertRaises(ValueError):c.read_cache(folder,s)
    def test_auth_error_stops_instead_of_becoming_nodata(self):
        with tempfile.TemporaryDirectory() as d:
            s=c.spec('viirs',2020,1,0,0)
            r=c.save_response(Path(d),s,b'{"error":{"code":499}}',200)
            self.assertEqual(r['status'],'AUTH_OR_RATE_LIMIT')
    def test_target_completeness_requires_three_valid_months_and_population(self):
        obs=[{'product':'viirs','status':'OK','value':i} for i in (1,2,3)]
        obs += [{'product':'viirs','status':'NO_DATA','value':None},{'product':'population_global1','status':'OK','value':0}]
        self.assertTrue(c.summary(obs)['complete']);self.assertEqual(c.summary(obs)['viirs_median'],2)
        obs[-1]['status']='NO_DATA';obs[-1]['value']=None
        self.assertFalse(c.summary(obs)['complete'])
    def test_attempted_record_is_not_success(self):
        self.assertFalse(c.summary([{'product':'viirs','status':'REQUEST_FAILED','value':None}])['complete'])

if __name__=='__main__':unittest.main()
