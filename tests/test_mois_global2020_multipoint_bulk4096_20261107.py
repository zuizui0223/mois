"""MOIS 2020 next4096 bulk extraction contracts; no false GIS years or zero fill."""
import unittest
from scripts.analysis import mois_global2020_multipoint_bulk4096_20261107 as b

def fake(n=1024):
    return [{'target_id':f'2020|{i}_0|1','land_class':'mainland',
             'query_longitude':str(-170+i*.01),'query_latitude':'30',
             'n_attempted':'1','n_strict9':'0','source_subarms':'historical_main'}
             for i in range(n)]

class Test4096Bulk(unittest.TestCase):
    def test_4096_make_64_provider_requests(self):
        plans=b.plan_requests(fake(4096))
        self.assertEqual(len(plans),64)
        self.assertEqual(sum(p['n'] for p in plans if p['product']=='viirs'),4096*4)
        self.assertEqual(sum(p['n'] for p in plans if p['product']=='population'),4096)
        self.assertEqual(max(p['n'] for p in plans if p['product']=='viirs'),512)
        self.assertEqual(max(p['n'] for p in plans if p['product']=='population'),128)
    def test_remaining30630_uses_bounded_512_viirs_and_128_population(self):
        plans=b.plan_requests(fake(30630))
        self.assertEqual(len(plans),480)
        self.assertEqual(sum(p['n'] for p in plans if p['product']=='viirs'),30630*4)
        self.assertEqual(sum(p['n'] for p in plans if p['product']=='population'),30630)
        self.assertEqual(max(p['n'] for p in plans if p['product']=='viirs'),512)
        self.assertEqual(max(p['n'] for p in plans if p['product']=='population'),128)
    def test_no_duplicate_target_and_coordinates(self):
        x=fake(4);x.append(x[1].copy())
        with self.assertRaisesRegex(ValueError,'duplicate'):
            b.validate_source_targets(x)
    def test_bad_period_fails_before_api(self):
        x=fake(4);x[0]['target_id']='2024|0_0|1'
        with self.assertRaisesRegex(ValueError,'2020'):
            b.validate_source_targets(x)
    def test_both_mainland_and_island_retained(self):
        x=fake(4);x[0]['land_class']='island'
        self.assertEqual(len(b.validate_source_targets(x)),4)
    def test_summary_success_after_3_months(self):
        statuses={1:('OK',1.),4:('OK',3.),7:('OK',2.),10:('LOW_CLOUD_COUNT',None)}
        result=b.summarize_one(statuses,('OK',0.))
        self.assertTrue(result['complete'])
        self.assertEqual(result['viirs_valid_months'],3)
        self.assertEqual(result['viirs_median'],2.)
        self.assertEqual(result['population_density'],0.)
    def test_missing_population_is_not_zero(self):
        statuses={m:('OK',5.) for m in (1,4,7,10)}
        result=b.summarize_one(statuses,('NO_DATA',None))
        self.assertFalse(result['complete'])
        self.assertIsNone(result['population_density'])
    def test_two_valid_months_not_eligible(self):
        statuses={1:('OK',1.),4:('LOW_CLOUD_COUNT',None),7:('NO_DATA',None),10:('OK',2.)}
        self.assertFalse(b.summarize_one(statuses,('OK',5.))['complete'])
    def test_wrong_source_year_and_metadata_fails_closed(self):
        spec=b.provider_params('viirs',1,fake(2),2020)
        self.assertIn('"points"',spec['geometry'])
        with self.assertRaisesRegex(ValueError,'2020'):
            b.provider_params('viirs',1,fake(2),2024)
    def test_raw_band_numeric_zero_not_missing(self):
        from scripts.analysis import mois_global_multipoint_numeric_20261107 as old
        z={'samples':[{'locationId':0,'location':{'x':-170.,'y':30.},'value':0}]}
        v=old.extract_samples('population',0,fake(2),z)
        self.assertEqual(v[0]['status'],'OK')
        self.assertEqual(v[0]['value'],0.)
if __name__=='__main__':
    unittest.main()
