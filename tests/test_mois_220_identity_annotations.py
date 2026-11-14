import unittest
from scripts.analysis import mois_220_identity_annotations as identity
class TestIdentity(unittest.TestCase):
    def test_separate_id_namespaces(self):
        self.assertEqual(identity.kind('12345'),'numeric_id')
        self.assertEqual(identity.kind('c5763598-2e2e-40cc-a912-a7eb3bf11e1a'),'uuid')
        with self.assertRaises(ValueError):identity.kind('nan')
    def test_aliases_canonicalized_and_photo_verified(self):
        uid='c5763598-2e2e-40cc-a912-a7eb3bf11e1a'
        rows=[dict(photo_id='77',observation_id='123',inat_taxon_id='5'),dict(photo_id='78',observation_id=uid,inat_taxon_id='5')]
        obs=[dict(id=123,uuid=uid,taxon={'id':5},photos=[{'id':77},{'id':78}],annotations=[{'controlled_attribute_id':12,'controlled_value_id':13}])]
        result=identity.link(rows,obs,{(12,13):('Plant Phenology','Flowering')})
        self.assertEqual([r['canonical_observation_uuid'] for r in result],[uid,uid])
        self.assertTrue(all(r['photo_membership']=='matched' for r in result))
        self.assertTrue(all(r['photo_flower_stage']=='unreviewed' for r in result))
        self.assertIn('Flowering',result[0]['observation_annotations'])
    def test_missing_record_not_absent_flower(self):
        result=identity.link([dict(photo_id='77',observation_id='123',inat_taxon_id='5')],[],{})
        self.assertEqual(result[0]['record_status'],'not_returned')
        self.assertEqual(result[0]['photo_flower_stage'],'unreviewed')
    def test_wrong_photo_never_accepted(self):
        obs=[dict(id=123,uuid='c5763598-2e2e-40cc-a912-a7eb3bf11e1a',taxon={'id':5},photos=[{'id':99}])]
        result=identity.link([dict(photo_id='77',observation_id='123',inat_taxon_id='5')],obs,{})
        self.assertEqual(result[0]['photo_membership'],'not_in_current_record')
if __name__=='__main__':unittest.main()
