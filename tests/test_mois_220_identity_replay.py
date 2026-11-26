import unittest
from scripts.analysis.mois_220_identity_annotations import link
from scripts.analysis.mois_220_identity_replay import canonicalize
class TestCanonicalize(unittest.TestCase):
    def test_unqueried_uuid_is_source_identity_not_live_annotation(self):
        uid='c5763598-2e2e-40cc-a912-a7eb3bf11e1a'
        rows=[dict(photo_id='77',observation_id=uid,inat_taxon_id='5',review_order='1')]
        result=canonicalize(link(rows,[],{}))
        self.assertEqual(result[0]['canonical_observation_uuid'],uid)
        self.assertEqual(result[0]['canonical_identity_evidence'],'source_uuid_only')
        self.assertEqual(result[0]['record_status'],'not_fetched_uuid_batch_rejected')
        self.assertEqual(result[0]['photo_flower_stage'],'unreviewed')
    def test_alias_pair_keeps_both_images_one_independent_observation(self):
        uid='c5763598-2e2e-40cc-a912-a7eb3bf11e1a'
        rows=[dict(photo_id='77',observation_id='123',inat_taxon_id='5',review_order='1'),dict(photo_id='78',observation_id=uid,inat_taxon_id='5',review_order='2')]
        obs=[dict(id=123,uuid=uid,taxon={'id':5},photos=[{'id':77},{'id':78}])]
        result=canonicalize(link(rows,obs,{}))
        self.assertEqual(len(result),2)
        self.assertEqual(sum(r['independent_representative'] for r in result),1)
        self.assertEqual(result[1]['represented_by_photo_id'],'77')
    def test_unresolved_numeric_id_never_claimed_canonical_uuid(self):
        row=dict(photo_id='77',observation_id='123',inat_taxon_id='5',review_order='1')
        with self.assertRaisesRegex(ValueError,'numeric'):canonicalize(link([row],[],{}))
if __name__=='__main__':unittest.main()
