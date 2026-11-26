import unittest
import pandas as pd
from scripts.analysis.mois_global_identity_consensus_audit import canonicalize,namespace_bounds
U1='00000000-0000-4000-8000-000000000001'
U2='00000000-0000-4000-8000-000000000002'
def source(pid='1',old='101',tax='9',uid=U1,newtax='9'):
    return {'original_photo_id':pid,'original_numeric_observation_id':old,'original_taxon_id':tax,'snapshot_observation_uuid':uid,'snapshot_taxon_id':newtax}
class TestGlobalIdentity(unittest.TestCase):
    def test_same_photo_source_collapses_id_alias_without_deleting_photo(self):
        p=pd.DataFrame({'photo_id':['1','2','3'],'observation_id':['101',U1,U2],'inat_taxon_id':['9']*3,'strict9':[True]*3})
        r=canonicalize(p,pd.DataFrame([source()]))
        self.assertEqual(r.observation_id.tolist(),[U1,U1,U2]);self.assertEqual(len(r),3)
        self.assertEqual(r.source_observation_id.tolist(),['101',U1,U2])
    def test_taxonomic_disagreement_held(self):
        p=pd.DataFrame({'photo_id':['1'],'observation_id':['101'],'inat_taxon_id':['9']})
        self.assertEqual(canonicalize(p,pd.DataFrame([source(newtax='10')])).observation_id.tolist(),['101'])
    def test_old_numeric_id_disagreement_held(self):
        p=pd.DataFrame({'photo_id':['1'],'observation_id':['999'],'inat_taxon_id':['9']})
        self.assertEqual(canonicalize(p,pd.DataFrame([source()])).observation_id.tolist(),['999'])
    def test_uuid_and_missing_not_reassigned(self):
        p=pd.DataFrame({'photo_id':['1','2'],'observation_id':[U2,''],'inat_taxon_id':['9','9']})
        self.assertEqual(canonicalize(p,pd.DataFrame([source()])).observation_id.tolist(),[U2,''])
    def test_old_numeric_id_with_two_snapshot_uuids_held(self):
        p=pd.DataFrame({'photo_id':['1','2'],'observation_id':['101','101'],'inat_taxon_id':['9','9']})
        self.assertEqual(canonicalize(p,pd.DataFrame([source(),source(pid='2',uid=U2)])).observation_id.tolist(),['101','101'])
    def test_duplicated_source_photo_rejected(self):
        p=pd.DataFrame({'photo_id':['1'],'observation_id':['101'],'inat_taxon_id':['9']})
        with self.assertRaises(ValueError):canonicalize(p,pd.DataFrame([source(),source()]))
    def test_two_plus_two_not_claimed_three_independent_observations(self):
        p=pd.DataFrame({'inat_taxon_id':['9']*4,'observation_id':['101','102',U1,U2]})
        r=namespace_bounds(p).loc['9']
        self.assertEqual(r.independent_lower_bound,2);self.assertEqual(r.independent_upper_bound,4)
    def test_three_in_one_namespace_survive_all_cross_namespace_aliases(self):
        p=pd.DataFrame({'inat_taxon_id':['9']*5,'observation_id':['101','102','103',U1,U2]})
        self.assertEqual(namespace_bounds(p).loc['9','independent_lower_bound'],3)
if __name__=='__main__':unittest.main()
