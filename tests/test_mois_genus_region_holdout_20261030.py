"""Tests for independent accepted-genus and heldout species provenance."""
import unittest
import numpy as np
import pandas as pd
from scripts.analysis import mois_genus_region_holdout_20261030 as m

class TaxonomyHoldoutTests(unittest.TestCase):
    def setUp(self):
        self.names=pd.DataFrame({
            'inat_taxon_id':['1','2','3','4','5','6','7','8'],
            'wcvp_match_state':['unique_exact_name_to_accepted_id']*6+['multiple_accepted_ids_ambiguous','unique_exact_name_to_accepted_id'],
            'wcvp_accepted_plant_name_id':['101','102','102','103','104','105','999','106']})
        self.accepted=pd.DataFrame({
            'wcvp_accepted_plant_name_id':['101','102','103','104','105','106'],
            'wcvp_accepted_name':['Genusa alpha','Genusa beta','Genusa gamma','Genusb delta','Genusb epsilon','Genusc zeta']})
    def test_taxonomy_drops_ambiguous_and_proves_accepted_ids(self):
        mapping,meta=m.accepted_taxonomy(self.names,self.accepted)
        self.assertEqual(len(mapping),7)
        self.assertNotIn('7',set(mapping.inat_taxon_id))
        self.assertEqual(mapping.set_index('inat_taxon_id').loc['8','genus'],'Genusc')
        self.assertEqual(meta['rejected_ambiguous_names'],1)
    def test_taxonomy_rejects_conflicting_accepted_id_names(self):
        bad=pd.concat([self.accepted,self.accepted.iloc[[0]].assign(wcvp_accepted_name='Wronggenus species')])
        with self.assertRaisesRegex(ValueError,'Conflicting accepted names'):
            m.accepted_taxonomy(self.names,bad)
    def test_synonym_original_taxa_collapse_with_observation_weights(self):
        mapping,_=m.accepted_taxonomy(self.names,self.accepted)
        t=pd.DataFrame({'sector':['arctic']*3,'inat_taxon_id':['1','2','3'],'n':[1,3,1],
                        'strict9_white':[.2,.8,.0],'strict9_bronze':[0.,.1,.2]})
        species,cov=m.collapse_to_accepted_species(t,mapping,['strict9_white','strict9_bronze'])
        self.assertEqual(len(species),2)
        self.assertAlmostEqual(species.set_index('accepted_id').loc['102','strict9_white'],.6)
        self.assertEqual(species.set_index('accepted_id').loc['102','n'],4)
        self.assertEqual(cov['ambiguous_or_unmapped_input_rows'],0)
    def test_missing_taxonomy_not_silently_mapped(self):
        mapping,_=m.accepted_taxonomy(self.names,self.accepted)
        t=pd.DataFrame({'sector':['arctic']*2,'inat_taxon_id':['1','7'],'n':[1,1],'strict9_white':[.2,.9]})
        species,cov=m.collapse_to_accepted_species(t,mapping,['strict9_white'])
        self.assertEqual(len(species),1)
        self.assertEqual(cov['ambiguous_or_unmapped_input_rows'],1)
    def test_holdout_requires_three_distinct_reference_species(self):
        species=pd.DataFrame({'sector':['focal','focal','reference','reference','reference','reference'],
            'accepted_id':['10','11','1','2','3','4'],'genus':['A','B','A','A','A','B'],
            'n':[2,1,1,1,1,1],'strict9_white':[.4,.8,.2,.3,.4,.8]})
        r,scored=m.score_holdout(species,'focal','reference',['strict9_white'],3)
        self.assertEqual(r['supported_novel_species'],1)
        self.assertEqual(r['novel_target_species'],2)
        self.assertEqual(scored.accepted_id.tolist(),['10'])
        self.assertAlmostEqual(scored.pred_strict9_white.iloc[0],.3)
    def test_same_species_in_train_never_becomes_novel(self):
        species=pd.DataFrame({'sector':['focal','focal','reference','reference','reference','reference'],
            'accepted_id':['1','10','1','2','3','4'],'genus':['A']*6,'n':[1]*6,'strict9_white':[.9,.4,.1,.2,.3,.4]})
        r,scored=m.score_holdout(species,'focal','reference',['strict9_white'],3)
        self.assertEqual(r['novel_target_species'],1)
        self.assertEqual(scored.accepted_id.tolist(),['10'])
    def test_no_support_has_null_skill_and_no_fake_predictions(self):
        species=pd.DataFrame({'sector':['focal','reference'],'accepted_id':['1','2'],'genus':['A','B'],'n':[1,1],'strict9_white':[.4,.5]})
        r,s=m.score_holdout(species,'focal','reference',['strict9_white'],3)
        self.assertEqual(r['supported_novel_species'],0)
        self.assertEqual(len(s),0)
        self.assertIsNone(r['metrics']['strict9_white']['mse_skill_vs_reference_global_mean'])
    def test_holdout_records_exact_reference_source(self):
        species=pd.DataFrame({'sector':['focal','reference','reference','reference'],
            'accepted_id':['11','1','2','3'],'genus':['A']*4,'n':[1]*4,'strict9_white':[.8,.2,.3,.4]})
        r,s=m.score_holdout(species,'focal','reference',['strict9_white'],3)
        self.assertAlmostEqual(s.pred_strict9_white.iloc[0],.3)
        self.assertTrue((s.train_species_in_genus>=3).all())
        self.assertNotIn('11',set(species.loc[species.sector.eq('reference'),'accepted_id']))
    def test_genus_decomposition_exact_closure(self):
        species=pd.DataFrame({'sector':['focal']*2+['reference']*2,
         'accepted_id':['1','2','3','4'],'genus':['A','B','A','C'],
         'n':[1,1,1,1],'strict9_white':[.8,.4,.6,.2]})
        r=m.genus_composition(species,'focal','reference','strict9_white')
        self.assertAlmostEqual(r['genus_equal_difference_pp'],20)
        self.assertAlmostEqual(sum(r[k] for k in ['common_genus_trait_pp','common_genus_weight_pp',
          'focal_only_genera_pp','reference_only_genera_pp']),20)
    def test_input_palette_all_values_finite_for_measured_taxa(self):
        mapping,_=m.accepted_taxonomy(self.names,self.accepted)
        t=pd.DataFrame({'sector':['arctic'],'inat_taxon_id':['1'],'n':[1],'strict9_white':[np.nan]})
        with self.assertRaisesRegex(ValueError,'Nonfinite'):
            m.collapse_to_accepted_species(t,mapping,['strict9_white'])

if __name__=='__main__':unittest.main()
