"""Kew accepted-family provenance and genus-vs-family region holdout tests."""
import tempfile,zipfile,unittest
from pathlib import Path
import pandas as pd
import numpy as np
from scripts.analysis import mois_family_palette_holdout_20261030 as m

class FamilyHoldoutTests(unittest.TestCase):
    def zip(self,lines):
        self.temp=tempfile.TemporaryDirectory();path=Path(self.temp.name)/'wcvp.zip'
        with zipfile.ZipFile(path,'w') as z:
            z.writestr('wcvp_names.csv','plant_name_id|taxon_name|taxon_status|taxon_rank|family|genus\n'+'\n'.join(lines)+'\n')
        self.addCleanup(self.temp.cleanup)
        return path
    def test_family_uses_accepted_species_only(self):
        path=self.zip(['1|Acer alpha|accepted|species|Sapindaceae|Acer',
                       '2|Acer beta|synonym|species|Sapindaceae|Acer',
                       '3|Betula gamma|accepted|species|Betulaceae|Betula',
                       '4|Acer delta|accepted|genus|Sapindaceae|Acer'])
        result,stats=m.read_accepted_families(path,{'1','2','3','4','missing'})
        self.assertEqual(set(result.accepted_id),{'1','3'})
        self.assertEqual(result.set_index('accepted_id').loc['1','family'],'Sapindaceae')
        self.assertEqual(stats['unresolved_accepted_ids'],3)
    def test_duplicate_family_conflict_rejected(self):
        path=self.zip(['1|Acer alpha|accepted|species|Sapindaceae|Acer',
                       '1|Acer alpha|accepted|species|Differentaceae|Acer'])
        with self.assertRaisesRegex(ValueError,'Conflicting official'):
            m.read_accepted_families(path,{'1'})
    def test_missing_family_remains_missing(self):
        path=self.zip(['1|Acer alpha|accepted|species||Acer'])
        result,stats=m.read_accepted_families(path,{'1'})
        self.assertTrue(result.empty)
        self.assertEqual(stats['unresolved_accepted_ids'],1)
    def test_taxonomic_contradiction_is_hard_error(self):
        src=pd.DataFrame({'sector':['north'],'accepted_id':['1'],'genus':['Acer'],'n':[1],'strict9_white':[.2]})
        fam=pd.DataFrame({'accepted_id':['1'],'accepted_name':['Wronggenus alpha'],'genus':['Wronggenus'],'family':['Sapindaceae']})
        with self.assertRaisesRegex(ValueError,'WCVP genus drift'):m.join_family(src,fam)
    def test_family_genus_conflict_is_hard_error(self):
        src=pd.DataFrame({'sector':['north','south'],'accepted_id':['1','2'],'genus':['Acer','Acer'],'n':[1,1],'strict9_white':[.2,.5]})
        fam=pd.DataFrame({'accepted_id':['1','2'],'accepted_name':['Acer alpha','Acer beta'],
                          'genus':['Acer','Acer'],'family':['Sapindaceae','Otheraceae']})
        with self.assertRaisesRegex(ValueError,'Cross-family genus'):m.join_family(src,fam)
    def test_family_holdout_no_target_species_leak(self):
        a=pd.DataFrame({'sector':['focal']*2+['reference']*6,
          'accepted_id':['10','1','1','2','3','4','5','6'],
          'family':['FamA','FamA','FamA','FamA','FamA','FamA','FamA','FamB'],
          'strict9_white':[.7,.9,.1,.2,.3,.4,.5,.6]})
        r,s=m.score_family_holdout(a,'focal','reference',['strict9_white'],5)
        self.assertEqual(r['novel_target_species'],1)
        self.assertEqual(r['supported_novel_species'],1)
        self.assertEqual(s.accepted_id.tolist(),['10'])
        self.assertAlmostEqual(s.pred_strict9_white.iloc[0],.3)
    def test_no_family_holdout_support_not_imputed(self):
        a=pd.DataFrame({'sector':['focal','reference'],'accepted_id':['10','1'],
                        'family':['FamA','FamB'],'strict9_white':[.8,.2]})
        r,s=m.score_family_holdout(a,'focal','reference',['strict9_white'],5)
        self.assertEqual(r['supported_novel_species'],0)
        self.assertEqual(len(s),0)
    def test_exact_genus_vs_family_comparison_common_target_only(self):
        fam=pd.DataFrame({'focal':['focal','focal'],'reference':['ref','ref'],
                          'accepted_id':['10','11'],'strict9_white':[.8,.5],
                          'pred_strict9_white':[.4,.6]})
        genus=pd.DataFrame({'focal':['focal'],'reference':['ref'],
                            'accepted_id':['10'],'pred_strict9_white':[.7]})
        metrics,matched=m.compare_family_genus(fam,genus,['strict9_white'])
        self.assertEqual(metrics['compared_novel_species'],1)
        self.assertEqual(matched.accepted_id.tolist(),['10'])
        self.assertAlmostEqual(metrics['metrics']['strict9_white']['family_mse'],.16)
        self.assertAlmostEqual(metrics['metrics']['strict9_white']['genus_mse'],.01)
    def test_family_rank_partition_reconciles(self):
        src=pd.DataFrame({'sector':['a','a','b','b'],
                          'accepted_id':['1','2','3','4'],'family':['F1','F2','F1','F3'],
                          'strict9_white':[.8,.4,.6,.2]})
        d=m.family_composition(src,'a','b','strict9_white')
        self.assertAlmostEqual(d['family_equal_difference_pp'],20)
        self.assertAlmostEqual(sum(d[k] for k in ['common_family_trait_pp','common_family_weight_pp',
           'focal_only_families_pp','reference_only_families_pp']),20)

if __name__=='__main__':unittest.main()
