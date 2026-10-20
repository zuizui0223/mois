"""Offline source-identity tests; fixtures are not research outcomes."""
import csv, importlib.util, json, tempfile, unittest
from pathlib import Path
MODULE = Path(__file__).resolve().parents[1]/'scripts/acquisition/mois_completed_pixel_checkpoint.py'

class TestCheckpoint(unittest.TestCase):
    def setUp(self):
        if not MODULE.exists(): self.fail('Completed-pixel checkpoint implementation missing')
        spec=importlib.util.spec_from_file_location('checkpoint',MODULE)
        self.m=importlib.util.module_from_spec(spec);spec.loader.exec_module(self.m)
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name);self.source=self.root/'source';self.workers=self.root/'workers'
        (self.source/'shards').mkdir(parents=True);self.workers.mkdir()
        self.method={'synthetic_fixture_only':True}
        self.rows=[dict(photo_id=str(i),measurement_id='m'+str(i),observation_id='o'+str(i),inat_taxon_id='7',species='Fixture plant',shard=str((i-1)//2)) for i in range(1,5)]
        for u in range(2):self.write(self.source/'shards'/f'shard-{u:03d}.csv',self.rows[u*2:u*2+2])
        self.plan=dict(schema='mois_full_source_observation_sampling_cohort_v1',colour_measured=False,full_frame_photo_ids=93209351,worker_ids=[0,1],work_units=2,new_photo_attempts_planned=4,photo_shard_sha256={str(i):self.m.sha(self.source/'shards'/f'shard-{i:03d}.csv') for i in range(2)})
        self.save_plan();self.seal(0)
    def write(self,path,rows):
        path.parent.mkdir(parents=True,exist_ok=True)
        with path.open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    def save_plan(self):(self.source/'manifest.json').write_text(json.dumps(self.plan))
    def seal(self,u):
        folder=self.workers/f'group-{u}'/f'shard_{u:03d}'
        a,b=self.rows[u*2:u*2+2]
        outcomes=[dict(measurement_id=a['measurement_id'],morph='white',measurement_status='classified_four_state_morph',image_sha256='a'*64),dict(measurement_id=b['measurement_id'],morph='mixed_uncertain',measurement_status='image_acquisition_failed',image_sha256='')]
        t=folder/'sealed'/f'terminal_p{u:03d}_b0000000.csv';self.write(t,outcomes)
        self.write(folder/'photos.csv',[dict(x,**{k:v for k,v in y.items() if k!='measurement_id'},classified_four_hue=str(y['morph']=='white').lower()) for x,y in zip((a,b),outcomes)])
        r=dict(schema='mois_full_original93m_real_RGB_photo_unit_v1',shard=u,attempted=2,downloaded=1,classifiable=1,acquisition_failed=1,unclassifiable_acquired=0,source_sha256=self.plan['photo_shard_sha256'][str(u)],photos_sha256=self.m.sha(folder/'photos.csv'),terminal_sha256=self.m.sha(t),photo_colours_measured=True,raw_licensed_image_bytes_persisted=False,method=self.method)
        (folder/'receipt.json').write_text(json.dumps(r));return folder
    def collect(self):return self.m.collect_wave(self.source,self.workers,self.method)
    def test_only_complete_units_count(self):
        rows,r=self.collect();self.assertEqual([x['photo_id'] for x in rows],['1','2']);self.assertEqual(r['attempted'],2);self.assertEqual(r['classifiable'],1);self.assertEqual(r['acquisition_failed'],1);self.assertEqual(r['pending_photo_attempts'],2);self.assertEqual(r['missing_units'],[1]);self.assertFalse(r['wave_complete'])
    def test_all_units_complete(self):
        self.seal(1);rows,r=self.collect();self.assertEqual(len(rows),4);self.assertEqual(r['classifiable'],2);self.assertTrue(r['wave_complete'])
    def test_no_finished_units_not_flower_absence(self):
        import shutil
        shutil.rmtree(self.workers);self.workers.mkdir();rows,r=self.collect();self.assertEqual(rows,[]);self.assertEqual(r['pending_photo_attempts'],4);self.assertFalse(r['wave_complete'])
    def test_changed_bytes(self):
        p=next(self.workers.rglob('photos.csv'));p.write_text(p.read_text().replace('white','red_pink'))
        with self.assertRaisesRegex(ValueError,'SHA'):self.collect()
    def change_photo(self,key,value):
        p=next(self.workers.rglob('photos.csv'));rows=self.m.read_rows(p);rows[0][key]=value;self.write(p,rows)
        rp=p.with_name('receipt.json');r=json.loads(rp.read_text());r['photos_sha256']=self.m.sha(p);rp.write_text(json.dumps(r))
    def test_changed_source_identity_with_new_hash(self):
        self.change_photo('observation_id','wrong')
        with self.assertRaisesRegex(ValueError,'identity'):self.collect()
    def test_changed_hue_with_new_hash(self):
        self.change_photo('morph','red_pink')
        with self.assertRaisesRegex(ValueError,'terminal'):self.collect()
    def test_wrong_method_or_count(self):
        p=next(self.workers.rglob('receipt.json'));r=json.loads(p.read_text());r['classifiable']=2;p.write_text(json.dumps(r))
        with self.assertRaisesRegex(ValueError,'count'):self.collect()
        r['classifiable']=1;r['method']={};p.write_text(json.dumps(r))
        with self.assertRaisesRegex(ValueError,'method'):self.collect()
    def test_duplicate_units(self):
        import shutil
        p=next(self.workers.rglob('receipt.json')).parent;shutil.copytree(p,self.workers/'duplicate'/p.name)
        with self.assertRaisesRegex(ValueError,'Duplicate'):self.collect()
    def test_source_count_reconciles(self):
        self.plan['new_photo_attempts_planned']=5;self.save_plan()
        with self.assertRaisesRegex(ValueError,'source count'):self.collect()
    def test_terminal_required(self):
        next(self.workers.rglob('terminal_*.csv')).unlink()
        with self.assertRaisesRegex(ValueError,'terminal'):self.collect()
    def test_union_preserves_old_hues_unknowns_and_rejects_overlap(self):
        import pandas as pd
        old=pd.DataFrame([dict(photo_id='90',morph='white',classified=True,white=True,observation_id='z',inat_taxon_id='5',species='Old plant',arm='legacy',measurement_status='classified_four_state_morph',technical_failure=False)])
        fresh,_=self.collect()
        for r in fresh:r['checkpoint_wave']='p0004'
        union,report=self.m.union_projection(old,pd.DataFrame(fresh))
        self.assertEqual(union.photo_id.tolist(),['90','1','2']);self.assertEqual(union.morph.tolist(),['white','white','mixed_uncertain']);self.assertEqual(union.classified.tolist(),[True,True,False]);self.assertEqual(report['new_actual_four_hue_classified'],1);self.assertEqual(report['combined_unique_photo_attempts'],3)
        fresh[0]['photo_id']='90'
        with self.assertRaisesRegex(ValueError,'overlap'):self.m.union_projection(old,pd.DataFrame(fresh))
    def test_legacy_unclassified_labels_preserved_without_inventing_hues(self):
        import pandas as pd
        old=pd.DataFrame([dict(photo_id=str(90+i),morph=label,classified=False,white=False,
            observation_id='legacy'+str(i),inat_taxon_id='5',species='Old plant',arm='legacy',
            measurement_status='not_evaluable',technical_failure=False)
            for i,label in enumerate(('', 'legacy_roi_unknown'))])
        fresh,_=self.collect()
        for r in fresh:r['checkpoint_wave']='p0004'
        union,report=self.m.union_projection(old,pd.DataFrame(fresh))
        self.assertEqual(union.morph.iloc[:2].tolist(),['','legacy_roi_unknown'])
        self.assertFalse(union.classified.iloc[:2].any())
        self.assertFalse(union.white.iloc[:2].any())
        self.assertEqual(report['combined_four_hue_classified'],1)

    def test_new_noncanonical_unknown_is_still_rejected(self):
        import pandas as pd
        old=pd.DataFrame([dict(photo_id='90',morph='white',classified=True,white=True,
                              arm='legacy',measurement_status='classified_four_state_morph')])
        fresh,_=self.collect()
        for r in fresh:r['checkpoint_wave']='p0004'
        fresh[1]['morph']='new_unknown_unsupported'
        with self.assertRaisesRegex(ValueError,'hue'):
            self.m.union_projection(old,pd.DataFrame(fresh))

    def test_combined_projection_preserves_geographic_selection_probability_and_arm(self):
        import pandas as pd
        old=pd.DataFrame([dict(photo_id='90',morph='white',classified=True,white=True,
            arm='original_legacy',observation_id='old_obs',
            measurement_status='classified_four_state_morph')])
        fresh,_=self.collect()
        for index,row in enumerate(fresh):
            row['checkpoint_wave']='geo05'
            row['source_arm']='global93m_geographic_enrichment_tier_05'
            row['selection_probability']=('0.0016','0.0008')[index]
        frozen=old.copy(deep=True)
        union,counts=self.m.union_projection(old,pd.DataFrame(fresh))
        self.assertTrue(old.equals(frozen))
        self.assertEqual(union.photo_id.tolist(),['90','1','2'])
        self.assertEqual(union.selection_probability.tolist(),['','0.0016','0.0008'])
        self.assertEqual(union.source_arm.tolist(),[
            'original_legacy','global93m_geographic_enrichment_tier_05',
            'global93m_geographic_enrichment_tier_05'])
        self.assertEqual(union.checkpoint_wave.tolist(),['baseline','geo05','geo05'])
        self.assertEqual(union.morph.tolist(),['white','white','mixed_uncertain'])
        self.assertEqual(counts['combined_four_hue_classified'],2)

    def test_unparseable_or_nonpositive_new_probability_cannot_enter_combined_checkpoint(self):
        import pandas as pd
        old=pd.DataFrame([dict(photo_id='90',morph='white',classified=True,white=True,
            arm='historical',observation_id='old_obs')])
        for invalid in ('', 'not-a-number', 'NaN', '0', '-0.1','1.2'):
            fresh,_=self.collect()
            for row in fresh:
                row['checkpoint_wave']='geo05'
                row['selection_probability']='0.0008'
            fresh[0]['selection_probability']=invalid
            with self.subTest(invalid=invalid):
                with self.assertRaisesRegex(ValueError, 'selection probability'):
                    self.m.union_projection(old,pd.DataFrame(fresh))

    def test_same_observation_cannot_have_two_inclusion_probabilities(self):
        import pandas as pd
        old=pd.DataFrame([dict(photo_id='90',morph='white',classified=True,white=True,
            arm='historical',observation_id='old_obs')])
        fresh,_=self.collect()
        for index,row in enumerate(fresh):
            row['checkpoint_wave']='geo05'
            row['observation_id']='same_observation'
            row['selection_probability']=('0.0016','0.0008')[index]
        with self.assertRaisesRegex(ValueError,'selection probability'):
            self.m.union_projection(old,pd.DataFrame(fresh))

if __name__=='__main__':unittest.main()
