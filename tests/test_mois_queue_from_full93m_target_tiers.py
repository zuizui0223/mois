"""Only new image IDs from one target tier enter full FCP ROI acquisition."""
import csv
import gzip
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch
from scripts.acquisition import mois_onepass_full93m_target_tiers as r
from scripts.acquisition import mois_queue_from_full93m_target_tiers as q


class TestFromTargetTier(unittest.TestCase):
    def setUp(self):
        self.tmp=TemporaryDirectory()
        self.root=Path(self.tmp.name)
        index=self.root/'index'
        (index/'photos').mkdir(parents=True)
        (index/'observations').mkdir(parents=True)
        shards=[]
        for part,photo,obs in [
            (0,[('1','A','1','cc-by','jpg','p1'),('2','A','1','cc-by','jpg','p1')],
               [('A','1','one','2025-07-01','20','30','100')]),
            (1,[('3','B','2','cc-by','gif','p2')],
               [('B','2','two','2025-08-01','15','45','100')]),
        ]:
            pf=index/'photos'/f'part-{part:03d}.csv.gz'
            of=index/'observations'/f'part-{part:03d}.csv.gz'
            for f,cols,rows in ((pf,r.SOURCE_PHOTO_COLUMNS,photo),
                                (of,r.SOURCE_OBS_COLUMNS,obs)):
                with gzip.open(f,'wt',newline='') as stream:
                    writer=csv.writer(stream);writer.writerow(cols);writer.writerows(rows)
            shards.append({'part':part,'photo_count':len(photo),
                           'observation_count':len(obs),
                           'photos_sha256':r.sha256(pf),
                           'observations_sha256':r.sha256(of)})
        self.index=index
        self.source=index/'manifest.json'
        self.source.write_text(json.dumps({
            'schema':r.SOURCE_SCHEMA,'photos_emitted':3,'taxa_in_frame':2,
            'photo_caps':None,'photo_colours_measured':False,
            'original_photos_are_globally_unique_by_source_SQLite_primary_key':True,
            'shards':shards}))
        with patch.object(r,'obs_draw',lambda _:0.0001):
            self.receipt=r.materialize_parts(index,self.root/'tiers',[0,1],
                          expected_total=3,expected_taxa=2,expected_shards=2)
        self.registry=self.root/'registry.json'
        self.registry.write_text(json.dumps({
            'schema':'mois_complete_full93m_original_photo_target_registry_v1',
            'exact_original_photo_IDs_in_all_disjoint_tiers':3,
            'original_named_taxa_from_full_source':2,
            'full_source_manifest_sha256':r.sha256(self.source),
            'new_target_counts_by_tier':self.receipt['counts_by_disjoint_tier'],
            'native_vs_introduced_not_hard_filtered':True,
            'source_only_NO_photo_pixels_or_hue_measured':True,
        }))
        self.prior=self.root/'previous.csv'
        self.prior.write_text('photo_id\n1\n')
        self.counts=self.root/'counts.csv'
        self.counts.write_text('inat_taxon_id,species\n1,Plant one\n2,Plant two\n')

    def tearDown(self):
        self.tmp.cleanup()

    def make(self):
        with patch.object(q,'EXPECTED_PHOTOS',3),patch.object(q,'EXPECTED_TAXA',2), \
             patch.object(q,'obs_draw',lambda _:0.0001):
            return q.queue_tier(self.registry,[self.root/'tiers'/'receipt.json'],
                     self.root/'tiers'/'targets',self.counts,[self.prior],
                     self.root/'queue',tier=0,compute_bucket=0,
                     compute_bucket_count=1,expected_taxa=2,expected_shards=2,
                     check_taxa_sha=False)

    def test_signed_source_tier_only_unseen_supported_original_images(self):
        d=self.make()
        self.assertEqual(d['source_tier_total_photo_candidates'],3)
        self.assertEqual(d['source_tier_compute_bucket_total_photos'],3)
        self.assertEqual(d['previously_measured_photo_ids_skipped'],1)
        self.assertEqual(d['selected_new_photos_unsupported_by_frozen_loader'],1)
        self.assertEqual(d['new_photo_attempts_planned'],1)
        self.assertTrue(d['native_photo_status_not_inferred_or_hard_filtered'])
        with (self.root/'queue'/'shards'/'shard-000.csv').open(newline='') as f:
            rows=list(csv.DictReader(f))
        self.assertEqual([x['photo_id'] for x in rows],['2'])
        self.assertEqual(rows[0]['source_arm'],'global93m_onepass_target_tier_00')
        self.assertTrue(rows[0]['photo_url_large'].endswith('/2/large.jpg'))
        self.assertTrue((self.root/'queue'/'unsupported_selected_original_photo_formats.jsonl').read_text().find('"3"')!=-1)

    def test_modified_photo_tier_fails_before_any_image_request(self):
        file=self.root/'tiers'/'targets'/'tier-00'/'part-000.csv.gz'
        file.write_bytes(file.read_bytes()+b'tamper')
        with self.assertRaisesRegex(ValueError,'source bytes changed'):
            self.make()

    def test_same_observation_remains_in_one_compute_bucket(self):
        self.assertEqual(q.observation_compute_bucket('A',64),
                         q.observation_compute_bucket('A',64))
        self.assertTrue(0<=q.observation_compute_bucket('A',64)<64)


if __name__=='__main__':unittest.main()