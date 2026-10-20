"""Actual source-only target tier accounting; synthetic test data contain no flowers."""
import csv
import gzip
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch
from scripts.acquisition import mois_onepass_full93m_target_tiers as m


class TestOnePassTargetTiers(unittest.TestCase):
    def setUp(self):
        self.work = TemporaryDirectory()
        self.root = Path(self.work.name)
        (self.root / 'index' / 'photos').mkdir(parents=True)
        (self.root / 'index' / 'observations').mkdir(parents=True)
        photo_records = [
            [('11','OA','1','cc-by','jpg','u1'),
             ('12','OA','1','cc-by','JPG','u1'),
             ('13','OB','1','all-rights-reserved','gif','u2')],
            [('14','OC','2','cc-by','png','u3'),
             ('15','OD','2','cc-by','jpeg','u4'),
             ('16','OE','2','cc-by','','u4')],
        ]
        obs_records = [
            [('OA','1','obs1','2025-07-11','33.1','130','100'),
             ('OB','1','obs2','2025-99-22','12','40','100')],
            [('OC','2','obs3','2025-01-01','44','133','25'),
             ('OD','2','obs4','2025-02-01','80','-40','30'),
             ('OE','2','obs5','2025-03-01','91','200','')],
        ]
        shards = []
        for n in range(2):
            a = self.root / 'index' / 'photos' / f'part-{n:03d}.csv.gz'
            b = self.root / 'index' / 'observations' / f'part-{n:03d}.csv.gz'
            for p,header,rows in (
                (a,m.SOURCE_PHOTO_COLUMNS,photo_records[n]),
                (b,m.SOURCE_OBS_COLUMNS,obs_records[n])):
                with gzip.open(p,'wt',newline='') as f:
                    w=csv.writer(f);w.writerow(header);w.writerows(rows)
            shards.append({'part':n,'photo_count':len(photo_records[n]),
                           'observation_count':len(obs_records[n]),
                           'photos_sha256':m.sha256(a),'observations_sha256':m.sha256(b)})
        self.manifest=self.root/'index'/'manifest.json'
        self.manifest.write_text(json.dumps({
            'schema':m.SOURCE_SCHEMA,'photos_emitted':6,'taxa_in_frame':2,
            'photo_caps':None,'photo_colours_measured':False,
            'original_photos_are_globally_unique_by_source_SQLite_primary_key':True,
            'shards':shards}))
        self.score={'OA':0.0001,'OB':0.0003,'OC':0.0007,'OD':0.5,'OE':0.9}

    def tearDown(self):
        self.work.cleanup()

    def build(self):
        with patch.object(m,'obs_draw',side_effect=lambda uuid:self.score[uuid]):
            report=m.materialize_parts(self.root/'index',self.root/'out',[0,1],
                 expected_total=6,expected_taxa=2,expected_shards=2)
        return report

    def rows(self,tier):
        records=[]
        for n in range(2):
            p=self.root/'out'/'targets'/f'tier-{tier:02d}'/f'part-{n:03d}.csv.gz'
            with gzip.open(p,'rt',newline='') as f:
                records.extend(csv.DictReader(f))
        return records

    def test_all_targets_are_present_once_even_gif_and_bad_dates(self):
        report=self.build()
        self.assertEqual(report['original_photo_IDs_in_this_group'],6)
        self.assertEqual([report['counts_by_disjoint_tier'][str(x)] for x in (0,1,2,12,13)],
                         [2,1,1,1,1])
        rows=[r for t in range(len(m.TIER_ENDS)) for r in self.rows(t)]
        self.assertEqual({r['photo_id'] for r in rows}, {'11','12','13','14','15','16'})
        self.assertEqual(len(rows),6)
        self.assertEqual({r['observation_inclusion_score'] for r in self.rows(0)}, {'0.0001'})
        row={x['photo_id']:x for x in rows}
        self.assertEqual(row['13']['frozen_loader_supported'],'false')
        self.assertEqual(row['13']['valid_observed_date'],'false')
        self.assertEqual(row['13']['photo_license'],'all-rights-reserved')
        self.assertEqual(row['16']['valid_geopoint'],'false')
        self.assertEqual(row['12']['image_extension'],'jpg')
        self.assertFalse(report['flower_anthesis_certified'])
        self.assertEqual(report['colour_pixels_measured'],0)
        self.assertFalse(report['native_TDWG3_individual_or_island_status_inferred'])

    def test_nesting_and_new_tier_releases_previous_source_photos_unchanged(self):
        self.build()
        first={r['photo_id'] for r in self.rows(0)}
        second={r['photo_id'] for r in self.rows(1)}
        third={r['photo_id'] for r in self.rows(2)}
        self.assertEqual(first,{'11','12'})
        self.assertEqual(second,{'13'})
        self.assertEqual(third,{'14'})
        self.assertFalse(first & second or second & third or first & third)
        # A later wave only downloads its newly admitted target tier, not 93m again.
        self.assertEqual(sum(len(self.rows(t)) for t in range(3)),4)

    def test_original_photo_source_byte_drift_rejected(self):
        p=self.root/'index'/'photos'/'part-000.csv.gz'
        p.write_bytes(p.read_bytes()+b'changed')
        with self.assertRaisesRegex(ValueError,'SHA'):
            self.build()

    def test_duplicate_photo_in_shard_rejected(self):
        p=self.root/'index'/'photos'/'part-000.csv.gz'
        with gzip.open(p,'rt') as f:
            rows=list(csv.reader(f))
        rows.append(rows[1][:])
        with gzip.open(p,'wt',newline='') as f:
            w=csv.writer(f);w.writerows(rows)
        source=json.loads(self.manifest.read_text())
        source['shards'][0]['photos_sha256']=m.sha256(p)
        source['shards'][0]['photo_count']+=1
        self.manifest.write_text(json.dumps(source))
        with self.assertRaisesRegex(ValueError,'duplicate source photo ID'):
            self.build()

    def test_complete_union_proven_only_after_all_partition_receipts(self):
        self.build()
        with patch.object(m,'EXPECTED_SHARDS',2),patch.object(m,'EXPECTED_PHOTOS',6):
            report=m.verify_all_groups(self.manifest,[self.root/'out'/'receipt.json'])
        self.assertEqual(report['exact_original_photo_IDs_in_all_disjoint_tiers'],6)
        self.assertTrue(report['native_vs_introduced_not_hard_filtered'])
        self.assertTrue(report['source_only_NO_photo_pixels_or_hue_measured'])


if __name__=='__main__':unittest.main()