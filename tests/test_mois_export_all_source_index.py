"""Small SQLite source verifies uncapped, no duplicate/colour invention."""
import csv
import gzip
import sqlite3
from pathlib import Path
import tempfile
import unittest
from scripts.acquisition.mois_export_all_source_index import export_index

class TestFullSource(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.root=Path(self.temp.name)
        self.db=self.root/'source.sqlite'
        with sqlite3.connect(self.db) as c:
            c.execute('CREATE TABLE observations(uuid TEXT PRIMARY KEY,taxon INTEGER,observer TEXT,observed_on TEXT,lat REAL,lon REAL,accuracy TEXT)')
            c.execute('CREATE TABLE photos(photo_id TEXT PRIMARY KEY,uuid TEXT,license TEXT,extension TEXT,photographer TEXT,taxon INTEGER,ambiguous INTEGER)')
            c.executemany('INSERT INTO observations VALUES (?,?,?,?,?,?,?)',[
                ('A',1,'x','2025-05-10',30.,135.,'100'),
                ('B',2,'y','2025-06-12',45.,138.,'100'),
                ('C',2,'y','2025-07-12',44.,137.,'500')])
            c.executemany('INSERT INTO photos VALUES (?,?,?,?,?,?,?)',[
                ('101','A','cc-by','jpg','x',1,0),
                ('102','A','cc-by','jpg','x',1,0),
                ('103','B','cc-by','jpg','y',2,0),
                ('104','C','cc-by-sa','jpeg','y',2,0),
                ('105','C','cc-by','jpg','y',2,1)])
        self.csv=self.root/'counts.csv'
        with self.csv.open('w',newline='') as f:
            w=csv.writer(f)
            w.writerow(['inat_taxon_id','species','n_unique_source_photo_ids','snapshot_date'])
            w.writerows([(1,'Plant alpha',2,'2026-09-27'),(2,'Plant beta',2,'2026-09-27'),
                         (3,'Plant zero',0,'2026-09-27')])
    def tearDown(self):self.temp.cleanup()
    def export(self):
        return export_index(self.db,self.csv,self.root/'index',parts=4,expected_photos=4,expected_taxa=3)
    def test_entire_original_photo_frame_and_observation_repeats_retained(self):
        d=self.export()
        self.assertEqual((d['photos_emitted'],d['taxa_in_frame'],d['taxa_with_original_photos']),
                         (4,3,2))
        photos=[]
        for p in (self.root/'index'/'photos').glob('*.csv.gz'):
            with gzip.open(p,'rt',newline='') as f:photos+=list(csv.DictReader(f))
        self.assertEqual({p['photo_id'] for p in photos},{'101','102','103','104'})
        self.assertEqual(sum(p['observation_uuid']=='A' for p in photos),2)
        self.assertFalse(d['photo_colours_measured'])
        self.assertIsNone(d['photo_caps'])
    def test_conflicting_photo_ids_quarantined(self):
        d=self.export()
        self.assertEqual(d['conflicting_photo_ids_excluded'],1)
        self.assertEqual(sum(p['photo_count'] for p in d['shards']),4)
    def test_source_count_mismatch_fails_closed(self):
        self.csv.write_text(self.csv.read_text().replace('Plant alpha,2','Plant alpha,3'))
        with self.assertRaisesRegex(ValueError,'count|mismatch'):
            self.export()
        self.assertFalse((self.root/'index'/'manifest.json').exists())
    def test_existing_complete_manifest_never_overwritten(self):
        self.export()
        with self.assertRaisesRegex(ValueError,'overwrite'):
            self.export()

if __name__=='__main__': unittest.main()
