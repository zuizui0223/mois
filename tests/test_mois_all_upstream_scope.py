import csv, importlib.util, json, sqlite3
from pathlib import Path
import pytest
FILE=Path(__file__).parents[1]/'scripts/acquisition/mois_all_upstream_scope.py'
spec=importlib.util.spec_from_file_location('scope', FILE)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def test_every_row_and_every_photo_count_no_species_quota(tmp_path):
 p=tmp_path/'counts.csv'
 with p.open('w',newline='') as f:
  w=csv.writer(f);w.writerow(['inat_taxon_id','species','n_unique_source_photo_ids','snapshot_date'])
  w.writerows([(1,'A a',10000,'2026-09-27'),(2,'B b',0,'2026-09-27'),(3,'C c',900001,'2026-09-27')])
 rows,n=m.read_counts(p)
 assert len(rows)==3 and n==910001
 assert rows[1]['n_unique_source_photo_ids']==0

def test_stage_accepts_multiple_photos_same_observation_and_taxon(tmp_path):
 db=tmp_path/'full.sqlite'
 c=sqlite3.connect(db)
 c.execute('CREATE TABLE photos(photo_id TEXT PRIMARY KEY, uuid TEXT, taxon INTEGER, ambiguous INTEGER)')
 c.executemany('INSERT INTO photos VALUES(?,?,?,?)', [('1','same',1,0),('2','same',1,0),('3','other',2,0),('4','bad',2,1)])
 c.commit();c.close()
 assert m.database_counts(db)==({1:2,2:1},1)

def test_source_counts_reject_duplicate_taxon_instead_of_silently_overwriting(tmp_path):
 p=tmp_path/'counts.csv';p.write_text('inat_taxon_id,species,n_unique_source_photo_ids,snapshot_date\n1,A a,2,2026-09-27\n1,A a,3,2026-09-27\n')
 with pytest.raises(ValueError, match='Duplicate'):m.read_counts(p)

def test_negative_counts_rejected(tmp_path):
 p=tmp_path/'counts.csv';p.write_text('inat_taxon_id,species,n_unique_source_photo_ids,snapshot_date\n1,A a,-1,2026-09-27\n')
 with pytest.raises(ValueError):m.read_counts(p)
