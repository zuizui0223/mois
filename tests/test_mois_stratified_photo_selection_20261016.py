"""Synthetic-only tests; never represent real 93M-source biological outcomes."""
import csv
import importlib.util
import json
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parents[1] / 'scripts/acquisition/mois_stratified_photo_selection_20261016.py'
spec = importlib.util.spec_from_file_location('mois_sampler', SRC)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def fixtures(tmp_path):
    counts = tmp_path / 'complete.csv'
    with counts.open('w', newline='') as f:
        w = csv.writer(f)
        w.writerow(m.COUNT_COLUMNS)
        for tax, sp, n in ((1,'Alpha alba',5),(2,'Beta rosea',4),(3,'Gamma picta',0)):
            w.writerow([tax,sp,n,'complete','a'*64])
    audit = tmp_path / 'audit.json'
    audit.write_text(json.dumps({'schema':'fcp_complete_photo_ID_source_vs_independent_provider_audit_v1',
        'status':'REAL_ALL_42111_ORIGINAL_PHOTO_ID_COUNTS_RECONCILED_NO_PIXELS',
        'snapshot_date':m.SNAPSHOT,'fully_frozen_count_ledger_sha256':m.sha256_file(counts),
        'original_photo_IDs_covered_without_count_loss':True,
        'full_study_original_pixels_classified':False,
        'flower_colour_H2_W_or_environment_results':False,
        'species':3,'distinct_original_photo_ids':9}))
    headers = ['inat_taxon_id','genus','photo_id','observation_id','latitude','longitude','bio12_mm','photo_license']
    data = [
        (1,'Alpha','11','o1','-12','10','250','CC-BY'),
        (1,'Alpha','12','o1','-12','10','250','CC-BY'),
        (1,'Alpha','13','o2','50','125','1800','CC-BY'),
        (1,'Alpha','14','o3','','','','CC-BY'),
        (1,'Alpha','15','o4','25','-135','1000','CC-BY'),
        (2,'Beta','21','o9','-15','80','500','CC0'),
        (2,'Beta','22','o10','-15','80','500','CC0'),
        (2,'Beta','23','o11','-15','80','500','CC0'),
        (2,'Beta','24','o12','-15','80','500','CC0'),
    ]
    manifest = tmp_path / 'candidate.csv'
    with manifest.open('w',newline='') as f:
        w=csv.writer(f)
        w.writerow(headers)
        w.writerows(data)
    return counts,audit,manifest,headers,data


def test_selection_complete_and_observer_single_photo(tmp_path):
    c,a,mf,_,_=fixtures(tmp_path)
    out = tmp_path/'selected.csv'; receipt = tmp_path/'receipt.json'
    info=m.select(mf,c,a,out,receipt,_synthetic=True)
    assert info['source_photo_opportunities_verified']==9
    assert info['distinct_source_observations']==8
    assert info['selected_original_photo_ids']==8
    assert info['ready_for_13_environment_association'] is False
    selected=list(csv.DictReader(out.open()))
    assert len({(r['inat_taxon_id'],r['observation_id']) for r in selected})==8
    assert any(int(r['source_photos_in_observation'])==2 for r in selected)
    assert any(abs(float(r['photo_inclusion_probability'])-0.5)<1e-12 for r in selected)
    with pytest.raises(ValueError, match='Refuse to overwrite'):
        m.select(mf,c,a,out,receipt,_synthetic=True)


def test_outcome_blind_and_source_fail_closed(tmp_path):
    c,a,mf,head,data=fixtures(tmp_path)
    a.write_text(json.dumps({**json.loads(a.read_text()),'status':'VERIFIED_COMPLETE_COUNTS_ONLY_NOT_COLOUR_MEASUREMENT'}))
    with pytest.raises(ValueError, match='count-only'):
        m.select(mf,c,a,tmp_path/'sel.csv',tmp_path/'result.json',_synthetic=True)
    a.write_text(json.dumps({**json.loads(a.read_text()),'status':'REAL_ALL_42111_ORIGINAL_PHOTO_ID_COUNTS_RECONCILED_NO_PIXELS'}))
    with mf.open('w',newline='') as f:
        w=csv.writer(f);w.writerow(head+['flower_colour']);w.writerows([list(r)+['white'] for r in data])
    with pytest.raises(ValueError,match='forbidden'):
        m.select(mf,c,a,tmp_path/'sel.csv',tmp_path/'result.json',_synthetic=True)
    assert not (tmp_path/'sel.csv').exists()


def test_missing_photo_or_duplicates_abort_without_output(tmp_path):
    c,a,mf,head,data=fixtures(tmp_path)
    with mf.open('w',newline='') as f:
        w=csv.writer(f);w.writerow(head);w.writerows(data[:-1])
    with pytest.raises(ValueError,match='source loss'):
        m.select(mf,c,a,tmp_path/'sel.csv',tmp_path/'result.json',_synthetic=True)
    assert not (tmp_path/'sel.csv').exists()
    with mf.open('w',newline='') as f:
        w=csv.writer(f);w.writerow(head);w.writerows(data[:-1]+[data[-2]])
    with pytest.raises(ValueError,match='duplicate'):
        m.select(mf,c,a,tmp_path/'sel.csv',tmp_path/'result.json',_synthetic=True)


def test_all_strata_positive_and_budget():
    # 40 mutually occupied cells can each receive a unit under the frozen cap.
    groups = {(f'{g}',f'{p}'): 200 for g in range(10) for p in range(4)}
    q = m.quotas(groups)
    assert sum(q.values())==40 and all(v>=1 for v in q.values())
    q=m.quotas({('A','dry'):100,('B','wet'):2})
    assert sum(q.values())==40 and q[('B','wet')]==2


def test_synthetic_cannot_pass_production_gate(tmp_path):
    c,a,_,_,_=fixtures(tmp_path)
    with pytest.raises(ValueError,match='population fingerprint changed'):
        m.read_verified_counts(c,a)
