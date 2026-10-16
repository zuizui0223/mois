"""Mois regression tests for source-verified seasonal opportunity (not colours)."""
import csv
import json
from pathlib import Path
import pytest
from scripts.analysis.fcp_source_month_season_priority import work, checksum, local_season

KEYS=['inat_taxon_id','species','year_month','lat_0p25_cell','lon_0p25_cell',
      'all_source_photo_ids_with_5km_accuracy','distinct_source_observations','distinct_source_observers']


def fixture(tmp_path):
    f=tmp_path/'source.csv';r=tmp_path/'receipt.json'
    def make(tax,name,ym,ilat,ilon,obs=4):return [tax,name,ym,ilat,ilon,obs,obs,obs]
    rows=[make('1','North','2024-01',600,100),make('1','North','2024-07',600,100),
          make('1','North','2024-11',600,100,1),
          make('2','South','2024-01',200,100),make('2','South','2024-07',200,100),
          make('2','South','2026-08',200,100)]
    with f.open('w',newline='') as o:
        w=csv.writer(o);w.writerow(KEYS);w.writerows(rows)
    r.write_text(json.dumps(dict(status='COMPLETE_ORIGINAL_SOURCE_CALENDAR_CENSUS_NO_FLOWER_PIXELS',
            original_source_site_month_file_sha256=checksum(f),
            source_photo_opportunity_is_not_flowering_detection=True,
            new_flower_colour_or_white_frequency_estimated=False,
            n_taxa=2,n_original_photo_ids=16)))
    return f,r


def test_hemisphere_seasons():
    assert local_season(1,600)=='winter'
    assert local_season(7,600)=='summer'
    assert local_season(1,200)=='summer'
    assert local_season(7,200)=='winter'
    assert local_season(4,360)=='tropical_no_four_season_label'


def test_signed_realized_candidate_ranking_without_pixel_outcomes(tmp_path):
    f,r=fixture(tmp_path)
    p=work(f,r,tmp_path/'out',min_obs=4,test_only=True)
    assert p['status']=='SYNTHETIC_NOT_REAL'
    assert p['n_eligible_geocell_years']==2 and p['n_eligible_taxa']==2
    assert p['candidate_photos_NOT_proven_flowering_or_white'] is True
    path=tmp_path/'out'/'season_depth_24_source_only_geocell_year_priority.csv'
    assert checksum(path)==p['output_priority_sha256']
    rows=list(csv.DictReader(path.open()))
    assert len(rows)==2
    north=next(q for q in rows if q['inat_taxon_id']=='1')
    south=next(q for q in rows if q['inat_taxon_id']=='2')
    assert {north['local_season_a'],north['local_season_b']}=={'winter','summer'}
    assert {south['local_season_a'],south['local_season_b']}=={'winter','summer'}
    assert north['independent_observers_ACROSS_seasons_verified']=='False'
    assert north['is_flowering_photo_or_colour']=='False'
    with pytest.raises(ValueError,match='overwrite'):
        work(f,r,tmp_path/'out',min_obs=4,test_only=True)


def test_changed_original_site_data_fails_signed_sha(tmp_path):
    f,r=fixture(tmp_path)
    f.write_text(f.read_text().replace('2024-01','2025-01',1))
    with pytest.raises(ValueError,match='SHA-signed'):
        work(f,r,tmp_path/'out',min_obs=4,test_only=True)


def test_synthetic_census_never_promoted_to_real(tmp_path):
    f,r=fixture(tmp_path)
    with pytest.raises(ValueError,match='SHA-signed'):
        work(f,r,tmp_path/'out',min_obs=4,test_only=False)


def test_require_at_least_two_seasons_with_enough_independent_observations(tmp_path):
    f,r=fixture(tmp_path)
    p=work(f,r,tmp_path/'out',min_obs=5,test_only=True)
    assert p['n_eligible_geocell_years']==0
    assert p['n_eligible_taxa']==0


def test_invalid_unverified_year_not_inferred(tmp_path):
    f,r=fixture(tmp_path)
    with pytest.raises(ValueError,match='Unsupported'):
        work(f,r,tmp_path/'out',min_obs=1,test_only=True)
    with pytest.raises(ValueError,match='Unsupported'):
        work(f,r,tmp_path/'out2',latest_climate_year=2026,test_only=True)
