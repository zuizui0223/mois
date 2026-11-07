import pandas as pd
from scripts.analysis import mois_viirs_2020_population_match as m

def obj(val='1.23',cloud='5',year=2020,month=1):
    return {'value':val,'catalogItems':{'features':[{'attributes':{'year':year,'month':month}}]},'properties':{'Values':[val+' '+cloud]}}

def test_raw_radiance_requires_correct_month_and_cloudfree():
    assert m._parse_viirs_response(obj(),2020,1)==(1.23,5.,'OK')
    assert m._parse_viirs_response(obj(month=7),2020,1)[2]=='SOURCE_MONTH_DRIFT'
    assert m._parse_viirs_response(obj(cloud='0'),2020,1)[2]=='INSUFFICIENT_CLOUD_FREE'
    assert m._parse_viirs_response(obj(val='0',cloud='0'),2020,1)[0] is None
    assert m._parse_viirs_response(obj(val='0',cloud='5'),2020,1)==(0.,5.,'OK')

def test_time_month_february_has_real_leap_day():
    from datetime import datetime,timezone
    start,end=m.month_range_ms(2020,2).split(',')
    assert datetime.fromtimestamp(int(start)/1000,timezone.utc).day==1
    assert datetime.fromtimestamp(int(end)/1000,timezone.utc).day==29

def test_same_population_and_climate_but_different_light_can_be_paired():
    entries=[
      dict(cell_id='a',block_id='0_0',source_complete=True,viirs_quarterly_median=1.,population_density=100.,
           bio5=25.,bio12=1000.,elev=150.,latitude=1.,longitude=1.,
           white_taxon_equal=.6,strict9_yield=.7,n_measured_taxa=5),
      dict(cell_id='b',block_id='0_0',source_complete=True,viirs_quarterly_median=20.,population_density=105.,
           bio5=25.5,bio12=1050.,elev=170.,latitude=1.1,longitude=1.1,
           white_taxon_equal=.2,strict9_yield=.8,n_measured_taxa=5),
      dict(cell_id='c',block_id='0_0',source_complete=True,viirs_quarterly_median=80.,population_density=1.,
           bio5=25.5,bio12=1050.,elev=170.,latitude=1.1,longitude=1.1,
           white_taxon_equal=.2,strict9_yield=.8,n_measured_taxa=5)]
    p,r=m.match_cells(pd.DataFrame(entries))
    assert len(p)==1 and p.iloc[0]['bright_cell_id']=='b'
    assert p.iloc[0]['white_diff_bright_minus_dark_pp']==-40.
    assert r['uses_outcome_for_matching'] is False

def test_without_overlap_no_artificial_effect_is_returned():
    e=pd.DataFrame([dict(cell_id='a',block_id='0_0',source_complete=True,viirs_quarterly_median=5.,
                         population_density=100.,bio5=20.,bio12=1000.,elev=300.,
                         latitude=1.,longitude=1.,white_taxon_equal=.3,strict9_yield=.6,n_measured_taxa=4)])
    p,r=m.match_cells(e)
    assert p.empty and r['matched_pairs']==0

def test_under_20_matches_does_not_promote_effect():
    d=pd.DataFrame([dict(source_complete=True,n_attempted=5,n_strict9=3)])
    p=pd.DataFrame([dict(white_diff_bright_minus_dark_pp=-40.,strict_yield_diff_pp=15.,
                         bright_population_per_km2=100.,dark_population_per_km2=110.)])
    report=m.outcome_report(d,p,{})
    assert report['matched_cell_white_difference_mean_pp'] is None
    assert report['pairs']==1
    assert report['photographed_taxon_colour_composition_NOT_species_absence']

def test_observation_unknown_is_not_nonwhite():
    assert m._parse_viirs_response({'value':'0'},2020,1)[2]=='SOURCE_MONTH_NOT_SINGLE'
