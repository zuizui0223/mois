"""Pre-event evidence guard: historical reporting index cannot use same-year/future first records."""
import numpy as np
import pandas as pd
from scripts.analysis.mois_firstrecords_v40_historical_effort_audit_20261124 import (
    canonical_v4_plant_records,source_opportunity_for_species,prepare_sample,score_sample
)


def sample_source():
    taxa=pd.DataFrame({
        'taxonID':['1','2','3','4'],
        'taxaGroup':['Vascular plants']*4,
        'GBIFstatus':['ACCEPTED']*4,
        'GBIFmatchType':['EXACT']*4,
        'GBIFtaxonRank':['SPECIES']*4,
        'species':['A alpha','A beta','A gamma','A delta'],
        'genus':['A']*4,
    })
    loc=pd.DataFrame({
        'locationID':['AA','BB','CC'],
        'location':['region1','region2','region3'],
        'continent':['Asia','Asia','Europe'],
    })
    table=[('1','AA','1900','1900','established'),
           ('1','BB','1900','1900','established'),
           ('2','AA','1890','1890','established'),
           ('3','CC','1899','1899','established'),
           ('4','BB','1910','1910','casual'),
           ('4','BB','1885','1885','established')]
    raw=pd.DataFrame(table,columns=['taxonID','locationID','firstRecordEvent','verbatimFirstRecordEvent','degreeOfEstablishment'])
    raw['confidenceFirstRecordEvent']='high confidence'
    raw['occurrenceStatus']='present'
    raw['establishmentMeans']='introduced'
    raw['accessRights']='Free'
    strict,anyfirst,audit=canonical_v4_plant_records(raw,taxa,loc)
    return strict,anyfirst,audit


def test_historical_opportunity_excludes_current_and_future_year():
    strict,_,audit=sample_source()
    assert strict.loc[strict.species.eq('A alpha')].shape[0]==2
    assert strict.loc[strict.species.eq('A delta')].year.tolist()==[1885]
    focal=pd.DataFrame({'key':['k1','k2','k3','k4'],
         'species':['A alpha','A beta','A gamma','A delta'],
         'first_year':[1900,1910,1900,1915],
         'first_record_continent':['Asia','Asia','Europe','MULTIPLE_CONTINENTS']})
    value,present=source_opportunity_for_species(focal,strict)
    # In Asia at 1900 the 1890 and 1885 source-first-record events qualify,
    # but the two 1900 region events do NOT qualify.
    assert value[0,0]==2
    assert value[0,1]==2
    # In Asia at 1910 two 1900 events AND the 1890 and 1885 events qualify.
    assert value[1,0]==4
    assert value[1,1]==2
    # Europe at 1900: previous 1899 event only.
    assert value[2,0]==1
    assert not present[3]


def test_first_record_is_not_always_first_year_in_all_source():
    strict,first_any,_=sample_source()
    assert first_any.loc['A alpha']==1900
    assert first_any.loc['A delta']==1885
    focal=pd.DataFrame({'key':['k1','k2'],'species':['A alpha','A delta'],
      'first_year':[1900,1915], 'first_record_continent':['Asia','Asia']})
    frame=prepare_sample(focal,strict,first_any)
    assert frame.source_has_earlier_any_year.tolist()==[False,True]
    assert frame.past25_source_records.notna().all()


def test_insufficient_sample_not_pseudopvalue():
    f=pd.DataFrame({'key':['a'],'species':['A alpha'],'genus':['A'],
                    'any_new_regions_50yr':[1]})
    result,score=score_sample(f)
    assert result['status'].startswith('insufficient') and score is None