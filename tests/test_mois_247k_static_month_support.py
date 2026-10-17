"""MOIS original measured photo static and calendar month source support."""
import pandas as pd
import numpy as np
import pytest
from scripts.analysis.mois_247k_static_month_support_analysis import local_season, classify_support

def test_hemispheres_and_tropics_not_equated():
    assert local_season(7,42.5)=="summer"
    assert local_season(7,-42.5)=="winter"
    assert local_season(1,-35)=="summer"
    assert local_season(7,0)=="tropical_month_only"
    assert local_season(0,42.5)=="unknown"
    assert local_season(7,float("nan"))=="unknown"

def test_no_missing_hour_is_fabricated():
    # Full original photo columns are required; the historical archive has
    # no time-of-day; date must not be converted to an invented local hour.
    a={
        "photo_id":["100","101"],"inat_taxon_id":["5","5"],
        "species":["Example sp.","Example sp."],
        "morph":["white","mixed_uncertain"],"old_mois_or_new_fcp":["historical","new_fcp"],
        "observation_id":["obs1","obs2"],
        "observed_on":["","2025-07-21"],"time_observed_at":["",""],
    }
    from scripts.analysis.mois_247k_static_month_support_analysis import STATIC,CLIMATOLOGY
    for k in (*STATIC,*CLIMATOLOGY): a[k]=[np.nan,np.nan]
    a["latitude"]=[-35,42.5]
    a["longitude"]=[148,-71]
    f,c=classify_support(pd.DataFrame(a))
    assert c["historical"]["observation_calendar_day_available"]==0
    assert c["new_fcp"]["observation_calendar_day_available"]==1
    assert all(x["observation_clock_time_validated"]==0 for x in c.values())
    assert f.loc[1,"local_season"]=="summer"
    assert f["is_classified"].sum()==1

def test_duplicate_source_photo_must_fail():
    from scripts.analysis.mois_247k_static_month_support_analysis import STATIC,CLIMATOLOGY
    r=pd.DataFrame({
      "photo_id":["10","10"],"morph":["white","blue_purple"],
      "observed_on":["2023-06-01","2023-06-01"],
      "time_observed_at":["",""],"latitude":[35,35],"longitude":[130,130],
      "old_mois_or_new_fcp":["historical","new_fcp"],
    })
    for k in STATIC+CLIMATOLOGY:
        if k not in r:r[k]=np.nan
    with pytest.raises(ValueError,match="Duplicate"):
        classify_support(r)
