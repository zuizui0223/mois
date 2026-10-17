"""Synthetic scientific guards for MOIS seasonal photo source and climates."""
from __future__ import annotations
import gzip
import hashlib
import io
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from scripts.analysis.mois_seasonal_267k_climate_normals import (
    ENERGY_PER_WATT_DAY_KJ, equal_group_within_estimate, labelled_calendar,
    matched_climate_panel, photo_season_table, read_signed,
)


def toy() -> pd.DataFrame:
    x=[]
    for tax in (1, 2):
        for season, mo, hue in [("winter", 1, "white"), ("winter", 2, "white"),
                                ("summer", 6, "red_pink"), ("summer", 7, "red_pink")]:
            x.append(dict(photo_id=str(tax*10+mo), observation_id=str(tax*10+mo),
                          inat_taxon_id=str(tax), morph=hue, latitude=38.3,
                          longitude=138.4, observed_on=f"2025-{mo:02d}-12",
                          wc2_month_tavg_normal_10m_C=1 if season=="winter" else 19,
                          wc2_month_prec_normal_10m_mm=30 if season=="winter" else 100,
                          # Raw 8640 kJ/m2/day is 100 W/m2, not 8640 W/m2.
                          wc2_month_srad_normal_10m_Wm2=8640 if season=="winter" else 17280))
    return labelled_calendar(pd.DataFrame(x))


def test_hemisphere_reversal_december_attaches_next_season_year():
    d=labelled_calendar(pd.DataFrame([
        dict(morph="white",observed_on="2024-12-15",latitude=40),
        dict(morph="white",observed_on="2025-01-15",latitude=40),
        dict(morph="white",observed_on="2024-12-15",latitude=-40),
        dict(morph="white",observed_on="2024-12-15",latitude=0),
    ]))
    assert d.season.tolist()==["winter","winter","summer","month_only"]
    assert d.season_year.tolist()==[2025,2025,2025,2025]


def test_min_two_observations_in_both_seasons():
    d=toy()
    panel=matched_climate_panel(d,2)
    assert panel.group.nunique()==2 and len(panel)==8
    assert matched_climate_panel(d,3).empty
    duplicates=d.copy()
    duplicates.loc[duplicates.observed_on.str.contains("-02-"),"observation_id"] = \
        duplicates.loc[duplicates.observed_on.str.contains("-01-"),"observation_id"].to_numpy()
    assert matched_climate_panel(duplicates,2).empty


def test_original_solar_energy_to_power_conversion_changes_units_not_data():
    q=matched_climate_panel(toy(),2)
    assert ENERGY_PER_WATT_DAY_KJ==86.4
    physical=equal_group_within_estimate(q,"wc2_month_srad_normal_10m_Wm2",50,86.4,bootstrap=20)
    raw=equal_group_within_estimate(q,"wc2_month_srad_normal_10m_Wm2",50,1,bootstrap=20)
    assert np.isclose(physical['point_white_pp_per_unit'],
                      raw['point_white_pp_per_unit']*86.4,atol=.0002)
    assert physical['point_white_pp_per_unit'] < 0


def test_roi_unknowns_not_black_white_or_flowering():
    d=toy()
    d.loc[0,"morph"]="mixed_uncertain"
    d["classified"] = d.morph.isin(("white", "yellow_orange", "red_pink", "blue_purple"))
    d["white"] = d.morph.eq("white")
    x=photo_season_table(d)
    w=x[x.season.eq("winter")].iloc[0]
    assert w.photo_attempts==4
    assert w.classified_hues==3
    assert w.white_hues==3


def test_signed_receipt_rejects_tampered_bytes(tmp_path:Path):
    fake=pd.DataFrame([{"photo_id":"9","morph":"white"}])
    csvbytes=fake.to_csv(index=False).encode()
    data=gzip.compress(csvbytes)
    rec={"test_sha":hashlib.sha256(data).hexdigest()}
    zpath=tmp_path/"source.zip"
    with zipfile.ZipFile(zpath,"w") as z:
        z.writestr("data.csv.gz",data)
        z.writestr("receipt.json",json.dumps(rec))
    good,_,_=read_signed(zpath,"data.csv.gz","receipt.json","test_sha",["photo_id","morph"])
    assert good.photo_id.tolist()==["9"]
    with zipfile.ZipFile(zpath,"w") as z:
        z.writestr("data.csv.gz",data + b"tampered")
        z.writestr("receipt.json",json.dumps(rec))
    with pytest.raises(ValueError,match="SHA"):
        read_signed(zpath,"data.csv.gz","receipt.json","test_sha",["photo_id","morph"])
