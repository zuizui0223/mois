"""MOIS original TerraClimate v1.1 CF-packed NCSS scientific unit regressions.

No fake source weather is promoted to real climate. No flower-organ outcomes.
"""
from __future__ import annotations
import csv,io,json,math
from pathlib import Path
import pytest

from scripts.analysis.mois_fetch_terraclimate_v11_2025_by_original_month import (
    metadata,unpack,csv_series,run,sha
)


def attributes(var="tmin",version="V1.1",scale=.01,offset=-99.0):
    return f'''Attributes {{
  {var} {{
    String long_name "air_temperature";
    String units "degC";
    Float64 add_offset {offset};
    Int16 _FillValue -32768;
    Float64 scale_factor {scale};
  }}
  NC_GLOBAL {{
    String version "{version}";
  }}
}}'''


def series_file():
    f=io.StringIO()
    w=csv.writer(f,lineterminator="\n")
    w.writerow(["time",'''latitude[unit="degrees_north"]''',
                '''longitude[unit="degrees_east"]''','''tmin[unit="degC"]'''])
    for year in range(1991,2026):
        for month in range(1,13):
            raw=11750 if (year,month)==(2025,7) else 11500
            w.writerow([f"{year}-{month:02d}-01T00:00:00Z",42.375,-71.125,raw])
    return f.getvalue().encode()


def test_2025_original_provider_ncss_raw_tmin_11750_is_18p5_not_11750():
    v=metadata("tmin",attributes())
    assert v['dataset_version']=="V1.1"
    assert v['scale']==.01
    assert v['offset']==-99.
    assert unpack("tmin",11750,v)==pytest.approx(18.5)
    vals=csv_series("tmin",series_file(),v)
    assert len(vals)==420
    assert vals["2025-07"]==pytest.approx(18.5)
    assert vals["1991-07"]==pytest.approx(16.0)


def test_new_v11_radiation_has_independent_factor():
    v=metadata("srad",attributes("srad",scale=.1,offset=0))
    assert unpack("srad",1930,v)==pytest.approx(193.0)
    with pytest.raises(ValueError,match="physically"):
        unpack("srad",999999,v)


def test_reject_v10_or_missing_original_scale():
    with pytest.raises(ValueError,match="v1.1"):
        metadata("tmin",attributes(version="V1.0"))
    with pytest.raises(ValueError,match="scale_factor"):
        metadata("tmin",attributes().replace("scale_factor","unknown_factor"))
    with pytest.raises(ValueError,match="CF metadata"):
        metadata("vpd",attributes())


def test_reject_fill_value_and_reject_unpacked_physical_impossibilities():
    v=metadata("tmin",attributes())
    with pytest.raises(ValueError,match="fill"):
        unpack("tmin",-32768,v)
    with pytest.raises(ValueError,match="physically"):
        unpack("tmin",30000,v)


def test_full_baseline_and_original_2025_month_required():
    v=metadata("tmin",attributes())
    raw=series_file().replace(b'2025-07-01T00:00:00Z',b'2024-07-01T00:00:00Z',1)
    with pytest.raises(ValueError,match="duplicate|Incomplete"):
        csv_series("tmin",raw,v)


def test_no_monthly_climate_with_unverified_botanical_or_source_SHA(tmp_path):
    rows=tmp_path/"fake.csv";rows.write_text("year_month,tc_lat_idx\n2025-07,1143\n")
    receipt=tmp_path/"receipt.json"
    receipt.write_text(json.dumps({
       "status":"REAL_DATE_AND_ROI_ATTEMPTS_NO_BOTANICAL_WHITE_CLIMATE_EFFECT",
       "original_month_climate_grid_csv_sha256":sha(rows.read_bytes()),
       "n_exact_2025_grid_and_month_climate_requests":40,
       "n_original_source_grid_cells":21,
       "open_flower_botanically_validated":True,
    }))
    with pytest.raises(ValueError,match="original"):
        run(rows,receipt,tmp_path/"climate")
