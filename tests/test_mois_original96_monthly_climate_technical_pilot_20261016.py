"""Synthetic original 96 photo-pilot site/month weather link, non-inference guard."""
from __future__ import annotations
import csv
import importlib.util
import io
import json
import zipfile
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
def mod(name,path):
    spec=importlib.util.spec_from_file_location(name,ROOT/path)
    x=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(x)
    return x
fake=mod("mois_dynamic_src_fixture", "tests/test_mois_join_monthly_dynamic_climate_20261016.py")
m=mod("mois_original96_jointest", "scripts/analysis/mois_original96_monthly_climate_technical_pilot_20261016.py")


def original(tmp):
    climate=fake.normal(tmp)
    rows=[]
    for i in range(96):
        if i<20:
            hue="white" if i in (0,1,2) else ("blue_purple" if i<14 else "yellow_orange")
            status="classified_four_state_morph"
        else:
            hue="mixed_uncertain"
            status=("not_evaluable_roi_or_flip_gate" if i<87
                    else "not_evaluable_ambiguous_palette_composition")
        rows.append({
           "inat_taxon_id":str(101+i//24),"species":"Species_"+str(101+i//24),
           "observation_id":"OBS_"+str(i),"photo_id":str(1000+i),
           "measurement_id":"FCPR-"+str(i),"latitude":"35.7","longitude":"139.7",
           "observed_on":"2020-04-09" if i<40 else "2020-07-09",
           "morph":hue,"measurement_status":status,
        })
    output=io.StringIO()
    writer=csv.DictWriter(output,fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
    raw=output.getvalue().encode()
    receipt={
        "schema":"fcp_real_original_roi_v4_case_sensitive_extension_repair_v1",
        "status":"REPAIRED_REAL_PHOTO_TECHNICAL_PILOT_NOT_42111_CENSUS",
        "source_photo_id_count":96,
        "n_actual_original_image_measurements_total":96,
        "n_classified_four_state_photo_rows":20,
        "terminal_csv_sha256":m.hashlib.sha256(raw).hexdigest(),
        "original_method":{"location_blind":True}
    }
    archive=tmp/"source96.zip"
    with zipfile.ZipFile(archive,"w") as z:
        z.writestr(m.SOURCE_CSV,raw)
        z.writestr(m.SOURCE_RECEIPT,json.dumps(receipt))
    return climate,archive


def test_exact_original_four_species_photo_source_classifiability_and_monthly_climate(tmp_path):
    climate,archive=original(tmp_path)
    out=m.audit(archive,climate,tmp_path/"aggregate.json",synthetic=True)
    assert out["status"]==m.STATUS
    assert out["n_original_source_photo_IDs"]==96
    assert out["n_original_distinct_observations"]==96
    assert out["n_original_taxa"]==4
    assert out["n_roi_original_four_colour_classifiable"]==20
    assert out["n_roi_original_classified_white"]==3
    assert out["n_monthly_source_site_climate_joined"]==96
    assert out["real_10m_source_climate_cells_with_any_original_photo"]==1
    assert out["real_10m_source_climate_cells_with_original_photos_in_two_or_more_months"]==1
    assert out["real_10m_source_climate_cells_with_classifiable_photos"]==1
    assert out["real_10m_source_climate_cells_with_classifiable_photos_in_two_or_more_months"]==0
    assert out["cell_FE_primary_four_dynamic_months_estimable_in_this_pilot"] is False
    assert out["prespecified_site_cell_season_FE_within_month_estimated"] is False
    assert out["original_source_population_42111_taxa_photo_IDs_joined"] is False
    rows=out["observed_season_rows"]
    assert len(rows)==2
    assert sum(x["sampled_original_photo_opportunities"] for x in rows)==96
    assert sum(x["classified_white"] for x in rows)==3
    spring=next(x for x in rows if x["local_season"]=="SPRING")
    assert spring["four_hue_classifiable"]==20
    assert spring["white_fraction_among_classifiable"]==.15
    assert spring["n_classifiable_with_known_monthly_climate"]==20
    assert spring["source_classified_monthly_climate_means"]["monthly_tavg_c"]>15
    assert "measurement_id" not in (tmp_path/"aggregate.json").read_text()
    with pytest.raises(ValueError,match="overwrite"):
        m.audit(archive,climate,tmp_path/"aggregate.json",synthetic=True)


def test_fake_source_cannot_promote_synthetic_monthlies_as_real_original(tmp_path):
    climate,archive=original(tmp_path)
    with pytest.raises(ValueError,match="Synthetic raster"):
        m.audit(archive,climate,tmp_path/"aggregate.json",synthetic=False)
    with zipfile.ZipFile(archive,"a") as z:
        z.writestr("unexpected_extra.txt","not from original archive")
    with pytest.raises(ValueError,match="source pilot ZIP"):
        m.audit(archive,climate,tmp_path/"aggregate.json",synthetic=True)
