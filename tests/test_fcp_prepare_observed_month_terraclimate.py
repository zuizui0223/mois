"""Exact-observation-date climate requests: synthetic provenance contracts only.

Never downloads real provider pixels or historical 149,900 colour outcomes.
"""
from __future__ import annotations

import csv
import json

import pytest

from scripts.analysis.fcp_prepare_observed_month_terraclimate import (
    REQ_COLS, CLIM_COLS, prepare, sha, grid, observed_month, precise,
)


def make_fixture(tmp_path):
    terminal=tmp_path/"genuine_only_in_production.csv"
    roidoc=tmp_path/"roi_receipt.json"
    monthdoc=tmp_path/"dated_site_receipt.json"
    specs=[
        # Photo 1/2 are ONE real observation with two matching image labels
        (1,101,"A", "o1","u1","2024-08-10",35.,135.,150,"white",1.,0.),
        (2,101,"A", "o1","u1","2024-08-10",35.,135.,150,"white",1.,0.),
        (3,101,"A", "o2","u2","2024-08-10",35.,135.,100,"yellow_orange",0.,1.),
        # Same inflorescence has white ray and yellow disc: no morph assignment
        (4,101,"A", "o3","u3","2024-08-10",35.,135.,150,"white",.4,.4),
        (5,101,"A", "o4","u4","2024-08-10",35.,135.,150,"mixed_uncertain",0.,0.),
        (6,101,"A", "o5","u5","2026-04-01",35.,135.,100,"white",1.,0.),
        (7,104,"B", "o6","u6","2024-07-10",-33.,180.,0,"blue_purple",0.,0.),
        (8,104,"B", "o7","u7","2024-07-10",-33.,-180.,200,"red_pink",0.,0.),
        (9,104,"B", "o8","u8","2024-08",-33.,20.,150,"white",1.,0.),
        (10,104,"B","o9","u9","2024-07-10",60.,120.,150,"white",1.,0.),
        (11,104,"B","o9","u9","2024-07-10",60.,120.,150,"yellow_orange",0.,1.),
        (12,104,"B","o10","u10","1957-10-01",60.,120.,150,"white",1.,0.),
    ]
    cols=["inat_taxon_id","species","photo_id","observation_id","observer_id",
          "observed_on","latitude","longitude","positional_accuracy_m",
          "measurement_status","morph","flower_fraction_white",
          "flower_fraction_yellow","source_observation_month"]
    with terminal.open("w",newline="",encoding="utf-8") as h:
        writer=csv.DictWriter(h,fieldnames=cols,lineterminator="\n")
        writer.writeheader()
        for pid,tid,name,obs,observer,day,lat,lon,accuracy,morph,white,yellow in specs:
            good=morph!="mixed_uncertain"
            writer.writerow({
                "inat_taxon_id":str(tid),"species":name,
                "photo_id":str(pid),"observation_id":obs,"observer_id":observer,
                "observed_on":day,"latitude":str(lat),"longitude":str(lon),
                "positional_accuracy_m":str(accuracy),
                "measurement_status":("classified_four_state_morph" if good
                                      else "not_evaluable_roi_or_flip_gate"),
                "morph":morph,"flower_fraction_white":str(white),
                "flower_fraction_yellow":str(yellow),
                "source_observation_month":observed_month(day),
            })
    roidoc.write_text(json.dumps({
        "status":"SYNTHETIC_ORIGINAL_GEOGRAPHIC_PAIR_TEST_NOT_BIOLOGY",
        "original_geographic_source_join_terminal_sha256":sha(terminal),
        "new_original_photo_pixels_processed":len(specs),
        "old_measured_colours_imported":False,
    }))
    monthdoc.write_text(json.dumps({
        "status":"SYNTHETIC_TEST_NO_BIOLOGICAL_FINDING",
        "new_roi_terminal_sha256":sha(terminal),
        "source_snapshot_date":"2026-09-27",
        "new_roi_terminal_photos":len(specs),
        "original_source_is_not_confirmed_flowers":True,
        "white_fraction_uses_independent_original_observations_not_photos":True,
        "calendar_date_was_never_used_to_select_or_classify_original_pixels":True,
        "all_93209351_photo_pixels_classified":False,
    }))
    return terminal,roidoc,monthdoc


def read(path):
    with path.open(newline="",encoding="utf-8") as h:
        return list(csv.DictReader(h))


def reseal(term,roi,monthly):
    for path,key in ((roi,"original_geographic_source_join_terminal_sha256"),
                     (monthly,"new_roi_terminal_sha256")):
        q=json.loads(path.read_text())
        q[key]=sha(term)
        path.write_text(json.dumps(q))


def test_terra_24th_degree_grid_antimeridian_and_source_dates():
    assert grid(-33,180)[:2]==grid(-33,-180)[:2]
    assert grid(90,180)[:2]==(0,0)
    assert grid(-90,-180)[:2]==(4319,0)
    assert grid(0,0)[:2]==(2160,4320)
    assert observed_month("2024-01-09")=="2024-01"
    for raw in ("2024-01","2026-10-01","2026-00-12",
                "invalid","",None):
        assert observed_month(raw)==""
    assert precise("5000") is True
    assert precise("5001") is False
    assert precise("0") is False
    assert precise("") is False


def test_precise_observation_month_and_colour_confounding_gates(tmp_path):
    args=make_fixture(tmp_path)
    out=tmp_path/"out"
    r=prepare(*args,out,_synthetic_test=True)
    assert r["status"]=="SYNTHETIC_OBSERVED_MONTH_CLIMATE_REQUEST_NOT_BIOLOGY"
    assert r["original_independent_observations"]==10
    assert r["real_original_photo_pixels"]==12
    assert r["n_exact_source_month_raster_cells_to_request"]==3
    assert r["last_independently_verified_monthly_climate_year"]==2025
    assert r["no_2026_climate_values_fabricated_from_2025_or_1991_2020_normals"] is True
    assert r["no_climate_raster_download_or_temp_radiation_water_effect_estimated"] is True
    assert r["source_based_date_and_position_join_happens_after_blind_image_roi"] is True
    q=r["observations_by_climate_request_eligibility"]
    assert q["original_observation_month_grid_request_eligible_not_fetched"]==6
    assert q["source_position_uncertainty_above_5km_or_unknown"]==1
    assert q["missing_or_invalid_original_observation_day"]==1
    assert q["year_not_in_independently_verified_monthly_climate_release"]==1
    assert q["source_observation_year_before_supported_climate_history"]==1
    states=r["observations_by_photo_colour_phenotype"]
    assert states["white"]==4
    assert states["chromatic"]==3
    assert states["same_flower_white_yellow_organ_mix_unresolved"]==1
    assert states["unclassifiable_or_not_flower_unknown"]==1
    assert states["discordant_photo_morphs_unresolved"]==1
    obs=read(out/"real_original_roi_observation_terraclimate_month_keys.csv")
    assert len(obs)==10
    assert list(obs[0])==list(REQ_COLS)
    o1=next(x for x in obs if x["observation_id"]=="o1")
    assert o1["n_new_roi_photo_attempts"]=="2"
    assert o1["phenotype_status"]=="white"
    assert o1["climate_join_status"]=="original_observation_month_grid_request_eligible_not_fetched"
    o3=next(x for x in obs if x["observation_id"]=="o3")
    assert o3["phenotype_status"]=="same_flower_white_yellow_organ_mix_unresolved"
    o5=next(x for x in obs if x["observation_id"]=="o5")
    assert o5["climate_join_status"]=="year_not_in_independently_verified_monthly_climate_release"
    request=read(out/"exact_observed_month_terraclimate_grid_requests.csv")
    assert len(request)==3
    assert list(request[0])==list(CLIM_COLS)
    japan=next(x for x in request if x["year_month"]=="2024-08")
    assert japan["n_independent_original_observations"]=="4"
    assert japan["n_photo_roi_attempts"]=="5"
    assert japan["n_photo_supported_white_observations"]=="1"
    assert japan["n_photo_supported_chromatic_observations"]=="1"
    assert japan["n_original_observations_without_confident_flower_colour"]=="2"
    assert r["unique_terra_cell_month_requests_csv_sha256"]==sha(
        out/"exact_observed_month_terraclimate_grid_requests.csv")
    with pytest.raises(ValueError,match="overwrite"):
        prepare(*args,out,_synthetic_test=True)


def test_production_cannot_run_on_synthetic_pixels(tmp_path):
    args=make_fixture(tmp_path)
    with pytest.raises(ValueError,match="authenticated"):
        prepare(*args,tmp_path/"prod",_synthetic_test=False)


def test_source_date_rewritten_after_roi_is_rejected_by_receipts(tmp_path):
    term,roi,monthly=make_fixture(tmp_path)
    text=term.read_text()
    term.write_text(text.replace("2024-08-10","2024-09-10",1))
    with pytest.raises(ValueError,match="authenticated"):
        prepare(term,roi,monthly,tmp_path/"bad",_synthetic_test=True)


def test_fake_observation_month_in_resealed_pixel_table_fails(tmp_path):
    term,roi,monthly=make_fixture(tmp_path)
    rows=read(term)
    rows[0]["source_observation_month"]="2024-02"
    with term.open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n")
        writer.writeheader();writer.writerows(rows)
    reseal(term,roi,monthly)
    with pytest.raises(ValueError,match="previously recovered"):
        prepare(term,roi,monthly,tmp_path/"bad",_synthetic_test=True)


def test_one_observation_with_changed_observer_is_not_independent(tmp_path):
    term,roi,monthly=make_fixture(tmp_path)
    text=term.read_text().replace("2,A,2,o1,u1,","2,A,2,o1,replaced,")
    # Editing literal CSV count or changing mock data schema must not create
    # a valid second source observer for the same original observation.
    if text==term.read_text():
        rows=read(term)
        rows[1]["observer_id"]="replaced"
        with term.open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]))
            w.writeheader();w.writerows(rows)
    else:
        term.write_text(text)
    reseal(term,roi,monthly)
    with pytest.raises(ValueError,match="disagree on date/site/observer"):
        prepare(term,roi,monthly,tmp_path/"bad",_synthetic_test=True)


def test_duplicate_photo_id_is_rejected_even_with_updated_sealed_receipts(tmp_path):
    term,roi,monthly=make_fixture(tmp_path)
    rows=read(term)
    rows[-1]["photo_id"]="1"
    with term.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n")
        w.writeheader();w.writerows(rows)
    reseal(term,roi,monthly)
    with pytest.raises(ValueError,match="duplicate"):
        prepare(term,roi,monthly,tmp_path/"bad",_synthetic_test=True)


def test_unverified_climate_year_is_never_replaced_by_later_or_earlier_year(tmp_path):
    args=make_fixture(tmp_path)
    r=prepare(*args,tmp_path/"out",last_verified_climate_year=2024,_synthetic_test=True)
    assert r["last_independently_verified_monthly_climate_year"]==2024
    assert r["observations_by_climate_request_eligibility"][
        "year_not_in_independently_verified_monthly_climate_release"]==1
    assert r["n_exact_source_month_raster_cells_to_request"]==3
    with pytest.raises(ValueError,match="Invalid independently verified"):
        prepare(*args,tmp_path/"bad",last_verified_climate_year=1900,_synthetic_test=True)
