"""Check publication figure data are byte-identical to original full256 FCP source.

The raw 309439-photo original-ID archive remains the 90-day signed GitHub
Actions artifact; public Git includes only its verified aggregate photo counts.
"""
from __future__ import annotations
import csv
import hashlib
import json
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
RESULT=ROOT/"results/mois_real309439_fullFCP256_original_photo_date_hue_season_decision_20261017.json"
SEASONS=ROOT/"results/mois_309439_fcp_full256_original_photo_season_figure_data_20261017.csv"
LOCALS=ROOT/"results/mois_309439_fcp_full256_locality_year_season_figure_data_20261017.csv"


def sample():
    j=json.loads(RESULT.read_text(encoding="utf-8"))
    with SEASONS.open(newline="",encoding="utf-8") as f:
        seasons=list(csv.DictReader(f))
    with LOCALS.open(newline="",encoding="utf-8") as f:
        spatial=list(csv.DictReader(f))
    return j,seasons,spatial


def test_source_announcement_is_exact_actual_finished_photo_pixel_measurement():
    a,s,l=sample()
    assert a["FCP_original_256_worker_partitions_completed"]==256
    assert a["FCP_new_wave4_original_49800_classifiable_ROIs"]==24585
    assert a["new_original_FCP_full_photo_IDs_after_exact_dedup"]==7197
    assert a["FCP_full_original_photo_ID_overlaps_with_MOIS_existing_measured"]==42603
    assert 7197+42603==49800
    assert a["new_automatic_hue_classifications_total"]==3601
    assert 1786+921+401+493==3601
    assert a["new_technical_unclassifiable_ROIs"]==3596
    assert 3601+3596==7197
    assert 302242+7197==a["full_MOIS_source_photo_attempts"]==309439
    assert 146061+3601==a["full_MOIS_source_four_hue_photo_labels"]==149662
    assert 146052+3601==a["full_MOIS_source_photo_ROIs_with_observation_day_and_four_hue"]==149653
    assert a["full_MOIS_source_photo_attempts_with_observation_day"]==309420
    assert a["full_MOIS_classifiable_hue_missing_observation_day"]==9
    assert a["same_photo_hue_taxon_identity_changes"]==0


def test_aggregated_figure_data_identical_to_verified_github_artifact_file_sha():
    a,s,l=sample()
    assert hashlib.sha256(SEASONS.read_bytes()).hexdigest()==a[
        "completed_MOIS_309439_hemisphere_season_csv_sha256"]
    assert hashlib.sha256(LOCALS.read_bytes()).hexdigest()==a[
        "completed_MOIS_309439_locality_season_year_source_csv_sha256"]
    assert len(s)==8 and len(l)==8
    assert len({(x["source_hemisphere"],x["meteorological_season"]) for x in s})==8
    for x in s:
        n=int(x["classified_original_photo_ROIs"])
        count=int(x["original_white_label_ROIs"])
        assert 0<=count<=n<=int(x["original_photo_attempts"])
        assert round(count/n,8)==pytest.approx(
            float(x["white_fraction_among_actual_classifiable_ROIs"]),abs=0.00000001)


@pytest.mark.parametrize("zone,raw_delta,local_delta,groups",[
    ("north_extratropical",2.841928,6.695455,27),
    ("south_extratropical",11.162867,4.492673,28),
])
def test_photo_census_global_seasonal_result_not_mistaken_for_local_weather(zone,raw_delta,local_delta,groups):
    j,seasons,geo=sample()
    q={(x["source_hemisphere"],x["meteorological_season"]):x for x in seasons}
    summer=q[(zone,"summer")]
    winter=q[(zone,"winter")]
    global_pp=100*(int(summer["original_white_label_ROIs"])/int(summer["classified_original_photo_ROIs"])-
                   int(winter["original_white_label_ROIs"])/int(winter["classified_original_photo_ROIs"]))
    assert global_pp==pytest.approx(raw_delta,abs=0.00001)
    x=[x for x in geo if x["hemisphere"]==zone and float(x["cell_degrees"])==0.25][0]
    assert float(x["equal_local_cell_year_summer_minus_winter_white_pp"])==pytest.approx(local_delta,abs=0.00001)
    assert int(x["same_cell_same_season_year_photo_community_groups"])==groups
    interval=json.loads(x["exploratory_999_geocell_bootstrap_percentile95_pp"])
    assert interval[0]<0<interval[1]
    assert j["local_original_photo_community_same_geocell_and_meteorological_season_year"][zone][
        "quarter_degree_site_year_groups"]==groups
