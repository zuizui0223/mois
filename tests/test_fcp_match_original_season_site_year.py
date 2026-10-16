"""Synthetic validation of source-month matched-season FCP post-ROI comparisons.

NOT measured plants, NOT proof of within-population seasonality.
"""
from __future__ import annotations

import csv
import json

import pytest

from scripts.analysis.fcp_match_original_season_site_year import (
    INCOLS, run, sha, local_season_year,
)


def source(tmp_path):
    qc=tmp_path/"monthly.csv"
    receipt=tmp_path/"monthly_receipt.json"
    def row(tax,species,ym,hem,season,lat,lon,white,other):
        n=white+other
        return dict(
            inat_taxon_id=str(tax),species=species,year_month=ym,
            calendar_month=str(int(ym[-2:])),hemisphere=hem,season=season,
            lat5=str(lat),lon5=str(lon),
            dated_original_photo_opportunities=str(n+2),
            new_ROI_photos_with_known_date=str(n),
            classifiable_flower_ROIs=str(n),
            same_ROI_white_yellow_organ_mixed="0",
            unclassifiable_ROIs="0",usable_flower_photo_labels=str(n),
            independent_usable_observations=str(n),
            white_observations=str(white),nonwhite_observations=str(other),
            contradictory_multiphoto_observations="0",
            distinct_usable_observers=str(min(3,n)),
            white_rate_per_observation_descriptive=f"{white/n:.8f}" if n else "")
    rows=[
        # Northern winter 2026 includes Dec 2025 and Jan 2026 at same cell.
        row(101,"north","2025-12","north","local_winter",25,63,2,8),
        row(101,"north","2026-01","north","local_winter",25,63,1,9),
        row(101,"north","2026-07","north","local_summer",25,63,8,2),
        # Second matched broad geography cell supports same seasons.
        row(101,"north","2026-01","north","local_winter",26,64,0,10),
        row(101,"north","2026-07","north","local_summer",26,64,4,6),
        # Unmatched, heavily oversampled summer geography MUST NOT enter
        # the matched standardized comparison; useful Simpson's paradox fixture.
        row(101,"north","2026-07","north","local_summer",27,65,90,10),
        # Southern Hemisphere December/January is local summer of next
        # meteorological season-year, June/July is winter.
        row(104,"south","2025-12","south","local_summer",11,34,1,9),
        row(104,"south","2026-01","south","local_summer",11,34,1,9),
        row(104,"south","2026-07","south","local_winter",11,34,7,3),
        # Tropics must remain monthly; arbitrary four-season comparison barred.
        row(105,"trop","2026-01","tropical","tropical_month_only",18,65,5,5),
    ]
    with qc.open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=INCOLS,lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    receipt.write_text(json.dumps({
        "status":"SYNTHETIC_TEST_NO_BIOLOGICAL_FINDING",
        "source_snapshot_date":"2026-09-27",
        "monthly_site_QC_csv_sha256":sha(qc),
        "original_source_is_not_confirmed_flowers":True,
        "roi_failed_is_not_called_nonflowering":True,
        "white_fraction_uses_independent_original_observations_not_photos":True,
        "calendar_date_was_never_used_to_select_or_classify_original_pixels":True,
        "all_93209351_photo_pixels_classified":False,
    }))
    return qc,receipt


def rows(path):
    with path.open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))


def test_months_in_meteorological_season_year_and_hemisphere():
    assert local_season_year("2025-12","north","local_winter")==2026
    assert local_season_year("2026-01","north","local_winter")==2026
    assert local_season_year("2025-12","south","local_summer")==2026
    assert local_season_year("2026-07","south","local_winter")==2026
    assert local_season_year("2026-01","tropical","tropical_month_only")==2026
    with pytest.raises(ValueError,match="Tropical"):
        local_season_year("2026-01","tropical","local_summer")
    with pytest.raises(ValueError,match="Hemisphere"):
        local_season_year("2026-07","south","local_summer")


def test_pooled_versus_matched_geography_year_descriptive_rates(tmp_path):
    qc,receipt=source(tmp_path)
    out=tmp_path/"results"
    r=run(qc,receipt,out,_test_only=True)
    assert r["status"]=="SYNTHETIC_ONLY_NOT_BIOLOGY"
    assert r["taxa_with_at_least_one_valid_original_source_month"]==3
    assert r["n_input_original_month_cells"]==10
    assert r["n_merged_season_year_cell_rows"]==8
    assert r["season_year_cell_csv_sha256"]==sha(
        out/"observed_month_to_local_season_cell_year.csv")
    assert r["matched_season_contrasts_csv_sha256"]==sha(
        out/"within_species_matched_coarse_cell_season_contrasts.csv")
    source_rows=rows(out/"observed_month_to_local_season_cell_year.csv")
    decjan=next(x for x in source_rows if x["inat_taxon_id"]=="101" and
                x["local_season"]=="local_winter" and x["lat5"]=="25")
    assert decjan["season_year"]=="2026"
    assert decjan["n_calendar_months"]=="2"
    assert decjan["independent_flower_observations"]=="20"
    assert decjan["white_observations"]=="3"
    matched=rows(out/"within_species_matched_coarse_cell_season_contrasts.csv")
    n=next(x for x in matched if x["inat_taxon_id"]=="101" and
           x["season_a"]=="local_winter" and x["season_b"]=="local_summer")
    assert n["comparison_status"]=="MATCHED_BROAD_CELL_YEAR_DESCRIPTIVE_ONLY"
    assert n["matched_5degree_cell_seasonyear_strata"]=="2"
    assert n["matched_independent_observations_a"]=="30"
    assert n["matched_independent_observations_b"]=="20"
    assert n["balanced_min_observation_weight"]=="20"
    assert abs(float(n["balanced_white_fraction_a"])-.075)<1e-8
    assert abs(float(n["balanced_white_fraction_b"])-.6)<1e-8
    assert abs(float(n["balanced_white_fraction_difference_a_minus_b"])+.525)<1e-8
    assert abs(float(n["unmatched_species_white_fraction_a"])-.1)<1e-8
    assert abs(float(n["unmatched_species_white_fraction_b"])-.85)<1e-8
    assert n["causal_or_genetic_inference"]=="False"
    assert n["source_observers_are_not_independent_across_seasons_verified"]=="False"
    missing=next(x for x in matched if x["inat_taxon_id"]=="101" and
                 x["season_a"]=="local_spring" and x["season_b"]=="local_autumn")
    assert missing["comparison_status"]=="NO_MATCHED_CELL_YEAR_NO_TEMPORAL_IDENTIFICATION"
    assert missing["balanced_white_fraction_difference_a_minus_b"]==""
    s=next(x for x in matched if x["inat_taxon_id"]=="104" and
           x["season_a"]=="local_winter" and x["season_b"]=="local_summer")
    assert s["matched_5degree_cell_seasonyear_strata"]=="1"
    assert abs(float(s["balanced_white_fraction_a"])-.7)<1e-8
    assert abs(float(s["balanced_white_fraction_b"])-.1)<1e-8
    trop=[x for x in matched if x["inat_taxon_id"]=="105"]
    assert all(x["comparison_status"]=="NO_MATCHED_CELL_YEAR_NO_TEMPORAL_IDENTIFICATION"
               for x in trop)
    with pytest.raises(ValueError,match="overwrite"):
        run(qc,receipt,out,_test_only=True)


def test_source_month_data_tampering_cannot_be_used_with_old_sha(tmp_path):
    qc,receipt=source(tmp_path)
    qc.write_text(qc.read_text().replace("local_winter","local_autumn",1))
    with pytest.raises(ValueError,match="not attested"):
        run(qc,receipt,tmp_path/"bad",_test_only=True)


def test_synthetic_source_cannot_be_passed_to_production(tmp_path):
    qc,receipt=source(tmp_path)
    with pytest.raises(ValueError,match="not attested"):
        run(qc,receipt,tmp_path/"bad",_test_only=False)


def test_invalid_month_hemisphere_guard_even_with_updated_sha(tmp_path):
    qc,receipt=source(tmp_path)
    qc.write_text(qc.read_text().replace("2025-12","2025-13",1))
    doc=json.loads(receipt.read_text())
    doc["monthly_site_QC_csv_sha256"]=sha(qc)
    receipt.write_text(json.dumps(doc))
    with pytest.raises(ValueError):
        run(qc,receipt,tmp_path/"bad",_test_only=True)


def test_impossible_observation_denominator_denied(tmp_path):
    qc,receipt=source(tmp_path)
    qc.write_text(qc.read_text().replace("10,2,10,10,0,0,10,10,2,8", "10,2,10,10,0,0,10,12,2,8"))
    # A modified source SHA must not hide the original receipt, even if
    # string replacement happens to not match due to CSV formatting.
    doc=json.loads(receipt.read_text())
    if sha(qc)==doc["monthly_site_QC_csv_sha256"]:
        with qc.open(newline="") as h:
            rows_list=list(csv.DictReader(h))
        rows_list[0]["independent_usable_observations"]="9000"
        with qc.open("w",newline="") as h:
            writer=csv.DictWriter(h,fieldnames=INCOLS)
            writer.writeheader()
            writer.writerows(rows_list)
    doc["monthly_site_QC_csv_sha256"]=sha(qc)
    receipt.write_text(json.dumps(doc))
    with pytest.raises(ValueError,match="denominators"):
        run(qc,receipt,tmp_path/"bad",_test_only=True)
