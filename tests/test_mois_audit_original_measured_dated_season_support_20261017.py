"""No false blooming or independent seasonal climate effect from MOIS 247k measured photo IDs."""
from __future__ import annotations
import csv,gzip,hashlib,io,json,zipfile
from pathlib import Path

import pytest

from scripts.analysis.mois_audit_original_measured_dated_season_support_20261017 import (
    PHOTO_SHA,UNION_IDS,MONTH_COLS,SUMMARY_COLS,MATCH_COLS,
    season,grid,audit,sha
)


def union(tmp_path,*,photo_replacement=None):
    zip_file=tmp_path/"original_measured.zip"
    columns=("photo_id","inat_taxon_id","species","morph","source_original_measured",
             "old_mois_or_new_fcp","observation_id","observed_on","latitude",
             "longitude","measurement_status")
    rows=[]
    for i in range(1,5):
        rows.append(dict(
            photo_id=str(500000+i),inat_taxon_id=str(101+i),species=f"Old taxon {i}",
            morph=("white" if i==1 else "yellow_orange" if i==2 else
                   "mixed_uncertain" if i==3 else ""),
            source_original_measured="historical_mois_20261014",
            old_mois_or_new_fcp="historical",observation_id="old"+str(i),
            observed_on="",latitude="35.0",longitude="-72.5",measurement_status="",
        ))
    for i in range(1,5):
        rows.append(dict(
            photo_id=str(600000+i),inat_taxon_id="100",species="New taxon test",
            morph=("white" if i in (1,3) else "red_pink"),
            source_original_measured="fcp_discovery",old_mois_or_new_fcp="new_fcp",
            observation_id="source-original-uuid-"+str(i),
            observed_on=("2025-06-01" if i<=2 else "2025-09-01"),
            latitude="40.0",longitude="-72.5",
            measurement_status="classified_four_state_morph",
        ))
    if photo_replacement:
        photo_replacement(rows)
    buff=io.StringIO()
    w=csv.DictWriter(buff,fieldnames=columns,lineterminator="\n")
    w.writeheader();w.writerows(rows)
    b=gzip.compress(buff.getvalue().encode("utf-8"),mtime=0)
    receipt={
        "status":"REAL_FROZEN_FOUR_HUE_PHOTO_UNION_NO_NEW_IMAGE_PIXELS",
        "combined_unique_photo_ids":8,
        "historical_global_genus_moisture_slopes_reestimated":False,
        "cross_source_measurement_independence":False,
        "new_fcp_photo_site_month_climate_join_complete":False,
        "photo_hue_union_csv_sha256":hashlib.sha256(b).hexdigest(),
    }
    with zipfile.ZipFile(zip_file,"w") as z:
        z.writestr("mois_fcp_combined_original_247032_photo_hues.csv.gz",b)
        z.writestr("mois_fcp_measured_photo_union_receipt.json",json.dumps(receipt))
    return zip_file


def csvrows(path):
    with path.open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))


def test_source_hemisphere_and_december_meteorological_season_year():
    assert season(2022,12,40)==("winter",2023)
    assert season(2022,12,-40)==("summer",2023)
    assert season(2023,7,-40)==("winter",2023)
    assert season(2025,7,5)==("tropical_month_only",2025)
    assert grid(40.,-72.5,.25)==(520,430)
    assert grid(90.,180.,.25)==(719,0)


def test_real_measured_photo_source_dated_four_state_support_is_photo_not_botany(tmp_path):
    z=union(tmp_path)
    a=audit(z,tmp_path/"out",test_only=True)
    assert a["status"]=="REAL_MEASURED_PHOTO_SEASONAL_SUPPORT_NOT_BOTANICAL_OPEN_FLOWER_SEASONAL_EFFECT"
    assert a["n_unique_original_photo_ids"]==8
    assert a["source_photo_and_automatic_colour_funnel"]["historical_photos"]==4
    assert a["source_photo_and_automatic_colour_funnel"]["new_fcp_photos"]==4
    assert a["source_photo_and_automatic_colour_funnel"]["new_fcp_fourstate_ROIs"]==4
    assert a["n_independent_original_colour_classified_observations"]==4
    assert a["fine_grid_within_species_siteyear_climate_scope_min2_pair_count"]==1
    assert a["historical_mois_100543_observed_on_missing_in_frozen_union"] is True
    assert a["source_photographer_observer_disjointness_not_proven_in_union"] is True
    assert a["source_photo_vs_true_flower_organ_anthesis_known"] is False
    assert a["no_genetic_polymorphism_or_true_blooming_white_flower_rate_estimated"] is True
    assert a["no_monthly_weather_link_or_causal_eco_effect_performed"] is True
    assert sha(tmp_path/"out"/"mois_new_fcp_measured_photo_original_observed_month_ROI_opportunities.csv")==a["dated_photo_month_ROI_opportunity_sha256"]
    summary=csvrows(tmp_path/"out"/"mois_measured_source_season_matched_observation_support_by_spatial_resolution.csv")
    assert len(summary)==4*2*5
    assert set(summary[0])==set(SUMMARY_COLS)
    fine=next(x for x in summary if float(x["grid_width_degrees"])==.25 and
           x["years_in_scope"]=="climate_available_1958_to_2025" and
           int(x["minimum_independent_source_observations_in_each_of_two_seasons"])==2)
    assert fine["matched_taxon_grid_years"]=="1"
    assert fine["matched_named_taxa"]=="1"
    assert fine["source_coordinate_accuracy_verified"]=="False"
    strata=csvrows(tmp_path/"out"/"mois_matched_0p25_original_observation_siteyear_two_season_ROI_support.csv")
    assert len(strata)==1
    assert set(strata[0])==set(MATCH_COLS)
    assert {strata[0]["season_a"],strata[0]["season_b"]}=={"summer","autumn"}
    assert {int(strata[0]["independent_observations_a"]),int(strata[0]["independent_observations_b"])}=={2}
    assert {int(strata[0]["automated_photo_white_labels_a"]),int(strata[0]["automated_photo_white_labels_b"])}=={1}
    months=csvrows(tmp_path/"out"/"mois_new_fcp_measured_photo_original_observed_month_ROI_opportunities.csv")
    assert len(months)==2
    assert months[0]["year_month"]=="2025-06" and months[1]["year_month"]=="2025-09"
    assert [int(x["source_automatic_four_state_ROIs"]) for x in months]==[2,2]
    with pytest.raises(ValueError,match="overwrite"):
        audit(z,tmp_path/"out",test_only=True)


def test_changed_historical_photo_identifier_fails_same_source_identity(tmp_path):
    z=union(tmp_path,photo_replacement=lambda rows:rows[3].update(photo_id=rows[0]["photo_id"]))
    with pytest.raises(ValueError,match="duplicate"):
        audit(z,tmp_path/"out",test_only=True)


def test_fcp_2026_out_of_snapshot_original_day_fails_source_date(tmp_path):
    z=union(tmp_path,photo_replacement=lambda rows:rows[4].update(observed_on="2026-10-30"))
    with pytest.raises(ValueError,match="beyond snapshot"):
        audit(z,tmp_path/"out",test_only=True)


def test_original_source_400_id_overlap_must_not_fake_fourstate_labels(tmp_path):
    z=union(tmp_path)
    source=tmp_path/"400_source.csv"
    with source.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=["photo_id"])
        w.writeheader();w.writerow({"photo_id":"600001"});w.writerow({"photo_id":"600002"})
    r=audit(z,tmp_path/"out",400_source=source,test_only=True)
    assert r["source_400_original_photo_ID_intersection"]["new_fcp"]["photo_ids"]==2
    assert r["source_400_original_photo_ID_intersection"]["new_fcp"]["fourstate_labels"]==2
    assert r["source_400_original_photo_ID_intersection"]["new_fcp"]["white_labels"]==1


def test_synthetic_unverified_union_cannot_claim_full_real_247032(tmp_path):
    z=union(tmp_path)
    with pytest.raises(ValueError,match="exact frozen original"):
        audit(z,tmp_path/"out",test_only=False)


def test_conflicting_images_of_same_observation_quarantined(tmp_path):
    def rewrite(rows):
        rows[7]["observation_id"]=rows[6]["observation_id"]
    z=union(tmp_path,photo_replacement=rewrite)
    r=audit(z,tmp_path/"out",test_only=True)
    assert r["source_photo_and_automatic_colour_funnel"]["repeated_photo_from_same_observation"]==1
    assert r["source_photo_and_automatic_colour_funnel"]["conflicting_multi_photo_observation"]==1
    assert r["n_independent_original_colour_classified_observations"]==2
    assert r["fine_grid_within_species_siteyear_climate_scope_min2_pair_count"]==0
