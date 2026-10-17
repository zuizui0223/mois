"""Synthetic source contract tests; no original photo/date claims from fixtures."""
import pandas as pd
import pytest

from scripts.analysis import mois_final_original_photo_s3_api_seasonal_composition as m
from scripts.analysis.mois_seasonal_267k_climate_normals import labelled_calendar


def tiny():
    photo=pd.DataFrame([
      dict(photo_id="101",observation_id="201",inat_taxon_id="301",
           morph="white",old_mois_or_new_fcp="historical",observed_on=None,
           latitude=40.01,longitude=139.1),
      dict(photo_id="102",observation_id="202",inat_taxon_id="302",
           morph="blue_purple",old_mois_or_new_fcp="historical",observed_on=None,
           latitude=40.01,longitude=139.1),
      dict(photo_id="103",observation_id="203",inat_taxon_id="303",
           morph="mixed_uncertain",old_mois_or_new_fcp="historical",observed_on=None,
           latitude=40.01,longitude=139.1),
      dict(photo_id="104",observation_id="204",inat_taxon_id="304",
           morph="white",old_mois_or_new_fcp="new_fcp",observed_on="2025-07-10",
           latitude=40.01,longitude=139.1),
      dict(photo_id="105",observation_id="205",inat_taxon_id="305",
           morph="red_pink",old_mois_or_new_fcp="new_fcp",observed_on="2025-01-10",
           latitude=40.01,longitude=139.1)])
    s3=pd.DataFrame([
      dict(original_photo_id="101",original_numeric_observation_id="201",
           original_taxon_id="301",original_frozen_ROI_morph="white",
           snapshot_observed_on="2025-07-10"),
      dict(original_photo_id="102",original_numeric_observation_id="202",
           original_taxon_id="302",original_frozen_ROI_morph="blue_purple",
           snapshot_observed_on=None),
      dict(original_photo_id="103",original_numeric_observation_id="203",
           original_taxon_id="303",original_frozen_ROI_morph="mixed_uncertain",
           snapshot_observed_on=None)])
    api=pd.DataFrame([
      dict(original_photo_id="101",original_observation_id="201",
           source_taxon_id="301",provenance_status=
           "source_photo_and_observation_verified_with_current_public_date",
           retrieved_observed_on="2025-07-11"),
      dict(original_photo_id="102",original_observation_id="202",
           source_taxon_id="302",provenance_status=
           "source_photo_and_observation_verified_with_current_public_date",
           retrieved_observed_on="2025-01-10")])
    return photo,s3,api


def patch_counts(monkeypatch):
    for key,value in {
      "HISTORICAL_N":3,"MOIS_ROWS":5,"MOIS_HUES":4,"API_N":2,"EXPECTED_S3_DATED":1,
      "EXPECTED_API_FALLBACK":1,"EXPECTED_DATE_BOTH":1,
      "EXPECTED_DATE_CONFLICT":1,"EXPECTED_DATED_HUES":4,
      "EXPECTED_DATED_PHOTOS":4}.items():
        monkeypatch.setattr(m,key,value)


def test_snapshot_date_prioritized_original_hue_unchanged_and_photo_api_fallback(monkeypatch):
    patch_counts(monkeypatch)
    photo,s3,api=tiny()
    h,p,receipt=m.rejoin(photo,s3,api)
    assert receipt["history_original_dates_total_with_both_source_precedence"]==2
    assert receipt["still_unresolved_original_historical_dates"]==1
    assert receipt["current_API_vs_official_S3_observed_on_disagreements"]==1
    assert receipt["newly_source_date_recovered_historical_original_four_hue_photos"]==2
    assert receipt["real_original_photo_hue_and_date_after_all_verified_rejoins"]==4
    assert h.loc[h.photo_id=="101","qualified_observed_on"].iloc[0]=="2025-07-10"
    assert h.loc[h.photo_id=="101","qualified_date_source"].iloc[0]=="s3_20260927"
    assert h.loc[h.photo_id=="102","qualified_date_source"].iloc[0]=="iNat_current_API_20261017"
    assert h.loc[h.photo_id=="103","qualified_date_source"].iloc[0]=="unresolved_original_photo_date"
    assert p.loc[p.photo_id=="103","date"].isna().all()
    assert int(p.classified.sum())==4


def test_wrong_historical_original_taxon_or_photo_hue_rejected(monkeypatch):
    patch_counts(monkeypatch)
    p,s,a=tiny()
    s.loc[0,"original_taxon_id"]="wrong"
    with pytest.raises(ValueError,match="taxon/observation"):
        m.rejoin(p,s,a)
    p,s,a=tiny()
    s.loc[0,"original_frozen_ROI_morph"]="red_pink"
    with pytest.raises(ValueError,match="hue predictions"):
        m.rejoin(p,s,a)


def test_api_date_rejected_unless_photo_and_original_numeric_observation_matched(monkeypatch):
    patch_counts(monkeypatch)
    p,s,a=tiny()
    a.loc[1,"original_observation_id"]="wrong"
    with pytest.raises(ValueError,match="identity"):
        m.rejoin(p,s,a)
    p,s,a=tiny()
    a.loc[1,"provenance_status"]="not_found_in_current_public_api"
    with pytest.raises(ValueError,match="overlap or conflict"):
        m.rejoin(p,s,a)


def test_full_photo_seasonal_taxon_accounting_is_additive_across_missing_taxa():
    rows=[]
    for tax,month,year,white,n in [
        ("a",1,2025,False,3),("b",1,2025,True,2),
        ("a",7,2025,True,5),("b",7,2025,True,2),
        ("c",7,2025,False,4),("d",1,2025,True,2)]:
        for i in range(n):
            rows.append(dict(photo_id=f"{tax}-{month}-{i}",observation_id=f"{tax}-{month}-{i}",
                inat_taxon_id=tax,morph="white" if white else "blue_purple",
                observed_on=f"{year}-{month:02d}-06",
                latitude=-35.0,longitude=150.0))
    d=labelled_calendar(pd.DataFrame(rows))
    result=m.seasonal_taxon_decomposition(d,"south_extratropical")
    assert result["common_taxa_both_seasons"]==2
    assert result["summer_photographed_named_taxa"]==3
    assert result["winter_photographed_named_taxa"]==3
    assert result["is_exact_additive_source_photo_algebra_NOT_causal_biological_mechanism"]


def test_community_locality_uses_distinct_observation_and_geographic_width():
    rows=[]
    for month,latitude,morph in [(1,40.01,"red_pink"),(7,40.30,"white")]:
        for i in range(5):
            rows.append(dict(photo_id=f"{month}-{i}",observation_id=f"{month}-{i}",
                inat_taxon_id=str(i%2),morph=morph,observed_on=f"2025-{month:02d}-07",
                latitude=latitude,longitude=138.1))
    d=labelled_calendar(pd.DataFrame(rows))
    x=m.local_support(d)
    def groups(w,n):
        return next(v["same_cell_same_meteorological_season_year_comparisons"]
                    for v in x if v["geographic_cell_degrees"]==w and
                    v["minimum_distinct_original_observation_ID_per_summer_and_winter"]==n
                    and v["hemisphere"]=="north_extratropical"
                    and v["minimum_named_original_taxa_per_season"]==2)
    assert groups(.25,5)==0
    assert groups(.5,5)==1
