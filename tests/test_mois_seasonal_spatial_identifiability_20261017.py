"""Synthetic source-independent mathematical guards for MOIS seasonal overlap."""
import pandas as pd
import pytest
from scripts.analysis.mois_seasonal_spatial_identifiability_20261017 import (
    full_photo_turnover, eligible_site_photo_support,
)


def example():
    d=[]
    for tax,season,year,lat,white,n in (
        ("common","summer",2024,40.1,True,3),
        ("common","winter",2025,40.1,False,3),
        ("summer_only","summer",2024,40.1,True,2),
        ("winter_only","winter",2025,40.1,False,2),
    ):
        for i in range(n):
            d.append({"photo_id":f"{tax}-{season}-{i}",
                "observation_id":f"{tax}-{season}-{i}","inat_taxon_id":tax,
                "latitude":lat,"longitude":140.1,"white":white,
                "season":season,"season_year":year})
    return pd.DataFrame(d)


def test_full_photo_summer_winter_includes_taxa_absent_other_season():
    summary=full_photo_turnover(example())
    assert summary["all_photo_summer_minus_winter_white_pp"]==100
    assert summary["n_taxa_common_both_seasons"]==1
    assert summary["n_taxa_observed_only_summer"]==1
    assert summary["n_taxa_observed_only_winter"]==1
    assert summary["photographic_only_season_summer_fraction"]==0.4
    assert summary["is_exact_full_source_photo_algebra_not_phenological_causation"]


def test_same_taxon_and_site_not_same_season_year():
    d=example()
    same_site=eligible_site_photo_support(d,.25,False,2)
    assert same_site=={"matched_groups":1,"matched_taxa":1,
                       "matched_photos":6,"matched_original_observations":6}
    same_year=eligible_site_photo_support(d,.25,True,2)
    assert same_year["matched_groups"]==0


def test_coarsening_site_cells_increases_apparent_season_overlap():
    d=example()
    d.loc[d.season=="summer","latitude"]=40.01
    d.loc[d.season=="winter","latitude"]=40.49
    assert eligible_site_photo_support(d,.25,False,2)["matched_groups"]==0
    assert eligible_site_photo_support(d,1.0,False,2)["matched_groups"]==1


def test_original_observation_count_is_not_photo_count():
    d=example()
    mask=(d.inat_taxon_id=="common") & d.season.eq("winter")
    # Three winter photos become a single observation; not enough.
    d.loc[mask,"observation_id"]="same_observer_record"
    assert eligible_site_photo_support(d,.25,False,2)["matched_groups"]==0
    assert eligible_site_photo_support(d,.25,False,1)["matched_groups"]==1
