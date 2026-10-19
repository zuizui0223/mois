"""Synthetic algebra, date semantics and observation-sensitive checks for MOIS 350399."""
import pandas as pd
from scripts.analysis.mois_real350399_taxon_season_decomposition_20261019 import (
    photographic_calendar, observation_equal_white, estimate,
)


def make(pid, taxon, lat, month, hue, obs=None, arm="core"):
    return dict(photo_id=str(pid), inat_taxon_id=str(taxon),
                observation_id=obs or f"o{pid}",
                species="Campanula punctata", latitude=lat, longitude=140.1,
                date=pd.Timestamp(2025,month,13),
                morph=hue, classified=hue!="mixed_uncertain",
                white=hue=="white", arm=arm)


def fixture():
    # Two shared summer/winter taxa per hemisphere; distinct exclusive taxa.
    x=[
        make(1,10,40,6,"white",obs="n1"),
        make(2,10,40,7,"white",obs="n1"),
        make(3,10,40,1,"blue_purple"),
        make(4,20,40,6,"red_pink"),
        make(5,20,40,12,"white"),
        make(6,30,40,7,"white",arm="tail16384"),
        make(7,10,-40,12,"white"),
        make(8,10,-40,1,"red_pink"),
        make(9,10,-40,7,"yellow_orange"),
        make(10,20,-40,1,"white"),
        make(11,20,-40,6,"red_pink"),
        make(12,30,-40,6,"white"),
        make(13,30,40,11,"white"),  # non-summer/winter: excluded
        make(14,20,40,7,"mixed_uncertain"),  # unknown: excluded
    ]
    return pd.DataFrame(x)


def test_calendar_hemisphere_and_djf_local_year():
    q=photographic_calendar(fixture())
    assert len(q)==12
    assert q.loc[q.photo_id.eq("7"),"season"].iloc[0]=="summer"
    assert q.loc[q.photo_id.eq("7"),"season_year"].iloc[0]==2026
    assert q.loc[q.photo_id.eq("5"),"season_year"].iloc[0]==2026
    assert q.loc[q.photo_id.eq("3"),"season"].iloc[0]=="winter"


def test_fourway_exact_decomposition_and_observation_collapse():
    d=fixture()
    r=estimate(d,draws=19)
    assert r["classified_extratropical_dated_summer_winter_photos"]==12
    for zone,x in r["hemispheres"].items():
        reconstructed=sum(x[k] for k in (
            "within_common_taxon_photographic_white_rate_component_pp",
            "changed_common_taxon_relative_photo_frequency_component_pp",
            "summer_photo_exclusive_taxon_support_component_pp",
            "winter_photo_exclusive_taxon_support_component_pp"))
        assert abs(reconstructed-x["raw_photo_white_summer_minus_winter_pp"])<1e-4, zone
        assert x["source_named_taxa_common_both_summer_winter"]==2
        assert x["exact_full_census_mixture_4terms_sum_to_raw"]
    s=r["observation_equal_sensitivity"]
    assert s["photos_with_observation_id"]==12
    assert s["observations"]==11
    assert s["max_images_per_observation"]==2
    assert r["source_arm_seasonal_photo_counts"]["tail16384"]["north_extratropical"]["summer"]["classifiable_dated_photo_ROIs"]==1


def test_missing_observation_ids_dont_create_fake_replication():
    d=fixture()
    d.loc[d.photo_id.eq("3"),"observation_id"]=""
    q=photographic_calendar(d)
    r=observation_equal_white(q)
    assert r["photos_with_observation_id"]==11
    assert r["observations"]==10
