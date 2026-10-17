"""Prove source-photo seasonal four-term identity without biological imputation."""
import pandas as pd
import pytest

from scripts.analysis.mois_full146052_original_photo_season_taxon_composition import (
    decompose_photographic_taxa
)


def sample_photographs():
    rows=[]
    for taxon,season,labels in (
        ("A","summer",[1,1]),("A","winter",[0,0]),
        ("B","summer",[0,0]),("B","winter",[1,1]),
        ("C","summer",[1]),("D","winter",[0])):
        for iswhite in labels:
            rows.append(dict(photo_id=str(len(rows)+1),observation_id=str(111+len(rows)),
                inat_taxon_id=taxon,season=season,white=iswhite))
    return pd.DataFrame(rows)


def test_all_photographic_source_taxa_including_exclusive_taxa_sum_to_raw():
    r=decompose_photographic_taxa(sample_photographs(),draws=27,seed=11)
    assert r["source_named_taxa_common_both_summer_winter"]==2
    assert r["original_named_taxa_union_both_seasons"]==4
    assert r["raw_photo_white_summer_minus_winter_pp"]==pytest.approx(20)
    assert r["within_common_taxon_photographic_white_rate_component_pp"]==pytest.approx(0)
    assert r["changed_common_taxon_relative_photo_frequency_component_pp"]==pytest.approx(0)
    assert r["summer_photo_exclusive_taxon_support_component_pp"]==pytest.approx(10)
    assert r["winter_photo_exclusive_taxon_support_component_pp"]==pytest.approx(10)
    assert r["exact_full_census_mixture_4terms_sum_to_raw"]


def test_taxon_not_photographed_in_winter_is_not_assigned_a_winter_white_colour():
    d=sample_photographs()
    one=d[d.inat_taxon_id.ne("D")]
    # winter photo-only correction zero if no exclusive winter photos;
    # this does not assign winter a nonexistent C photo.
    result=decompose_photographic_taxa(one,draws=12)
    assert result["winter_photo_exclusive_taxon_support_component_pp"]==pytest.approx(0)
    assert result["source_named_taxa_common_both_summer_winter"]==2


def test_only_one_shared_taxon_does_not_support_broad_comparative_decomposition():
    d=sample_photographs()
    q=d[d.inat_taxon_id.ne("B")]
    with pytest.raises(ValueError,match="common photographed taxa"):
        decompose_photographic_taxa(q,draws=5)


def test_sample_all_photo_ids_unique_required():
    d=sample_photographs()
    d.loc[1,"photo_id"]=d.loc[0,"photo_id"]
    with pytest.raises(ValueError,match="photo"):
        decompose_photographic_taxa(d,draws=5)
