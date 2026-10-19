"""Named taxon and matched genus-geography original-photo abiotic slope contracts.

The source is synthetic: no ecological effect or global colour frequency
is inferred from these examples.
"""
import numpy as np
import pandas as pd
from scripts.analysis.mois_real350399_annual_abiotic_geogenus_white_20261019 import (
    PRED, taxon_geocells, local_slopes, within_predictor_correlations,
)


def synthetic():
    # Twelve independent named genera, each with 3 species in EACH photo
    # recruitment arm in one cell; all abiotic variables are correlated as
    # required to demonstrate this is univariate NOT causal inference.
    rows=[]
    letters="ABCDEFGHIJKL"
    for g,letter in enumerate(letters):
        for arm in ("core","new_species_first"):
            for tax in range(3):
                rain=10*(tax+1)+g
                rows.append({
                    "photo_id":f"{g}-{arm}-{tax}",
                    "inat_taxon_id":f"{g}-{arm}-{tax}",
                    "genus":f"{letter}aa",
                    "cohort":arm,
                    "latitude":35.1, "longitude":132.1,
                    "white_float":float(tax/2 if arm=="core" else 1-tax/2),
                    "log1p_annual_rainfall":np.log1p(rain),
                    "mean_atmospheric_vapour":float(rain/10),
                    "annual_solar":float(200+2*rain),
                    "warmest_month_temperature":float(10+rain/10),
                    "elevation":float(5*rain),
                })
    return pd.DataFrame(rows)


def find(results,cohort,pred):
    return next(r for r in results if r["cohort"]==cohort
                and r["predictor"]==pred)


def test_source_cohort_stratification_does_not_create_spurious_common_slope():
    t=taxon_geocells(synthetic(),10)
    assert len(t)==72
    assert t.n_photos.eq(1).all()
    slopes,s=local_slopes(t,degrees=10,min_taxa=2,draws=99)
    old=find(slopes,"core","log1p_annual_rainfall")
    new=find(slopes,"new_species_first","log1p_annual_rainfall")
    both=find(slopes,"combined_arm_conditioned","log1p_annual_rainfall")
    assert old["white_pp_per_one_full_taxon_cell_SD"]>0
    assert new["white_pp_per_one_full_taxon_cell_SD"]<0
    assert abs(both["white_pp_per_one_full_taxon_cell_SD"])<1e-8
    assert old["named_genera_with_within_cell_variation"]==12
    assert len(old["exploratory_nominal_genus_bootstrap_95"])==2
    assert old["selected_named_taxon_cells"]==36
    assert both["selected_named_taxon_cells"]==72
    assert len(within_predictor_correlations(s))==10


def test_single_taxon_genus_cells_fail_matched_support_not_zero_effect():
    x=synthetic()
    x=x[x.inat_taxon_id.str.endswith("-0")]
    t=taxon_geocells(x,10)
    result,_=local_slopes(t,degrees=10,min_taxa=2,draws=9)
    assert all(r["white_pp_per_one_full_taxon_cell_SD"] is None for r in result)
    assert all(r["selected_named_taxon_cells"]==0 for r in result)
    assert all(r["exploratory_nominal_genus_bootstrap_95"] is None for r in result)


def test_photo_repeats_within_named_taxon_do_not_change_taxon_cell_weighting():
    x=synthetic()
    dup=x.iloc[[0]].copy()
    dup.loc[:,"photo_id"]=["additional-independent-photo-same-taxon"]
    dup.loc[:,"white_float"]=[1.0]
    t=taxon_geocells(pd.concat([x,dup],ignore_index=True),10)
    assert len(t)==72
    assert int(t.n_photos.sum())==73
    actual=t[t.inat_taxon_id.eq(x.iloc[0].inat_taxon_id)].iloc[0]
    assert int(actual.n_photos)==2
    assert actual.p_white==.5


def test_spatial_grid_does_not_pool_same_genus_distant_islands():
    x=synthetic()
    special=x.iloc[[0]].copy()
    special.loc[:,"latitude"]=[-30.]
    special.loc[:,"inat_taxon_id"]=["moved-species"]
    special.loc[:,"photo_id"]=["moved-real-photo"]
    t=taxon_geocells(pd.concat([x,special],ignore_index=True),10)
    assert len(t)==73
    slopes,_=local_slopes(t,degrees=10,min_taxa=2,draws=9)
    result=find(slopes,"core","log1p_annual_rainfall")
    assert result["selected_named_taxon_cells"]==36


def test_all_five_predictors_predeclared():
    assert set(PRED)=={
        "log1p_annual_rainfall","mean_atmospheric_vapour",
        "annual_solar","warmest_month_temperature","elevation"
    }
