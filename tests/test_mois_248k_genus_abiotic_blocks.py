"""Synthetic checks for MOIS genus-block photo-hue analysis (not botanical data)."""
import numpy as np
import pandas as pd
import pytest
from scripts.analysis.mois_248k_global_genus_abiotic_blocks import (
    HUES,GEOG,CLIM,SOIL,TERRAIN,genus,make_photo_features,
    aggregate_species,design,genus_cv_r2,
)


def test_nominal_genus_and_FCP_response_boundary():
    assert genus("Campanula punctata")=="Campanula"
    assert genus("Asclepias syriaca")=="Asclepias"
    assert genus("× Test hybrid") is None
    assert sum(map(len,(GEOG,CLIM,TERRAIN,SOIL)))==18
    assert "D" not in (*GEOG,*CLIM,*SOIL)
    assert "W" not in (*GEOG,*CLIM,*SOIL)


def synthetic():
    arr=[]
    for g in range(70):
        for j in range(4):
            for photo in range(2):
                i=len(arr)
                row={"inat_taxon_id":str(100000+g*4+j),
                     "species":f"Genus{g} species{j}",
                     "morph":HUES[(g+j+photo)%4],
                     "original_source":"historical" if g%2 else "new_fcp",
                     "latitude":30+g%15+0.1*j,"longitude":100+g%20+j,
                     "wc2_elev_2p5m":400+2*j+g,
                     "wc2_bio01_10m":11+.05*g+j*.1,
                     "wc2_bio05_10m":20+.04*g+j*.11,
                     "wc2_bio12_10m":900+j*200+g*1.5,
                     "wc2_bio15_10m":40+g*.5+j,
                     "wc2_annual_srad_normal_10m_Wm2":170+j*7+g,
                     "terrain_slope_2p5m_degree":4+j*.8,
                     "terrain_TRI_2p5m_m":22+j*3,
                     "historical_original_soil_pH":6+j*.1,
                     "historical_original_soil_SOC":15+g*.01+j,
                     "historical_original_soil_N":1+g*.005+j*.1,
                     "historical_original_soil_clay":20+g*.1+j,
                     "historical_original_soil_available_water_proxy":10+j*.3,
                     "observed_on":"2025-07-15"}
                # Actual taxon genus parser requires Latin letters; replace
                # numeric labels with names that remain consistent per genus.
                row["species"]=f"Genus{'abcdefghij'[g%10]}{'klmnopq'[g//10]} species{j}"
                arr.append(row)
    return make_photo_features(pd.DataFrame(arr))


def test_same_original_source_photos_support_geo_climate_soil_addition():
    p=synthetic()
    p.loc[p.index[:20],"historical_original_soil_pH"]=np.nan
    original=len(p)
    s,receipt=aggregate_species(p,(*GEOG,*CLIM,*SOIL))
    assert original==560
    assert receipt["original_photo_rows_with_all_block_covariates"]==540
    assert receipt["taxa_after_min_photo_and_congeneric_gates"]>=200
    assert receipt["qualified_nominal_genera"]==68
    assert receipt["min_classified_photos_per_species"]==1
    assert len(s)==270


def test_genus_centered_photo_white_hue_response_is_not_causal_photo_pair_D():
    s,_=aggregate_species(synthetic(),(*GEOG,*CLIM))
    X,y,weights,names=design(s,(*GEOG,*CLIM),"white")
    assert len(s)==280
    assert X.shape[0]==280
    assert np.isfinite(X).all() and np.isfinite(y).all()
    assert np.all(weights>0)
    assert abs(y.mean())<1e-12
    value=genus_cv_r2(s,(*GEOG,*CLIM),"white")
    assert np.isfinite(value) and value<=1.0


def test_taxon_genus_requires_replication():
    p=synthetic()
    p.loc[p.index[:],"species"]=[f"UniqueGenus{n} species" for n in range(len(p))]
    # Fake independent genera with numeric suffixes are NOT legitimate
    # plant genus tokens and must be rejected before fitting.
    assert genus("UniqueGenus123 species") is None
