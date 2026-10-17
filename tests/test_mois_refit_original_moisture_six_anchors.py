"""Original MOIS macro 100543-method bootstrap and six frozen water/climate anchors."""
import numpy as np
import pandas as pd
import pytest

from scripts.analysis import fcp_abiotic_hue_direction_macro_20261013 as original
from scripts.analysis.mois_refit_original_100k_water_six_anchors_248k import (
    ANCHORS, GEO, HUES, BOOTSTRAPS, EXPECTED_CLASSIFIED,
    EXPECTED_PHOTOS, extract_water
)


def test_exact_original_six_anchor_family_and_genus_method_unchanged():
    assert ANCHORS == original.ANCHOR_CLIMATE
    assert GEO == original.GEO
    assert HUES == original.CLASSES
    assert BOOTSTRAPS == 999
    assert len(ANCHORS)==6
    assert "wc_bio12" in ANCHORS
    assert "wc_vapr_annual_mean" in ANCHORS
    assert EXPECTED_PHOTOS==248312
    assert EXPECTED_CLASSIFIED==120807


def test_source_restricted_photo_model_can_hold_insufficient_species():
    f=pd.DataFrame({
      "inat_taxon_id":[1,2],
      "genus":["Testus","Testus"],
      "photo_id":["1","2"],
      "latitude":[10.,20.],"longitude":[12.,13.],
      "abs_lat":[10.,20.],"lon_sin":[0.1,0.2],"lon_cos":[.9,.8],
      "morph":["white","red_pink"],
      **{v:[1.,2.] for v in ANCHORS}
    })
    frame=original.species_frame(f,ANCHORS)
    groups,meta=original.cluster_suff(
        frame,ANCHORS,
        f[list((*GEO,*ANCHORS))].mean(),
        f[list((*GEO,*ANCHORS))].std().replace(0.,1.),
        geos=True)
    assert groups is None
    assert meta["n_genera"]==0


def test_water_summary_refuses_fake_positive_when_model_unestimable():
    fake={
      "original_method_FCP_macro_six_anchor_genus_equal":{"status":"OK","estimates":[
        {"anchor":"wc_bio12","hue":"white","point_pp_per_original_photo_env_sd":2.0,
         "scope_anchor_hue_maxT_ci95_pp":[-1,5],"maxT_excludes_zero":False},
        {"anchor":"wc_vapr_annual_mean","hue":"white","point_pp_per_original_photo_env_sd":1.2,
         "scope_anchor_hue_maxT_ci95_pp":[-1,3],"maxT_excludes_zero":False},
      ]},
      "original_method_FCP_macro_six_anchor_species_equal":{"status":"HOLD","estimates":[]},
      "same_genus_coarse_cell_original_method_sensitivity":{"status":"HOLD"}
    }
    d=extract_water(fake)
    assert not d["genus_equal__wc_bio12"]["maxT_excludes_zero"]
    assert d["coarse_cell__wc_bio12"]["not_reported_positive"]


def test_no_source_order_signal_in_sampled_environment_coefficients():
    assert original.GEO==("abs_lat","lon_sin","lon_cos")
    assert len(original.CLASSES)==4
    assert original.BOOT_SEED==2026101348
