"""The six anchors and original 4-hue 999-genus maxT family cannot drift."""
import numpy as np
import pandas as pd
from scripts.analysis.fcp_abiotic_hue_direction_macro_20261013 import (
    ANCHOR_CLIMATE, CLASSES, GEO
)
from scripts.analysis.mois_original_6anchor_water_maxT_expanded import (
    CLIMATE, one_scope
)

def test_original_six_anchor_fixed_moisture_and_source_solar():
    assert ANCHOR_CLIMATE == ("wc_bio5","wc_bio12","wc_bio14","wc_bio15",
                              "wc_srad_annual_mean","wc_vapr_annual_mean")
    assert len(CLASSES)*len(ANCHOR_CLIMATE)==24
    assert len(GEO)==3 and set(CLIMATE)==set(ANCHOR_CLIMATE)-{"wc_vapr_annual_mean"}

def test_original_null_one_fixed_six_anchor_four_hue_synthetic():
    # 45 replicated test genera, 4 artificial species each, no real plant data.
    fake=[]
    genera=[f"Genus{'abcde'[i%5]}{'fghijklmn'[i//5]}" for i in range(45)]
    rng=np.random.default_rng(62)
    for g,name in enumerate(genera):
        for species in range(4):
            row={"inat_taxon_id":str(10000+g*4+species),
                 "species":f"{name} testing",
                 "genus":name,"photo_id":str(10000+g*4+species),
                 "latitude":40+g*.1+species*.05,
                 "longitude":-80+g*.2+species*.1,
                 "morph":CLASSES[(g+species)%4],
                 "abs_lat":40+g*.1+species*.05,
                 "lon_sin":np.sin(np.deg2rad(-80+g*.2+species*.1)),
                 "lon_cos":np.cos(np.deg2rad(-80+g*.2+species*.1))}
            for j,k in enumerate(ANCHOR_CLIMATE):
                row[k]=10+g*.3+species*.7+rng.normal(0,0.1)+(j+1)*g*0.07
            fake.append(row)
    result=one_scope(pd.DataFrame(fake),"world_genus_geo",999,44)
    assert result["status"]=="OK"
    assert result["n_bootstrap"]==999
    assert len(result["estimates"])==24
    assert len({v["anchor"] for v in result["estimates"]})==6
    assert result["n_named_genus"]==45
    for line in result["estimates"]:
        assert len(line["scope_anchor_hue_maxT_ci95_pp"])==2
        assert np.all(np.isfinite(line["scope_anchor_hue_maxT_ci95_pp"]))
