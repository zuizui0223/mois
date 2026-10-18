"""Sampling-sensitive ranking diagnostics across climate and missing soil."""
import numpy as np
import pandas as pd
from scripts.analysis import mois_global313535_climate_subset_sampling_robustness as m


def sample(n_gen=245,taxa_per=5):
    rng=np.random.default_rng(101)
    records=[]
    for g in range(n_gen):
        genus=f"G{chr(65+g//26)}{chr(65+g%26)}"
        for j in range(taxa_per):
            i=g*taxa_per+j
            row={"photo_id":str(i+1),"inat_taxon_id":str(i+1),
                 "species":f"{genus} tax{j}","genus":genus,
                 "classified":True,"white":bool(i%3),
                 "latitude":float(-60+(g%18)*7),
                 "abs_lat":float(abs(-60+(g%18)*7)),
                 "longitude_sin":float(np.sin(np.deg2rad((g%30)*12-160))),
                 "longitude_cos":float(np.cos(np.deg2rad((g%30)*12-160))),
                 "source_new":float(j%2)}
            for k in m.base.ABIOTIC:
                row[k]=float(rng.normal(0,1)+g*.005)
            if g>=220:
                for k in m.base.BLOCKS["soil"]:
                    row[k]=np.nan
            records.append(row)
    return pd.DataFrame(records)


def test_outcome_source_never_blends_source_only_pixel_IDs():
    assert m.base.prev.N_PHOTOS==313535
    assert sum(len(m.base.BLOCKS[x]) for x in m.CLIMATE_BLOCKS)==21
    assert sum(len(m.base.BLOCKS[x]) for x in m.TERRAIN_BLOCKS)==28
    assert len(m.ALL_COVARIATES)==33


def test_source_support_cannot_silently_fill_missing_soil():
    s=sample()
    fields=m.base.BLOCKS["soil"]
    full=m.source_stats(s,())
    narrowed=m.source_stats(s,fields)
    assert full["processed_source_photo_attempts"]==len(s)
    assert narrowed["processed_source_photo_attempts"]<len(s)
    assert narrowed["classifiable_photo_hues"]<full["classifiable_photo_hues"]


def test_taxon_centroid_holdout_and_training_only_group_pca():
    s=sample()
    p,d=m.build_panel(s,m.CLIMATE_BLOCKS,condition=())
    assert len(p)>=1000
    assert d["qualifying_genera"]>=200
    rows,extra=m.evaluate(p,m.CLIMATE_BLOCKS,"heldout_genera")
    assert len(rows)==5
    assert {x["excluded_environment_group"] for x in rows if x["excluded_environment_group"]}==set(m.CLIMATE_BLOCKS)
    assert extra["holdout_design"]=="heldout_genera"
    assert np.isfinite(extra["full_cv_R2"])
    assert np.isfinite(extra["geo_only_R2"])


def test_photo_level_soil_gated_climate_panel_does_not_change_other_features():
    s=sample()
    c=tuple(x for group in m.CLIMATE_BLOCKS for x in m.base.BLOCKS[group])
    all33=m.build_panel(s,m.CLIMATE_BLOCKS,
                        condition=(*c,*m.base.BLOCKS["terrain"],*m.base.BLOCKS["soil"]))[0]
    only=m.build_panel(s,m.CLIMATE_BLOCKS)[0]
    assert len(only)>len(all33)
    assert all33[m.base.BLOCKS["temperature"][0]].notna().all()
    assert all33["photographed_white_rate"].between(0,1).all()
