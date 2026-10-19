#!/usr/bin/env python3
"""Full 350399 real-source-photo macroecological white × annual abiotic local-gradient test.

Responses are automatic photographic four-hue ROI labels, NOT verified petals,
metabolic costs or within-individual season effects. The test operates on
*named taxa within photographed genus × coarse geographic cells* and removes
genus-cell source mean before fitting univariate climatic differences. Source
cohort is an explicit part of every stratum, because taxon-first and deep
iNaturalist photography use different unknown sampling probabilities.

This is exploratory co-occurrence, not isolation of causal climate predictors:
water vapour, rain, heat, altitude and solar are correlated, and a 5°/10°/20°
cell is not a local experimental population.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis.mois_actual_photo_cohort_season_audit_20261019 import (
    assemble, require, sha256,
)

FULL_PHOTOS=350399
ENVFILE="mois350399_all19_annual_vapr_solar_terrain_source_photoids.parquet"
ENVRECEIPT="mois350399_immutable_annual_climate_source_receipt.json"
OUT_JSON="mois350399_real_abiotic_within_geogenus_regression.json"
OUT_CSV="mois350399_real_abiotic_within_geogenus_slopes.csv"
OUT_CORR="mois350399_abiotic_taxon_cell_within_environment_correlation.csv"
GRID=(5,10,20)
MIN_GROUP_TAXA=(2,3)
PRED={
    "log1p_annual_rainfall": "wc2_bio12_10m",
    "mean_atmospheric_vapour": "wc2_annual_vapr_10m",
    "annual_solar": "annual_solar_Wm2",
    "warmest_month_temperature": "wc2_bio05_10m",
    "elevation": "wc2_elev_2p5m"
}
COHORTS=("core","new_species_first","combined_arm_conditioned")
GENUS_RE=re.compile(r"^[A-Z][A-Za-z-]+$")


def source_data(photos, env, receipt):
    info=json.loads(Path(receipt).read_text())
    require(info["schema"]=="mois_full350399_source_photograph_original_worldclim_all19_annual_v1"
            and info["source_photo_attempts"]==FULL_PHOTOS
            and info["previous_unchanged_photo_climate"]==317631
            and info["new_only_original_image_photo_sites_sampled"]==32768
            and sha256(env)==info["sha256"][ENVFILE],
            "Not signed actual full350399 original WorldClim photo source")
    abiotic=pd.read_parquet(env)
    require(len(abiotic)==FULL_PHOTOS and abiotic.photo_id.is_unique
            and photos.photo_id.tolist()==abiotic.photo_id.astype(str).tolist(),
            "Actual same original photo ID and climate site is not aligned")
    x=photos.loc[photos.classified].copy()
    e=abiotic.loc[photos.classified.values]
    require(len(e)==len(x),"ROI original photo source and climate row mask mismatch")
    for key,value in PRED.items():
        require(value in e,"Missing original WorldClim measurement: "+value)
        x[key]=pd.to_numeric(e[value],errors="coerce").to_numpy(float)
    x["genus"]=x.species.str.strip().str.split().str[0]
    x=x[x.genus.map(lambda s:isinstance(s,str) and bool(GENUS_RE.fullmatch(s)))
          &x.inat_taxon_id.ne("")&
          x.latitude.between(-90,90)&x.longitude.between(-180,180)].copy()
    for key in PRED:
        x=x[np.isfinite(x[key])].copy()
    require(len(x)>10000, "Too few genuine 350399 source photographic hues with complete climate")
    x["log1p_annual_rainfall"]=np.log1p(x["log1p_annual_rainfall"].clip(lower=0))
    x["cohort"]=np.where(x.arm.eq("core"),"core","new_species_first")
    x["white_float"]=x.white.astype(float)
    return x, info


def taxon_geocells(x, degrees):
    """Classifiable original photographic taxa are the units, not photo count."""
    require(degrees in GRID, "Predeclared geography resolution only")
    v=x.copy()
    v["latcell"]=np.floor((v.latitude+90)/degrees).astype(int)
    v["loncell"]=np.floor((v.longitude+180)/degrees).astype(int)
    keys=["cohort","genus","latcell","loncell","inat_taxon_id"]
    aggs={"p_white":("white_float","mean"),"n_photos":("photo_id","size")}
    aggs.update({key:(key,"mean") for key in PRED})
    t=v.groupby(keys,observed=True).agg(**aggs).reset_index()
    require(len(t)>0 and t[keys].notna().all().all(),
            "Photographic taxon group missing identity")
    require(int(t.n_photos.sum())==len(v),
            "Classified photo IDs not fully partitioned into exact taxon × geography")
    return t


def local_slopes(t, *, degrees, min_taxa, seed=20261019, draws=999):
    """Predeclared univariate taxon- and source-arm-conditioned local slopes.

    β is white-label percentage points per one SD of a source-taxon-cell
    environmental predictor, using within-stratum demeaned responses.
    A nominal genus is the bootstrap cluster (not a phylogeny). No compound
    climatic causal effects or p values are inferred.
    """
    group=["cohort","genus","latcell","loncell"]
    require(degrees in GRID and min_taxa in MIN_GROUP_TAXA, "Unsupported spatial/overlap grid")
    frame=t.copy()
    n=frame.groupby(group,observed=True).inat_taxon_id.transform("size")
    frame=frame.loc[n.ge(min_taxa)].copy()
    raw_taxon_cells={k:len(t[t.cohort.eq(k)]) for k in ("core","new_species_first")}
    retained={k:len(frame[frame.cohort.eq(k)]) for k in ("core","new_species_first")}
    group_n=int(frame.groupby(group,observed=True).ngroups)
    # Predictor SD from ALL available photographed taxa before selecting
    # matched within-cell support. Never set sample SD by WHITE outcomes.
    scales={k:float(t[k].std(ddof=0)) for k in PRED}
    require(all(np.isfinite(s) and s>0 for s in scales.values()),
            "Original environmental predictor has no taxon geographic variance")
    for col in ("p_white",*PRED):
        frame["within_"+col]=frame[col]-frame.groupby(
            group,observed=True)[col].transform("mean")
    out=[]
    for cohort in COHORTS:
        q=(frame if cohort=="combined_arm_conditioned"
           else frame.loc[frame.cohort.eq(cohort)])
        for col in PRED:
            dx=q["within_"+col].to_numpy(float)/scales[col]
            dy=q["within_p_white"].to_numpy(float)
            sd=float(np.dot(dx,dx))
            if not len(q) or sd<1e-14:
                slope=None; bounds=None; count=0
            else:
                slope=float(100*np.dot(dx,dy)/sd)
                terms=pd.DataFrame({
                    "genus":q.genus.values,
                    "xy":dx*dy,
                    "xx":dx*dx
                }).groupby("genus",observed=True)[["xy","xx"]].sum()
                terms=terms.loc[terms.xx.gt(1e-14)]
                count=len(terms)
                # Reassemble all within-genus strata for every resampled
                # nominal genus. Null if fewer than 10 independent named
                # genera are represented with real predictor variation.
                if count>=10:
                    z=terms.to_numpy(float)
                    rng=np.random.default_rng(seed+degrees*100+min_taxa+
                                              list(PRED).index(col)*13+
                                              COHORTS.index(cohort)*107)
                    ix=rng.integers(0,len(z),(draws,len(z)))
                    sums=z[ix].sum(axis=1)
                    beta=100*sums[:,0]/sums[:,1]
                    bounds=[float(v) for v in np.quantile(beta,[.025,.975])]
                else: bounds=None
            out.append(dict(
                cohort=cohort, grid_degrees=degrees, min_named_taxa_per_genus_cell_arm=min_taxa,
                predictor=col,
                white_pp_per_one_full_taxon_cell_SD=slope,
                exploratory_nominal_genus_bootstrap_95=bounds,
                predictor_source_taxon_cell_SD=scales[col],
                named_genera_with_within_cell_variation=count,
                selected_named_taxon_cells=int(len(q)),
                selected_genus_grid_arm_groups=int(q.groupby(group,observed=True).ngroups),
                all_available_source_taxon_cells_by_arm=raw_taxon_cells,
                retained_source_taxon_cells_by_arm=retained,
                source_arm_fixed_within_cell=True,
                photo_classifiability_selected=True,
                causal_or_phylogeny_adjusted=False,
                group_n_in_both_arms=group_n))
    return out,frame


def within_predictor_correlations(frame):
    keys=list(PRED)
    if frame.empty:return []
    y=frame[["within_"+k for k in keys]]
    corr=y.corr()
    return [{"predictor_1":a,"predictor_2":b,
             "within_genus_geocell_cor":float(corr.iloc[i,j])}
            for i,a in enumerate(keys) for j,b in enumerate(keys)
            if i<j and np.isfinite(corr.iloc[i,j])]


def run(core,broad,wave2,cached8192,tail16384,environment,env_receipt,output):
    output=Path(output)
    require(not output.exists() or not any(output.iterdir()),
            "No previous source regression archive may be overwritten")
    d,arms=assemble(core,broad,wave2,cached8192,tail16384,strict=True)
    x,src=source_data(d,environment,env_receipt)
    data=[]
    correlations=[]
    reports={}
    for deg in GRID:
        t=taxon_geocells(x,deg)
        reports[str(deg)]={"actual_photo_classified_climate_complete":len(x),
                           "source_taxon_cell_rows":len(t),
                           "photo_ROIs_in_taxon_cells":int(t.n_photos.sum())}
        for gate in MIN_GROUP_TAXA:
            slopes,eligible=local_slopes(t,degrees=deg,min_taxa=gate)
            data.extend(slopes)
            correlations.extend({"degrees":deg,"min_taxa":gate,**r}
                                for r in within_predictor_correlations(eligible))
    output.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(data).to_csv(output/OUT_CSV,index=False)
    pd.DataFrame(correlations).to_csv(output/OUT_CORR,index=False)
    result={
        "schema":"mois_real350399_annual_climate_taxon_geogenus_white_within_locality_v1",
        "source_photo_attempts":len(d),
        "source_four_hue_classified":int(d.classified.sum()),
        "source_classified_and_all_five_annual_predictors_plus_geocoords":len(x),
        "source_annual_climate_SHA256":sha256(environment),
        "photographic_pixel_immutable_source_SHA256":{
            n:sha256(p) for n,p in zip(
                ("core","broad","wave2","cached8192","tail16384"),
                (core,broad,wave2,cached8192,tail16384))
        },
        "predetermined_environment_predictors":PRED,
        "rain_unit":"log1p(WorldClim annual precipitation mm)",
        "other_units":"WorldClim annual mean vapour, solar Wm2, BIO5 Celsius, DEM elevation m",
        "grid_degrees":GRID,
        "min_distinct_named_taxa_per_source_arm_genus_cell":MIN_GROUP_TAXA,
        "source_design":"Original 350399 photo IDs and original frozen 19BIO/VAPR/SRAD/DEM, no pixels relabelled; photo-species first cohort explicitly conditioned within genus x geographic stratum.",
        "photographic_environment_estimand":"Univariate within nominal genus × geographic cell × source arm taxon-cell slopes in WHITE ROIs, pp per pooled source-taxon cell predictor SD",
        "comparison_overviews":reports,
        "slopes":data,
        "limitations":[
            "Exploratory local within-genus among-species photographic association, not within-individual pigment change or causal selection.",
            "Source photo counts and available named taxa are not biological flowering abundances.",
            "No official observed-year local temperature, light, water deficit or precipitation anomalies fitted.",
            "No botanist verified open-flower petal organs in all 350399 original photo ROIs.",
            "Slope of each annual climatic factor is univariate; atmospheric vapour, rainfall, solar, heat, altitude are correlated, and simultaneous causal isolation is unsupported.",
            "Nominal genera and 5-20 degree regions are not independent clades/local populations, bootstrap is exploratory, no multiplicity correction.",
            "Original source photo capture/classifiability and species-first selection probabilities remain nonrandom or unknown; even within genus-region means are not a natural flora sample."
        ],
        "output_sha256":{
            OUT_CSV:sha256(output/OUT_CSV),
            OUT_CORR:sha256(output/OUT_CORR)
        }
    }
    (output/OUT_JSON).write_text(json.dumps(result,sort_keys=True,
                                           ensure_ascii=False,indent=2,
                                           allow_nan=False)+"\n")
    print("MOIS_FULL_REAL350399_ANNUAL_ABIOTIC_WITHIN_GENUS_REGION_SLOPES",
          json.dumps({"photo_hues_with_complete_abiotic":len(x),
            "selected_main_models":[{
                "grid":r["grid_degrees"],"min_taxa":r["min_named_taxa_per_genus_cell_arm"],
                "cohort":r["cohort"],"predictor":r["predictor"],
                "white_pp_per_SD":r["white_pp_per_one_full_taxon_cell_SD"],
                "genus_bootstrap95":r["exploratory_nominal_genus_bootstrap_95"],
                "genus_clusters":r["named_genera_with_within_cell_variation"],
                "taxon_cell_n":r["selected_named_taxon_cells"]
              } for r in data if r["grid_degrees"]==10 and
                    r["min_named_taxa_per_genus_cell_arm"]==2]},
          ensure_ascii=False,sort_keys=True),flush=True)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ("core","broad","wave2","cached8192","tail16384",
                "environment","env-receipt","output"):
        p.add_argument("--"+arg,type=Path,required=True)
    a=p.parse_args()
    run(a.core,a.broad,a.wave2,a.cached8192,a.tail16384,
        a.environment,a.env_receipt,a.output)


if __name__=="__main__":
    main()
