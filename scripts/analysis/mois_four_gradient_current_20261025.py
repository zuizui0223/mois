#!/usr/bin/env python3
"""Re-fit MOIS moisture, elevation, island and seasonal gradients from real signed photos.

350,399 *actual* original-image ROI-v4 attempts, 164,297 real automatic four-hue
labels, the SHA-matched 350,399-point original WorldClim annual covariates,
and the independently signed island GSHHG classification on their unchanged
313,535-photo overlapping baseline. Not a new 93m pixel-census result.
Outcome is photographic automatic WHITE among 4 classificable hues; all missing
ROIs stay missing. Neither phylogeny nor individual flowering are measured.
"""
from __future__ import annotations
import argparse, json, hashlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
from scipy.stats import t as student_t
import statsmodels.api as sm
import statsmodels.formula.api as smf

from scripts.analysis.mois_actual_photo_cohort_season_audit_20261019 import (
    assemble, sha256, summary, require,
)
from scripts.analysis.mois_real350399_annual_abiotic_geogenus_white_20261019 import (
    source_data, taxon_geocells, local_slopes, PRED,
)
from scripts.analysis.mois_real350399_taxon_season_decomposition_20261019 import (
    estimate as seasonal_estimate,
)

SCHEMA="mois_20261025_actual350399_four_macroecological_axes_refit_v1"
HUES=("white","yellow_orange","red_pink","blue_purple")
GEO_SRC=313535
GEO_CLASSES=151239
GEO_WHITE=71427
MAIN_N=350399
CLASS_N=164297

def safe_float(v):
    if v is None:return None
    f=float(v)
    return f if np.isfinite(f) else None

def check_geography(root):
    root=Path(root)
    meta=json.loads((root/"result.json").read_text(encoding="utf-8"))
    require(meta["schema"]=="mois_island_mainland_macro_photo_geography_exploratory_v1"
            and meta["source_photo_ids"]==GEO_SRC,
            "Geographic coastline/photograph source changed")
    for file,digest in meta["output_sha256"].items():
        require(sha256(root/file)==digest,
                "GSHHG source-derived photo/cell SHA changed: "+file)
    g=pd.read_csv(root/"photo_landmass_cache.csv.gz",
                  dtype={"photo_id":"string","landmass":"string"},
                  keep_default_na=False)
    require(len(g)==GEO_SRC and g.photo_id.is_unique and
            g.land_type.isin(("mainland","island","unassigned")).all(),
            "Signed island/mainland source cache is not original-photo unique")
    return meta,g

def island_comparisons(d,root):
    meta,geo=check_geography(root)
    photo=d.merge(geo,on="photo_id",how="inner",validate="one_to_one")
    require(len(photo)==GEO_SRC and photo.photo_id.is_unique,
            "The global 350399 source no longer includes ALL 313535 mapped IDs")
    require(int(photo.classified.sum())==GEO_CLASSES and
            int(photo.loc[photo.classified,"white"].sum())==GEO_WHITE,
            "The 313535 original photo/ROI classification was changed")
    raw=[]
    for land in ("island","mainland","unassigned"):
        v=photo[photo.land_type.eq(land)]
        known=v[v.classified]
        raw.append({"land_type":land,"photo_attempts":len(v),
                    "four_hue_classified":len(known),
                    "white":int(known.white.sum()),
                    "white_percent":100*float(known.white.mean()) if len(known) else None})
    island=next(q for q in raw if q["land_type"]=="island")
    mainland=next(q for q in raw if q["land_type"]=="mainland")
    require(island["four_hue_classified"]==18864 and
            mainland["four_hue_classified"]==129440,
            "Published original coastline classification counts changed")
    cells=pd.read_csv(root/"matched_island_mainland_cells.csv")
    expected=pd.read_csv(root/"island_mainland_adjusted.csv")
    formula="p_white ~ island + z_heat + z_rain + z_solar + z_elevation + C(region)"
    model=smf.ols(formula,cells).fit(
        cov_type="cluster",cov_kwds={"groups":cells.block},use_t=True)
    ci=model.conf_int().loc["island"]
    published=expected.loc[expected.model.eq("geography_climate_altitude")].iloc[0]
    effect=100*float(model.params["island"])
    require(abs(effect-published.island_minus_mainland_pp)<1e-8,
            "Spatial-adjusted source island model fails exact independent replay")

    # Reproduce the original caliper-1 optimum 1-to-1 (never use hue for matching).
    predictors=["z_heat","z_rain","z_solar","z_elevation"]
    pairs=[]
    for region,group in cells.groupby("region"):
        a=group[group.island.eq(1)]
        b=group[group.island.eq(0)]
        if a.empty or b.empty:continue
        dist=cdist(a[predictors].to_numpy(float),b[predictors].to_numpy(float))
        cost=np.concatenate(
            (np.where(dist<=1,dist,10000.0),
             np.full((len(a),len(a)),100.0)),axis=1)
        iis,jjs=linear_sum_assignment(cost)
        for i,j in zip(iis,jjs):
            if j<len(b) and dist[i,j]<=1:
                ia=a.iloc[i];mb=b.iloc[j]
                pairs.append({"region":region,"landmass":str(ia.landmass),
                              "island_white":float(ia.p_white),
                              "mainland_white":float(mb.p_white),
                              "diff":float(ia.p_white-mb.p_white),
                              "distance":float(dist[i,j])})
    pair=pd.DataFrame(pairs)
    require(len(pair)==102 and pair.region.nunique()==28,
            "Original environmental caliper support/pair count changed")
    paired=sm.OLS(pair["diff"],np.ones((len(pair),1))).fit(
        cov_type="cluster",cov_kwds={"groups":pair.region},use_t=True)
    paired_ci=paired.conf_int().iloc[0]
    # Equal matched pairs is NOT the same population/weighting as 1394-cell OLS.
    return {
        "source_photo_attempts":GEO_SRC,"classifiable":GEO_CLASSES,
        "raw_by_land_type":raw,
        "raw_island_minus_mainland_white_pp":
            island["white_percent"]-mainland["white_percent"],
        "spatial_climate_adjusted_1394_cell_white_pp":effect,
        "spatial_climate_adjusted_nominal_cluster_95_pp":
            [100*float(ci.iloc[0]),100*float(ci.iloc[1])],
        "spatial_climate_adjusted_nominal_p":float(model.pvalues["island"]),
        "spatial_climate_adjusted_source_cells":int(len(cells)),
        "original_coloured_outcomes_in_1394_cells":int(cells.n_photos.sum()),
        "matched_pairs":len(pair),
        "matched_island_landmasses":int(pair.landmass.nunique()),
        "matched_regions":int(pair.region.nunique()),
        "matched_pair_white_pp":100*float(paired.params.iloc[0]),
        "matched_pair_nominal_region_cluster_95_pp":
            [100*float(paired_ci.iloc[0]),100*float(paired_ci.iloc[1])],
        "matched_pair_nominal_region_p":float(paired.pvalues.iloc[0]),
        "matched_pair_island_white_percent":100*float(pair.island_white.mean()),
        "matched_pair_mainland_white_percent":100*float(pair.mainland_white.mean()),
        "no_actual_pollinator_availability_inferred":True,
    },pair

def paired_elevation(cells,*,width=5,min_taxa=5,min_photos=10):
    c=cells.loc[cells.n_taxa.ge(min_taxa)&cells.n_photos.ge(min_photos)].copy()
    c=c[np.isfinite(c.elevation)&np.isfinite(c.latitude)&np.isfinite(c.longitude)]
    c["latblock"]=np.floor((c.latitude+90)/width).astype(int)
    c["lonblock"]=np.floor((c.longitude+180)/width).astype(int)
    keys=["landmass","latblock","lonblock"]
    x=c.loc[c.elevation.ge(0)&c.elevation.lt(500)]
    y=c.loc[c.elevation.ge(1500)]
    def avg(frame):
        return frame.groupby(keys,dropna=False).agg(
            white=("p_white","mean"),
            n_cells=("p_white","size"),photos=("n_photos","sum"),
            n_taxa_cells=("n_taxa","sum")).reset_index()
    paired=avg(y).merge(avg(x),on=keys,validate="one_to_one",
                        suffixes=("_high","_low"))
    require(len(paired)>=15,"No geographically matched high/low altitude support")
    differences=100*(paired.white_high-paired.white_low)
    low=float(differences.mean())
    if len(differences)>1:
        sd=float(differences.std(ddof=1)/np.sqrt(len(differences)))
        interval=student_t.interval(.95,df=len(differences)-1,loc=low,scale=sd)
    else:interval=[None,None]
    result={
        "cell_min_source_named_taxa":min_taxa,
        "cell_min_classifiable_photos":min_photos,
        "local_neighbourhood_width_degrees":width,
        "paired_landmass_geographic_neighbourhoods":int(len(paired)),
        "mean_low_0_500m_white_percent":100*float(paired.white_low.mean()),
        "mean_high_ge1500m_white_percent":100*float(paired.white_high.mean()),
        "high_minus_low_white_pp":low,
        "nominal_neighbourhood_t_95_pp":[safe_float(z) for z in interval],
        "neighbourhoods_with_lower_high_elevation_white":
            int((differences<0).sum()),
        "classified_low_photo_records":int(paired.photos_low.sum()),
        "classified_high_photo_records":int(paired.photos_high.sum()),
        "low_source_cells":int(paired.n_cells_low.sum()),
        "high_source_cells":int(paired.n_cells_high.sum()),
    }
    return result,paired

def abiotic_slopes(x):
    outputs=[]
    for degrees,minimum,draws in ((5,2,399),(10,3,399)):
        t=taxon_geocells(x,degrees)
        fits,retained=local_slopes(t,degrees=degrees,min_taxa=minimum,draws=draws)
        outputs.extend(fits)
    selected=[
        dict(grid_degrees=r["grid_degrees"],
             min_taxa_per_genus_region_arm=r["min_named_taxa_per_genus_cell_arm"],
             cohort=r["cohort"],predictor=r["predictor"],
             pp_per_full_source_taxon_cell_sd=r["white_pp_per_one_full_taxon_cell_SD"],
             genus_bootstrap_nominal_95_pp=r["exploratory_nominal_genus_bootstrap_95"],
             named_genera_support=r["named_genera_with_within_cell_variation"],
             taxon_geocell_support=r["selected_named_taxon_cells"])
        for r in outputs
        if r["predictor"] in (
            "log1p_annual_rainfall","mean_atmospheric_vapour","elevation")
    ]
    # Descriptive wet-to-dry and low-to-high ranked source taxon cells.
    t=taxon_geocells(x,5)
    quintiles={}
    for variable in ("log1p_annual_rainfall","elevation"):
        valid=t[np.isfinite(t[variable])].copy()
        valid["bin"]=pd.qcut(valid[variable],q=5,labels=False,duplicates="drop")
        quintiles[variable]=[
            {"quintile_low_to_high":int(k)+1,
             "taxon_geocells":len(g),
             "white_percent_taxon_geocell_equal":100*float(g.p_white.mean()),
             "median_original_environment":float(g[variable].median())}
            for k,g in valid.groupby("bin",observed=True)
        ]
    return selected,quintiles

def run(core,broad,wave2,cached8192,tail16384,environment,env_receipt,
        geography,output):
    output=Path(output)
    require(not output.exists(),"Source-backed result must not overwrite an earlier analysis")
    d,arms=assemble(core,broad,wave2,cached8192,tail16384,strict=True)
    require(len(d)==MAIN_N and d.photo_id.is_unique
            and int(d.classified.sum())==CLASS_N,
            "Not the full 350399 observed-photo four-hue source")
    x,env_proof=source_data(d,environment,env_receipt)
    water_slopes,quintiles=abiotic_slopes(x)

    geography_report,pairs=island_comparisons(d,geography)
    cell=pd.read_csv(Path(geography)/"species_equal_1degree_landmass_cells.csv.gz")
    # Frozen 313535 GSHHG site assignments and taxon-cell means are aligned to
    # a true subset of the same 350399 source, not the full 93m metadata list.
    elev_primary,nearby=paired_elevation(cell,width=5)
    elev_sensitivity,_=paired_elevation(cell,width=10)
    season=seasonal_estimate(d,draws=199,seed=20261025)

    report={
        "schema":SCHEMA,
        "source_actual_photo_attempts":len(d),
        "source_four_hue_classified":int(d.classified.sum()),
        "technical_unknowns_not_reassigned":int((~d.classified).sum()),
        "source_photo_hues":{h:int(d.morph.eq(h).sum()) for h in HUES},
        "source_classified_complete_annual_five_predictors":len(x),
        "source_photo_and_climate_original_identity_sha256":{
            n:sha256(f) for n,f in zip(
                ("core","cumulative_broad","wave2","cached8192","tail16384","annual_climate"),
                (core,broad,wave2,cached8192,tail16384,environment))},
        "estimand":"Photography-derived original frozen four-hue white within classifiable ROIs",
        "axes":{
            "wet_dry":{
                "measured_photos_with_all_annual_env":len(x),
                "univariate_within_genus_region_and_source_arm_slopes":water_slopes,
                "source_taxon_geocell_environment_quintiles":quintiles,
                "causal_environmental_driver_isolated":False,
            },
            "elevation":{
                "landmass_5deg_313535_baseline":elev_primary,
                "landmass_10deg_sensitivity":elev_sensitivity,
                "same_species_elevation_plasticity_proven":False,
            },
            "islands":geography_report,
            "season":{
                "classifiable_dated_summer_winter_photos":
                    season["classified_extratropical_dated_summer_winter_photos"],
                "hemispheres":season["hemispheres"],
                "observation_equal_sensitivity":season["observation_equal_sensitivity"],
                "source_arm_seasonal_photo_counts":season["source_arm_seasonal_photo_counts"],
                "phenotypic_within_individual_seasonal_change_proven":False,
            }
        },
        "limits":[
            "350399 real photo attempts, 164297 four-hue machine ROI labels: the latest 1019272/449741 checkpoint lacks a source-complete covariate/island/season join.",
            "Source first-taxon enrichment and opportunistic photography do not identify global blooming flora prevalence.",
            "Water, heat and altitude are correlated; reported within-genus locality slopes are univariate not causal.",
            "Island GSHHG coverage is only the original 313535-site/151239-colour nested baseline, not every 350399 photo.",
            "Seasonal decomposition is exact algebra over photograph source-taxonomy/effort; not pigment plasticity or weather anomalies.",
            "Unknown ROI, fruit/leaf ambiguity, visible-only colour and open-flower status are not independently botanically validated.",
            "Nominal genus/spatial CIs and multiple exploratory axes lack familywise correction.",
        ],
    }
    output.mkdir(parents=True)
    pd.DataFrame(water_slopes).to_csv(output/"mois_four_axis_water_elevation_local_taxon_slopes.csv",index=False)
    pd.DataFrame([elev_primary,elev_sensitivity]).to_csv(
        output/"mois_four_axis_paired_elevation_by_geographic_scale.csv",index=False)
    pairs.to_csv(output/"mois_four_axis_island_matched_pairs_source.csv",index=False)
    nearby.to_csv(output/"mois_four_axis_elevation_matched_cells_source.csv",index=False)
    with (output/"mois_four_axis_20261025_source_reanalysis.json").open("w") as f:
        json.dump(report,f,indent=2,sort_keys=True,ensure_ascii=False,allow_nan=False)
        f.write("\n")
    report["output_file_sha256"]={
        p.name:sha256(p) for p in output.iterdir() if p.is_file()
    }
    (output/"mois_four_axis_result_receipt.json").write_text(
        json.dumps({"schema":SCHEMA,
                    "source_photo_attempts":MAIN_N,
                    "classified_photo_ROIs":CLASS_N,
                    "output_sha256":report["output_file_sha256"],
                    "new_image_downloads_or_reinferences":0},
                   indent=2,sort_keys=True)+"\n")
    print("MOIS_20261025_REAL_350399_FOUR_ECOLOGICAL_GRADIENTS",
          json.dumps({
            "attempted_photos":len(d),
            "classified":int(d.classified.sum()),
            "annual_complete":len(x),
            "moisture_local_slopes_main_5deg":[v for v in water_slopes
                if v["grid_degrees"]==5
                and v["cohort"]=="combined_arm_conditioned"
                and v["predictor"] in ("log1p_annual_rainfall","mean_atmospheric_vapour","elevation")],
            "nearby_elevation":elev_primary,
            "matched_island":{k:geography_report[k] for k in (
                "raw_island_minus_mainland_white_pp",
                "spatial_climate_adjusted_1394_cell_white_pp",
                "matched_pairs","matched_pair_white_pp",
                "matched_pair_nominal_region_cluster_95_pp")},
            "season_photo_summer_minus_winter":{
                hemi:zone["raw_photo_white_summer_minus_winter_pp"]
                for hemi,zone in season["hemispheres"].items()},
        },sort_keys=True,allow_nan=False),flush=True)
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for field in ("core","broad","wave2","cached8192","tail16384",
                  "environment","env-receipt","geography","output"):
        p.add_argument("--"+field,type=Path,required=True)
    a=p.parse_args()
    run(a.core,a.broad,a.wave2,a.cached8192,a.tail16384,
        a.environment,a.env_receipt,a.geography,a.output)

if __name__=="__main__":main()
