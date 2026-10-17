#!/usr/bin/env python3
"""MOIS main ecological result: near-half photographic WHITE composition along abiotic gradients.

TRUE 313535 original ROI pixel attempts; 151239 classifiable, 71427 white.
The source is photographic, NOT a 50:50 evolutionary equilibrium, pigment
synthesis-cost measurement, or sample of verified open flowers. Season is SI.
"""
from __future__ import annotations
import argparse, hashlib, json, re
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

RELEASE="mois-primary-313535-fullfourth4096-featurelake-v2"
PHOTO="primary313535_fullfourth4096_original_photo_hues.parquet"
ANNUAL="primary313535_fullfourth4096_annual_worldclim.parquet"
STATIC="primary313535_original_partial_soil_terrain.parquet"
MANIFEST="primary313535_fullfourth4096_featurelake_manifest.json"
N_PHOTOS,N_HUES,N_DATED,WHITE_N=313535,151239,151230,71427
HUES=("white","yellow_orange","red_pink","blue_purple")
RAIN="wc2_bio12_10m"
DRY="wc2_bio14_10m"
SEASON="wc2_bio15_10m"
HEAT="wc2_bio05_10m"
SOLAR="annual_solar_Wm2"
CLIMATE=(RAIN,DRY,SEASON,HEAT,SOLAR)
GEO=("latitude","abs_lat","lon_sin","lon_cos","new_source_fraction")
TERRAIN=("wc2_elev_2p5m","terrain_slope_2p5m_degree","terrain_TRI_2p5m_m")
SOIL=("soil_pH","soil_SOC","soil_N","soil_clay","soil_available_water_proxy")
OUT_TWOD="mois_313535_white_fractions_water_light_heat_tertiles.csv"
OUT_LAT="mois_313535_latitude_solar_stratified_wet_dry_contrast.csv"
OUT_FITS="mois_313535_white_genus_conditioned_abiotic_terms.csv"
OUT_BLOCK="mois_313535_equal_genus_nested_habitat_blocks.csv"
OUT_DECISION="mois_313535_white_abiotic_environment_partition_decision.json"

def require(ok,why):
    if not ok: raise ValueError(why)

def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):h.update(b)
    return h.hexdigest()

def load_source(root):
    rec=json.loads((root/MANIFEST).read_text())
    require(rec["schema"]=="mois_313535_genuine_original_ROI_feature_lake_v2"
        and rec["actual_original_photo_attempts"]==N_PHOTOS
        and rec["four_hue_labels"]==N_HUES
        and rec["dated_four_hue_labels"]==N_DATED
        and rec["new_photo_pixels"]==4096
        and rec["new_classifiable_visible_hues"]==1577
        and rec["old512_1792_photos_not_double_counted"] is True
        and rec["season_targeted_256_source_excluded_from_primary_hue_balance"] is True,
        "Source-verified all-original RGB/photo Parquet release v2 missing")
    for name in (PHOTO,ANNUAL,STATIC):
        require((root/name).is_file() and digest(root/name)==rec["files_sha256"][name],
                "Immutable original site, colour or topographic Parquet SHA drift: "+name)
    p=pd.read_parquet(root/PHOTO)
    a=pd.read_parquet(root/ANNUAL)
    s=pd.read_parquet(root/STATIC)
    require(len(p)==len(a)==len(s)==N_PHOTOS and p.photo_id.is_unique
            and p.photo_id.tolist()==a.photo_id.tolist()==s.photo_id.tolist(),
            "Original photo identifiers must align across the three immutable sources")
    d=p.copy()
    for v in CLIMATE:d[v]=pd.to_numeric(a[v],errors="coerce")
    for v in (*SOIL,*TERRAIN):d[v]=pd.to_numeric(s[v],errors="coerce")
    for v in ("latitude","longitude"):d[v]=pd.to_numeric(d[v],errors="coerce")
    require(int(d.classified.sum())==N_HUES and int(d.white.sum())==WHITE_N
            and int((d.classified&d.qualified_date.astype(str).str.fullmatch(
                 r"\d{4}-\d{2}-\d{2}")).sum())==N_DATED
            and d.classified.equals(d.morph.isin(HUES))
            and d.white.equals(d.morph.eq("white")),
            "Photographic hue, date and unknown-state labels changed")
    return d,rec

def genus(name):
    if not isinstance(name,str) or not name.strip():return None
    word=name.strip().split()[0]
    return word if re.fullmatch(r"[A-Z][A-Za-z]+",word) else None

def rates(x):
    c=x[x.classified]
    taxa=c.groupby("inat_taxon_id",observed=True).agg(
        white=("white","mean"),genus=("genus","first"))
    genera=taxa.dropna(subset=["genus"]).groupby("genus",observed=True).white.mean()
    return {
        "photo_attempts":int(len(x)),
        "classified_original_photo_ROIs":int(len(c)),
        "photo_white_ROIs":int(c.white.sum()),
        "photo_chromatic_ROIs":int(len(c)-c.white.sum()),
        "unclassified_original_photo_ROIs":int(len(x)-len(c)),
        "classification_yield":float(len(c)/len(x)) if len(x) else None,
        "photo_white_fraction":float(c.white.mean()) if len(c) else None,
        "equal_named_taxon_white_fraction":float(taxa.white.mean()) if len(taxa) else None,
        "equal_genus_white_fraction":float(genera.mean()) if len(genera) else None,
        "classifiable_named_taxa":int(c.inat_taxon_id.nunique()),
        "classifiable_genera":int(len(genera)),
        **{"hue_count_"+k:int(c.morph.eq(k).sum()) for k in HUES},
    }

def mark_tertiles(d):
    x=d.copy();cuts={}
    for var,name in ((RAIN,"rain"),(SOLAR,"solar"),(HEAT,"heat")):
        z=pd.to_numeric(x[var],errors="coerce")
        valid=z.notna()&np.isfinite(z)
        require(valid.sum()>100 and z[valid].nunique()>3,
                "Missing actually sampled climate range: "+var)
        q1,q2=np.quantile(z[valid].to_numpy(float),[1/3,2/3])
        require(q1<q2,"Environmental climate thirds are degenerate")
        x[name+"_band"]=pd.cut(z,[-np.inf,q1,q2,np.inf],
            labels=("low","middle","high"),include_lowest=True)
        cuts[var]={"low_to_middle":float(q1),"middle_to_high":float(q2),
                   "cutoffs_derived_from_ALL_original_processed_photos":True}
    return x,cuts

def two_way_tables(d):
    rows=[]
    for other in ("solar_band","heat_band"):
        q=d.dropna(subset=["rain_band",other])
        for (rb,ob),g in q.groupby(["rain_band",other],observed=True):
            rows.append({
                "conditioning_axis":other,
                "rain_tertile":str(rb),"conditioning_tertile":str(ob),
                "median_annual_rainfall_mm":float(g[RAIN].median()),
                "median_annual_solar_W_m2":float(g[SOLAR].median()),
                "median_warmest_month_temp_C":float(g[HEAT].median()),
                **rates(g),"descriptive_not_measured_floral_fitness":True})
    return pd.DataFrame(rows)

def latitude_strata_wet_dry(d,*,min_photo=80,min_taxa=20,draws=999):
    q=d[d.latitude.between(-90,90)&d.longitude.between(-180,180)
        &d.solar_band.notna()&d.rain_band.isin(("low","high"))].copy()
    q["abs_lat_10deg"]=np.floor(q.latitude.abs()/10).astype(int)
    q["hemisphere"]=np.where(q.latitude.ge(0),"north","south")
    rows=[]
    for (hemi,lat,solar),v in q.groupby(
            ["hemisphere","abs_lat_10deg","solar_band"],observed=True):
        d0,w0=rates(v[v.rain_band.eq("low")]),rates(v[v.rain_band.eq("high")])
        if min(d0["classified_original_photo_ROIs"],
               w0["classified_original_photo_ROIs"])<min_photo:continue
        if min(d0["classifiable_named_taxa"],w0["classifiable_named_taxa"])<min_taxa:continue
        rows.append({
          "hemisphere":str(hemi),"abs_latitude_band_10deg":int(lat),
          "solar_tertile":str(solar),
          "dry_classified_original_photos":d0["classified_original_photo_ROIs"],
          "wet_classified_original_photos":w0["classified_original_photo_ROIs"],
          "dry_named_taxa":d0["classifiable_named_taxa"],
          "wet_named_taxa":w0["classifiable_named_taxa"],
          "dry_photo_white_fraction":d0["photo_white_fraction"],
          "wet_photo_white_fraction":w0["photo_white_fraction"],
          "wet_minus_dry_white_percentage_points":100*(
              w0["photo_white_fraction"]-d0["photo_white_fraction"]),
          "within_latitude_solar_only_not_full_geographic_matching":True,
        })
    t=pd.DataFrame(rows)
    if t.empty:return t,{"matched_strata":0,"bootstrap95_pp":None,
                           "point_equal_stratum_pp":None}
    val=t.wet_minus_dry_white_percentage_points.to_numpy(float)
    rng=np.random.default_rng(20261018)
    boot=val[rng.integers(0,len(val),(draws,len(val)))].mean(axis=1)
    return t,{
      "matched_strata":int(len(t)),
      "photo_observations_in_these_strata":int(
          (t.dry_classified_original_photos+t.wet_classified_original_photos).sum()),
      "point_equal_stratum_pp":float(val.mean()),
      "bootstrap95_pp":[float(a) for a in np.percentile(boot,[2.5,97.5])],
      "bootstrap_over_latitude_solar_strata_not_geographically_independent":True,
      "not_causal_moisture_response":True,
    }

def species_panel(d,required,min_genus_taxa=2):
    p=d[d.classified&d.inat_taxon_id.notna()].dropna(subset=list(required)).copy()
    p=p[p.genus.notna()&p.inat_taxon_id.astype(str).ne("")]
    cols=[k for k in set((*GEO,*CLIMATE,*TERRAIN,*SOIL)) if k in p.columns]
    agg={k:"mean" for k in cols}
    agg.update({"white":"mean","photo_id":"size","genus":"first"})
    s=p.groupby("inat_taxon_id",observed=True).agg(agg).rename(
        columns={"white":"photo_white_rate","photo_id":"n_photo_ROIs"})
    counts=s.groupby("genus",observed=True).size()
    s=s[s.genus.isin(counts[counts>=min_genus_taxa].index)].copy()
    require(len(s)>150 and s.genus.nunique()>=50,
            "Insufficient taxa and genera on matched photographic source cohort")
    return s,{"classified_photo_rows_mapped_to_species":int(len(p)),
              "taxa_with_congeneric_support":int(len(s)),
              "qualifying_genera":int(s.genus.nunique()),
              "nominal_taxon_not_independent_genomic_lineage":True}

def model_fit(s,features,*,interaction_features=frozenset()):
    names=[];cols=[]
    for field in features:
        v=s[field].to_numpy(float)
        require(np.isfinite(v).all(),"Photographic environmental covariate missing")
        if field in interaction_features:z=v
        else:
            sd=float(v.std())
            if sd<1e-12:continue
            z=(v-float(v.mean()))/sd
        names.append(field);cols.append(z)
    require(len(cols)>0,"No environmental model terms")
    x=pd.DataFrame(np.column_stack(cols),index=s.index,columns=names)
    x=x-x.groupby(s.genus).transform("mean")
    y=(s.photo_white_rate -
       s.photo_white_rate.groupby(s.genus).transform("mean")).to_numpy(float)
    x=x.loc[:,x.std(axis=0).to_numpy(float)>1e-10]
    names=x.columns.tolist()
    from scipy.linalg import qr
    matrix=x.to_numpy(dtype=float)
    _,r,pivot=qr(matrix,mode="economic",pivoting=True)
    diagonal=np.abs(np.diag(r))
    rank=int((diagonal>diagonal.max()*1e-9).sum())
    keep=np.sort(pivot[:rank])
    dropped=[v for i,v in enumerate(names) if i not in keep]
    matrix=matrix[:,keep]
    names=[names[i] for i in keep]
    w=1/s.groupby("genus").genus.transform("size").to_numpy(float)
    fit=sm.WLS(y,matrix,weights=w).fit(cov_type="cluster",
                    cov_kwds={"groups":s.genus.to_numpy()},use_t=False)
    require(np.isfinite(fit.params).all(),"Genus-cluster photo model invalid")
    r2=1-float(np.sum(w*(y-fit.fittedvalues)**2)/
               max(1e-20,np.sum(w*y*y)))
    ci=np.asarray(fit.conf_int())
    terms=[]
    for i,field in enumerate(names):
        terms.append({
            "variable":field,
            "coefficient_photo_white_pp_per_standardised_environment_unit":
                float(100*fit.params[i]),
            "nominal_genus_cluster_95CI_low_pp":float(100*ci[i,0]),
            "nominal_genus_cluster_95CI_high_pp":float(100*ci[i,1]),
            "nominal_genus_cluster_p":float(fit.pvalues[i]),
            "interaction_per_one_SD_rain_times_one_SD_other_axis":
                bool(field in interaction_features),
            "not_original_multiplicity_protected_maxT_24_family":True,
        })
    return {"within_genus_in_sample_R2":r2,"n_genus_taxa":int(len(s)),
            "n_genera":int(s.genus.nunique()),
            "dropped_collinear_predictors":dropped,"effects":terms}

def fit_habitat_blocks(d):
    b=d.copy()
    b["abs_lat"]=b.latitude.abs()
    b["lon_sin"]=np.sin(np.deg2rad(pd.to_numeric(b.longitude,errors="coerce")))
    b["lon_cos"]=np.cos(np.deg2rad(pd.to_numeric(b.longitude,errors="coerce")))
    b["new_source_fraction"]=(b.old_mois_or_new_fcp.astype(str)
                               .ne("historical")).astype(float)
    b["genus"]=b.species.map(genus)
    main=(*GEO,*CLIMATE)
    s,info=species_panel(b,main)
    for v,name in ((RAIN,"rain_z"),(SOLAR,"solar_z"),(HEAT,"heat_z")):
        a=s[v].to_numpy(float)
        require(a.std()>0,"No mapped climate variation")
        s[name]=(a-a.mean())/a.std()
    s["rain_x_solar"]=s.rain_z*s.solar_z
    s["rain_x_heat"]=s.rain_z*s.heat_z
    before=model_fit(s,main)
    after=model_fit(s,(*main,"rain_x_solar","rain_x_heat"),
                    interaction_features={"rain_x_solar","rain_x_heat"})
    summary={"all_annual_climate_source":{
        **info,
        "geography_and_annual_climate_R2":before["within_genus_in_sample_R2"],
        "geography_climate_plus_interactions_R2":after["within_genus_in_sample_R2"],
        "interaction_delta_in_sample_R2":after["within_genus_in_sample_R2"]-
            before["within_genus_in_sample_R2"],
        "not_out_of_geography_or_genus_predictive_improvement":True,
    }}
    terms=[{"model_block":label,**term}
           for label,res in (("annual_all_base",before),
                             ("annual_all_plus_rain_solar_heat",after))
           for term in res["effects"]]
    selected=[v for v in terms
              if v["model_block"]=="annual_all_plus_rain_solar_heat" and
              v["variable"] in ("rain_x_solar","rain_x_heat")]
    if selected:
        q=multipletests([v["nominal_genus_cluster_p"] for v in selected],
                        method="fdr_bh")[1]
        for row,val in zip(selected,q):
            row["BH_q_exploratory_two_interaction_tests"]=float(val)
    blocks=[]
    for label,extra in (("landscape_terrain",TERRAIN),
                        ("partially_covered_soil",SOIL)):
        sub,matched=species_panel(b,(*main,*extra))
        old=model_fit(sub,main);new=model_fit(sub,(*main,*extra))
        summary[label]={
            **matched,
            "matched_original_photo_support_within_block":True,
            "baseline_geography_climate_R2":old["within_genus_in_sample_R2"],
            "added_static_block_R2":new["within_genus_in_sample_R2"],
            "matched_photo_source_delta_in_sample_R2":
                new["within_genus_in_sample_R2"]-old["within_genus_in_sample_R2"],
            "cannot_compare_between_differently_covered_soil_terrain_cohorts":True,
        }
        terms.extend({"model_block":label,**v} for v in new["effects"])
        blocks.append({"block":label,**summary[label]})
    return pd.DataFrame(terms),pd.DataFrame(blocks),summary

def run(source,output):
    require(not output.exists() or not any(output.iterdir()),
            "Never replace earlier signed ecological source data")
    photo,receipt=load_source(source)
    photo["genus"]=photo.species.map(genus)
    b,cuts=mark_tertiles(photo)
    composition=two_way_tables(b)
    lat_strata,lat_receipt=latitude_strata_wet_dry(b)
    fit,blocks,models=fit_habitat_blocks(photo)
    require(not composition.empty and not fit.empty and not blocks.empty,
            "Genuine 313535 photographed environmental sample unavailable")
    output.mkdir(parents=True,exist_ok=True)
    outputs=((OUT_TWOD,composition),(OUT_LAT,lat_strata),
             (OUT_FITS,fit),(OUT_BLOCK,blocks))
    for name,table in outputs:
        table.to_csv(output/name,index=False,lineterminator="\n")
    totals=rates(photo)
    result={
        "schema":"mois_real313535_white_water_solar_heat_spatial_environment_partition_v1",
        "source_photo_release":RELEASE,
        "actual_RGB_original_photo_attempts":N_PHOTOS,
        "classifiable_original_four_hue_photo_ROIs":N_HUES,
        "original_white_ROI_photos":WHITE_N,
        "original_photo_source_frozen_by_SHA256":
            receipt["full_16_workers_original_pixel_source_SHA256"],
        "primary_white_chromatic_photo_counts":totals,
        "annual_climate_tertile_cutoffs":cuts,
        "water_light_and_water_heat_3x3_source_cells":int(len(composition)),
        "latitude_solar_stratified_wet_dry":lat_receipt,
        "genus_cluster_observational_abiotic_model_receipts":models,
        "four_hue_not_equal_to_anthocyanin_chemistry":True,
        "photo_hue_not_botanically_confirmed_open_petals":True,
        "true_global_flowering_population_white_fraction_not_identified":True,
        "pigment_synthesis_cost_fitness_NOT_measured":True,
        "no_50_50_adaptive_equilibrium_is_assumed":True,
        "site_year_monthly_climate_seasonal_SI_not_fitted":True,
        "1970_2000_annual_climate_normals_not_observed_year_month_weather":True,
        "season_enriched_photo_cohort_EXCLUDED_from_main_sample":True,
        "all_four_hue_unknowns_remain_missing_not_nonwhite":True,
        "all_model_tests_exploratory_beyond_original_maxT_24_family":True,
        "output_SHA256":{name:digest(output/name) for name,_ in outputs},
    }
    (output/OUT_DECISION).write_text(json.dumps(
        result,indent=2,sort_keys=True,ensure_ascii=False,
        allow_nan=False)+"\n")
    print("MOIS_GLOBAL_WHITE_WATER_SOLAR_HEAT_HABITAT",json.dumps({
        "n_classified":totals["classified_original_photo_ROIs"],
        "photo_white_share":totals["photo_white_fraction"],
        "latitude_solar_stratified_wet_dry":lat_receipt,
        "annual_interaction_model":models["all_annual_climate_source"],
        "matched_terrain":models["landscape_terrain"],
        "matched_soil":models["partially_covered_soil"],
        "annual_climate_and_interaction_effects": [
            {k:row[k] for k in (
                "variable",
                "coefficient_photo_white_pp_per_standardised_environment_unit",
                "nominal_genus_cluster_95CI_low_pp",
                "nominal_genus_cluster_95CI_high_pp",
                "nominal_genus_cluster_p",
            )}
            for row in fit.loc[
                fit.model_block.eq("annual_all_plus_rain_solar_heat") &
                fit.variable.isin((RAIN, SOLAR, HEAT, "rain_x_solar", "rain_x_heat"))
            ].to_dict("records")
        ],
        "twod_rain_solar":composition[composition.conditioning_axis=="solar_band"][
            ["rain_tertile","conditioning_tertile",
             "classified_original_photo_ROIs","photo_white_fraction",
             "equal_genus_white_fraction"]].to_dict("records")
    },sort_keys=True),flush=True)
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--photo-feature-release",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    x=p.parse_args()
    run(x.photo_feature_release,x.output)
