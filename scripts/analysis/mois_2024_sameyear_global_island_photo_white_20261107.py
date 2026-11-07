#!/usr/bin/env python3
"""Recompute 2024 original MOIS colour × VIIRS/WorldPop Global2 same-year association.

ONLY historical_main source photos, 2024 original strict ROI9 white-mass
measurements and uniquely WCVP-accepted taxa. A taxon within a cell gets
equal weight independent of photo number; cells are also equally weighted.
Photo cells overlapping more than one GSHHG original landmass are excluded
rather than joined to an incorrect island/mainland. NO genuine disappearance,
moth-mediated selection, or unknown-colour imputation is inferred.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np,pandas as pd,statsmodels.formula.api as sf

SOURCE_SHA={
 "canonical":"f033661d76e07898a5c13199a55c5d0ce6e55353daa48b0058c456979e90644c",
 "wcvp":"dc5bb0e3686fde0cfd327e4336245d0ceae3c7657b3bd12841bb814cc94d1eb2",
 "photo_landmass":"f7f65e6bc167089b3034555c655e777f528ecdd29952ab7ee9833d8722f39267",
 "sameyear_2024_environment":"6e0a077d496e90848d05d043bffae57e4772519f4572a457221c5ef4da0fb813"}
COVARS="bio5 + log_rain + elevation + log_attempts"
def check(x,msg):
    if not bool(x):raise ValueError("MOIS_2024_WHITE_SOURCE: "+msg)
def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(4<<20),b""):h.update(b)
    return h.hexdigest()
def cell_id(lon,lat):
    x=np.minimum(np.floor((lon+180)/.05).astype(int),7199)
    y=np.minimum(np.floor((lat+90)/.05).astype(int),3599)
    return x.astype(str)+"_"+y.astype(str)
def make_frames(canonical,wcvp,geo,env):
    check(len(canonical)==1032882 and canonical.photo_id.is_unique and
          canonical.observation_id.is_unique,"One-photo-per-observation ID drift")
    check(len(geo)==69196 and geo.target_id.is_unique,"Original photo landmass source drift")
    check(len(env)==66456 and env.target_id.is_unique and int(env.complete.sum())==43595,
          "Same-year VIIRS+Global2 original source drift")
    check(set(env.target_id)==set(geo.loc[geo.land_class.isin(("island","mainland")),"target_id"]),
          "2024 photo/environment accepted landmass IDs not the same cohort")
    # Do not guess which original landmass the image belonged to if the
    # coordinate cell has more than one L1 landmass member.
    ambiguous=geo.cell_id.duplicated(keep=False)
    pure=geo.loc[~ambiguous & geo.land_class.isin(("island","mainland"))].copy()
    check(pure.cell_id.is_unique,"Source cell still ambiguously shared across landmasses")
    d=canonical.loc[
       canonical.observed_on.astype(str).str[:4].eq("2024") &
       canonical.source_subarm.eq("historical_main") &
       canonical.latitude.between(-90,90)&canonical.longitude.between(-180,180)].copy()
    d["cell_id"]=cell_id(d.longitude,d.latitude)
    keep=pure[["cell_id","target_id","land_id","land_class"]]
    d=d.merge(keep,on="cell_id",how="inner",validate="many_to_one")
    check(d.photo_id.is_unique and d.observation_id.is_unique,
          "Landmass photo source join duplicated original observation")
    lookup=wcvp.loc[wcvp.wcvp_match_state.eq("unique_exact_name_to_accepted_id"),
       ["inat_taxon_id","wcvp_accepted_plant_name_id"]].copy()
    check(lookup.inat_taxon_id.is_unique,"WCVP exact accepted source duplicated")
    d=d.merge(lookup,on="inat_taxon_id",how="left",validate="many_to_one")
    climatic=d.groupby("target_id").agg(
       photo_attempts=("photo_id","size"),strict9_success=("strict9","sum"),
       bio5=("wc2_bio05_10m","mean"),bio12=("wc2_bio12_10m","mean"),
       elevation=("wc2_elev_2p5m","mean")).reset_index()
    win=d.loc[d.strict9 & d.wcvp_accepted_plant_name_id.notna() &
              d.strict9_white.notna(),["target_id","wcvp_accepted_plant_name_id","strict9_white"]].copy()
    check(win.strict9_white.between(0,1).all(),"ROI9 white-mass values invalid")
    per_species=win.groupby(["target_id","wcvp_accepted_plant_name_id"],
                             observed=True).strict9_white.mean().reset_index()
    per_cell=per_species.groupby("target_id",observed=True).agg(
       white_equal=("strict9_white","mean"),
       accepted_species_n=("wcvp_accepted_plant_name_id","nunique")).reset_index()
    data=env.merge(per_cell,on="target_id",how="left",validate="one_to_one")
    data=data.merge(climatic,on="target_id",how="left",validate="one_to_one")
    data=data.merge(pure[["target_id","land_id"]],on="target_id",how="left",validate="one_to_one")
    data["photo_eligible"]=data.photo_attempts.ge(5)&data.strict9_success.ge(3)&data.accepted_species_n.ge(3)
    data["viirs_eligible"]=data.viirs_valid_months.ge(3)&data.viirs_median.notna()
    data["pop_eligible"]=data.population_source_status.eq("OK")&data.population_density.notna()
    data["white_pp"]=100*data.white_equal
    data["log_light"]=np.log1p(data.viirs_median.clip(lower=0))
    data["log_pop"]=np.log1p(data.population_density.clip(lower=0))
    data["log_rain"]=np.log1p(data.bio12.clip(lower=0))
    data["log_attempts"]=np.log1p(data.photo_attempts.clip(lower=0))
    data["region5"]=(np.floor((data.query_longitude+180)/5).astype("int64").astype(str)+"_"+
                     np.floor((data.query_latitude+90)/5).astype("int64").astype(str))
    b=data.loc[data.photo_eligible&data.viirs_eligible].dropna(
       subset=["white_pp","log_light","log_rain","bio5","elevation","log_attempts","region5"]).copy()
    m=b.loc[b.pop_eligible].dropna(subset=["log_pop"]).copy()
    b["is_island"]=b.land_class.eq("island").astype(int)
    m["is_island"]=m.land_class.eq("island").astype(int)
    return d,data,b,m,per_species,int(ambiguous.sum())
def cluster_fit(data,formula,cluster):
    check(data[cluster].nunique()>=5,"Too few source geographic clusters")
    fit=sf.ols(formula,data).fit(cov_type="cluster",
        cov_kwds={"groups":data[cluster]},use_t=True)
    out={}
    for term in ("log_light","is_island","log_light:is_island"):
        if term not in fit.params:continue
        lo,hi=fit.conf_int().loc[term]
        out[term]={"beta":float(fit.params[term]),
           "ci_low":float(lo),"ci_high":float(hi),
           "SE":float(fit.bse[term]),"p":float(fit.pvalues[term])}
    return out
def island_bootstrap(island,B=1500,seed=20261107):
    xs=["log_light","log_pop","bio5","log_rain","elevation","log_attempts"]
    arrays=[]
    for _,part in island.groupby("land_id"):
        x=part[xs].to_numpy(float);y=part.white_pp.to_numpy(float)
        arrays.append((x-x.mean(axis=0),y-y.mean()))
    check(len(arrays)==17,"Original island number changed")
    rg=np.random.default_rng(seed)
    draws=[]
    for i in range(B):
        inds=rg.choice(len(arrays),len(arrays),replace=True)
        X=np.concatenate([arrays[k][0] for k in inds])
        Y=np.concatenate([arrays[k][1] for k in inds])
        beta=np.linalg.lstsq(X,Y,rcond=None)[0][0]
        check(np.isfinite(beta),"Nonfinite bootstrap source coefficient")
        draws.append(float(beta))
    return {"bootstrap_samples":B,"seed":seed,
       "CI_low":float(np.quantile(draws,.025)),
       "CI_high":float(np.quantile(draws,.975)),
       "fraction_negative":float(np.mean(np.array(draws)<0))}
def run(paths,out):
    out=Path(out)
    check(not out.exists(),"Refuse to overwrite previous actual 2024 ecological results")
    for key,p in paths.items():
        check(digest(p)==SOURCE_SHA[key],"Frozen input SHA mismatch: "+key)
    canonical=pd.read_csv(paths["canonical"],compression="gzip",low_memory=False,
      dtype={"photo_id":str,"observation_id":str,"inat_taxon_id":str},
      usecols=["photo_id","observation_id","inat_taxon_id","source_subarm","observed_on",
               "latitude","longitude","strict9","strict9_white",
               "wc2_bio05_10m","wc2_bio12_10m","wc2_elev_2p5m"])
    wcvp=pd.read_csv(paths["wcvp"],compression="gzip",dtype=str,
      usecols=["inat_taxon_id","wcvp_match_state","wcvp_accepted_plant_name_id"])
    geo=pd.read_csv(paths["photo_landmass"],low_memory=False,dtype={
       "target_id":str,"cell_id":str,"land_id":str})
    env=pd.read_csv(paths["sameyear_2024_environment"],low_memory=False,
      dtype={"target_id":str})
    original,data,b,m,per_species,ambiguous=make_frames(canonical,wcvp,geo,env)
    check(len(original)==105441 and int(original.strict9.sum())==42522,"Original historical_main target source changed")
    check(len(b)==2091 and len(m)==1617,"Original 2024 photographed colour matched support changed")
    fits=[]
    models=[
      ("VIIRS_all_photo_eligible",b,
        "white_pp ~ log_light + "+COVARS+" + C(region5)","region5"),
      ("VIIRS_pop2024_region_adjusted",m,
        "white_pp ~ log_light + log_pop + "+COVARS+" + C(region5)","region5"),
      ("VIIRS_pop2024_island_indicator",m,
        "white_pp ~ is_island + log_light + log_pop + "+COVARS+" + C(region5)","region5"),
      ("VIIRS_pop2024_island_light_interaction",m,
        "white_pp ~ log_light * is_island + log_pop + "+COVARS+" + C(region5)","region5")]
    for name,d,form,cluster in models:
        for term,stat in cluster_fit(d,form,cluster).items():
            fits.append({"model":name,"term":term,"cells":len(d),
                "clusters":int(d[cluster].nunique()),**stat})
    # Island-fixed comparisons require at least TWO local exposure-distinct
    # source cells. No single-island photo is used to estimate intra-island slope.
    island=m.loc[m.land_class.eq("island")].copy()
    supported=island.groupby("land_id").agg(
        n_cells=("target_id","size"),n_distinct_light=("viirs_median","nunique"))
    ids=set(supported.loc[supported.n_cells.ge(2)&supported.n_distinct_light.ge(2)].index)
    island=island.loc[island.land_id.isin(ids)].copy()
    check(len(island)==142 and island.land_id.nunique()==17,
          "2024 true island fixed-effect source support changed")
    main_form="white_pp ~ log_light + log_pop + "+COVARS+" + C(land_id)"
    core=cluster_fit(island,main_form,"land_id")["log_light"]
    fits.append({"model":"within_same_island_pop2024_islandFE","term":"log_light",
                 "cells":len(island),"clusters":17,**core})
    check(abs(core["beta"]-(-4.884245016310668))<1e-5,
          "Independent original source island coefficient drift")
    loo=[]
    for id in sorted(ids):
        sample=island.loc[island.land_id.ne(id)]
        r=cluster_fit(sample,main_form,"land_id")["log_light"]
        loo.append({"omitted_landmass_id":str(id),"omitted_cell_count":int(
           island.land_id.eq(id).sum()),"cells":len(sample),
           "islands":sample.land_id.nunique(),**r})
    loo=pd.DataFrame(loo)
    boot=island_bootstrap(island)
    no_pop=cluster_fit(island,
       "white_pp ~ log_light + "+COVARS+" + C(land_id)","land_id")["log_light"]
    island_viirs=b.loc[b.land_class.eq("island")].copy()
    info=island_viirs.groupby("land_id").agg(
        n=("target_id","size"),n_light=("viirs_median","nunique"))
    retain=set(info.loc[info.n.ge(2)&info.n_light.ge(2)].index)
    island_viirs=island_viirs.loc[island_viirs.land_id.isin(retain)]
    no_pop_full=cluster_fit(island_viirs,
        "white_pp ~ log_light + "+COVARS+" + C(land_id)","land_id")["log_light"]
    # Accepted species within the same island at multiple source cells:
    rep=per_species.merge(island[["target_id","land_id","log_light","log_pop","bio5",
                                  "log_rain","elevation","log_attempts"]],on="target_id",
                          how="inner",validate="many_to_one")
    rep["species_island"]=rep.land_id.astype(str)+"|"+rep.wcvp_accepted_plant_name_id.astype(str)
    support=rep.groupby("species_island").agg(
        n=("target_id","nunique"),light_span=("log_light",lambda x:x.max()-x.min()))
    repeated=rep.loc[rep.species_island.isin(
       set(support.loc[support.n.ge(2)&support.light_span.gt(.1)].index))].copy()
    repeated["species_photo_white_pp"]=100*repeated.strict9_white
    # With only 7 independent source islands, the species FE p value is
    # diagnostic only and must NOT be elevated to confirmatory evidence.
    species_fit=cluster_fit(repeated,
        "species_photo_white_pp ~ log_light + log_pop + "+COVARS+" + C(species_island)","land_id")["log_light"]
    groups=[]
    for name,t in repeated.groupby("species_island"):
        if t.log_light.max()-t.log_light.min()<.25:continue
        beta=np.polyfit(t.log_light,t.species_photo_white_pp,1)[0]
        groups.append(float(beta))
    population=b.groupby("pop_eligible").agg(
        n=("target_id","size"),median_viirs=("viirs_median","median"),
        mean_white_pp=("white_pp","mean"),
        island_cells=("is_island","sum")).reset_index()
    descriptive=b.groupby("land_class").agg(
        n_cells=("target_id","size"),
        mean_white_component_pp=("white_pp","mean"),
        median_viirs=("viirs_median","median"),
        accepted_species_cell_sum=("accepted_species_n","sum")).reset_index()
    counts=island.groupby("land_id").agg(
        n_cells=("target_id","size"),source_photo_attempts=("photo_attempts","sum"),
        accepted_species_cell_sum=("accepted_species_n","sum"),
        mean_white_pp=("white_pp","mean"),
        median_viirs=("viirs_median","median")).reset_index()
    out.mkdir(parents=True)
    pd.DataFrame(fits).to_csv(out/"2024_sameyear_original_white_vs_nightlight_models.csv",index=False)
    population.to_csv(out/"2024_population_missingness_source_exposure.csv",index=False)
    descriptive.to_csv(out/"2024_island_vs_mainland_photo_white_descriptive.csv",index=False)
    loo.to_csv(out/"2024_same_island_leave_one_out.csv",index=False)
    counts.to_csv(out/"2024_same_island_cell_counts.csv",index=False)
    receipt={
      "schema":"MOIS_2024_original_photo_white_vs_sameyear_VIIRS_G2_source_20261107",
      "source_SHA256":SOURCE_SHA,
      "photographic_source_arm":"historical_main",
      "original_2024_photo_rows_in_unambiguous_landmass_cells":len(original),
      "original_existing_strict9_success_in_unambiguous_cells":int(original.strict9.sum()),
      "ambiguous_GSHHG_cell_landmass_rows_omitted":ambiguous,
      "photo_eligible_VIIRS_cells":len(b),
      "photo_eligible_with_sameyear_2024_population_cells":len(m),
      "within_island_num_cells":len(island),
      "within_island_num_independent_islands":int(island.land_id.nunique()),
      "original_within_island_sameyear_pop_adj_loglight_pp":core,
      "island_without_population_adjustment_on_population_valid_rows":no_pop,
      "island_without_population_adjustment_including_population_missing_cells":no_pop_full,
      "leave_one_island_out_count":len(loo),
      "leave_one_island_out_all_negative":bool(loo.beta.lt(0).all()),
      "leave_one_island_out_significant_nominal_95_negative":int(loo.ci_high.lt(0).sum()),
      "leave_one_island_out_including_zero":int(((loo.ci_low<=0)&(loo.ci_high>=0)).sum()),
      "island_replacement_cluster_bootstrap":boot,
      "repeated_same_accepted_species_island_groups":int(repeated.species_island.nunique()),
      "repeated_same_species_islands":int(repeated.land_id.nunique()),
      "repeated_same_species_cell_rows":len(repeated),
      "repeated_same_species_fixed_effect_exploratory_nominal":species_fit,
      "repeated_same_species_sign_groups_n":len(groups),
      "repeated_same_species_negative_signs":int(np.sum(np.array(groups)<0)),
      "island_mainland_cross_section_light_interaction":next(
          x for x in fits if x["model"]=="VIIRS_pop2024_island_light_interaction"
          and x["term"]=="log_light:is_island"),
      "island_light_interaction_not_a_causal_or_botanical_species_rate_effect":True,
      "multiple_model_search_exploratory":True,
      "new_photos_downloaded_or_ROI_inferred":0,
      "unknown_photo_colours_NOT_imputed":True,
      "moth_effectiveness_and_botanical_white_species_fraction_NOT_measured":True,
      "no_ecological_causal_or_extirpation_claim":True}
    receipt["output_SHA256"]={p.name:digest(p) for p in sorted(out.iterdir()) if p.suffix==".csv"}
    (out/"source_receipt.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print("MOIS_2024_ALAN_WHITE_SPATIAL_MODELS_CONFIRMED",json.dumps({
       "photo_cells":len(b),"same_year_population_cells":len(m),
       "islands":len(ids),"island_cells":len(island),
       "negative_island_slope_pp":core["beta"],
       "island_mainland_light_interaction":receipt["island_mainland_cross_section_light_interaction"],
       "LOO_negative":int(loo.beta.lt(0).sum()),
       "LOO_significant_negative":int(loo.ci_high.lt(0).sum()),
       "bootstrap_CI":[boot["CI_low"],boot["CI_high"]]},sort_keys=True),flush=True)
    return receipt
def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in SOURCE_SHA:p.add_argument("--"+name.replace("_","-"),dest=name,type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    run({k:getattr(a,k) for k in SOURCE_SHA},a.out)
if __name__=="__main__":
    main()
