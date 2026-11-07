#!/usr/bin/env python3
"""GloBI evidence exclusivity × existing 2020 VIIRS: supported photo-colour contrast only.

Moth recorded as a visitor != nocturnal effective pollination. Butterfly
record != exclusion of moth pollination. No source-GloBI hit != diurnal.
The outcome is existing continuous ROI9 white mass, NOT white-flowered species
presence, species extinction, native flora or a three-way causal difference.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

PINS={
 "canonical":"f033661d76e07898a5c13199a55c5d0ce6e55353daa48b0058c456979e90644c",
 "names":"dc5bb0e3686fde0cfd327e4336245d0ceae3c7657b3bd12841bb814cc94d1eb2",
 "accepted":"9cd6aa7ca29f00c64283bcc08dcb778e212b644b837e64070051339b11d8b7b3",
 "direct":"7e86170218814abe2845da9a0282aecfdc930819b9d0c439bae8df4a0679c577",
 "globi":"6116acab80f2277af7ec109715e5815bf6c4461767ae61580da8d842c9e230e4",
 "viirs":"3562db1e5783e77d3bb408bd33817b4d7cccf31ad4d4ad42586469708fd407e3",
}
CATEGORIES=("moth_only","butterfly_only","both")
WHITE={"white":"white","yellow_orange":"nonwhite","red_pink":"nonwhite",
       "blue_purple":"nonwhite","green_brown_inconspicuous":"nonwhite"}

def check(ok,msg):
    if not bool(ok):raise ValueError("MOIS_ALAN_GLOBI_EVIDENCE: "+msg)

def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(4<<20),b""):h.update(block)
    return h.hexdigest()

def evidence_groups(globi):
    check({"plant_name","guild"} <= set(globi),"GloBI source schema changed")
    check(globi.guild.isin(("moth_family_visit","butterfly_family_visit")).all(),
          "Unspecified GloBI channel")
    names=globi[["plant_name","guild"]].drop_duplicates()
    # Group plants by GloBI evidence, NOT by missing evidence.
    flags=names.assign(evidence=True).pivot(index="plant_name",columns="guild",
                                          values="evidence").fillna(False)
    m=flags.get("moth_family_visit",pd.Series(False,index=flags.index)).astype(bool)
    b=flags.get("butterfly_family_visit",pd.Series(False,index=flags.index)).astype(bool)
    flags["evidence_group"]=np.select([m&~b,b&~m,m&b],CATEGORIES,default="none")
    return flags["evidence_group"]

def cell_colour_support(species_rows):
    x=species_rows.copy()
    check({"cell_id","accepted_species","evidence_group","strict9_white","direct_colour"} <= set(x),
          "Required same-cell accepted-species evidence missing")
    check(not x.duplicated(["cell_id","accepted_species"]).any(),
          "Repeated same accepted species; collapse photos by species before analysis")
    check(x.strict9_white.between(0,1).all(),"Invalid original ROI9 white mass")
    exclusive=x.loc[x.evidence_group.isin(("moth_only","butterfly_only"))].copy()
    groups=exclusive.groupby(["cell_id","evidence_group"],observed=True).agg(
        n_species=("accepted_species","nunique"),mean_white=("strict9_white","mean")
    ).reset_index()
    cnt=groups.pivot(index="cell_id",columns="evidence_group",values="n_species").fillna(0)
    mean=groups.pivot(index="cell_id",columns="evidence_group",values="mean_white")
    for label in ("moth_only","butterfly_only"):
        if label not in cnt:cnt[label]=0
        if label not in mean:mean[label]=np.nan
    eligible=(cnt.moth_only.ge(1)&cnt.butterfly_only.ge(1))
    matched=pd.DataFrame({"cell_id":cnt.index[eligible]})
    matched["moth_species"]=cnt.loc[eligible,"moth_only"].astype(int).to_numpy()
    matched["butterfly_species"]=cnt.loc[eligible,"butterfly_only"].astype(int).to_numpy()
    matched["white_moth"]=mean.loc[eligible,"moth_only"].to_numpy()
    matched["white_butterfly"]=mean.loc[eligible,"butterfly_only"].to_numpy()
    matched["moth_minus_butterfly_white_pp"]=100*(matched.white_moth-matched.white_butterfly)
    d=exclusive.loc[exclusive.direct_colour.eq("white")]
    tw=d.groupby(["cell_id","evidence_group"]).accepted_species.nunique().unstack(fill_value=0)
    for label in ("moth_only","butterfly_only"):
        if label not in tw:tw[label]=0
    direct_both=int(((tw.moth_only>=1)&(tw.butterfly_only>=1)).sum())
    known=exclusive.loc[exclusive.direct_colour.isin(("white","nonwhite"))]
    f=known.groupby(["cell_id","evidence_group","direct_colour"]).accepted_species.nunique().unstack(
        ["evidence_group","direct_colour"],fill_value=0)
    def n(group,state):return f.get((group,state),pd.Series(0,index=f.index))
    full=(n("moth_only","white").ge(1)&n("butterfly_only","white").ge(1)&
          n("moth_only","nonwhite").ge(1)&n("butterfly_only","nonwhite").ge(1))
    return matched,{
        "distinct_cells_with_both_exclusive_evidence_groups_ROI9":len(matched),
        "cells_with_at_least_2_distinct_species_in_each_group":int(
            (matched.moth_species.ge(2)&matched.butterfly_species.ge(2)).sum()),
        "cells_with_at_least_3_distinct_species_in_each_group":int(
            (matched.moth_species.ge(3)&matched.butterfly_species.ge(3)).sum()),
        "cells_with_direct_species_white_both_evidence_groups":direct_both,
        "cells_with_direct_white_nonwhite_by_both_evidence_groups":int(full.sum())
    }

def local_models(matched,exposure):
    x=matched.merge(exposure,on="cell_id",how="left",validate="one_to_one")
    check(x.cell_id.is_unique and x.source_complete.all(),"Missing source-locked VIIRS")
    x["night"]=np.log1p(x.viirs_quarterly_median.clip(lower=0))
    x["population"]=np.log1p(x.population_density.clip(lower=0))
    x["rain"]=np.log1p(x.bio12.clip(lower=0))
    x["photo_effort"]=np.log1p(x.n_attempted)
    specs=[
      ("unadjusted","moth_minus_butterfly_white_pp ~ night"),
      ("adjusted_without_geographic_FE","moth_minus_butterfly_white_pp ~ night+population+bio5+rain+elev+photo_effort"),
      ("adjusted_with_5deg_FE","moth_minus_butterfly_white_pp ~ night+population+bio5+rain+elev+photo_effort+C(block_id)")
    ]
    fits=[]
    check(len(x)>50 and x.block_id.nunique()>25,"Insufficient source support")
    for name,formula in specs:
        fit=smf.ols(formula,data=x).fit(cov_type="cluster",
                    cov_kwds={"groups":x.block_id},use_t=True)
        interval=fit.conf_int().loc["night"]
        fits.append({"model":name,
                     "coefficient_pp_per_log1p_VIIRS":float(fit.params["night"]),
                     "nominal_95ci_low_pp":float(interval.iloc[0]),
                     "nominal_95ci_high_pp":float(interval.iloc[1]),
                     "nominal_cluster_p":float(fit.pvalues["night"]),
                     "number_cells":len(x),"blocks":int(x.block_id.nunique()),
                     "p_value_confirmatory":False})
    return x,fits

def run(paths,out):
    out=Path(out)
    check(not out.exists(),"Do not overwrite prior source audit")
    for k,p in paths.items():check(digest(p)==PINS[k],"Frozen source SHA differs: "+k)
    a=pd.read_csv(paths["canonical"],compression="gzip",
        usecols=["photo_id","observation_id","inat_taxon_id","source_subarm",
                 "strict9","strict9_white","latitude","longitude","observed_on"],
        dtype={"photo_id":str,"observation_id":str,"inat_taxon_id":str},low_memory=False)
    check(len(a)==1032882 and a.photo_id.is_unique and a.observation_id.is_unique,
          "Canonical photo/observation identity drift")
    a=a.loc[a.source_subarm.eq("historical_main")&
          a.strict9&a.latitude.between(-60,60)&a.longitude.between(-180,180)&
          a.observed_on.astype(str).str.startswith("2020")].copy()
    lookup=pd.read_csv(paths["names"],compression="gzip",
        usecols=["inat_taxon_id","wcvp_match_state","wcvp_accepted_plant_name_id"],dtype=str)
    lookup=lookup.loc[lookup.wcvp_match_state.eq("unique_exact_name_to_accepted_id"),
           ["inat_taxon_id","wcvp_accepted_plant_name_id"]]
    check(lookup.inat_taxon_id.is_unique,"Ambiguous WCVP original taxon identity")
    accepted=pd.read_csv(paths["accepted"],compression="gzip",
        usecols=["wcvp_accepted_plant_name_id","wcvp_accepted_name"],dtype=str).drop_duplicates()
    check(accepted.groupby("wcvp_accepted_plant_name_id").wcvp_accepted_name.nunique().le(1).all(),
          "Conflicting accepted species names")
    lookup=lookup.merge(accepted.drop_duplicates("wcvp_accepted_plant_name_id"),
                       on="wcvp_accepted_plant_name_id",how="left",validate="many_to_one")
    source=pd.read_csv(paths["globi"],compression="gzip")
    guild=evidence_groups(source)
    direct=pd.read_csv(paths["direct"],compression="gzip",dtype=str,
          usecols=["accepted_species","trait_name","normalized_value",
                   "quality","resolution_status","classification"])
    direct=direct.loc[direct.trait_name.eq("flower_primary_color")&
        direct.quality.isin(("high","medium"))&
        direct.resolution_status.eq("resolved")&
        ~direct.classification.isin(("cultivar_contamination","partial_self_compatibility")),
        ["accepted_species","normalized_value"]].copy()
    check(direct.accepted_species.is_unique,"Direct colour source duplicated")
    direct["direct_colour"]=direct.normalized_value.map(WHITE)
    exposure=pd.read_csv(paths["viirs"],dtype={"cell_id":str})
    check(len(exposure)==1059 and int(exposure.source_complete.sum())==1044,
          "Prior 2020 population and VIIRS source changed")
    exposure=exposure.loc[exposure.source_complete].copy()
    a["cell_id"]=(np.floor((a.longitude+180)/.05).astype(int).astype(str)+"_"+
                  np.floor((a.latitude+90)/.05).astype(int).astype(str))
    a=a.loc[a.cell_id.isin(set(exposure.cell_id))]
    a=a.merge(lookup[["inat_taxon_id","wcvp_accepted_name"]],on="inat_taxon_id",
              how="left",validate="many_to_one")
    a=a.merge(direct[["accepted_species","direct_colour"]],
              left_on="wcvp_accepted_name",right_on="accepted_species",
              how="left",validate="many_to_one")
    a["evidence_group"]=a.wcvp_accepted_name.map(guild)
    a=a.loc[a.evidence_group.isin(("moth_only","butterfly_only","both"))]
    by_species=a.groupby(["cell_id","wcvp_accepted_name","evidence_group"],as_index=False).agg(
        strict9_white=("strict9_white","mean"),
        direct_colour=("direct_colour","first"))
    by_species=by_species.rename(columns={"wcvp_accepted_name":"accepted_species"})
    check(by_species.accepted_species.notna().all(),"Missing matched plant name")
    matched,support=cell_colour_support(by_species)
    x,models=local_models(matched,exposure)
    # This observational contrast is NOT the user's documented white-species fraction
    # or a causal triple interaction: only 1 cell has both direct white classes.
    out.mkdir(parents=True)
    aggregate=pd.DataFrame(models)
    aggregate.to_csv(out/"ROI9_guild_evidence_spatial_exploratory_models.csv",index=False)
    receipt={
       "schema":"mois_ALAN_GloBI_exclusive_moth_butterfly_2020_exploratory_v1",
       "source_sha256":PINS,"original_VIIRS_2020_source_cells":1044,
       "strict9_success_photo_records_guild_evidence":len(a),
       "accepted_species_cell_evidence_rows":len(by_species),
       "genuine_exclusive_evidence_moth_is_not_nocturnal_truth":True,
       "butterfly_only_means_no_recorded_moth_not_proven_moth_absence":True,
       "GloBI_source_unversioned_and_some_pages_truncated":True,
       "moth_and_butterfly_both_evidence_excluded_from_exclusive_comparison":True,
       **support,
       "descriptive_equal_cell_group_white_ROI_difference_pp":float(x.moth_minus_butterfly_white_pp.mean()),
       "candidate_nightlight_regression_on_photographic_white_contrast_not_DDD":True,
       "third_difference_model_fit":False,
       "original_white_species_loss_measured":False,
       "strong_independent_plant_moth_effectiveness_evidence":False,
       "photographic_white_is_not_direct_botanical_white":True,
       "old_missing_photo_outcomes_unchanged":True,
       "all_fits_exploratory_and_unadjusted_multiple_comparisons":True,
       "new_original_image_downloads_or_inferences":0,
       "MNAR_sensitivity_scenarios_run":False,
       "results":models,
       "output_sha256":digest(out/"ROI9_guild_evidence_spatial_exploratory_models.csv")
    }
    (out/"receipt.json").write_text(json.dumps(receipt,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("MOIS_GLOBI_EXCLUSIVE_2020_SOURCE_GATES",
          json.dumps({k:receipt[k] for k in
            ("distinct_cells_with_both_exclusive_evidence_groups_ROI9",
             "cells_with_at_least_2_distinct_species_in_each_group",
             "cells_with_direct_species_white_both_evidence_groups",
             "cells_with_direct_white_nonwhite_by_both_evidence_groups")},sort_keys=True),
          flush=True)
    print("MOIS_GLOBI_ROI9_EXPLORATORY_ASSOCIATION",
          json.dumps(models,sort_keys=True),flush=True)
    return receipt

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in PINS:p.add_argument("--"+k,required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    run({k:getattr(a,k) for k in PINS},a.output)

if __name__=="__main__":main()
