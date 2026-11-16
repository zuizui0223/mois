#!/usr/bin/env python3
"""GBIF + WCVP native/introduction genus overlap, conditioning on FAMILY.

For each present-day GBIF 100-km cell, compare observed overlap of WCVP
native- and introduced-range species' genera with TWO exact expectations:
  A: shuffle native/introduced labels across all included species in the cell
  B: shuffle labels only WITHIN official WCVP family in that cell, holding
     family-specific native/introduced species totals fixed.
In both nulls also keep every genus' recorded species count unchanged.

No photographic hue is used in genus/Jaccard outcomes or these null models.
This does not reconstruct LGM plant occupancy or infer a causal introduction
filter. WCVP TDWG3 is BROAD REGION native status, not wild individual status.
"""
from __future__ import annotations
import argparse, io, json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd

from scripts.analysis.mois_lgm_gbif_genus_overlap_conditional_null_20261116 import (
    p_all_members_same_status, block_mean
)
from scripts.analysis.mois_lgm_gbif_native_intro_genus_turnover_20261116 import (
    signed_context, model
)
from scripts.analysis.mois_lgm_native_intro_common_genus_20261116 import (
    regional_species_cell_rows
)
from scripts.analysis.mois_lgm_native_introduced_matched_cells_20261116 import (
    same_cell_cohort
)
from scripts.analysis.mois_lgm_wcvp_native_regional_block_20261116 import (
    get_wcvp_source, match_gbif_to_wcvp, assign_cells_broad_region
)
from scripts.acquisition.mois_wcvp_TDWG3_native_source_20261019 import (
    wcvp_reader, _member_name
)
from scripts.analysis.mois_lgm_climate_history_species_20261116 import (
    sha, require
)

OFFICIAL_WCVP_SHA="d32ea2b3a85e489b14e83bcc9eae7274532e1d113753f7be290d4b2dfde573fa"
TARGET={"strict":3106,"moderate":3602}

def accepted_wcvp_families(wcvp_zip, want):
    """Read official 2026 version; select exact accepted family by accepted PID."""
    require(sha(wcvp_zip)==OFFICIAL_WCVP_SHA,
            "WCVP original taxonomic archive is not the version used by regional status")
    selected={}
    conflicts=set()
    nrows=0
    with ZipFile(wcvp_zip) as z:
        m=_member_name(z,{"wcvp_names.csv"})
        for row in wcvp_reader(z,m):
            nrows+=1
            accepted_id=row["plant_name_id"].strip()
            if accepted_id not in want:
                continue
            if row["taxon_rank"].strip().lower()!="species" or \
               row["taxon_status"].strip().lower()!="accepted":
                continue
            family=row.get("family","").strip()
            if not family:
                continue
            if accepted_id in selected and selected[accepted_id]!=family:
                conflicts.add(accepted_id)
            else:
                selected[accepted_id]=family
    for x in conflicts:
        selected.pop(x,None)
    require(nrows>100000 and len(selected)>30000
            and not conflicts, "unexpected accepted-ID family coverage or conflicting family")
    return pd.DataFrame({
        "wcvp_accepted_plant_name_id":list(selected),
        "WCVP_family":[selected[k] for k in selected]
    }),{"scanned_WCVP_name_rows":int(nrows),
        "accepted_IDs_with_unique_family":len(selected),
        "family_count":len(set(selected.values())),
        "conflicting_accepted_IDs":len(conflicts),
        "source_SHA256":OFFICIAL_WCVP_SHA}

def genus_null_both_species_counts(p):
    """Conditional exact E[shared genera], species totals and family totals fixed.

    Family-genus identity used as structural key to guard homonymous genus
    strings. Each genus contributes 1 to union independent of shuffled label.
    """
    required={"grid_x","grid_y","specieskey","genus","WCVP_family","status"}
    require(required.issubset(p.columns),"required exact species-by-cell keys unavailable")
    require(not p[["specieskey","grid_x","grid_y"]].duplicated().any(),
            "species repeated within one GBIF grid")
    require(p.WCVP_family.notna().all()
            and p.status.isin(["native_at_TDWG3","introduced_at_TDWG3"]).all(),
            "missing family or non-binary status in null comparison")
    q=p[["grid_x","grid_y","genus","WCVP_family","status"]].copy()
    q["native"]=q.status.eq("native_at_TDWG3").astype("int8")
    unit=["grid_x","grid_y","WCVP_family","genus"]
    gen=q.groupby(unit,sort=False,observed=True).agg(
        genus_species=("native","size"),genus_native=("native","sum")
    ).reset_index()
    fam=q.groupby(["grid_x","grid_y","WCVP_family"],sort=False,observed=True).agg(
        family_species=("native","size"),family_native=("native","sum")
    ).reset_index()
    cell=q.groupby(["grid_x","grid_y"],sort=False).agg(
        species_total=("native","size"),native_total=("native","sum")
    ).reset_index()
    x=gen.merge(fam,on=["grid_x","grid_y","WCVP_family"],validate="many_to_one")
    x=x.merge(cell,on=["grid_x","grid_y"],validate="many_to_one")
    x["family_introduced"]=x.family_species-x.family_native
    x["introduced_total"]=x.species_total-x.native_total
    require((x.genus_native<=x.genus_species).all(),"genus class counts invalid")
    genus_n=x.genus_species.to_numpy(int)
    fam_n=x.family_species.to_numpy(int)
    fam_native=x.family_native.to_numpy(int)
    fam_intro=x.family_introduced.to_numpy(int)
    global_n=x.species_total.to_numpy(int)
    global_native=x.native_total.to_numpy(int)
    global_intro=x.introduced_total.to_numpy(int)
    x["expected_both_within_family"]=1-(
        p_all_members_same_status(fam_n,fam_native,genus_n)
        +p_all_members_same_status(fam_n,fam_intro,genus_n))
    x["expected_both_no_family_condition"]=1-(
        p_all_members_same_status(global_n,global_native,genus_n)
        +p_all_members_same_status(global_n,global_intro,genus_n))
    x["observed_both"]=(
        (x.genus_native>0)&(x.genus_native<x.genus_species)).astype("int8")
    for c in ("expected_both_within_family","expected_both_no_family_condition"):
        x[c]=x[c].clip(0,1)
    summary=x.groupby(["grid_x","grid_y"],sort=False).agg(
        genus_union=("genus","size"),
        observed_shared_genera=("observed_both","sum"),
        family_conditioned_expected_shared_genera=("expected_both_within_family","sum"),
        unstratified_expected_shared_genera=("expected_both_no_family_condition","sum")
    ).reset_index()
    families=fam.groupby(["grid_x","grid_y"],sort=False).agg(
        family_union=("WCVP_family","size"),
        family_shared=("family_native",lambda a: 0)
    ).reset_index().drop(columns=["family_shared"])
    fam["present_on_both"]=(fam.family_native>0)&(
        fam.family_native<fam.family_species)
    share=fam.groupby(["grid_x","grid_y"],sort=False).agg(
        shared_families=("present_on_both","sum"),
        native_families=("family_native",lambda a: int((a>0).sum())),
        introduced_families=("family_native",
                             lambda a: int((a<fam.loc[a.index,"family_species"]).sum()))
    ).reset_index()
    summary=summary.merge(families,on=["grid_x","grid_y"],validate="one_to_one")
    summary=summary.merge(share,on=["grid_x","grid_y"],validate="one_to_one")
    summary=summary.merge(cell,on=["grid_x","grid_y"],validate="one_to_one")
    summary["observed_jaccard"]=summary.observed_shared_genera/summary.genus_union
    summary["family_conditioned_expected_jaccard"]=(
        summary.family_conditioned_expected_shared_genera/summary.genus_union)
    summary["unstratified_expected_jaccard"]=(
        summary.unstratified_expected_shared_genera/summary.genus_union)
    summary["residual_after_family_conditioning"]=(
        summary.observed_jaccard-summary.family_conditioned_expected_jaccard)
    summary["residual_after_unstratified_conditioning"]=(
        summary.observed_jaccard-summary.unstratified_expected_jaccard)
    for c in ("observed_jaccard","family_conditioned_expected_jaccard",
              "unstratified_expected_jaccard"):
        require(summary[c].between(-1e-9,1+1e-9).all(),
                f"invalid conditional Jaccard {c}")
    return summary

def run(gbif,env,lgm,native,spatial,common,wcvp,tdwg,official_wcvp,out):
    out=Path(out)
    require(not out.exists(),"refuse overwriting confirmed source")
    pairs,species,climate,native_s,spatial_s,prior=signed_context(
        gbif,env,lgm,native,spatial,common)
    names,distribution,source_receipt=get_wcvp_source(wcvp)
    m=match_gbif_to_wcvp(species,names)
    ids=m.wcvp_accepted_plant_name_id.dropna().astype(str).unique().tolist()
    family,faudit=accepted_wcvp_families(official_wcvp,set(ids))
    mapframe=m[["key","wcvp_accepted_plant_name_id"]].merge(
        family,on="wcvp_accepted_plant_name_id",how="left",validate="many_to_one")
    require(mapframe.key.is_unique and len(mapframe)==37104,
            "GBIF accepted species-family join lost or duplicated species")
    region=assign_cells_broad_region(climate,tdwg)
    groups={}; id_frames=[]
    for mode in TARGET:
        x=same_cell_cohort(native_s,spatial_s,mode,min_each=30)
        require(len(x)==TARGET[mode],"signed original matched sample changed")
        groups[mode]=x
        id_frames.append(x[["grid_x","grid_y"]])
    ids=pd.concat(id_frames,ignore_index=True).drop_duplicates()
    p=regional_species_cell_rows(
        pairs,species,names,distribution,region,ids)
    p=p.merge(mapframe[["key","WCVP_family"]],left_on="specieskey",
              right_on="key",how="left",validate="many_to_one")
    unknown=int(p.WCVP_family.isna().sum())
    known=p.loc[p.WCVP_family.notna()].copy()
    require(len(known)>500000 and known.WCVP_family.nunique()>100,
            "insufficient original WCVP family coverage")
    family_null=genus_null_both_species_counts(known)
    result={
      "GBIF_accepted_species":37104,
      "GBIF_species_100km_pairs":5493931,
      "null_taxonomic_outcome_uses_no_photo_hue":True,
      "WCVP_accepted_ID_family_same_pinned_release_as_native_range":True,
      "TDWG3_region_native_not_wild_individual":True,
      "historical_refugia_reconstructed":False,
      "WCVP_family_acquisition":faudit,
      "WCVP_native_introduced_species_cells":len(p),
      "pairs_family_unknown_excluded":unknown,
      "family_known_coverage_fraction":len(known)/len(p),
      "regional_models":{}
    }
    tables=[]
    for mode,d in groups.items():
        z=d.merge(family_null,on=["grid_x","grid_y"],how="left",
                  validate="one_to_one")
        # only current matched cells with >=30 accepted-family species each
        z=z.loc[(z.native_total>=30)
                &(z.species_total-z.native_total>=30)
                &z.residual_after_family_conditioning.notna()].copy()
        require(len(z)>2000,"too many familyless species or missing paired cells")
        diff=z.residual_after_family_conditioning
        result["regional_models"][mode]={
            "original_matched_cells":len(d),
            "family_known_cells":len(z),
            "mean_native_families":float(z.native_families.mean()),
            "mean_introduced_families":float(z.introduced_families.mean()),
            "mean_shared_families":float(z.shared_families.mean()),
            "mean_observed_jaccard":float(z.observed_jaccard.mean()),
            "mean_unstratified_expected_jaccard":float(
                z.unstratified_expected_jaccard.mean()),
            "mean_family_conditioned_expected_jaccard":float(
                z.family_conditioned_expected_jaccard.mean()),
            "fraction_cells_below_family_expected":float((diff<0).mean()),
            "mean_excess_vs_family_null_spatial_CI":block_mean(z,
                    "residual_after_family_conditioning"),
            "models":{}}
        for col in ("observed_jaccard",
                    "unstratified_expected_jaccard",
                    "family_conditioned_expected_jaccard",
                    "residual_after_unstratified_conditioning",
                    "residual_after_family_conditioning"):
            result["regional_models"][mode]["models"][col]={
                "500km":model(z,col,"500"),
                "1000km":model(z,col,"1000")}
        z["region_mode"]=mode
        columns=["grid_x","grid_y","region_mode","family_union",
                 "shared_families","native_families","introduced_families",
                 "species_total","native_total","genus_union",
                 "observed_shared_genera","observed_jaccard",
                 "unstratified_expected_jaccard",
                 "family_conditioned_expected_jaccard",
                 "residual_after_unstratified_conditioning",
                 "residual_after_family_conditioning",
                 "threeGCM_delta_temperature_C_median","lat"]
        tables.append(z[columns].copy())
        rr=result["regional_models"][mode]
        print("MOIS_GBIF_FAMILY_CONDITIONAL_GENUS_NULL",
              json.dumps({"mode":mode,"n":len(z),
                   "observed":rr["mean_observed_jaccard"],
                   "unstratified_expected":rr["mean_unstratified_expected_jaccard"],
                   "family_expected":rr["mean_family_conditioned_expected_jaccard"],
                   "fraction_below_family_expected":rr["fraction_cells_below_family_expected"],
                   "family_residual_CI":rr["mean_excess_vs_family_null_spatial_CI"]["95CI_cluster"],
                   "LGM_beta_500km":rr["models"]["residual_after_family_conditioning"]["500km"]["effects"]["LGM_warming"]["beta_per_SD"],
                   "LGM_p_500km":rr["models"]["residual_after_family_conditioning"]["500km"]["effects"]["LGM_warming"]["p"]},
                   sort_keys=True),flush=True)
    out.mkdir(parents=True)
    table=out/"GBIF100km_WCVP_family_conditioned_native_intro_genus_null.parquet"
    pd.concat(tables,ignore_index=True).to_parquet(table,index=False)
    models=out/"models.json"
    models.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    audit={
      "source_SHA256":{"GBIF":sha(gbif),"ENV":sha(env),
         "LGM":sha(lgm),"native":sha(native),"spatial":sha(spatial),
         "common":sha(common),"original_WCVP_family":sha(official_wcvp),
         "WCVP_source_names":sha(Path(wcvp)/"mois350399_wcvp_source_taxon_exact_name_resolution.csv.gz"),
         "TDWG3":sha(tdwg)},
      "family_native_intro_counts_fixed_within_each_present_GBIF_cell":True,
      "genus_species_counts_fixed_within_each_present_GBIF_cell":True,
      "no_flower_colour_used_for_taxonomic_overlap":True,
      "no_historical_species_range_reconstructed":True,
      "results":result,
      "output_SHA256":{table.name:sha(table),models.name:sha(models)}
    }
    (out/"receipt.json").write_text(json.dumps(audit,indent=2,sort_keys=True)+"\n")
    return audit

if __name__=="__main__":
    p=argparse.ArgumentParser()
    for name in ("gbif","env","lgm","native","spatial","common",
                 "wcvp","tdwg","official-wcvp","output"):
        p.add_argument("--"+name,required=True)
    a=p.parse_args()
    run(a.gbif,a.env,a.lgm,a.native,a.spatial,a.common,
        a.wcvp,a.tdwg,a.official_wcvp,a.output)
