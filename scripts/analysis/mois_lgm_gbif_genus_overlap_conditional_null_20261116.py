#!/usr/bin/env python3
"""True GBIF/WCVP 100km cell genus overlap conditional on native/introduction counts.

Exact hypergeometric expectation under random assignment of WCVP status to
recorded species, conditional on each cell's native/introduced species totals
and species counts per genus. No photo hue involved in the taxonomic outcome.
No causal conclusion about immigration, climate, or true glacial refugia.
"""
from __future__ import annotations
import argparse,io,json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.special import gammaln

from scripts.analysis.mois_lgm_gbif_native_intro_genus_turnover_20261116 import (
    signed_context,genus_overlap,model
)
from scripts.analysis.mois_lgm_wcvp_native_regional_block_20261116 import (
    assign_cells_broad_region,get_wcvp_source
)
from scripts.analysis.mois_lgm_native_intro_common_genus_20261116 import (
    regional_species_cell_rows
)
from scripts.analysis.mois_lgm_native_introduced_matched_cells_20261116 import (
    same_cell_cohort
)
from scripts.analysis.mois_lgm_climate_history_species_20261116 import (
    require,sha
)

TARGET={"strict":3106,"moderate":3602}

def p_all_members_same_status(popsize, status_count, genus_size):
    """Choose(status_count,genus_size) / Choose(popsize,genus_size)."""
    M=np.asarray(popsize,dtype=np.float64)
    K=np.asarray(status_count,dtype=np.float64)
    m=np.asarray(genus_size,dtype=np.float64)
    require((M>=1).all() and (K>=0).all() and (m>=1).all()
            and (K<=M).all() and (m<=M).all(),
            "hypergeometric conditioning counts impossible")
    out=np.zeros(np.shape(M),dtype=np.float64)
    valid=m<=K
    if np.any(valid):
        q=(gammaln(K[valid]+1)-gammaln(m[valid]+1)
           -gammaln(K[valid]-m[valid]+1)
           -gammaln(M[valid]+1)+gammaln(m[valid]+1)
           +gammaln(M[valid]-m[valid]+1))
        out[valid]=np.exp(q)
    return np.clip(out,0,1)

def expected_genus_overlap(source):
    """Expected shared genera given observed genus sizes and each group size."""
    source=source.copy()
    require(not source[["specieskey","grid_x","grid_y"]].duplicated().any(),
            "species counted more than once per GBIF grid")
    source["native"]=(source.status=="native_at_TDWG3").astype("int8")
    require(source.status.isin(["native_at_TDWG3","introduced_at_TDWG3"]).all(),
            "non-WCVP status cannot be randomly relabelled as introduced/native")
    genus=source.groupby(["grid_x","grid_y","genus"],sort=False).agg(
        genus_species=("specieskey","size"),
        native_species=("native","sum")).reset_index()
    totals=genus.groupby(["grid_x","grid_y"],sort=False).agg(
        species_total=("genus_species","sum"),
        native_total=("native_species","sum"),
        genus_union=("genus","size")).reset_index()
    totals["introduced_total"]=totals.species_total-totals.native_total
    genus=genus.merge(totals,on=["grid_x","grid_y"],how="left",
                      validate="many_to_one")
    M=genus.species_total.to_numpy()
    Nn=genus.native_total.to_numpy()
    Ni=genus.introduced_total.to_numpy()
    k=genus.genus_species.to_numpy()
    p_native=p_all_members_same_status(M,Nn,k)
    p_intro=p_all_members_same_status(M,Ni,k)
    p_both=np.clip(1-p_native-p_intro,0,1)
    genus["expected_shared_if_status_random"]=p_both
    genus["observed_shared"]=((genus.native_species>0)
                    &(genus.native_species<genus.genus_species)).astype(int)
    grouped=genus.groupby(["grid_x","grid_y"],sort=False).agg(
        expected_shared_genera=("expected_shared_if_status_random","sum"),
        observed_shared_genera=("observed_shared","sum")).reset_index()
    r=totals.merge(grouped,on=["grid_x","grid_y"],
                   validate="one_to_one")
    r["expected_jaccard_if_status_random"]=(
        r.expected_shared_genera/r.genus_union)
    r["observed_jaccard"]=r.observed_shared_genera/r.genus_union
    r["excess_jaccard_vs_status_random"]=(
        r.observed_jaccard-r.expected_jaccard_if_status_random)
    require(r.expected_jaccard_if_status_random.between(0,1).all()
            and r.observed_jaccard.between(0,1).all()
            and r.excess_jaccard_vs_status_random.between(-1,1).all(),
            "invalid conditional overlap")
    return r

def block_mean(d,col):
    y=d[col].astype(float)
    grp=d.block_500km.astype(str)
    z=sm.OLS(y,np.ones((len(y),1))).fit(
        cov_type="cluster",cov_kwds={"groups":grp},use_t=False)
    a,b=np.asarray(z.conf_int())[0].tolist()
    return {"n_cells":len(y),"n_500km_blocks":int(grp.nunique()),
            "mean":float(y.mean()),"95CI_cluster":[float(a),float(b)],
            "p_vs_zero":float(np.asarray(z.pvalues)[0])}

def run(gbif,env,lgm,native,spatial,common,wcvp,tdwg,out):
    out=Path(out)
    require(not out.exists(),"refuse overwriting actual null output")
    pairs,species,climate,native_s,spatial_s,prior=signed_context(
        gbif,env,lgm,native,spatial,common)
    nm,dist,_=get_wcvp_source(wcvp)
    regions=assign_cells_broad_region(climate,tdwg)
    matched={}; ids=[]
    for mode in TARGET:
        x=same_cell_cohort(native_s,spatial_s,mode,30)
        require(len(x)==TARGET[mode],"signed same GBIF cell universe changed")
        matched[mode]=x
        ids.append(x[["grid_x","grid_y"]])
    ids=pd.concat(ids,ignore_index=True).drop_duplicates()
    source=regional_species_cell_rows(pairs,species,nm,dist,regions,ids)
    null=expected_genus_overlap(source)
    old=genus_overlap(source)
    validate=null.merge(old[["grid_x","grid_y","jaccard_native_introduced_genus",
                          "shared_genus_n"]],on=["grid_x","grid_y"],
                        validate="one_to_one")
    require(np.allclose(validate.observed_jaccard,
                        validate.jaccard_native_introduced_genus,atol=1e-12)
            and np.array_equal(validate.observed_shared_genera.to_numpy(),
                               validate.shared_genus_n.to_numpy()),
            "exact genus set Jaccard does not match independent signed engine")
    result={
        "GBIF_current_source_taxa":37104,
        "source_species_x_GBIF_grid":5493931,
        "genus_overlap_and_null_do_not_use_any_photo_hue":True,
        "conditional_random_status_labels_are_not_historical_shuffling":True,
        "WCVP_native_is_broad_TDWG3_not_verified_local_wild":True,
        "does_not_observe_LGM_plant_ranges_or_refugia":True,
        "regional_models":{}}
    tables=[]
    for mode,d in matched.items():
        z=d.merge(null,on=["grid_x","grid_y"],how="left",validate="one_to_one")
        require(len(z)==TARGET[mode]
                and z.excess_jaccard_vs_status_random.notna().all(),
                "conditional null missing same GBIF cells")
        results={
            "n_matched_cells":len(z),
            "mean_observed_jaccard":float(z.observed_jaccard.mean()),
            "mean_conditional_expected_jaccard":float(
                z.expected_jaccard_if_status_random.mean()),
            "fraction_cells_below_conditional_null":float(
                (z.excess_jaccard_vs_status_random<0).mean()),
            "spatial_mean_excess":block_mean(z,"excess_jaccard_vs_status_random"),
            "models":{}}
        for response in ("observed_jaccard",
                         "expected_jaccard_if_status_random",
                         "excess_jaccard_vs_status_random",
                         "native_minus_intro_white"):
            results["models"][response]={
                "500km":model(z,response,"500"),
                "1000km":model(z,response,"1000")}
        # Test whether the residual genus nonrandomness, rather than raw
        # Jaccard overlap, changes the original photography-based association.
        # This is a sensitivity control, NOT causal mediation.
        d=z.rename(columns={
            "genus_turnover_one_minus_jaccard":"old_genus_turnover_rename"})
        # The imported regression helper's covariate-name interface is
        # re-used here solely to insert conditional-null EXCESS Jaccard.
        # Its term labelled 'genus_turnover' is NOT raw 1 - Jaccard in
        # this particular sensitivity.
        d["genus_turnover_one_minus_jaccard"]=(
            d.excess_jaccard_vs_status_random)
        results["models"]["white_hue_adjusted_for_excess_genus_overlap"]={
            "500km":model(d,"native_minus_intro_white","500",True),
            "1000km":model(d,"native_minus_intro_white","1000",True)}
        result["regional_models"][mode]=results
        z["tdwg_region_mode"]=mode
        tables.append(z[["grid_x","grid_y","tdwg_region_mode",
                "native_total","introduced_total","genus_union",
                "observed_shared_genera","expected_shared_genera",
                "observed_jaccard","expected_jaccard_if_status_random",
                "excess_jaccard_vs_status_random",
                "threeGCM_delta_temperature_C_median","lat"]].copy())
        e=results["models"]["excess_jaccard_vs_status_random"]["500km"]
        print("MOIS_GBIF_NATIVE_INTRO_GENUS_CONDITIONAL_NULL",
              json.dumps({"mode":mode,"n":len(z),
                    "mean_observed":results["mean_observed_jaccard"],
                    "mean_null":results["mean_conditional_expected_jaccard"],
                    "fraction_observed_below_null":results["fraction_cells_below_conditional_null"],
                    "mean_excess_CI":results["spatial_mean_excess"]["95CI_cluster"],
                    "LGM_warming_beta":e["effects"]["LGM_warming"]["beta_per_SD"],
                    "LGM_warming_p":e["effects"]["LGM_warming"]["p"]},
                    sort_keys=True),flush=True)
    out.mkdir(parents=True)
    f=out/"GBIF100km_native_intro_genus_overlap_random_status_null.parquet"
    pd.concat(tables,ignore_index=True).to_parquet(f,index=False)
    m=out/"models.json"
    m.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    receipt={
     "source_full_GBIF_accepted_species":37104,
     "source_full_GBIF_species_cell_pairs":5493931,
     "sources_sha256":{"GBIF":sha(gbif),"ENV":sha(env),"LGM":sha(lgm),
        "native":sha(native),"spatial":sha(spatial),"common":sha(common),
        "WCVP_names":sha(Path(wcvp)/"mois350399_wcvp_source_taxon_exact_name_resolution.csv.gz"),
        "TDWG3":sha(tdwg)},
     "source_native_vs_introduced_is_present_regional_not_individual":True,
     "no_photo_hue_in_null_taxonomic_outcome":True,
     "no_LGM_refugia_reconstructed":True,
     "model":result,"output_sha256":{f.name:sha(f),m.name:sha(m)}}
    (out/"receipt.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    return receipt

if __name__=="__main__":
    p=argparse.ArgumentParser()
    for n in ("gbif","env","lgm","native","spatial","common","wcvp","tdwg","output"):
        p.add_argument("--"+n,required=True)
    a=p.parse_args()
    run(a.gbif,a.env,a.lgm,a.native,a.spatial,a.common,a.wcvp,a.tdwg,a.output)
