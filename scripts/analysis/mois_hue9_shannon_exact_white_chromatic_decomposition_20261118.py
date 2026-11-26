#!/usr/bin/env python3
"""MOIS 9-hue global community entropy: exact Shannon chain-rule partition.

For each true 100-km GBIF species assemblage with 9 mutually exclusive
photographic hue classes, the Shannon identity is

  H9 = H(white vs nonwhite) + P(nonwhite)*H8(nonwhite colour classes).

Apply the SAME 199 whole-catalogue genus-by-occupancy-breadth species-label
permutations as the signed genus+breadth null. Contrast each observed
component with its permuted expectation, so total observed-minus-expected
H9 is exactly the sum of binary and conditional-chromatic components.

Does NOT use a different nonwhite-only null: this is a strict additive
decomposition of the EXISTING nine-hue null, not a causal mediator analysis
or independently confirmed botanically correct flower-petal colours.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd

from scripts.analysis.mois_global_hue9_community_assembly_genus_null_20261118 import (
    signed_original,baseline,counts_by_cell,empirical_entropy,SEED,
    GBIF_SHA,ENV_SHA,check,sha
)
from scripts.analysis.mois_nonwhite8_assembly_white_fixed_null_20261118 import (
    load_strict_nine_null,NULL_SHA
)

N_PERM=199

def entropy_parts(counts):
    n=counts.sum(axis=1).astype(float)
    nw=counts[:,0].astype(float)
    nk=n-nw
    w=np.divide(nw,np.maximum(n,1))
    q=np.divide(nk,np.maximum(n,1))
    binary=-(np.where(w>0,w*np.log(np.maximum(w,1e-12)),0)+
             np.where(q>0,q*np.log(np.maximum(q,1e-12)),0))
    h8=empirical_entropy(counts[:,1:])
    chromatic=q*h8
    h9=empirical_entropy(counts)
    check(np.allclose(binary+chromatic,h9,atol=1e-12,rtol=0),
          "Exact Shannon conditional-decomposition broken")
    return np.column_stack((binary,chromatic,h9))

def analyze(x,source):
    observed_counts=counts_by_cell(x,x["hue"])
    observed=entropy_parts(observed_counts)
    check(np.allclose(observed[:,2],source.observed_hue9_Shannon,
                      atol=1e-12,rtol=0),
          "Actual GBIF 9-hue Shannon differs from previous signed report")
    masks={
      "richness_ge30":x["cell_n"]>=30,
      "richness_ge100":x["cell_n"]>=100
    }
    rng=np.random.default_rng(SEED)
    accum=np.zeros_like(observed)
    draws={name:np.zeros((N_PERM,3)) for name in masks}
    for i in range(N_PERM):
        p=x["hue"].copy()
        for g in x["variable_groups"]:
            p[g]=x["hue"][g][rng.permutation(len(g))]
        counts=counts_by_cell(x,p)
        terms=entropy_parts(counts)
        accum+=terms
        for name,m in masks.items():
            draws[name][i,:]=terms[m].mean(axis=0)
        if (i+1)%50==0:
            print("MOIS_NINE_HUE_CHAIN_RULE_SIGNED_PERM",i+1,flush=True)
    expected=accum/N_PERM
    require_same=np.allclose(expected[:,2],
       source.expected_hue9_Shannon_given_genus.to_numpy(float),
       atol=1e-10,rtol=0)
    check(require_same,
          "The original 199 signed source permutations did not reproduce")
    results={}
    for name,m in masks.items():
        obs=observed[m].mean(axis=0)
        null=draws[name].mean(axis=0)
        delta=obs-null
        check(abs(delta[0]+delta[1]-delta[2])<1e-10,
              "Whole-community deficit does not equal components")
        r={}
        for k,j in (("white_vs_nonwhite",0),
                    ("nonwhite_weighted_chromatic8",1),
                    ("complete_nine_hue",2)):
            vals=draws[name][:,j]
            p_low=(1+np.count_nonzero(vals<=obs[j]+1e-12))/(N_PERM+1)
            p_two=(1+np.count_nonzero(
                abs(vals-null[j])>=abs(delta[j])-1e-12))/(N_PERM+1)
            r[k]={
              "observed_mean_Shannon_component":float(obs[j]),
              "null_mean_Shannon_component":float(null[j]),
              "observed_minus_null_Shannon":float(delta[j]),
              "null_distribution_95pct":[float(v) for v in
                  np.percentile(vals,[2.5,97.5])],
              "one_sided_more_clustering_p":float(p_low),
              "two_sided_permutation_p":float(p_two)}
        r["n_cells"]=int(m.sum())
        r["white_nonwhite_deficit_fraction_of_total"]=(
           float(delta[0]/delta[2]) if abs(delta[2])>1e-9 else None)
        r["chain_rule_total_residual"]=float(delta[0]+delta[1]-delta[2])
        results[name]=r
    coords=np.asarray(x["cell_key"].tolist(),dtype=np.int64)
    check(coords.shape==(x["ngrids"],2),"2D global GBIF coordinates changed")
    cols={}
    for k,j in (("white_nonwhite_binary",0),
                ("nonwhite_weighted_colour8",1),
                ("total_ninecolour",2)):
        cols[k+"_observed"]=observed[:,j]
        cols[k+"_null"]=expected[:,j]
        cols[k+"_anomaly"]=observed[:,j]-expected[:,j]
    gr=pd.DataFrame({"grid_x":coords[:,0],
                     "grid_y":coords[:,1],
                     "n_recorded_source_species":x["cell_n"],
                     "n_original_white_species":observed_counts[:,0],
                     "original_white_fraction":np.divide(
                         observed_counts[:,0],x["cell_n"],
                         out=np.zeros(x["ngrids"]),where=x["cell_n"]>0),
                     **cols})
    check(np.allclose(gr.white_nonwhite_binary_anomaly+
                      gr.nonwhite_weighted_colour8_anomaly,
                      gr.total_ninecolour_anomaly,atol=1e-11,rtol=0),
          "Observed-minus-permutation-mean entropy partition failed")
    return results,gr

def main(args):
    out=Path(args.out)
    check(not out.exists(),"Refuse overwriting numerical source analysis")
    pair,spp=signed_original(args.GBIF,args.ENV)
    x=baseline(pair,spp,null_model="genus_breadth")
    nine=load_strict_nine_null(args.NINE)
    coords=np.asarray(x["cell_key"].tolist(),dtype=np.int64)
    check(np.array_equal(coords,nine[["grid_x","grid_y"]].to_numpy()),
          "Same global grid order and accepted species source not reproduced")
    results,grid=analyze(x,nine)
    out.mkdir(parents=True)
    f=out/"same_cell_exact_white_nonwhite_chromatic8_entropy_decomposition.csv.gz"
    grid.to_csv(f,index=False,compression="gzip")
    report={
       "analysis":"Shannon chain rule exact components using original genus-by-breadth nine-colour label null",
       "n_accepted_species":37104,
       "n_source_species_grid_pairs":5493931,
       "n_spatial_100km_cells":18577,
       "n_global_shared_permutations":N_PERM,
       "same_seed_and_same_shuffled_species_groups_as_previous_signed_run":True,
       "previous_9hue_expected_entropy_exactly_reproduced":True,
       "components_add_up_to_total_to_machine_precision":True,
       "source_zip_sha256":{"GBIF":GBIF_SHA,"ENV":ENV_SHA,"original_genus_breadth9":NULL_SHA},
       "results":results,
       "not_field_verified_flower_petal_colours":True,
       "descriptive_identity_not_causal_white_mediation":True
    }
    summary=out/"ninecolour_shannon_component_null_summary.json"
    summary.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    report["output_sha256"]={p.name:sha(p) for p in (f,summary)}
    (out/"receipt.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print("MOIS_WHITE_COLOURED_SHANNON_CHAINRULE_REAL_SUCCESS",
          json.dumps(results,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    for name in ("GBIF","ENV","NINE","out"):p.add_argument("--"+name,required=True)
    main(p.parse_args())
