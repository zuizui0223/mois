#!/usr/bin/env python3
"""MOIS nine-hue assembly decomposition: does the extra clustering persist among
the EIGHT nonwhite photo hues after fixing every white species label in place?

Source: same signed 37,104 accepted-species photo9, 5,493,931 GBIF
species-cell records and 18,577 100km cells. One globally identical species
label is shuffled per genus × geographic occupancy-breadth bin, but only
among species already labelled nonwhite. Thus actual white label geography
AND all nonwhite species occupancy cells are unchanged in every permutation.

Null is descriptive; it neither validates petal pixels nor establishes floral
pollinator competition, photoreceptor preference, or evolutionary selection.
"""
from __future__ import annotations
import argparse,hashlib,io,json
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd

from scripts.analysis.mois_global_hue9_community_assembly_genus_null_20261118 import (
    signed_original,baseline,empirical_entropy,SEED,GBIF_SHA,ENV_SHA,check,sha
)

NULL_SHA="596166fa2fe8bfb28653c325a3a8f1a1495c0792688a7fdb781014cfedaa88c7"
ROOT="mois_global_9hue_community_assembly_breadth_20261118"
NAME="genuine_global_100km_photo9_diversity_vs_genus_null.csv.gz"
N_PERM=199

def load_strict_nine_null(zippath):
    check(sha(zippath)==NULL_SHA,"Previous strict 9-hue source ZIP digest changed")
    with ZipFile(zippath) as z:
        r=json.loads(z.read(ROOT+"/receipt.json"))
        b=z.read(ROOT+"/"+NAME)
    check(r["null_model"]=="genus_breadth"
          and r["source_species"]==37104
          and r["genus_preserving_global_species_hue_label_permutations"]==199
          and sha(b)==r["output_sha256"][NAME],
          "Original strict nine-colour source null SHA/metadata drift")
    d=pd.read_csv(io.BytesIO(b),compression="gzip")
    check(len(d)==18577 and not d[["grid_x","grid_y"]].duplicated().any(),
          "Previous strict nine-hue grid changed")
    return d

def count_chromatic(x,labels,chromatic_cell,chromatic_species):
    hue=labels[chromatic_species]-1
    check(hue.min()>=0 and hue.max()<8,
          "White entered fixed chromatic-species set in null permutation")
    ids=chromatic_cell.astype(np.int64)*8+hue
    return np.bincount(ids,minlength=x["ngrids"]*8).reshape((-1,8))

def test(x,full_nine):
    hue=x["hue"]
    pair_hue=hue[x["species"]]
    mask=pair_hue>0
    cell=x["cell"][mask]
    species=x["species"][mask]
    chromatic_n=np.bincount(cell,minlength=x["ngrids"])
    check((chromatic_n<=x["cell_n"]).all()
          and int(chromatic_n.sum())>1500000,
          "Source chromatic GBIF presence unexpectedly sparse")
    groups=[]
    for g in x["variable_groups"]:
        subset=g[hue[g]>0]
        if len(subset)>1 and len(set(hue[subset]))>1:
            groups.append(subset)
    check(len(groups)>300 and sum(map(len,groups))>1000,
          "Insufficient conditioned nonwhite interchangeability")
    counts=count_chromatic(x,hue,cell,species)
    check(np.array_equal(counts.sum(axis=1),chromatic_n),
          "Chromatic colour counts inconsistent with original GBIF species records")
    obs=empirical_entropy(counts)
    masks={"nonwhite_richness_ge30":chromatic_n>=30,
           "nonwhite_richness_ge100":chromatic_n>=100}
    check(masks["nonwhite_richness_ge30"].sum()>7000
          and masks["nonwhite_richness_ge100"].sum()>3000,
          "Chromatic species geographic coverage unexpected")
    rng=np.random.default_rng(SEED+803)
    null_acc=np.zeros(x["ngrids"])
    null_sums={k:[] for k in masks}
    for i in range(N_PERM):
        perm=hue.copy()
        for g in groups:
            perm[g]=hue[g][rng.permutation(len(g))]
        check(np.array_equal(perm==0,hue==0),
              "White photo hue must remain fixed in chromatic-only null")
        h=empirical_entropy(count_chromatic(x,perm,cell,species))
        null_acc+=h
        for name,condition in masks.items():
            null_sums[name].append(float(h[condition].mean()))
        if (i+1)%50==0:
            print("MOIS_CHROMATIC8_WHITE_FIXED_SOURCE_NULL_PERMUTATIONS",i+1,flush=True)
    mu=null_acc/N_PERM
    results={}
    for name,condition in masks.items():
        vals=np.asarray(null_sums[name])
        actual=float(obs[condition].mean())
        expected=float(vals.mean())
        delta=actual-expected
        p_one_low=(1+int((vals<=actual+1e-12).sum()))/(N_PERM+1)
        p_two=(1+int((np.abs(vals-expected)>=abs(delta)-1e-12).sum()))/(N_PERM+1)
        results[name]={
          "n_100km_cells":int(condition.sum()),
          "n_nonwhite_photo_species_mean":float(chromatic_n[condition].mean()),
          "observed_eight_chromatic_Shannon":actual,
          "genus_breadth_white_fixed_null_mean":expected,
          "observed_minus_null_eight_chromatic_Shannon":delta,
          "empirical_one_sided_more_chromatic_clustering_p":p_one_low,
          "empirical_two_sided_p":p_two,
          "null_distribution_95pct":
             [float(v) for v in np.percentile(vals,[2.5,97.5])],
          "mean_effective_number_chromatic_classes":
             float(np.exp(obs[condition]).mean()),
          "mean_original_full9_anomaly_same_cells":float(
              full_nine.loc[condition,"observed_minus_genus_null_Shannon"].mean()),
          "not_photo_color_truth":True
        }
    cells=x["cell_key"]
    coord=np.asarray(cells.tolist(),dtype=np.int64)
    check(coord.shape==(x["ngrids"],2),"Grid coordinates not 2-dimensional")
    pergrid=pd.DataFrame({
        "grid_x":coord[:,0],"grid_y":coord[:,1],
        "n_nonwhite_recorded_species":chromatic_n,
        "chromatic8_Shannon_observed":obs,
        "chromatic8_Shannon_null_white_fixed":mu,
        "chromatic8_Shannon_anomaly":obs-mu,
        "ninecolour_Shannon_anomaly":full_nine.observed_minus_genus_null_Shannon,
        "original_n_white_photo_species":full_nine.n_white_photo_species,
        "original_n_all_photo_species":full_nine.n_accepted_species,
    })
    return results,pergrid,len(groups),int(sum(map(len,groups)))

def run(args):
    root=Path(args.out)
    check(not root.exists(),"Refuse source-backed ecological result overwrite")
    gbif,spp=signed_original(args.GBIF,args.ENV)
    x=baseline(gbif,spp,null_model="genus_breadth")
    full=load_strict_nine_null(args.NINE)
    coords=np.asarray(x["cell_key"].tolist(),dtype=np.int64)
    check(np.array_equal(coords,full[["grid_x","grid_y"]].to_numpy()),
          "Strict original nine-hue source grid coordinates differ from fresh GBIF grid")
    results,grid,groups,n_species=test(x,full)
    root.mkdir(parents=True)
    csv=root/"global_nonwhite_eightcolour_genus_breadth_conditional_whitefixed.csv.gz"
    grid.to_csv(csv,index=False,compression="gzip")
    report={
       "study":"conditional_nonwhite_colour_assembly_after_fixing_global_white_species",
       "n_accepted_species":37104,"n_source_species_grid_pairs":5493931,
       "n_global_source_100km_cells":18577,
       "n_variable_chromatic_shuffle_strata":groups,
       "n_species_in_conditional_shuffle":n_species,
       "null_model":"within_genus_occupancy5_bins_permute_nonwhite_labels_only",
       "number_of_whole_species_label_permutations":N_PERM,
       "white_species_identity_and_geography_fixed_in_every_null":True,
       "nonwhite_species_footprints_fixed_in_every_null":True,
       "same_original_global_9hue_grid_checked":True,
       "source_sha256":{"GBIF":GBIF_SHA,"species_ENV":ENV_SHA,
                        "strict9_null":NULL_SHA},
       "results":results,
       "not_a_true_flower_petal_validation":True,
       "does_not_infer_pollinator_selection_or_character_displacement":True
    }
    summary=root/"chromatic8_conditional_null_summary.json"
    summary.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    report["output_sha256"]={p.name:sha(p) for p in (csv,summary)}
    (root/"receipt.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print("MOIS_CHROMATIC8_WHITE_FIXED_REAL_RESULT",json.dumps(report["results"]),flush=True)

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    for key in ("GBIF","ENV","NINE","out"):parser.add_argument("--"+key,required=True)
    run(parser.parse_args())
