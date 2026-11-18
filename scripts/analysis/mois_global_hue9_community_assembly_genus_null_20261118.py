#!/usr/bin/env python3
"""MOIS Global 9-hue assembly: genus-conditioned species-label permutation null.

Scientific question: across present-day GBIF 100-km plant assemblages, is
species-level photographic 9-hue diversity more concentrated/even than expected
from species richness, observed genus membership and each genus's global
photographic hue counts?

Species-label hue is shuffled globally *within each genus*, not per cell.
Species occurrence footprints, local species richness, genus counts and global
per-genus color frequencies are unchanged. All cells see the SAME shuffled
global species label per permutation. The null thus tests whether geography
within a genus is related to color, not whether plant flowers converge via
pollinator-mediated character displacement.

Outcome: Shannon entropy of counts of species in nine recorded photo hue
classes, one species once per present-day 100-km GBIF grid. No actual flowers
were collected or independently annotated; no local flower color observed.
"""
from __future__ import annotations
import argparse, hashlib, io, json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd

GBIF_SHA="422f6d6489e61c6e57c5d51db614e47f613039781a901e3fd3481c13b2b085c2"
GBIF_MEMBER="all_37104_species_x_EPSG6933_100km_cells.parquet"
GBIF_MEMBER_SHA="e00b0408055984abfe4ec0ff0a24833fc148888afd1ea0e0c9b270a077d08b11"
ENV_SHA="d8b97528188d9f8d63cb802d5b3f6e8f7f251937fd30af42aeae47e43ebad1fc"
ENV_MEMBER="BIO19_weather_NDVI_ESACCIv7_AGB_37104_accepted_species.parquet"
HUES=("white","yellow","orange","red","pink","magenta","purple","blue","bronze")
SEED=20261118

def check(x,msg):
    if not bool(x):
        raise ValueError("MOIS_COLOR_ASSEMBLY: "+msg)

def sha(x):
    if isinstance(x,(Path,str)):x=Path(x).read_bytes()
    return hashlib.sha256(x).hexdigest()

def signed_original(gbif_zip,env_zip):
    import pyarrow.parquet as pq
    check(sha(gbif_zip)==GBIF_SHA and sha(env_zip)==ENV_SHA,
          "Source immutable ZIP digest differs")
    with ZipFile(gbif_zip) as z:
        b=z.read(GBIF_MEMBER)
        rec=json.loads(z.read("receipt.json"))
    check(sha(b)==GBIF_MEMBER_SHA and
          rec["all_16_shards_source_verified"] and
          rec["GBIF_occurrences_passed_all_filters"]==396142825,
          "Original global GBIF 16-shard lineage invalid")
    pairs=pq.read_table(io.BytesIO(b),
                        columns=["specieskey","grid_x","grid_y","occurrence_count"]).to_pandas()
    check(len(pairs)==5493931 and int(pairs.occurrence_count.sum())==396142825,
          "Actual original accepted-species x grid-pair source changed")
    with ZipFile(env_zip) as z:
        b=z.read(ENV_MEMBER)
        rec=json.loads(z.read("receipt.json"))
    check(sha(b)==rec["output_SHA256"][ENV_MEMBER]
          and rec["source_global_accepted_species"]==37104,
          "Original accepted species and per-species hue table changed")
    species=pq.read_table(io.BytesIO(b)).to_pandas()
    species["key"]=species.key.astype(str)
    check(len(species)==37104 and species.key.is_unique
          and set(species.source_hue9)==set(HUES)
          and species.genus.notna().all(),"Original species-hue catalog bad")
    return pairs,species

def empirical_entropy(counts):
    cnt=counts.astype(float)
    n=cnt.sum(axis=1)
    p=cnt/np.maximum(n[:,None],1)
    h=-(np.where(p>0,p*np.log(np.maximum(p,1e-12)),0)).sum(axis=1)
    return h

def baseline(pairs,species,null_model="genus"):
    check(null_model in ("genus","genus_breadth"),
          "Unexpected conditional community null model")
    spp=pd.Index(species.key.astype(str))
    ix=spp.get_indexer(pairs.specieskey.astype(str))
    check((ix>=0).all(),"GBIF species not in signed 37,104 accepted color catalog")
    grid_index=pd.MultiIndex.from_frame(pairs[["grid_x","grid_y"]])
    codes,unique_grid=pd.factorize(grid_index,sort=True)
    ngrid=len(unique_grid)
    check(ngrid==18577,"Original 18,577 actual GBIF spatial cells changed")
    check(not pairs[["specieskey","grid_x","grid_y"]].duplicated().any(),
          "Duplicate species recorded in single 100km cell")
    cell_n=np.bincount(codes,minlength=ngrid)
    hue=np.array([HUES.index(k) for k in species.source_hue9],dtype=np.int8)
    genus=species.genus.to_numpy(str)
    gcodes,_=pd.factorize(genus,sort=True)
    # More restrictive alternative fixes genus AND species occupied-grid
    # footprint size bands. The colour-richness association with occupancy
    # breadth is thus not arbitrarily erased in every genus permutation.
    species_grid_n=np.bincount(ix,minlength=len(species)).astype(np.int32)
    available=np.flatnonzero(species_grid_n>0)
    check(len(available)>30000,
          "Too few genuine observed species have GBIF grid footprints")
    # DO NOT shuffle absent-source species into any real assemblage.
    # Their photo hues are outside the documented 100km GBIF species pool.
    if null_model=="genus":
        group_codes=gcodes[available].astype(np.int64)
    else:
        breadth=np.digitize(species_grid_n,[2,6,21,101]).astype(np.int32)
        group_codes=gcodes[available].astype(np.int64)*5+breadth[available]
    indexed_sort=np.argsort(group_codes,kind="stable")
    order=available[indexed_sort]
    sorted_group_codes=group_codes[indexed_sort]
    start=np.r_[0,1+np.flatnonzero(np.diff(sorted_group_codes))]
    end=np.r_[start[1:],len(order)]
    groups=[order[a:b] for a,b in zip(start,end)]
    variable_groups=[g for g in groups if len(g)>1 and len(set(hue[g]))>1]
    varying_species=int(sum(len(g) for g in variable_groups))
    check(len(variable_groups)>250 and varying_species>2500,
          "Insufficient multi-colour within-condition null permutations")
    return {
      "cell":codes.astype(np.int32),
      "species":ix.astype(np.int32),
      "cell_n":cell_n,"ngrids":ngrid,
      "cell_key":unique_grid,"hue":hue,
      "variable_groups":variable_groups,
      "n_genera":len(set(gcodes)),
      "n_variable_null_strata":len(variable_groups),
      "n_genera_with_variable_null_strata":
          len(set(int(gcodes[g[0]]) for g in variable_groups)),
      "n_species_in_variable_strata":varying_species,
      "n_accepted_species_with_at_least_one_GBIF_cell":int(len(available)),
      "n_accepted_species_without_GBIF_cells":int(len(species)-len(available)),
      "null_model":null_model,
    }

def counts_by_cell(x,labels):
    encoded=x["cell"].astype(np.int64)*9+labels[x["species"]]
    return np.bincount(encoded,minlength=x["ngrids"]*9).reshape((-1,9))

def test(x,permutations):
    check(permutations>=99 and permutations<=2000,
          "Permutation count out of specified safe range")
    observed=counts_by_cell(x,x["hue"])
    check(np.array_equal(observed.sum(axis=1),x["cell_n"]),
          "Source species counts no longer correspond to 9 colour classes")
    H_obs=empirical_entropy(observed)
    masks={"richness_ge30":x["cell_n"]>=30,
           "richness_ge100":x["cell_n"]>=100}
    check(masks["richness_ge30"].sum()>6000 and
          masks["richness_ge100"].sum()>1000,
          "Unexpected GBIF recorded plant species richness coverage")
    rng=np.random.default_rng(SEED)
    ngrid=x["ngrids"]
    acc=np.zeros(ngrid)
    acc2=np.zeros(ngrid)
    high=np.zeros(ngrid,dtype=np.int32)
    low=np.zeros(ngrid,dtype=np.int32)
    test_values={name:[] for name in masks}
    # Controls: exactly ONE original global hue permutation per species per run,
    # with all native/introduction and photographed source IDs unchanged.
    for p in range(permutations):
        perm=x["hue"].copy()
        for group in x["variable_groups"]:
            perm[group]=x["hue"][group][rng.permutation(len(group))]
        H=empirical_entropy(counts_by_cell(x,perm))
        acc+=H
        acc2+=H*H
        high+=(H>=H_obs-1e-12).astype(np.int32)
        low+=(H<=H_obs+1e-12).astype(np.int32)
        for name,m in masks.items():
            test_values[name].append(float(H[m].mean()))
        if (p+1)%50==0:
            print("MOIS_COLOR_ASSEMBLY_FULL_GBIF_GENUS_LABEL_SHUFFLES",p+1,flush=True)
    n=float(permutations)
    mu=acc/n
    sd=np.sqrt(np.maximum((acc2-n*mu**2)/(n-1),0))
    ses=np.divide(H_obs-mu,sd,out=np.full(ngrid,np.nan),
                  where=sd>1e-10)
    results={}
    for name,mask in masks.items():
        null_dist=np.asarray(test_values[name])
        obs_mean=float(H_obs[mask].mean())
        mean_null=float(null_dist.mean())
        offset=obs_mean-mean_null
        # Proper single global permutation statistic; cells are NOT treated
        # as independent samples or independent null replicates.
        p2=float((1+np.count_nonzero(
            np.abs(null_dist-mean_null)>=abs(offset)-1e-13
        ))/(permutations+1))
        p_high=float((1+np.count_nonzero(null_dist>=obs_mean-1e-13))
                     /(permutations+1))
        p_low=float((1+np.count_nonzero(null_dist<=obs_mean+1e-13))
                    /(permutations+1))
        results[name]={
          "n_global_grid_cells":int(mask.sum()),
          "mean_recorded_species_richness":float(x["cell_n"][mask].mean()),
          "observed_mean_Shannon_natural_log":obs_mean,
          "genus_preserving_null_mean_Shannon":mean_null,
          "observed_minus_null_Shannon":float(offset),
          "null_95pct_permutation_interval":[float(v) for v in
              np.percentile(null_dist,[2.5,97.5])],
          "mean_exp_Shannon_effective_photo_colours":
              float(np.exp(H_obs[mask]).mean()),
          "empirical_two_sided_permutation_p":p2,
          "empirical_one_sided_overdispersion_p":p_high,
          "empirical_one_sided_clustering_p":p_low,
          "note":"All ecological species-cell correlations preserved in shared whole-catalogue permutation"
        }
    z=x["cell_key"]
    # pandas.factorize(MultiIndex) may return an Index of tuples without
    # source level names. Preserve coordinate position, never infer labels.
    coords=np.asarray(z.tolist(),dtype=np.int64)
    check(coords.shape==(ngrid,2),"GBIF grid x/y key dimensions changed")
    cell=pd.DataFrame({
        "grid_x":coords[:,0],
        "grid_y":coords[:,1],
        "n_accepted_species":x["cell_n"],
        "observed_hue9_Shannon":H_obs,
        "expected_hue9_Shannon_given_genus":mu,
        "null_SD_Shannon_given_genus":sd,
        "observed_minus_genus_null_Shannon":H_obs-mu,
        "SES_genus_conditioned":ses,
        "n_unique_observed_photo_hues":(observed>0).sum(axis=1),
        "n_white_photo_species":observed[:,0]
    })
    for name,m in masks.items():
        cell[name]=m
    counts={k:int(np.sum(x["hue"]==i)) for i,k in enumerate(HUES)}
    return results,cell,counts

def synthetic():
    # Invariant test: species hue permuted globally, local species richness
    # and genus composition stay exactly unchanged.
    source={
        "cell":np.array([0,0,1,1,2,2,3,3],dtype=np.int32),
        "species":np.array([0,2,1,3,0,1,2,3],dtype=np.int32),
        "ngrids":4,
        "cell_n":np.array([2,2,2,2]),
        "hue":np.array([0,1,2,3],dtype=np.int8),
        "variable_groups":[np.array([0,1]),np.array([2,3])]
    }
    r=counts_by_cell(source,source["hue"])
    assert np.array_equal(r.sum(axis=1),source["cell_n"])
    assert np.all(empirical_entropy(r)>=0)
    p=source["hue"].copy()
    p[np.array([0,1])]=source["hue"][[1,0]]
    assert np.array_equal(counts_by_cell(source,p).sum(axis=1),source["cell_n"])
    print("MOIS_GENUS_CONDITIONAL_9HUE_SYNTHETIC_IDENTITY_PASS",flush=True)

def main(gbif,env,out,n,null_model):
    root=Path(out)
    check(not root.exists(),"Refuse overwriting previous genuine global analysis")
    pairs,species=signed_original(gbif,env)
    x=baseline(pairs,species,null_model=null_model)
    results,cell,huecounts=test(x,n)
    root.mkdir(parents=True)
    p=root/"genuine_global_100km_photo9_diversity_vs_genus_null.csv.gz"
    cell.to_csv(p,index=False,compression="gzip")
    summary={
        "source_species":37104,"source_species_grid_pairs":5493931,
        "source_photo_hue_counts":huecounts,
        "num_global_grid_cells":x["ngrids"],
        "n_original_genera":x["n_genera"],
        "null_model":x["null_model"],
        "n_variable_hue_genera":x["n_genera_with_variable_null_strata"],
        "n_variable_shuffle_strata":x["n_variable_null_strata"],
        "n_species_in_exchangeable_colour_strata":x["n_species_in_variable_strata"],
        "n_accepted_species_with_GBIF_cells":x["n_accepted_species_with_at_least_one_GBIF_cell"],
        "n_accepted_species_without_GBIF_cells":x["n_accepted_species_without_GBIF_cells"],
        "genus_preserving_global_species_hue_label_permutations":n,
        "random_seed":SEED,
        "source_GBIF_sha256":GBIF_SHA,
        "source_ENV_sha256":ENV_SHA,
        "results":results,
        "not_observed_local_flower_petal_colour":True,
        "GBIF_presence_is_not_plant_abundance":True,
        "introduced_and_native_species_not_filtered":True,
        "null_preserves_genus_palette_and_species_occupancy_footprints":True,
        "null_breaks_species_hue_geography_within_genus":True,
        "null_does_not_establish_pollinator_selection_or_trait_evolution":True,
    }
    s=root/"ecological_assembly_null_summary.json"
    s.write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    receipt={**summary,"output_sha256":{p.name:sha(p),s.name:sha(s)}}
    (root/"receipt.json").write_text(json.dumps(receipt,indent=2,allow_nan=False)+"\n")
    print("MOIS_GLOBAL_FLORAL_COLOUR_ASSEMBLY_REAL_RESULT",
          json.dumps({**results,"null_model":null_model,
                     "n_genera":x["n_genera"],
                     "variable_genera":x["n_genera_with_variable_null_strata"],
                     "variable_shuffle_strata":x["n_variable_null_strata"],
                     "exchangeable_species":x["n_species_in_variable_strata"],
                     "hue_species":huecounts},allow_nan=False),flush=True)

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--gbif");parser.add_argument("--env")
    parser.add_argument("--out")
    parser.add_argument("--permutations",type=int,default=199)
    parser.add_argument("--null-model",choices=("genus","genus_breadth"),
                        default="genus")
    parser.add_argument("--selftest",action="store_true")
    args=parser.parse_args()
    if args.selftest:synthetic()
    else:
        check(args.gbif and args.env and args.out,"GBIF+ENV genuine signed zip paths needed")
        main(args.gbif,args.env,args.out,args.permutations,args.null_model)
