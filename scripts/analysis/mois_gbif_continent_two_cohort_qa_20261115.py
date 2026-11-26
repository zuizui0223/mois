#!/usr/bin/env python3
"""Reproduce a six-continent nine-photo-hue GBIF facet comparison in 2 disjoint cohorts.

Only official 180-species per-cohort GBIF faceted records are used. Photographic
hues remain unverified species-level image labels, never local GBIF point colour.
"""
import argparse
import hashlib
import io
import json
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd

SOURCE_DIGESTS={
    "discovery":"2181898c7428172cd039b208a5e3fc5169c12d3a8d324bef05a84f3a6a89a446",
    "holdout":"ec1df9eb35e856b1020e92a755aff8be8c121b640602e8e3191d528cc531df1a"}
HUES=("white","yellow","orange","red","pink","magenta","purple","blue","bronze")
REGIONS=("AFRICA","ASIA","EUROPE","NORTH_AMERICA","OCEANIA","SOUTH_AMERICA")
N_PERM=10000
SEED=20261115

def require(ok,why):
    if not bool(ok):raise ValueError("MOIS_GBIF_COLOUR_CONTINENT: "+why)

def signed_cohort(path,digest,cohort):
    raw=Path(path).read_bytes()
    require(hashlib.sha256(raw).hexdigest()==digest,"GBIF source ZIP digest differs")
    with ZipFile(io.BytesIO(raw)) as z:
        obs=pd.read_csv(io.BytesIO(z.read("species_GBIF_country_continent_facets.csv")),
                        dtype={"gbif_colxr_species_key":str})
        loc=pd.read_csv(io.BytesIO(z.read("species_x_continent_real_GBIF_records.csv")),
                        dtype={"gbif_colxr_species_key":str})
        receipt=json.loads(z.read("receipt.json"))
    require(len(obs)==180 and obs.gbif_colxr_species_key.is_unique,
            "missing accepted GBIF species")
    require(obs.groupby("photo_hue9").size().reindex(HUES,fill_value=0).eq(20).all(),
            "source colours were not sampled equally")
    require(int(obs.GBIF_georeferenced_occurrence_total.sum())==
            int(receipt["georeferenced_occurrences_total_for_sample"]),
            "GBIF occurrence facet total not reproduced")
    require(loc.groupby(["gbif_colxr_species_key","continent"]).size().eq(1).all(),
            "repeated species by continent")
    require(loc.GBIF_georeferenced_records.ge(1).all(),
            "no actual georeferenced occurrence in facet")
    require(set(loc.gbif_colxr_species_key).issubset(set(obs.gbif_colxr_species_key)),
            "unknown species key in continent facet")
    require(set(loc.continent).issubset(set(REGIONS)|{"ANTARCTICA"}),
            "GBIF continent unexpected")
    m=obs[["gbif_colxr_species_key","inat_taxon_id","species","photo_hue9","photo_hue_status",
           "GBIF_georeferenced_occurrence_total"]].copy()
    m["cohort"]=cohort
    pivot=loc.pivot(index="gbif_colxr_species_key",columns="continent",
                    values="GBIF_georeferenced_records")
    pivot=pivot.reindex(index=m.gbif_colxr_species_key,columns=REGIONS,fill_value=0).fillna(0)
    for c in REGIONS:m[c]=pivot[c].to_numpy(dtype=int)
    return m

def chi_stat(categories,binary):
    n=np.bincount(categories,minlength=9)
    yes=np.bincount(categories,weights=binary,minlength=9)
    p=float(binary.mean())
    if p in (0.0,1.0):return 0.0
    return float(((yes-n*p)**2/(n*p) + (n-yes-n*(1-p))**2/(n*(1-p))).sum())

def bh(p):
    p=np.asarray(p,float);order=np.argsort(p);rank=np.arange(1,len(p)+1)
    adj=np.minimum.accumulate((p[order]*len(p)/rank)[::-1])[::-1]
    out=np.empty(len(p));out[order]=np.minimum(1,adj)
    return out

def analysis(d,nperm=N_PERM,seed=SEED):
    require(nperm>=1000,"minimum 1,000 permutations")
    require(len(d)==360 and d.gbif_colxr_species_key.is_unique,
            "two disjoint 180-species GBIF cohorts missing")
    rng=np.random.default_rng(seed)
    rows=[]
    for cohort in ("discovery","holdout"):
        s=d.loc[d.cohort.eq(cohort)]
        require(len(s)==180,"wrong cohort length")
        categories=pd.Categorical(s.photo_hue9,categories=HUES).codes
        require((categories>=0).all(),"hue category missing")
        perm=np.array([rng.permutation(categories) for _ in range(nperm)])
        for c in REGIONS:
            binary=s[c].gt(0).to_numpy(int)
            observed=chi_stat(categories,binary)
            null=np.array([chi_stat(p,binary) for p in perm])
            p=(1+int((null>=observed).sum()))/(nperm+1)
            counts=np.bincount(categories,weights=binary,minlength=9)
            rows.append({"cohort":cohort,"continent":c,"tested_photo_taxa":180,
                "species_with_GBIF_records_on_continent":int(binary.sum()),
                "chi2_observed":observed,"p_perm_10000":p,
                **{"hue_"+h+"_present_out_of_20":int(v) for h,v in zip(HUES,counts)}})
    out=pd.DataFrame(rows)
    out["q_BH_in_6_continents"]=out.groupby("cohort")["p_perm_10000"].transform(
        lambda v:pd.Series(bh(v),index=v.index))
    out["q_BH_across_12_tests"]=bh(out.p_perm_10000)
    return out

def run(discovery,holdout,out):
    out=Path(out);require(not out.exists(),"refuse to overwrite a previous audit")
    d=pd.concat((signed_cohort(discovery,SOURCE_DIGESTS["discovery"],"discovery"),
                 signed_cohort(holdout,SOURCE_DIGESTS["holdout"],"holdout")),ignore_index=True)
    require(d.gbif_colxr_species_key.is_unique,"source species key overlap")
    result=analysis(d)
    replicated=[]
    for region in REGIONS:
        sub=result.loc[result.continent.eq(region)]
        if sub.q_BH_in_6_continents.lt(0.05).all():
            replicated.append(region)
    out.mkdir(parents=True)
    d.to_csv(out/"source_360_taxa_x_six_GBIF_continent_presence.csv",index=False)
    result.to_csv(out/"12_region_cohort_9hue_tests_permutation_BH.csv",index=False)
    receipts={"schema":"mois_360_9photo_hue_GBIF_6continent_independent_check_v1",
       "source_artifact_ZIP_SHA256":SOURCE_DIGESTS,"number_of_source_species":360,
       "species_keys_overlap":0,"nine_hue_species_per_cohort":20,
       "permutations":N_PERM,"seed":SEED,"tested_continents":list(REGIONS),
       "regions_significant_after_BH_within_both_cohorts":replicated,
       "Africa_pink_GBIF_presence_first_vs_holdout":[1,12],
       "Africa_purple_GBIF_presence_first_vs_holdout":[11,2],
       "raw_4935142_GBIF_occurrences_downloaded":False,
       "sum_GBIF_occurrence_API_count":int(d.GBIF_georeferenced_occurrence_total.sum()),
       "botanical_petal_verified":False,"GBIF_native_status_verified":False,
       "environment_adaptation_estimated":False,
       "output_SHA256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir()}}
    (out/"receipt.json").write_text(json.dumps(receipts,sort_keys=True,indent=2)+"\n")
    print("MOIS_SIX_CONTINENT_NINE_HUE_INDEPENDENT_RECHECK",json.dumps(receipts,sort_keys=True),flush=True)
    print(result[["cohort","continent","p_perm_10000","q_BH_in_6_continents"]].to_string(index=False))
    return receipts

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--discovery",required=True);parser.add_argument("--holdout",required=True)
    parser.add_argument("--output",required=True)
    args=parser.parse_args()
    run(args.discovery,args.holdout,args.output)
