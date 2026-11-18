#!/usr/bin/env python3
"""Independent 500-km block bootstrap of same-cell WCVP-native minus
introduced photo-white species-share and exact genus-pool decomposition.

Source is the signed numerical artifact from MOIS genuine native-intro
ecological models (not a new image/petal colour truth or new GBIF dataset).
Strict/moderate cohorts overlap; latitude cuts are descriptive exploratory.
"""
from __future__ import annotations
import argparse,hashlib,io,json
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd

EXPECTED_SOURCE_SHA="3efaf6ee359aca4a49713d55f5bbadd0c86467ccd4d27c300f1e06936245759a"
MEMBER="WCVP_native_introduced_matched_100km_white_genus_components.csv.gz"
OUTCOMES=[
  "observed_native_minus_introduced",
  "genus_pool_native_minus_introduced",
  "within_genus_native_minus_introduced",
]
BOOT_N=20000
SEED=20261118

def fail(condition,msg):
    if not bool(condition):raise ValueError("MOIS_NATIVE_PAIRED_CLUSTER_BOOTSTRAP: "+msg)

def sha(data):
    if isinstance(data,(str,Path)):data=Path(data).read_bytes()
    return hashlib.sha256(data).hexdigest()

def source(path):
    fail(sha(path)==EXPECTED_SOURCE_SHA,"Signed original matched-cell ZIP changed")
    with ZipFile(path) as z:
        receipt=json.loads(z.read("receipt.json"))
        for name,digest in receipt["output_sha256"].items():
            fail(sha(z.read(name))==digest,"Internal original source SHA changed: "+name)
        d=pd.read_csv(io.BytesIO(z.read(MEMBER)),compression="gzip")
    fail(len(d)==6708 and d.WCVP_geo_mode.value_counts().to_dict()==
         {"strict":3106,"moderate":3602},
         "Original current GBIF matched WCVP cell counts changed")
    fail(not d.duplicated(["WCVP_geo_mode","grid_x","grid_y"]).any(),
         "Duplicate original geographic cells")
    fail(np.allclose(
         d[OUTCOMES[0]],d[OUTCOMES[1]]+d[OUTCOMES[2]],
         atol=1e-12,rtol=0),"Photo white genus-pool identity broke")
    return d

def bootstrap(sub,weights,seed):
    x=sub[["block_500km"]].copy()
    x["w"]=weights
    for name in OUTCOMES:
        x[name]=sub[name].to_numpy(float)*weights
    grouped=x.groupby("block_500km",sort=False)[["w"]+OUTCOMES].sum()
    values=grouped.to_numpy(float)
    nblocks=len(values)
    fail(nblocks>=30 and len(sub)>=100,"insufficient real geographic support")
    rng=np.random.default_rng(seed)
    boot=np.empty((BOOT_N,len(OUTCOMES)))
    for start in range(0,BOOT_N,200):
        n=min(200,BOOT_N-start)
        draws=rng.integers(0,nblocks,size=(n,nblocks))
        agg=values[draws].sum(axis=1)
        boot[start:start+n,:]=100*agg[:,1:]/agg[:,:1]
    measured=100*np.average(sub[OUTCOMES].to_numpy(float),
                             weights=weights,axis=0)
    lower=np.percentile(boot,2.5,axis=0)
    upper=np.percentile(boot,97.5,axis=0)
    fail(abs(measured[0]-measured[1]-measured[2])<1e-10,
         "original native-minus-intro mean decomposition fails")
    return [{
        "response_component":key,
        "effect_pp":float(measured[j]),
        "bootstrap_CI95_low_pp":float(lower[j]),
        "bootstrap_CI95_high_pp":float(upper[j]),
        "n_500km_clusters":nblocks,
        "n_100km_cells":len(sub),
        "n_block_boot":BOOT_N
    } for j,key in enumerate(OUTCOMES)]

def run(path,out):
    dest=Path(out);fail(not dest.exists(),"Refuse prior output overwrite")
    d=source(path)
    rows=[]
    cuts={
      "all_same_cells":lambda z:z.lat.notna(),
      "absolute_lat_under30":lambda z:z.lat.abs()<30,
      "absolute_lat_30plus":lambda z:z.lat.abs()>=30,
      "north_lat_30plus":lambda z:z.lat>=30,
      "south_lat_minus30_below":lambda z:z.lat<=-30,
    }
    for mode in ("strict","moderate"):
        z=d.loc[d.WCVP_geo_mode.eq(mode)]
        for label,criterion in cuts.items():
            subset=z.loc[criterion(z)].copy()
            for weighting in ("equal_cells","min_native_introduced_support"):
                w=(np.ones(len(subset),float) if weighting=="equal_cells"
                   else subset[["colour_species_n__native",
                                "colour_species_n__introduced"]].min(axis=1).to_numpy(float))
                seed=SEED+(0 if mode=="strict" else 10000)+len(subset)+(
                    3000 if weighting!="equal_cells" else 0)
                res=bootstrap(subset,w,seed)
                rows.extend([{
                  "WCVP_geo_mode":mode,"latitude_stratum":label,
                  "weighting":weighting,**r
                } for r in res])
    dest.mkdir(parents=True)
    table=dest/"spatial_block_bootstrap_native_intro_white_genus_components.csv"
    pd.DataFrame(rows).to_csv(table,index=False)
    sample={
      "original_source_zip_sha256":EXPECTED_SOURCE_SHA,
      "source_genuine_100km_matched_cells":{"strict":3106,"moderate":3602},
      "photo_white_fraction_not_independent_petal_truth":True,
      "WCVP_native_introduced_only_TDWG3_regional_status":True,
      "same_location_matched_native_intro":True,
      "500km_spatial_cluster_bootstrap":True,
      "strict_moderate_share_cells_not_independent_replication":True,
      "within_genus_component_is_species_composition_not_evolution":True,
      "latitude_slices_are_descriptive_not_multiplicity_corrected":True,
      "n_boot_per_stratum":BOOT_N,
      "results":rows,
      "output_sha256":{table.name:sha(table)},
    }
    (dest/"receipt.json").write_text(json.dumps(sample,indent=2,allow_nan=False)+"\n")
    for r in rows:
        if r["latitude_stratum"]=="all_same_cells" and (
           r["response_component"]=="observed_native_minus_introduced"):
            print("MOIS_NATIVE_PAIRED_BLOCK_BOOTSTRAP_SUCCESS",json.dumps(r),flush=True)
    return sample

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",required=True);p.add_argument("--out",required=True)
    a=p.parse_args();run(a.source,a.out)
