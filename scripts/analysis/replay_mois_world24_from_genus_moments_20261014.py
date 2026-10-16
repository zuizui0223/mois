#!/usr/bin/env python3
"""Reconstruct the world-climate26 genus-equal 24-comparison maxT intervals.

Input is a RESTRICTED-INTERNAL, genus-anonymous set of covariance moments:
1934 9x9 crossproducts, 1934 9x4 crossproducts and genus species counts.
No photo IDs, genus strings, species names or point coordinates are required.

This archive has not been rights/privacy cleared for open publication; do not
treat an internal file handoff as a stable anonymous peer-review data deposit.
The original source 999-resample *within-scope* maxT protects six axes × four
photographed hues, NOT the wider post-outcome search over locations or designs.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np

AXES=("wc_bio5","wc_bio12","wc_bio14","wc_bio15","wc_srad_annual_mean","wc_vapr_annual_mean")
HUES=("white","yellow_orange","red_pink","blue_purple")
SEED=2026101430

def replay(moments:Path, ledger:Path, manifest:Path|None=None) -> dict:
    if manifest is not None:
        meta=json.loads(manifest.read_text(encoding="utf8"))
        assert meta["n_genera"]==1934 and meta["bootstrap_seed"]==SEED
        for k in ("identifiers_exported","geocoordinates_exported","genus_names_exported"):
            assert meta[k] is False
    with np.load(moments,allow_pickle=False) as arr:
        x=np.asarray(arr["xtx_genus"],dtype=float)
        y=np.asarray(arr["xty_genus"],dtype=float)
        ns=np.asarray(arr["n_species_per_genus"],dtype=int)
    assert x.shape==(1934,9,9) and y.shape==(1934,9,4) and ns.shape==(1934,)
    assert np.isfinite(x).all() and np.isfinite(y).all()
    assert ns.min()>=3 and int(ns.sum())==20694
    n=len(ns)
    rng=np.random.default_rng(SEED)
    draws=rng.multinomial(n,np.repeat(1/n,n),size=999)
    bx=np.einsum("bg,gij->bij",draws,x,optimize=True)
    by=np.einsum("bg,gik->bik",draws,y,optimize=True)
    totalx=x.sum(axis=0);totaly=y.sum(axis=0)
    point=np.empty((6,4));boot=np.empty((999,6,4))
    for j in range(6):
        sub=np.array((0,1,2,j+3))
        point[j]=np.linalg.solve(totalx[np.ix_(sub,sub)],totaly[sub])[-1]
        boot[:,j]=np.linalg.solve(bx[:,sub,:][:,:,sub],by[:,sub,:])[:,-1]
    dev=boot-point[None,:,:]
    sd=dev.std(axis=0,ddof=1)
    assert (sd>1e-12).all()
    crit=float(np.quantile(np.max(np.abs(dev/sd),axis=(1,2)),0.95))
    ll=point-crit*sd;hh=point+crit*sd
    qlo=np.quantile(boot,0.025,axis=0);qhi=np.quantile(boot,0.975,axis=0)
    with ledger.open(newline="",encoding="utf8") as f:
        rows=list(csv.DictReader(f))
    assert len(rows)==1280
    maxerr=0.0;out=[]
    for j,a in enumerate(AXES):
        for k,h in enumerate(HUES):
            matching=[r for r in rows if all((r["mode"]=="climate",r["quality"]=="all",
                r["region"]=="world",r["adjustment"]=="genus_geo",r["weight"]=="genus_equal",
                r["anchor"]==a,r["hue"]==h))]
            assert len(matching)==1
            ref=matching[0]
            for observed,column in ((point[j,k]*100,"direction_pp_per_env_global_sd"),
                (ll[j,k]*100,"maxT_ci_lo_pp"),(hh[j,k]*100,"maxT_ci_hi_pp"),
                (qlo[j,k]*100,"cluster_ci_lo_pp"),(qhi[j,k]*100,"cluster_ci_hi_pp")):
                err=abs(observed-float(ref[column]))
                maxerr=max(maxerr,err)
                if err>1e-9:
                    raise AssertionError((a,h,column,err))
            supported=bool(ll[j,k]>0 or hh[j,k]<0)
            assert supported==(ref["maxT_excludes_zero"]=="True")
            out.append({"axis":a,"hue":h,"point_pp":float(point[j,k]*100),
                "maxT95_pp":[float(ll[j,k]*100),float(hh[j,k]*100)],
                "maxT_excludes_zero":supported})
    assert len(out)==24 and np.allclose(point.sum(axis=1),0,atol=1e-12)
    return {"schema":"mois_genus_anonymous_world24_bootstrap_replay_v1",
        "status":"PASS_ALL_24_FROZEN_POINT_NOMINAL_AND_MAXT_CI",
        "n_genera":n,"n_species":int(ns.sum()),"bootstrap_replicates":999,
        "seed":SEED,"maxT_critical_value":crit,"max_abs_difference_pp":maxerr,
        "moments_sha256":hashlib.sha256(moments.read_bytes()).hexdigest(),
        "not_source_independent_biology":True,"rights_clearance_required":True,
        "rows":out}

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--suffstats",type=Path,required=True)
    p.add_argument("--ledger",type=Path,required=True)
    p.add_argument("--manifest",type=Path)
    p.add_argument("--output",type=Path)
    args=p.parse_args()
    receipt=replay(args.suffstats,args.ledger,args.manifest)
    if args.output:
        args.output.write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf8")
    print(json.dumps({k:v for k,v in receipt.items() if k!="rows"},indent=2))
