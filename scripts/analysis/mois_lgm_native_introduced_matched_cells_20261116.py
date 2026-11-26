#!/usr/bin/env python3
"""WCVP TDWG3 broad-native vs broad-introduced photo-hue fractions at EXACT
same recorded GBIF 100km cell. Tests history associations descriptively, not
species-ancestral history, wild individual nativeness or floral adaptation.
"""
from __future__ import annotations
import argparse, hashlib, io, json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd
import statsmodels.api as sm

NATIVE_SHA="b72b7358ed5382d625583ad61b90952e4e77377296fe90ae677bf9481aa80a7f"
SPATIAL_SHA="e209257adbc185b4855f3b8d76aeaa50d2e2afe8a1feeb0f76013366e57c4d16"
NATIVE="WCVP_TDWG3_native_introduced_photo9_LGM_100km_cell_summary.parquet"
SPATIAL="LGM_white_photo9_species_assemblages_100km_cells.parquet"
EXPECTED_NATIVE_REGIONAL_MATCHES=36762

def sha(data):
    if isinstance(data,(Path,str)):
        data=Path(data).read_bytes()
    return hashlib.sha256(data).hexdigest()

def require(cond,msg):
    if not bool(cond):
        raise ValueError("MOIS_LGM_MATCHED: "+str(msg))

def signed_tables(nativezip,spatialzip):
    import pyarrow.parquet as pq
    require(sha(nativezip)==NATIVE_SHA and sha(spatialzip)==SPATIAL_SHA,
            "immutable WCVP native/spatial source ZIP differs")
    with ZipFile(nativezip) as z:
        nbytes=z.read(NATIVE)
        rec=json.loads(z.read("receipt.json"))
    require(sha(nbytes)==rec["output_sha256"][NATIVE]
            and rec["source_GBIF_species_cells"]==5493931
            and rec["source_taxa"]==37104
            and not rec["individual_is_wild_native_inferred"],
            "regional native-only artifact source signature / boundaries invalid")
    for key in ("strict","moderate"):
        require(rec["results"][key]["region_and_match_audit"][
            "matched_WCVP_unique_exact_species"]==EXPECTED_NATIVE_REGIONAL_MATCHES,
            "native/source species identity changed")
    native=pq.read_table(io.BytesIO(nbytes)).to_pandas()
    with ZipFile(spatialzip) as z:
        b=z.read(SPATIAL)
        rec2=json.loads(z.read("receipt.json"))
    require(sha(b)==rec2["output_sha256"][SPATIAL]
            and rec2["results"]["n_GBIF_species_x_cell_pairs"]==5493931,
            "original real spatial LGM source changed")
    spatial=pq.read_table(io.BytesIO(b)).to_pandas()
    require(native.tdwg3_assignment_mode.isin(["strict","moderate"]).all()
            and native.groupby(["tdwg3_assignment_mode","grid_x","grid_y"]).size().max()==1,
            "native TDWG cell duplication")
    require(len(spatial)==12560 and not spatial[["grid_x","grid_y"]].duplicated().any(),
            "original LGM 100km spatial cell table changed")
    return native,spatial

def same_cell_cohort(native,spatial,mode,min_each=30):
    s=native[native.tdwg3_assignment_mode.eq(mode)].copy()
    d=s.merge(spatial,on=["grid_x","grid_y"],how="inner",validate="one_to_one",
              suffixes=("_regional","_global"))
    require(not d[["grid_x","grid_y"]].duplicated().any(),
            "same native/intro cell duplicated")
    # Only cells with genuine both-group species-guild proportions are paired.
    d=d[(d.native_species_n>=min_each)&(d.introduced_species_n>=min_each)].copy()
    require(not d.empty,"no WCVP both-group cells")
    d["native_white_fraction"]=d.native_white_n/d.native_species_n
    d["intro_white_fraction"]=d.introduced_white_n/d.introduced_species_n
    d["native_minus_intro_white"]=d.native_white_fraction-d.intro_white_fraction
    require(d.native_white_fraction.between(0,1).all()
            and d.intro_white_fraction.between(0,1).all()
            and d.native_minus_intro_white.between(-1,1).all(),
            "native/introduced white fractions not valid")
    return d

def matrix(d):
    x=pd.DataFrame(index=d.index)
    x["LGM_warming"]=d.threeGCM_delta_temperature_C_median.astype(float)
    x["LGM_precip_delta"]=d.threeGCM_delta_log_precip_median.astype(float)
    x["present_temp"]=d.present_temperature_C.astype(float)
    x["log_present_precip"]=np.log1p(d.present_precip_mm.astype(float))
    x["log_woody_AGB"]=np.log1p(d.ESA_AGB_woody_Mg_ha_2020.astype(float))
    x["abs_lat"]=d.lat.abs().astype(float)
    x["log_native_n"]=np.log1p(d.native_species_n.astype(float))
    x["log_intro_n"]=np.log1p(d.introduced_species_n.astype(float))
    x["log_GBIF_occ"]=np.log1p(d.GBIF_occurrence_count.astype(float))
    sd=x.std()
    require(sd.gt(0).all() and np.isfinite(sd).all(),
            "zero variance/nonfinite model feature")
    x=(x-x.mean())/sd
    return x,float(sd["LGM_warming"])

def fit(d,kind="difference",block="500",species_weight=False):
    if len(d)<300:
        return {"status":"insufficient_same_cell_support","n_cells":len(d)}
    x,sd=matrix(d)
    y=d[{"difference":"native_minus_intro_white",
         "native":"native_white_fraction",
         "introduced":"intro_white_fraction"}[kind]].astype(float)
    groups=d.block_500km.astype(str)
    b=d.block_500km if block=="500" else (d.block_1000km if block=="1000" else None)
    if b is None:
        x["lat2"]=(d.lat.astype(float)**2)
        x["lon2"]=(d.lon.astype(float)**2)
        x["lat2"]=(x.lat2-x.lat2.mean())/x.lat2.std()
        x["lon2"]=(x.lon2-x.lon2.mean())/x.lon2.std()
        x=sm.add_constant(x)
    else:
        x=x-x.groupby(b).transform("mean")
        y=y-y.groupby(b).transform("mean")
        useful=x.abs().sum(axis=1).gt(0)
        x=x.loc[useful];y=y.loc[useful];groups=groups.loc[useful]
    require(len(y)>200 and groups.nunique()>30,
            "too few geographic clusters for grouped uncertainty")
    if species_weight:
        # Amount of taxonomy support, NOT photo effort and NOT individual abundance.
        # Weighted is a second sensitivity only.
        counts=d.loc[y.index,["native_species_n","introduced_species_n"]].min(axis=1).astype(float)
        weights=counts/counts.median()
        m=sm.WLS(y,x,weights=weights).fit(
            cov_type="cluster",cov_kwds={"groups":groups},use_t=False)
    else:
        m=sm.OLS(y,x).fit(cov_type="cluster",cov_kwds={"groups":groups},use_t=False)
    a,b=m.conf_int().loc["LGM_warming"].tolist()
    return {"status":"fitted","n_cells":len(y),"n_spatial_500km_clusters":int(groups.nunique()),
            "effect_white_fraction_per_SD_LGM_warming":float(m.params["LGM_warming"]),
            "95CI":[float(a),float(b)],"p":float(m.pvalues["LGM_warming"]),
            "R2":float(m.rsquared),"LGM_warming_SD_C":sd,
            "response":kind,"block_fixed_effect":block,
            "weights":"min_native_intro_species_support" if species_weight else "equal_cells",
            "historical_causality_demonstrated":False}

def analyze(native,spatial):
    results={}
    summaries=[]
    for mode in ("strict","moderate"):
        results[mode]={}
        for min_each in (20,30,50):
            d=same_cell_cohort(native,spatial,mode,min_each)
            label=f"min_{min_each}_each"
            fits={}
            for block in ("global","1000","500"):
                for kind in ("native","introduced","difference"):
                    fits[f"{kind}_{block}"]=fit(d,kind,block)
            # Core paired effect robust to weighting native and intro n.
            fits["difference_500_species_support_weighted"]=fit(
                d,"difference","500",True)
            results[mode][label]={"n_matched_cells":len(d),
                "n_500km_blocks":int(d.block_500km.nunique()),
                "mean_native_white":float(d.native_white_fraction.mean()),
                "mean_intro_white":float(d.intro_white_fraction.mean()),
                "models":fits}
            summaries.append({"region_mode":mode,"min_each":min_each,
                "matched_cells":len(d),"blocks":int(d.block_500km.nunique()),
                "native_mean_white":float(d.native_white_fraction.mean()),
                "introduced_mean_white":float(d.intro_white_fraction.mean()),
                "paired_LGM_warming_beta_500km":fits["difference_500"].get(
                     "effect_white_fraction_per_SD_LGM_warming"),
                "paired_p_500km":fits["difference_500"].get("p")})
            print("MOIS_WCVP_MATCHED_CELL",json.dumps(summaries[-1],sort_keys=True),flush=True)
    return results,pd.DataFrame(summaries)

def run(nzip,szip,out):
    out=Path(out)
    require(not out.exists(),"no overwrite")
    n,s=signed_tables(nzip,szip)
    results,table=analyze(n,s)
    out.mkdir(parents=True)
    table.to_csv(out/"same_GBIF_cell_native_intro_models.csv",index=False)
    (out/"models.json").write_text(json.dumps(results,indent=2,sort_keys=True)+"\n")
    receipt={
        "schema":"mois_photo9_LGM_climate_sameGBIFcell_WCVP_native_introduced_difference",
        "original_sources_SHA256":{"native_regional":NATIVE_SHA,"LGM_spatial":SPATIAL_SHA},
        "n_source_species":37104,"n_original_source_species_cell_pairs":5493931,
        "native_at_TDWG3_does_not_validate_wild_individual":True,
        "native_at_TDWG3_does_not_validate_island_native":True,
        "present_occurrences_do_not_reconstruct_LGM_refugia":True,
        "species_photo_modal9_is_not_verified_petal_hue":True,
        "results":results,
        "output_SHA256":{p.name:sha(p) for p in out.iterdir()}}
    (out/"receipt.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    return receipt

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--native",required=True)
    p.add_argument("--spatial",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    run(a.native,a.spatial,a.output)
