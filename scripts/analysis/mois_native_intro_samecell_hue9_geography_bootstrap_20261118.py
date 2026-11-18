#!/usr/bin/env python3
"""MOIS 2026-11-18: WCVP native vs introduced photo9 Shannon anomalies
in the SAME authentic GBIF 100-km grid cells, across island, altitude,
Arctic and hemisphere geographic strata.

Input is an immutable signed WCVP+TDWG3 source (199 genus×occupied-range
label permutations already executed globally) with 18,577 grid rows.
For each same grid with >=30 species of EACH WCVP status, define:

delta = [H9(native)-E_null(H9_native)] -
        [H9(introduced)-E_null(H9_introduced)].

Positive delta means introduced communities are MORE hue-concentrated
relative to the global species-label null. The null is a randomization of
PHOTO-derived nominal hue labels, not a null of ecological introduction.

500km cluster bootstrap resamples original GBIF geographic 500km blocks,
keeping both statuses within each 100km cell together. For geographical
high/low difference, also restrict to 500km blocks containing both strata.
No spatial bootstrap can turn an unbalanced geography into a controlled
biogeographic intervention. WCVP native is TDWG3 regional distribution,
not direct native individual or validated flower petal identity.
"""
from __future__ import annotations
import argparse,hashlib,io,json
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd

SOURCE_SHA="cbe1eef381e3e1a0062dbf5134b047f947338b7c499d2ecc0cd22ee8855909fd"
SOURCE_MEMBER="global_WCVP_TDWG3_native_only_photo9_island_alpine_polar_south_cells.csv.gz"
N_PERM=199
N_BOOT=20000
SEED=20261118
MODES=("strict","moderate")
OUTCOMES=(
    "native_minus_intro_9hue_null_anomaly",
    "native_minus_intro_observed_Shannon",
    "native_minus_intro_expected_null_Shannon",
    "native_minus_intro_white_fraction_pp",
)
REGION_CONTRASTS=(
    ("high_minus_low_elevation","high_elevation_1500plus","low_elevation_below500"),
    ("island_minus_mainland","island_7of9","mainland_7of9"),
    ("arctic_minus_boreal","arctic_66p562plus","boreal_50to66p562"),
    ("south_minus_north","southern_hemisphere","northern_hemisphere"),
    ("southern_mid_minus_northern_mid","southern_mid_30to55","northern_mid_30to55"),
)

def check(ok,msg):
    if not bool(ok):raise ValueError("MOIS_NATIVE_INTRO_9HUE_GEOGRAPHY: "+msg)

def sha(data):
    if isinstance(data,(str,Path)):data=Path(data).read_bytes()
    return hashlib.sha256(data).hexdigest()

def read_signed(path):
    check(sha(path)==SOURCE_SHA,"original WCVP native/introduced 9colour source ZIP SHA differs")
    with ZipFile(path) as z:
        receipt=json.loads(z.read("receipt.json"))
        for name,digest in receipt["output_sha256"].items():
            check(sha(z.read(name))==digest,
                  "original WCVP photo9 inner output SHA differs: "+name)
        d=pd.read_csv(io.BytesIO(z.read(SOURCE_MEMBER)),
                      compression="gzip",low_memory=False)
    check(len(d)==18577 and not d[["grid_x","grid_y"]].duplicated().any(),
          "18,577 geographic accepted-species source rows changed")
    check(receipt["source_37104_accepted_photo_hue_taxa"]==37104
          and receipt["source_full_GBIF_species_cell_pairs"]==5493931
          and receipt["source_null_199_same_global_genus_by_width_colour_label_permutations"] is True,
          "source original global photo9 accepted denominator or null drift")
    r=receipt["source_WCVP_exact_match_and_TDWG3_status"]
    check(r["cohorts_source_species_grid_pairs"]=={
          "strict_native":2271852,"moderate_native":2604870,
          "strict_introduced":436422,"moderate_introduced":522378},
          "genuine WCVP native+introduced per-TDWG3 source pair counts drift")
    check("center_latitude" in d and "sampled_land_median_elevation_m" in d,
          "original GSHHG WorldClim geography absent")
    return d,receipt

def sample(d,mode,min_species=30):
    check(mode in MODES and min_species in (30,100),
          "nonfrozen WCVP status or richness cutoff")
    a=f"{mode}_native";b=f"{mode}_introduced"
    nc=f"{a}_n_photo9_species";ic=f"{b}_n_photo9_species"
    use=d.loc[d[nc].ge(min_species) & d[ic].ge(min_species)].copy()
    check(not use[["grid_x","grid_y"]].duplicated().any(),
          "duplicate same-cell pair")
    if min_species==30:
        expected={"strict":3106,"moderate":3602}[mode]
        check(len(use)==expected,
              "original signed native/introduced same-grid cohort changed")
    else:
        check(len(use)>500,"No adequate native+intro >=100-species source support")
    use=use.reset_index(drop=True)
    use["block_500km"]=(use.grid_x.floordiv(5).astype(str)+"|"+
                        use.grid_y.floordiv(5).astype(str))
    use["block_1000km"]=(use.grid_x.floordiv(10).astype(str)+"|"+
                         use.grid_y.floordiv(10).astype(str))
    use["min_status_species"]=use[[nc,ic]].min(axis=1)
    use["native_minus_intro_observed_Shannon"]=(
        use[f"{a}_Shannon_observed"]-use[f"{b}_Shannon_observed"])
    use["native_minus_intro_expected_null_Shannon"]=(
        use[f"{a}_Shannon_null"]-use[f"{b}_Shannon_null"])
    use["native_minus_intro_9hue_null_anomaly"]=(
        use[f"{a}_Shannon_anomaly"]-use[f"{b}_Shannon_anomaly"])
    use["native_minus_intro_white_fraction_pp"]=100*(
        use[f"{a}_white_species_fraction"]-use[f"{b}_white_species_fraction"])
    check(np.allclose(
        use.native_minus_intro_observed_Shannon -
        use.native_minus_intro_expected_null_Shannon,
        use.native_minus_intro_9hue_null_anomaly,
        atol=1e-10,rtol=0),
        "source null-anomaly same-cell response does not exactly add")
    check(use[OUTCOMES].notna().all().all(),"missing paired ecological differences")
    check(use["min_status_species"].min()>=min_species,
          "native or introduced side lacked enough photographed source species")
    return use

def regions(d):
    lat=d.center_latitude.astype(float)
    elev=d.sampled_land_median_elevation_m.astype(float)
    land=d.gshhg_land_class
    return {
      "all_same_cells":np.ones(len(d),dtype=bool),
      "island_7of9":land.eq("island").to_numpy(),
      "mainland_7of9":land.eq("mainland").to_numpy(),
      "island_all9":(land.eq("island") & d.strict_all9_same_landmass).to_numpy(),
      "high_elevation_1500plus":elev.ge(1500).to_numpy(),
      "low_elevation_below500":elev.ge(-430).to_numpy() & elev.lt(500).to_numpy(),
      "very_high_elevation_2500plus":elev.ge(2500).to_numpy(),
      "arctic_66p562plus":lat.ge(66.562).to_numpy(),
      "boreal_50to66p562":lat.ge(50).to_numpy() & lat.lt(66.562).to_numpy(),
      "northern_hemisphere":lat.ge(0).to_numpy(),
      "southern_hemisphere":lat.lt(0).to_numpy(),
      "northern_mid_30to55":lat.ge(30).to_numpy() & lat.lt(55).to_numpy(),
      "southern_mid_30to55":lat.le(-30).to_numpy() & lat.gt(-55).to_numpy(),
      "north_tropical":lat.ge(0).to_numpy() & lat.lt(23.44).to_numpy(),
      "south_tropical":lat.lt(0).to_numpy() & lat.ge(-23.44).to_numpy(),
    }

def summarize_by_cluster(d,weighting,block_size):
    check(weighting in ("equal_cells","min_native_intro_species"),
          "unexpected WCVP status support weighting")
    key="block_500km" if block_size==500 else "block_1000km"
    w=(np.ones(len(d),float) if weighting=="equal_cells"
       else d.min_status_species.to_numpy(float))
    arr=d[list(OUTCOMES)].to_numpy(float)
    tmp=pd.DataFrame({"block":d[key].to_numpy(str),"weight":w})
    for k,col in enumerate(OUTCOMES):
        tmp[col]=arr[:,k]*w
    result=tmp.groupby("block",sort=False)[["weight",*OUTCOMES]].sum()
    return result

def bootstrap_groups(d,groupmask,weighting="equal_cells",km=500,
                     seed_offset=0,n_boot=N_BOOT):
    subset=d.loc[groupmask]
    if len(subset)<25:
        return {"status":"insufficient_25_source_cells","n_cells":int(len(subset))}
    block=summarize_by_cluster(subset,weighting,km)
    n=len(block)
    if n<8:
        return {"status":"insufficient_8_spatial_blocks",
                "n_cells":int(len(subset)),"n_blocks":n}
    arr=block.to_numpy(float)
    base=arr[:,1:].sum(axis=0)/arr[:,0].sum()
    rng=np.random.default_rng(SEED+seed_offset+len(subset))
    draws=np.empty((n_boot,len(OUTCOMES)))
    for start in range(0,n_boot,250):
        end=min(n_boot,start+250)
        idx=rng.integers(0,n,size=(end-start,n))
        sums=arr[idx].sum(axis=1)
        draws[start:end]=sums[:,1:]/sums[:,[0]]
    ci=np.percentile(draws,[2.5,97.5],axis=0)
    return {
        "status":"ok","n_cells":int(len(subset)),"n_blocks":n,
        "weighting":weighting,"spatial_block_km":km,
        "n_bootstrap":n_boot,
        "metrics":{name:{
            "mean":float(base[k]),
            "bootstrap_CI95":[float(ci[0,k]),float(ci[1,k])]}
                   for k,name in enumerate(OUTCOMES)}
    }

def bootstrap_geographic_contrast(d,mA,mB,weighting="equal_cells",
                                  km=500,seed_offset=0,n_boot=N_BOOT):
    a=d.loc[mA];b=d.loc[mB]
    if min(len(a),len(b))<25:
        return {"status":"insufficient_25_source_cells",
                "n_cells_A":len(a),"n_cells_B":len(b)}
    key="block_500km" if km==500 else "block_1000km"
    blockA=summarize_by_cluster(a,weighting,km)
    blockB=summarize_by_cluster(b,weighting,km)
    ids=blockA.index.union(blockB.index)
    x=blockA.reindex(ids,fill_value=0).to_numpy(float)
    y=blockB.reindex(ids,fill_value=0).to_numpy(float)
    nb=len(ids)
    if min(len(blockA),len(blockB))<8:
        return {"status":"insufficient_8_spatial_blocks"}
    means=x[:,1:].sum(axis=0)/x[:,0].sum()-y[:,1:].sum(axis=0)/y[:,0].sum()
    sims=np.full((n_boot,len(OUTCOMES)),np.nan)
    rng=np.random.default_rng(SEED+seed_offset+nb)
    for start in range(0,n_boot,250):
        end=min(n_boot,start+250)
        idx=rng.integers(0,nb,size=(end-start,nb))
        a_sum=x[idx].sum(axis=1);b_sum=y[idx].sum(axis=1)
        valid=a_sum[:,0].gt(0) if isinstance(a_sum,pd.DataFrame) else a_sum[:,0]>0
        valid &= b_sum[:,0]>0
        piece=np.full((end-start,len(OUTCOMES)),np.nan)
        piece[valid]=(a_sum[valid,1:]/a_sum[valid,[0]]
                     -b_sum[valid,1:]/b_sum[valid,[0]])
        sims[start:end]=piece
    check(np.isfinite(sims).mean()>.98,
          "Bootstrap regularly leaves out a geographic stratum")
    low=np.nanpercentile(sims,2.5,axis=0)
    high=np.nanpercentile(sims,97.5,axis=0)
    return {
       "status":"ok",
       "n_cells_A":len(a),"n_cells_B":len(b),
       "n_clusters_A":len(blockA),"n_clusters_B":len(blockB),
       "n_common_clusters":len(blockA.index.intersection(blockB.index)),
       "n_union_clusters":nb,"weighting":weighting,"spatial_block_km":km,
       "n_bootstrap":n_boot,
       "metrics":{name:{"A_minus_B":float(means[k]),
                        "bootstrap_CI95":[float(low[k]),float(high[k])]}
                  for k,name in enumerate(OUTCOMES)}
    }

def same_block_matched_contrast(d,mA,mB,weighting="equal_cells",km=500,
                                seed_offset=0,n_boot=N_BOOT):
    a=summarize_by_cluster(d.loc[mA],weighting,km)
    b=summarize_by_cluster(d.loc[mB],weighting,km)
    both=a.index.intersection(b.index)
    if len(both)<12:
        return {"status":"insufficient_within_block_matched_pairs",
                "n_shared_geographic_blocks":len(both)}
    a=a.loc[both].to_numpy(float);b=b.loc[both].to_numpy(float)
    # Regional average gives each genuinely shared 500km geographical block
    # equal weight, rather than favouring densely-recorded floras.
    diff=(a[:,1:]/a[:,:1]-b[:,1:]/b[:,:1])
    mean=diff.mean(axis=0)
    rng=np.random.default_rng(SEED+seed_offset+len(both))
    n=len(both)
    sims=np.empty((n_boot,len(OUTCOMES)))
    for start in range(0,n_boot,250):
        end=min(n_boot,start+250)
        idx=rng.integers(0,n,size=(end-start,n))
        sims[start:end]=diff[idx].mean(axis=1)
    ci=np.percentile(sims,[2.5,97.5],axis=0)
    return {
        "status":"ok","n_common_blocks":n,
        "weighting":weighting,"spatial_block_km":km,
        "n_bootstrap":n_boot,
        "within_shared_block_metrics":{name:{
          "mean_A_minus_B":float(mean[k]),
          "bootstrap_CI95":[float(ci[0,k]),float(ci[1,k])]}
          for k,name in enumerate(OUTCOMES)}
    }

def run(source,out):
    root=Path(out)
    check(not root.exists(),"Refuse overwriting source-verified result")
    src,rec=read_signed(source)
    analyses={}
    for mode in MODES:
        analyses[mode]={}
        for min_species in (30,100):
            d=sample(src,mode,min_species)
            mask=regions(d)
            name=f"min_species_each_status_{min_species}"
            output={
                "n_same_100km_cells":len(d),
                "n_source_500km_blocks":int(d.block_500km.nunique()),
                "n_source_1000km_blocks":int(d.block_1000km.nunique()),
                "weighting":{},
                "spatial_geographic_contrasts":{},
                "within_same_500km_block_contrasts":{}
            }
            for w in ("equal_cells","min_native_intro_species"):
                gr={}
                for n,m in mask.items():
                    res=bootstrap_groups(d,m,weighting=w,
                         seed_offset=(11 if w=="equal_cells" else 220))
                    gr[n]=res
                output["weighting"][w]=gr
                c={}
                for n,A,B in REGION_CONTRASTS:
                    c[n]=bootstrap_geographic_contrast(
                       d,mask[A],mask[B],weighting=w,
                       seed_offset=(23 if w=="equal_cells" else 300))
                output["spatial_geographic_contrasts"][w]=c
                if min_species==30:
                    within={}
                    for n,A,B in REGION_CONTRASTS:
                        within[n]=same_block_matched_contrast(
                            d,mask[A],mask[B],weighting=w,seed_offset=400)
                    output["within_same_500km_block_contrasts"][w]=within
            analyses[mode][name]=output
            first=output["weighting"]["equal_cells"]["all_same_cells"]
            check(first["status"]=="ok" and first["n_cells"]==len(d),
                  "Missing full-geographical paired source analysis")
            print("MOIS_NATIVE_INTRO_SAMECELL_PHOTO9_BOOTSTRAP",json.dumps({
               "mode":mode,"min_species_each":min_species,
               "n_cells":len(d),"n_500km_blocks":d.block_500km.nunique(),
               "mean_native_minus_intro_Shannon_anomaly":
                   first["metrics"]["native_minus_intro_9hue_null_anomaly"],
               "high_vs_low":output["spatial_geographic_contrasts"][
                   "equal_cells"]["high_minus_low_elevation"],
               "high_vs_low_within_same_block":
                   output["within_same_500km_block_contrasts"].get(
                   "equal_cells",{}).get("high_minus_low_elevation")
            },allow_nan=False),flush=True)
    root.mkdir(parents=True)
    a=root/"samecell_native_introduced_geographic_entropy_bootstrap.json"
    content={
       "analysis":"WCVP regional native vs introduced photo9 Shannon anomaly same real GBIF100km cells",
       "source_zip_sha256":SOURCE_SHA,
       "source_cell_rows":18577,
       "source_accepted_species":37104,
       "source_species_x_100km_grid":5493931,
       "n_original_global_genus_occupancywidth_shuffle":199,
       "bootstrap_method":"500km block resample keeping both WCVP statuses of same 100km cell paired",
       "n_bootstrap":N_BOOT,
       "strict_moderate_are_overlapping_sensitivities_not_replication":True,
       "geographic_strata_are_unmatched_except_within_shared_block_altitude_comparison":True,
       "null_expected_Shannon_fixed_from_original_199_shuffles_bootstrap_does_not_propagate_null_MC_error":True,
       "photo_hue_not_organ_verified_petal_colour":True,
       "WCVP_TDWG3_native_not_wild_individual_or_small_island_native_status":True,
       "areas_and_altitude_are_100km_grid_majority_or_median":True,
       "results":analyses,
    }
    a.write_text(json.dumps(content,indent=2,allow_nan=False)+"\n")
    rows=[]
    for mode,cohorts in analyses.items():
        for min_n,info in cohorts.items():
            for w,regions in info["weighting"].items():
                for geo,res in regions.items():
                    if res["status"]!="ok":continue
                    row={
                       "WCVP_mode":mode,"cohort":min_n,"weight":w,
                       "region":geo,"n_cells":res["n_cells"],
                       "n_500km_blocks":res["n_blocks"]
                    }
                    for key,v in res["metrics"].items():
                        row[key]=v["mean"]
                        row[key+"_ci_low"]=v["bootstrap_CI95"][0]
                        row[key+"_ci_high"]=v["bootstrap_CI95"][1]
                    rows.append(row)
    csv=root/"samecell_native_intro_9hue_geographic_status_strata.csv"
    pd.DataFrame(rows).to_csv(csv,index=False)
    content["numeric_output_sha256"]={p.name:sha(p) for p in (a,csv)}
    (root/"receipt.json").write_text(json.dumps(content,indent=2,allow_nan=False)+"\n")
    print("MOIS_NATIVE_INTRO_9HUE_GEOGRAPHIC_BLOCK_BOOTSTRAP_REAL_SUCCESS",
          json.dumps({
            "strict_30_global":analyses["strict"]["min_species_each_status_30"][
                 "weighting"]["equal_cells"]["all_same_cells"],
            "strict_30_high_low_500km":analyses["strict"][
                 "min_species_each_status_30"]["within_same_500km_block_contrasts"][
                 "equal_cells"]["high_minus_low_elevation"]
          },allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",required=True)
    p.add_argument("--out",required=True)
    a=p.parse_args()
    run(a.source,a.out)
