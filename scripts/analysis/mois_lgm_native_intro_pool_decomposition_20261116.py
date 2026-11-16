#!/usr/bin/env python3
"""GBIF-first decomposition of WCVP native-introduced hue differences by genus pool.

Given the SAME present-day GBIF 100km cell, separate the photographic
white-label difference into: hue contrast inside genera present on both
sides, differential representation of those shared genera, and unmatched-
genus (exclusive species pool) contribution. Sum must equal original
all-genus difference per cell AND its spatially adjusted coefficient.

This is a descriptive Kitagawa-style arithmetic decomposition, not a
causal mediation / glacial refugium / local botanical petal-colour test.
The colour trait is species-modal photographic ROI, inherited into each
current GBIF species-by-grid presence. WCVP nativeness is TDWG3 level only.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd
import statsmodels.api as sm

from scripts.analysis.mois_lgm_native_introduced_matched_cells_20261116 import (
    signed_tables, same_cell_cohort, matrix, sha, require
)

COMMON_SHA = "ba725090afbd732c9dacc8e098161aad2bf88042994f5c683db465519969356a"
COMMON_FILE = "samecell_native_intro_shared_genera_100km_summary.parquet"
EXPECTED_FITTED = {
    "strict": 0.003001374579903353,
    "moderate": 0.001192127107729269,
}
REQUIRED = ("n_shared_genera", "shared_native_species_n",
            "shared_intro_species_n", "shared_native_white_n",
            "shared_intro_white_n", "genus_equal_native_minus_intro")

def signed_common(archive):
    import pyarrow.parquet as pq
    require(sha(archive) == COMMON_SHA, "signed common-genus ZIP changed")
    with ZipFile(archive) as z:
        raw = z.read(COMMON_FILE)
        receipt = json.loads(z.read("receipt.json"))
    require(sha(raw) == receipt["output_sha256"][COMMON_FILE],
            "common-genus source parquet SHA differs")
    result = receipt["results"]
    require(result["source_GBIF_distinct_species"] == 37104
            and result["source_GBIF_species_cell_pairs"] == 5493931
            and result["LGM_warming_is_current_cell_paleoclimate_difference_not_refugia"],
            "common-genus provenance or scientific limitations changed")
    frame = pq.read_table(io.BytesIO(raw)).to_pandas()
    require(not frame[["grid_x", "grid_y"]].duplicated().any(),
            "duplicate signed common genus GBIF grid")
    require(set(REQUIRED).issubset(frame.columns),
            "needed common genus provenance columns absent")
    return frame, receipt

def decompose_cells(d):
    """Three terms exactly reconstruct all-genera native-minus-introduced hue.

    For each cell, let sn and si be fractions of native and introduced species
    belonging to genera present on BOTH sides. Let hn and hi be white fractions
    among those species. Then:
      overall delta = (sn+si)/2*(hn-hi)
                    + (sn-si)/2*(hn+hi)
                    + unique_native_white/total_native
                    - unique_intro_white/total_intro
    All genus-unmatched species (whether from genuinely exclusive genera
    or genus matches not resolved) are excluded from the common-genus pool.
    """
    require(d.native_species_n.ge(1).all()
            and d.introduced_species_n.ge(1).all(),
            "requires both WCVP status groups")
    n = d.native_species_n.astype(float)
    i = d.introduced_species_n.astype(float)
    sn = d.shared_native_species_n.astype(float)
    si = d.shared_intro_species_n.astype(float)
    snw = d.shared_native_white_n.astype(float)
    siw = d.shared_intro_white_n.astype(float)
    nw = d.native_white_n.astype(float)
    iw = d.introduced_white_n.astype(float)
    require((sn > 0).all() and (si > 0).all()
            and (sn <= n).all() and (si <= i).all()
            and (snw <= sn).all() and (siw <= si).all()
            and (nw >= snw).all() and (iw >= siw).all()
            and (nw-snw <= n-sn+1e-9).all()
            and (iw-siw <= i-si+1e-9).all(),
            "negative or impossible species-pool provenance")
    alpha_n, alpha_i = sn/n, si/i
    hn, hi = snw/sn, siw/si
    z = d.copy()
    z["native_species_share_in_common_genera"] = alpha_n
    z["introduced_species_share_in_common_genera"] = alpha_i
    z["unmatched_native_species_n"] = n-sn
    z["unmatched_introduced_species_n"] = i-si
    z["contrast_inside_common_genus_pool"] = (alpha_n+alpha_i)/2*(hn-hi)
    z["representation_of_common_genera"] = (alpha_n-alpha_i)/2*(hn+hi)
    z["exclusive_genus_pool_white_contribution"] = (nw-snw)/n-(iw-siw)/i
    z["reconstructed_white_difference"] = (
        z["contrast_inside_common_genus_pool"]
        +z["representation_of_common_genera"]
        +z["exclusive_genus_pool_white_contribution"])
    require(np.allclose(z.reconstructed_white_difference,
                        z.native_minus_intro_white, atol=1e-12, rtol=0),
            "cell-level white native minus introduced did not reconstruct")
    return z

COMPONENTS = ("native_minus_intro_white",
              "contrast_inside_common_genus_pool",
              "representation_of_common_genera",
              "exclusive_genus_pool_white_contribution",
              "genus_equal_native_minus_intro",
              "native_species_share_in_common_genera",
              "introduced_species_share_in_common_genera")

def fit_shared_design(d, block="500"):
    """One X and one fixed eligible cell set for ALL components (additivity)."""
    x, sd = matrix(d)
    b = d.block_500km.astype(str) if block=="500" else d.block_1000km.astype(str)
    x = x - x.groupby(b).transform("mean")
    usable = x.abs().sum(axis=1).gt(0)
    x = x.loc[usable]
    group = d.loc[usable, "block_500km"].astype(str)
    require(len(x)>300 and group.nunique()>40,
            "insufficient independent geography after fixed block")
    r = {}
    for name in COMPONENTS:
        y = d.loc[usable, name].astype(float).copy()
        y = y - y.groupby(b.loc[usable]).transform("mean")
        model = sm.OLS(y, x).fit(cov_type="cluster",
                                 cov_kwds={"groups":group}, use_t=False)
        low,high = model.conf_int().loc["LGM_warming"].tolist()
        r[name] = {
            "beta_fraction_per_SD": float(model.params["LGM_warming"]),
            "95CI": [float(low),float(high)],
            "p": float(model.pvalues["LGM_warming"]),
            "n_cells": len(x), "n_500km_clusters": int(group.nunique()),
            "R2": float(model.rsquared)
        }
    reconstructed=sum(r[name]["beta_fraction_per_SD"] for name in
                      COMPONENTS[1:4])
    actual=r["native_minus_intro_white"]["beta_fraction_per_SD"]
    require(abs(reconstructed-actual)<1e-10,
            "model-level decomposition coefficient not additive")
    return {
        "component_models":r,
        "predictor_LGM_warming_SD_C":sd,
        "block_fixed_effect_km":int(block),
        "standard_errors_clustered_500km":True,
        "additivity_check_max_abs":abs(reconstructed-actual),
        "mean_genus_shared_native_species_share":float(
            d.native_species_share_in_common_genera.mean()),
        "mean_genus_shared_intro_species_share":float(
            d.introduced_species_share_in_common_genera.mean()),
        "median_shared_genera":float(d.n_shared_genera.median()),
        "min_10_each_shared_species":True,
        "status":"fitted"
    }

def run(native_zip,spatial_zip,common_zip,output):
    output=Path(output)
    require(not output.exists(),"refuse overwriting signed analysis")
    native,spatial=signed_tables(native_zip,spatial_zip)
    common,receipt=signed_common(common_zip)
    results={}
    out_tables=[]
    for mode in ("strict","moderate"):
        # Exact pre-existing GBIF same-cell population, not a fresh cohort.
        raw=same_cell_cohort(native,spatial,mode,min_each=30)
        d=raw.merge(common,on=["grid_x","grid_y"],
                    how="left",validate="one_to_one")
        require(not d[["grid_x","grid_y"]].duplicated().any(),
                "duplicated same-cell accepted species universe")
        d=d.loc[(d.n_shared_genera>=5)
                &(d.shared_native_species_n>=10)
                &(d.shared_intro_species_n>=10)].copy()
        require(len(d)==({"strict":2516,"moderate":2926}[mode]),
                "species-genus sample changed from signed common-genus run")
        d=decompose_cells(d)
        tests={str(k): fit_shared_design(d,str(k)) for k in (500,1000)}
        # Check same cohort, design and outcome against independent source run.
        observed=tests["500"]["component_models"][
            "genus_equal_native_minus_intro"]["beta_fraction_per_SD"]
        require(abs(observed-EXPECTED_FITTED[mode])<1e-8,
                "current genus-equal coefficient differs from prior independent result")
        results[mode]={
            "source_matched_cells":int(len(raw)),
            "common_genus_eligible_cells":int(len(d)),
            "n_cells_shared_common_genera":int(d.n_shared_genera.ge(5).sum()),
            "n_unique_shared_genera_cell_combinations":
                int(d.n_shared_genera.sum()),
            "summary":{
               "native_shared_species_fraction_mean":float(
                   d.native_species_share_in_common_genera.mean()),
               "intro_shared_species_fraction_mean":float(
                   d.introduced_species_share_in_common_genera.mean()),
               "mean_original_white_fraction_diff":float(
                   d.native_minus_intro_white.mean()),
               "mean_genus_equal_white_fraction_diff":float(
                   d.genus_equal_native_minus_intro.mean()),
               "mean_unique_genus_white_component":float(
                   d.exclusive_genus_pool_white_contribution.mean())
            },
            "500km_models":tests["500"],
            "1000km_models":tests["1000"]
        }
        tab=d[["grid_x","grid_y","n_shared_genera","native_species_n",
               "introduced_species_n","native_white_n","introduced_white_n",
               "native_species_share_in_common_genera",
               "introduced_species_share_in_common_genera",
               *COMPONENTS]].copy()
        tab["region_assignment_mode"]=mode
        out_tables.append(tab)
        core=tests["500"]["component_models"]
        print("MOIS_GBIF_LGM_POOL_DECOMPOSITION",json.dumps({
            "mode":mode,"n":len(d),"n_shared_genus_cell_pairs":int(
                d.n_shared_genera.sum()),
            "all_beta":core["native_minus_intro_white"]["beta_fraction_per_SD"],
            "within_shared_pool_beta":core["contrast_inside_common_genus_pool"]["beta_fraction_per_SD"],
            "genus_representation_beta":core["representation_of_common_genera"]["beta_fraction_per_SD"],
            "exclusive_genus_beta":core["exclusive_genus_pool_white_contribution"]["beta_fraction_per_SD"],
            "common_genus_equal_beta":core["genus_equal_native_minus_intro"]["beta_fraction_per_SD"]
        },sort_keys=True),flush=True)
    output.mkdir(parents=True)
    p=output/"GBIF100km_WCVP_native_intro_shared_genus_exact_decomposition.parquet"
    pd.concat(out_tables,ignore_index=True).to_parquet(p,index=False)
    q=output/"models.json"
    q.write_text(json.dumps(results,indent=2,sort_keys=True)+"\n")
    report={
      "source_GBIF_accepted_species":37104,
      "source_GBIF_species_100km_pairs":5493931,
      "source_current_WCVP_status_is_broad_TDWG3_only":True,
      "source_flower_colour_is_species_photo_modal_proxy":True,
      "source_past_LGM_range_refugia_identified":False,
      "source_zip_SHA256":{
        "regional_native":sha(native_zip),
        "LGM_100km_spatial":sha(spatial_zip),
        "common_genus_paired":sha(common_zip)},
      "decomposition_algebra_exact":True,
      "paired_white_difference_explained_by_causal_taxonomic_history":False,
      "results":results,"output_sha256":{k.name:sha(k) for k in (p,q)}
    }
    (output/"receipt.json").write_text(
        json.dumps(report,indent=2,sort_keys=True)+"\n")
    return report

if __name__=="__main__":
    p=argparse.ArgumentParser()
    for k in ("native","spatial","common","output"):
        p.add_argument("--"+k,required=True)
    a=p.parse_args()
    run(a.native,a.spatial,a.common,a.output)
