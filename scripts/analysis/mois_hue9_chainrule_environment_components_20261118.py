#!/usr/bin/env python3
"""MOIS: trace ecological rainfall relationships into mathematically exact
white/nonwhite and nonwhite-colour Shannon community components.

Use two separate SHA-verified real-data CI artifacts (2026-11-18) for:
  1. 12,295 GBIF100km environmental complete-case grid observations
  2. Exact 18,577-grid Shannon H9 chain-rule decomposition from the
     SAME 199 global genus-by-range-size hue-label permutations.
Use the signed original 5,493,931 accepted species x GBIF grid table to
recover the authentic 396,142,825 record-effort denominator, ensuring the
fitted environmental X and 500km FE are IDENTICAL to the earlier real model.

This is an arithmetic decomposition of recorded species-colour photo labels,
not a causal mediation, true flower petal measurement, or palaeorefugium test.
"""
from __future__ import annotations
import argparse,hashlib,io,json
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd

from scripts.analysis.mois_hue9_assembly_environment_spatial_20261118 import (
    original_GBIF_record_effort,spatial_fit,ECOLOGY,require,sha
)

ENV_ZIP_SHA="4f5d8fc3a308482fc3645d796a779d1a4725c51ed2afbf6a06e61fd6213f3ee6"
CHAIN_ZIP_SHA="4b76d95df58d794a70338b71d61ee91a0e49a6c2c67bb721e806822000e0f804"
ENV_CSV="all_current_GBIF_100km_cells_genus_breadth_entropy_ecology.csv.gz"
CHAIN_CSV="same_cell_exact_white_nonwhite_chromatic8_entropy_decomposition.csv.gz"
RESPONSE={
    "total_ninecolour":"entropy_anomaly",
    "white_nonwhite_binary":"white_nonwhite_binary_anomaly",
    "weighted_chromatic8":"nonwhite_weighted_colour8_anomaly"
}

def artifact_csv(path,expected,source_name,checksum_field):
    require(sha(path)==expected,source_name+" signed source ZIP SHA changed")
    with ZipFile(path) as z:
        rec=json.loads(z.read("receipt.json"))
        for fname,digest in rec[checksum_field].items():
            require(sha(z.read(fname))==digest,
                    source_name+" internal numeric output SHA differs: "+fname)
        frame=pd.read_csv(io.BytesIO(z.read(
            ENV_CSV if source_name=="environment" else CHAIN_CSV)),
            compression="gzip")
    require(not frame[["grid_x","grid_y"]].duplicated().any(),
            source_name+" original grid has duplicate cells")
    return frame,rec

def assemble(environment,chain,gbif):
    require(len(environment)==12295 and len(chain)==18577,
            "Signed original complete-case or chainrule geographic sample drift")
    fields=["grid_x","grid_y","n_recorded_source_species",
            "original_white_fraction","total_ninecolour_anomaly",
            "white_nonwhite_binary_anomaly","nonwhite_weighted_colour8_anomaly"]
    require(set(fields).issubset(chain),"Chain-rule source missing components")
    d=environment.merge(chain[fields],on=["grid_x","grid_y"],how="left",
                        validate="one_to_one")
    require(len(d)==len(environment),"Original environmental cohort altered")
    require(d.total_ninecolour_anomaly.notna().all(),
            "Environmental cells missing Shannon components")
    require((d.n_accepted_species==d.n_recorded_source_species).all(),
            "Genus null source species numbers differ")
    require(np.allclose(d.white_species_fraction,d.original_white_fraction,
                        atol=1e-12,rtol=0),
            "Photo white source fraction changed in later arithmetic audit")
    require(np.allclose(d.entropy_anomaly,d.total_ninecolour_anomaly,
                        atol=1e-10,rtol=0),
            "Source nine-class environment response != exact Shannon decomposition")
    require(np.allclose(
        d.entropy_anomaly,
        d.white_nonwhite_binary_anomaly+d.nonwhite_weighted_colour8_anomaly,
        atol=1e-10,rtol=0),"GBIF cell Shannon components fail exact additive identity")
    effort=original_GBIF_record_effort(gbif)
    d=d.merge(effort,on=["grid_x","grid_y"],validate="one_to_one")
    require(len(d)==len(environment) and d.GBIF_occurrence_count.gt(0).all(),
            "True GBIF observation-effort merge does not preserve cells")
    return d

def check_total_reconstruction(previous,fit):
    old=previous["models"]["primary_500km_richness30plus"]
    require(old["n_fit_cells"]==fit["n_fit_cells"]
            and old["n_spatial_clusters"]==fit["n_spatial_clusters"],
            "Refitted model has different cells or 500km blocks")
    for term,p in old["coefficients"].items():
        q=fit["coefficients"][term]
        require(abs(p["beta_Shannon_per_global_SD"]-
                    q["beta_Shannon_per_global_SD"])<1e-9,
                "Old spatial ecosystem coefficient not reproduced: "+term)
    return True

def run(args):
    root=Path(args.out)
    require(not root.exists(),"Refuse overwriting verified arithmetic ecology results")
    original,env_receipt=artifact_csv(args.environment,ENV_ZIP_SHA,
                                      "environment","numeric_outputs_sha256")
    decomposition,chain_receipt=artifact_csv(args.chain,CHAIN_ZIP_SHA,
                                             "chainrule","output_sha256")
    require(chain_receipt["n_global_shared_permutations"]==199
            and chain_receipt["previous_9hue_expected_entropy_exactly_reproduced"]
            and env_receipt["source_cohort"]["strict_complete_30plus_cells"]==12295,
            "Original global permutation and 500km cohort audit fails")
    d=assemble(original,decomposition,args.GBIF)
    models={
        "primary_500km_richness30plus":{},
        "richness_ge100_500km":{},
        "1000km_spatial_FE":{}
    }
    for target,response in RESPONSE.items():
        models["primary_500km_richness30plus"][target]=spatial_fit(d,response=response)
        d100=d[d.n_accepted_species.ge(100)].reset_index(drop=True)
        models["richness_ge100_500km"][target]=spatial_fit(d100,response=response)
        models["1000km_spatial_FE"][target]=spatial_fit(
            d,response=response,block_km=1000)
    require(check_total_reconstruction(
          env_receipt,models["primary_500km_richness30plus"]["total_ninecolour"]),
          "Old total Shannon environmental association cannot be replicated")
    deviations={}
    for label,sets in models.items():
        a=sets["total_ninecolour"]["coefficients"]
        b=sets["white_nonwhite_binary"]["coefficients"]
        c=sets["weighted_chromatic8"]["coefficients"]
        devs={term:(a[term]["beta_Shannon_per_global_SD"]-
                    b[term]["beta_Shannon_per_global_SD"]-
                    c[term]["beta_Shannon_per_global_SD"]) for term in a}
        require(max(abs(v) for v in devs.values())<1e-10,
                "Environmental beta must exactly decompose at same GBIF cells")
        deviations[label]=float(max(map(abs,devs.values())))
    root.mkdir(parents=True)
    records=[]
    for cohort,targets in models.items():
        for label,r in targets.items():
            for term,e in r["coefficients"].items():
                records.append({
                    "cohort":cohort,"Shannon_component":label,
                    "term":term,"n_cells":r["n_fit_cells"],
                    "n_spatial_blocks":r["n_spatial_clusters"],
                    "beta_per_SD":e["beta_Shannon_per_global_SD"],
                    "CI95_low":e["CI95"][0],"CI95_high":e["CI95"][1],
                    "p_cluster":e["p_cluster"],
                    "q7_BH":e.get("q_BH_seven_ecological_exposures"),
                })
    numeric=root/"same_cell_exact_Shannon_component_environment_effects.csv"
    pd.DataFrame(records).to_csv(numeric,index=False)
    summary={
        "source":{
           "signed_environment_zip_SHA256":ENV_ZIP_SHA,
           "signed_chainrule_zip_SHA256":CHAIN_ZIP_SHA,
           "n_same_GBIF_cells":12295,
           "from_37104_accepted_photo_hue_species":True,
           "original_100km_species_grid_pairs":5493931,
           "original_199_genus_by_breadth_global_permutations":True,
           "original_GBIF_396M_record_effort_recovered":True,
           "original_500km_total_model_coefficients_exactly_reproduced":True,
        },
        "models":models,"additivity_max_abs_beta_error":deviations,
        "Shannon_identity":"H9=H_white_nonwhite + P_nonwhite*H8_nonwhite",
        "component_associations_not_causal_mediation":True,
        "source_flower_photo_hue_not_true_local_petal_color":True,
        "past_LGM_precip_at_present_GBIF_sites_not_refugia":True,
    }
    jf=root/"ecological_component_environment_model.json"
    jf.write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    summary["output_sha256"]={v.name:sha(v) for v in (numeric,jf)}
    (root/"receipt.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    print("MOIS_HUE9_EXACT_COMPONENT_ECOLOGY_REAL_SUCCESS",
          json.dumps({
           "n":models["primary_500km_richness30plus"][
                 "total_ninecolour"]["n_fit_cells"],
           "rainfall_effects":{
              term:{target:
                   models["primary_500km_richness30plus"][target][
                       "coefficients"][term]
                   for target in RESPONSE}
              for term in ("log_precip","LGM_precip_delta","NPP")},
           "reproduction":deviations,
          },allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--environment",required=True)
    p.add_argument("--chain",required=True)
    p.add_argument("--GBIF",required=True)
    p.add_argument("--out",required=True)
    run(p.parse_args())
