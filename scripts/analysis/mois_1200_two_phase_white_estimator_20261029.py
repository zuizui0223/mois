#!/usr/bin/env python3
"""Use *completed* independent 1200 human corolla labels for design-based E1.

No classification outcome is imputed from failed ROI9. An incomplete review
returns a BLOCKED receipt without white-colour estimates. Fully reviewed but
ambiguous flower/colour records produce honest partial-identification bounds.
The target is the *frozen measured four-arm represented-observation frame*,
not the 93m source, native flora or abundance-weighted flowering community.
"""
from __future__ import annotations
import argparse
from collections import Counter
import json
import math
from pathlib import Path

import pandas as pd

from scripts.analysis.mois_observation_probability_audit_20261029 import (
    SCHEMA as AUDIT_SCHEMA, require, sha_file
)

SCHEMA="mois_human1200_two_phase_observation_white_frequency_estimator_v1"
COLOURED=("yellow_orange","red_pink","blue_purple","other")
AMBIGUOUS_COLOUR=("mixed","uncertain")
FLOWER=("yes","no","uncertain","inaccessible")

def compute(roster_path, source_receipt, annotation_path, output):
    roster_path,source_receipt,annotation_path,output=map(
        Path,(roster_path,source_receipt,annotation_path,output))
    require(not output.exists(), "Do not overwrite existing estimator result")
    freeze=json.loads(source_receipt.read_text())
    require(freeze.get("schema")==AUDIT_SCHEMA and freeze.get("sample_total")==1200,
            "Wrong frozen observation-level audit provenance")
    require(freeze["output_sha256"]["audit1200_source_observation_roster.csv"]==
            sha_file(roster_path), "Frozen original 1200 roster SHA drifted")
    d=pd.read_csv(roster_path,dtype={"photo_id":str,"observation_id":str})
    r=pd.read_csv(annotation_path,dtype={"photo_id":str,"observation_id":str},
                  keep_default_na=False)
    needed=("photo_id","observation_id","human_review_completed",
            "human_open_corolla_visible","human_corolla_colour")
    require(all(c in d.columns for c in ("photo_id","observation_id","stratum",
        "stratum_N","stratum_n","audit_inclusion_q","review_inverse_q")) and
            all(c in r.columns for c in needed),
            "Missing required observation roster or human annotation columns")
    require(len(d)==1200 and d.photo_id.is_unique and d.observation_id.is_unique,
            "Original audit frame must be 1200 distinct observations")
    require(r.photo_id.is_unique and r.observation_id.is_unique,
            "Human annotations have duplicate photo/observation IDs")
    if "review_role" in r.columns:
        require(r.review_role.eq("R1").all(),
                "Only primary blinded reviewer R1 may enter population estimator")
    require(set(r.photo_id)==set(d.photo_id),"Annotation set must be exactly 1200 IDs")
    require((d.audit_inclusion_q>0).all() and (d.audit_inclusion_q<=1).all() and
            ((d.stratum_n/d.stratum_N-d.audit_inclusion_q).abs()<1e-11).all(),
            "Invalid frozen review selection probabilities")
    require(((d.review_inverse_q*d.audit_inclusion_q-1).abs()<1e-9).all(),
            "Review weights do not invert frozen selection probabilities")
    for stratum,g in d.groupby("stratum"):
        require(len(g)==g.stratum_n.iloc[0] and
                g.stratum_n.nunique()==1 and g.stratum_N.nunique()==1,
                "Stratum sample size mismatched")
    a=d.merge(r,on=["photo_id","observation_id"],
              how="left",validate="one_to_one",sort=False)
    marked=a.human_review_completed.astype(str).str.strip().str.lower().isin(
        ("true","yes","1"))
    pending=int((~marked).sum())
    result={
        "schema":SCHEMA,
        "frozen_1200_source_roster_SHA256":sha_file(roster_path),
        "human_annotation_input_SHA256":sha_file(annotation_path),
        "n_source_observations_in_audit":1200,
        "n_review_complete":int(marked.sum()),"n_review_pending":pending,
        "reference_population":"four-arm measured-photo represented observations only",
        "natural_flora_white_frequency_identified":False,
        "not_corrected_to_93m_or_global_native_plants":True,
    }
    if pending:
        result["status"]="BLOCKED_PENDING_HUMAN_REVIEW"
        result["design_corrected_white_fraction"]=None
        result["note"]="Do not use complete-case reviewed subset; it can bias the weighted result."
    else:
        require(a.human_open_corolla_visible.isin(FLOWER).all(),
                "Open-flower human review contains invalid/blank category")
        yes=a.human_open_corolla_visible.eq("yes")
        no=a.human_open_corolla_visible.eq("no")
        c=a.human_corolla_colour
        require(c[yes].isin(("white",)+COLOURED+AMBIGUOUS_COLOUR).all() and
                c[no].isin(("","not_applicable")).all(),
                "Corolla colour inconsistent with the independently labeled flower")
        w=a.review_inverse_q.astype(float)
        white=yes & c.eq("white")
        coloured=yes & c.isin(COLOURED)
        unknown_colour=yes & c.isin(AMBIGUOUS_COLOUR)
        unknown_flower=~yes & ~no
        W=float(w[white].sum());C=float(w[coloured].sum())
        U=float(w[unknown_colour].sum());V=float(w[unknown_flower].sum())
        require(W+C+U+V>0,"No reviewed or possible flowers to estimate")
        result["weighted_open_white"]=W
        result["weighted_open_nonwhite"]=C
        result["weighted_open_ambiguous_colour"]=U
        result["weighted_unknown_flower_status"]=V
        result["flower_review_counts"]={
            key:int(n) for key,n in
            a.human_open_corolla_visible.value_counts().items()
        }
        result["white_frequency_partial_bounds"]={
            "lower":W/(W+C+U+V),
            "upper":(W+U+V)/(W+C+U+V),
        }
        if U+V>0:
            result["status"]="PARTIALLY_IDENTIFIED_AFTER_COMPLETE_REVIEW"
            result["design_corrected_white_fraction"]=None
            result["note"]="Uncertain corolla/colour; bounds condition on all missing potentially open flowers. No point prevalence imputed."
        else:
            p=W/(W+C)
            z=(white.astype(int)-p*yes.astype(int)).to_numpy(float)
            a=a.copy()
            a["_linearized_z"]=z
            total_var=0.0
            for _,g in a.groupby("stratum",sort=True):
                n=int(g.stratum_n.iloc[0])
                N=int(g.stratum_N.iloc[0])
                if n==N:
                    continue
                require(n>=2,"Non-census singletons cannot estimate sampling variance")
                s2=float(g._linearized_z.var(ddof=1))
                total_var+=(N*N*(1-n/N)*s2/n)
            se=math.sqrt(max(0.0,total_var))/(W+C)
            result["status"]="DESIGN_BASED_FULLY_OBSERVED_HUMAN_E1"
            result["design_corrected_white_fraction"]=p
            result["taylor_linearization_standard_error"]=se
            result["approximate_normal_95_interval"]=[
                max(0.,p-1.96*se),min(1.,p+1.96*se)
            ]
            result["note"]="First-stage photo inclusion is not corrected; approximate normal interval conditions on measured observation frame and independently coded open corollas."
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("MOIS_HUMAN1200_BIAS_AUDIT_ESTIMATOR_GATE",
          json.dumps({"status":result["status"],
                      "pending":pending,
                      "design_corrected_white_fraction":
                        result.get("design_corrected_white_fraction")},sort_keys=True),
          flush=True)
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--roster",type=Path,required=True)
    p.add_argument("--source-receipt",type=Path,required=True)
    p.add_argument("--human-annotations",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    compute(a.roster,a.source_receipt,a.human_annotations,a.output)

if __name__=="__main__":
    main()
