#!/usr/bin/env python3
"""Independent numerical audit of original source-verified European DoPI v2
flower visitor guild × MOIS photo nine-hue results.

The parent workflow already SHA-verified 676,900 source raw rows and the
accepted MOIS photo9 source, and exported 90 anonymous positive study-guild
comparisons. This deliberately re-reads only the SHA-pinned aggregate
artifact, NOT a new source of positive records or a false host-absence table.

For every guild pair × 9 photo hues × two source scopes, compare:
(1) each publication's POSITIVE guild-A plant set vs POSITIVE guild-B plant
    set, with per-publication species frequency contrast; and
(2) the study-local accepted genera that have at least one POSITIVE
    species in both guild sets, equal-genus differences (n>=2 common genera).
The latter has different eligibility and cannot be subtracted and called
causal genus mediation. BH multiplicity across all supported tests, not the
most favourable hue. Source grade 3 is visitation ONLY, not confirmed
pollination. Missing source guild is unknown, not negative.
"""
from __future__ import annotations
import argparse,hashlib,io,json
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

PARENT_SHA="bf7dd9b41e8b71813a430f0a917ddfc949ed12156dfe833c91862c13f6dc7595"
TABLE="within_publication_positive_guild_pair_nine_hue_tests.csv"
PARENT_SUMMARY="source_and_inferential_limits.json"
FAMILY=90
SUPPORTED=63
N_SIG=10

def sha(b):
    if isinstance(b,(str,Path)):b=Path(b).read_bytes()
    return hashlib.sha256(b).hexdigest()
def require(ok,msg):
    if not bool(ok):raise ValueError("MOIS_POSITIVE_DOPI_GENUS_AUDIT: "+str(msg))

def read(p):
    require(sha(p)==PARENT_SHA,"Original 676900-row DoPI paired study test numeric ZIP changed")
    with ZipFile(p) as z:
        receipt=json.loads(z.read("receipt.json"))
        for name,checksum in receipt["output_SHA256"].items():
            require(sha(z.read(name))==checksum,
                    "Original DoPI numerical output unexpectedly changed "+name)
        df=pd.read_csv(io.BytesIO(z.read(TABLE)))
        summary=json.loads(z.read(PARENT_SUMMARY))
    require(len(df)==FAMILY and
            int(df.eligible_for_exploratory_test.sum())==SUPPORTED,
            "Original supported study×colour universe changed")
    require(summary["all_source_row_screen_counts"]["raw_source_rows"]==676900
            and summary["all_source_row_screen_counts"][
             "MOIS_matched_valid_positive_rows"]==428871,
            "Source-authenticated DoPI positive guild coverage changed")
    require(df["n_same_study_shared_genus_panels"].ge(0).all(),
            "Shared genus source study counts absent")
    return df,summary

def audit(df):
    valid=df.source_publication_t_p.notna()
    gen=df.within_publication_shared_genus_t_p.notna()
    require(valid.sum()==SUPPORTED and gen.sum()==SUPPORTED,
            "Positive-guild direct and common-genus study support not aligned")
    qs=false_discovery_control(df.loc[valid,"source_publication_t_p"].to_numpy(float),
                               method="bh")
    qg=false_discovery_control(
       df.loc[gen,"within_publication_shared_genus_t_p"].to_numpy(float),
       method="bh")
    df.loc[valid,"independent_BH_q_source_published_study"]=qs
    df.loc[gen,"independent_BH_q_shared_genus_study"]=qg
    require(np.max(np.abs(df.loc[valid,"BH_q_among_supported_90_screen"].to_numpy(float)-
                              qs))<1e-12,
            "Original multiple testing q did not independently reproduce")
    significant=df.independent_BH_q_source_published_study.lt(.05)
    common_sig=df.independent_BH_q_shared_genus_study.lt(.05)
    require(int(significant.sum())==N_SIG,
            "Original nominal colour exploratory positives changed")
    require(int(common_sig.sum())==0,
            "Common-genus FDR signal changed; recalibrate interpretation")
    require((df.loc[significant,"evidence_scope"]==
             "positive_visit_or_pollen_grade123").all(),
            "Any source-level positives now documented as pollen evidence alone")
    selected=df.loc[significant,[
       "evidence_scope","guild_A_positive_plants","guild_B_positive_plants",
       "photo9_hue","n_same_published_studies_both_positive_guilds",
       "n_same_study_shared_genus_panels",
       "positive_guildA_minus_guildB_photo_hue_fraction",
       "bootstrap95_by_publication_low","bootstrap95_by_publication_high",
       "source_publication_t_p","independent_BH_q_source_published_study",
       "within_publication_shared_genus_equal_effect",
       "within_publication_shared_genus_t_p","independent_BH_q_shared_genus_study"
    ]].copy().sort_values(["guild_A_positive_plants",
                            "guild_B_positive_plants","photo9_hue"])
    summary={
      "signed_source_zip_SHA256":PARENT_SHA,
      "prior_original_DoPI_source_rows":676900,
      "prior_true_positive_MOIS_matched_guild_rows":428871,
      "all_planned_test_panels":FAMILY,
      "both_scope_supported_test_panels":SUPPORTED,
      "photo9_BH_supported_publication_positives":N_SIG,
      "all_10_publication_positives_only_grade123_including_visit_only":True,
      "same_study_shared_genus_BH_supported":int(common_sig.sum()),
      "common_genus_min_unadjusted_p":float(
         df.within_publication_shared_genus_t_p.min()),
      "common_genus_min_BH_q":float(df.independent_BH_q_shared_genus_study.min()),
      "qualified_positive_pollen_grade12_BH_supported":int((
          df.loc[df.evidence_scope.eq("positive_pollen_evidence_grade12"),
                 "independent_BH_q_source_published_study"]<.05).sum()),
      "positive_source_study_gamma_not_observed_local_flower_preference":True,
      "source_unrecorded_visitors_not_negatives":True,
      "effect_changes_due_to_common_genus_restriction_are_not_causal_mediation":True,
      "not_an_independent_new_validation_dataset":True,
      "qualifying_highlighted_colour_patterns":selected.to_dict("records"),
    }
    return df,selected,summary

def run(src,out):
    dest=Path(out)
    require(not dest.exists(),"Do not overwrite signed ecological re-analysis")
    d,parent=read(src)
    full,chosen,summary=audit(d)
    dest.mkdir(parents=True)
    a=dest/"all_90_study_colour_panels_sharedgenus_BH_audit.csv"
    b=dest/"ten_source_visitor_colour_associations_and_genus_controls.csv"
    full.to_csv(a,index=False)
    chosen.to_csv(b,index=False)
    c=dest/"standalone_positive_guild_evidence_inferential_gate.json"
    c.write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    summary["output_SHA256"]={p.name:sha(p) for p in (a,b,c)}
    (dest/"receipt.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    print("MOIS_INDEPENDENT_DOPI_SAME_STUDY_GENUS_FDR_REAL_SUCCESS",
          json.dumps({
           "source_positives":summary["photo9_BH_supported_publication_positives"],
           "same_genus_supported":summary["same_study_shared_genus_BH_supported"],
           "pollen_evidence_only_supported":
              summary["qualified_positive_pollen_grade12_BH_supported"],
           "min_genus_q":summary["common_genus_min_BH_q"],
           "comparisons":[{
             "guild_a":r.guild_A_positive_plants,
             "guild_b":r.guild_B_positive_plants,"hue":r.photo9_hue,
             "source_photo_hue_percent_points":
                float(r.positive_guildA_minus_guildB_photo_hue_fraction*100),
             "source_BH_q":float(r.independent_BH_q_source_published_study),
             "shared_genus_percent_points":float(
                r.within_publication_shared_genus_equal_effect*100),
             "shared_genus_BH_q":float(r.independent_BH_q_shared_genus_study),
             "n_studies":int(r.n_same_published_studies_both_positive_guilds),
             "n_shared_genus_studies":int(r.n_same_study_shared_genus_panels)}
             for r in chosen.itertuples(index=False)]
          },allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",required=True)
    p.add_argument("--out",required=True)
    args=p.parse_args()
    run(args.source,args.out)
